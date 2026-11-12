using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Quantification;
using Readers;

namespace Test.FileReadingTests.ExternalFileReading
{
    /// <summary>
    /// Spectronaut_TRICReport.tsv is cut from the 595 MB TRIC export of MSV000090777 (SCIEX ZenoTOF DIA, 120 min). It holds
    /// whole transition groups, every rank row of each, from three runs: the file's head (run 0), and two byte ranges at
    /// 3/16 and 5/16 of the file (runs 1 and 2, where decoys are interleaved). That gives 63 candidate rows: 16 run and
    /// precursor keys, each with one rank-1 peak group, 4 of them decoys.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestSpectronautTricReportFile
    {
        private static string FixturePath =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes", "Spectronaut_TRICReport.tsv");

        private static List<SpectraFileInfo> SpectraFiles => new()
        {
            new SpectraFileInfo(@"D:\Data\22_0505_002_SSL2_ND_Non_23_01.wiff", "Non-injured", 0, 0, 0),
            new SpectraFileInfo(@"D:\Data\22_0505_003_SSL2_ND_Inj_27_01.wiff", "Injured", 0, 0, 0),
            new SpectraFileInfo(@"D:\Data\22_0505_004_SSL2_ND_Non_24_01.wiff", "Non-injured", 1, 0, 0),
        };

        [Test]
        public void KeepsTheRankOnePeakGroupOfEachRunAndPrecursor()
        {
            var file = new SpectronautTricReportFile(FixturePath);
            Assert.That(file.Count(), Is.EqualTo(16));
            Assert.That(file.KeysWithoutRankOne, Is.EqualTo(0));
            Assert.That(file.All(p => p.Rank == 1));
            Assert.That(file.GroupBy(p => p.RunId).ToDictionary(g => g.Key, g => g.Count()),
                Is.EqualTo(new Dictionary<string, int> { ["0"] = 6, ["1"] = 5, ["2"] = 5 }));
            Assert.That(file.Count(p => p.IsDecoy), Is.EqualTo(4));
            Assert.That(file.FileType, Is.EqualTo(SupportedFileType.SpectronautTricReport));
            Assert.That(file.Software, Is.EqualTo(Software.Spectronaut));
        }

        [Test]
        public void FirstRecordFieldByField()
        {
            SpectronautTricPrecursor first = new SpectronautTricReportFile(FixturePath).First(p => p.PrecursorId == "_YTGVPDR_.2");
            Assert.That(first.TransitionGroupId, Is.EqualTo("1409187921_0"));
            Assert.That(first.FullSequence, Is.EqualTo("_YTGVPDR_"));
            Assert.That(first.BaseSequence, Is.EqualTo("YTGVPDR"));
            Assert.That(first.PrecursorMz, Is.EqualTo(404.2034));
            Assert.That(first.ChargeState, Is.EqualTo(2));
            Assert.That(first.MonoisotopicMass, Is.EqualTo(404.2034 * 2 - 2 * 1.007276466812).Within(1e-6));
            Assert.That(first.ProteinGroupInfos.Single().proteinAccessions, Is.EqualTo("A0A075B5N7"));
            Assert.That(first.SpectraFileName, Is.EqualTo("22_0505_002_SSL2_ND_Non_23_01.wiff"));
            Assert.That(first.FileName, Is.EqualTo("22_0505_002_SSL2_ND_Non_23_01"));
            Assert.That(first.IsDecoy, Is.False);
            Assert.That(first.RetentionTime, Is.EqualTo(8.457559));
            Assert.That(first.MScore, Is.EqualTo(0.028300075526097448));
            Assert.That(first.DScore, Is.EqualTo(-23.467194553840834));
            Assert.That(first.Intensity, Is.EqualTo(3.6124958992004395));
        }

        [Test]
        public void ShuffledRowsGiveTheSameRecords()
        {
            string shuffled = Path.Combine(TestContext.CurrentContext.TestDirectory, $"shuffled_{Guid.NewGuid():N}_TRICReport.tsv");
            var lines = File.ReadAllLines(FixturePath);
            var random = new Random(42);
            File.WriteAllLines(shuffled, new[] { lines[0] }.Concat(lines.Skip(1).OrderBy(_ => random.Next())));
            try
            {
                var expected = new SpectronautTricReportFile(FixturePath).Select(p => (p.RunId, p.PrecursorId, p.ApexRetentionTime)).OrderBy(k => k).ToList();
                var actual = new SpectronautTricReportFile(shuffled).Select(p => (p.RunId, p.PrecursorId, p.ApexRetentionTime)).OrderBy(k => k).ToList();
                Assert.That(actual, Is.EqualTo(expected));
            }
            finally
            {
                File.Delete(shuffled);
            }
        }

        [Test]
        public void TwoRankOneRowsForOneKeyAreRefused()
        {
            string doubled = Path.Combine(TestContext.CurrentContext.TestDirectory, $"doubled_{Guid.NewGuid():N}_TRICReport.tsv");
            var lines = File.ReadAllLines(FixturePath).ToList();
            lines.Add(lines.First(l => l.Contains("\t_YTGVPDR_.2\t") && l.Split('\t')[19] == "1"));
            File.WriteAllLines(doubled, lines);
            try
            {
                var e = Assert.Throws<MzLibException>(() => new SpectronautTricReportFile(doubled).LoadResults());
                Assert.That(e!.Message, Does.Contain("_YTGVPDR_.2"));
            }
            finally
            {
                File.Delete(doubled);
            }
        }

        [Test]
        public void MakeIdentificationsCarriesEveryRankOneRecordToFlashLfq()
        {
            var file = new SpectronautTricReportFile(FixturePath);
            List<Identification> identifications = MzLibExtensions.MakeIdentifications(file, SpectraFiles);

            Assert.That(identifications, Has.Count.EqualTo(16));
            Assert.That(identifications.Count(id => id.IsDecoy), Is.EqualTo(4));

            Identification first = identifications.First(id => id.ModifiedSequence == "_YTGVPDR_");
            Assert.That(first.BaseSequence, Is.EqualTo("YTGVPDR"));
            Assert.That(first.Ms2RetentionTimeInMinutes, Is.EqualTo(8.457559));
            Assert.That(first.PrecursorChargeState, Is.EqualTo(2));
            Assert.That(first.FileInfo.FullFilePathWithExtension, Is.EqualTo(@"D:\Data\22_0505_002_SSL2_ND_Non_23_01.wiff"));
            Assert.That(first.QValue, Is.EqualTo(0.028300075526097448)); // m_score, not left at 0
            Assert.That(first.PsmScore, Is.EqualTo(-23.467194553840834)); // d_score
        }

        /// <summary>
        /// The TRIC export is not filtered on m_score, so FlashLFQ's confidence rule does the filtering: 10 targets and 1
        /// decoy pass below 1%.
        /// </summary>
        [Test]
        public void QuantifiedIdentificationsAreFilteredOnMScore()
        {
            var file = new SpectronautTricReportFile(FixturePath);
            var kept = MzLibExtensions.MakeQuantifiedIdentifications(file, SpectraFiles, out QuantifiedPsmTier? tier);
            Assert.That(tier, Is.EqualTo(QuantifiedPsmTier.QValue));
            Assert.That(kept.Count(id => !id.IsDecoy), Is.EqualTo(10));
            Assert.That(kept.Count(id => id.IsDecoy), Is.EqualTo(1));
            Assert.That(kept.All(id => id.QValue < 0.01));
        }

        [Test]
        public void RunsMapToSpectraFilesByStem()
        {
            var map = new SpectronautTricReportFile(FixturePath).FileNameToFilePath(SpectraFiles.Select(f => f.FullFilePathWithExtension).ToList());
            Assert.That(map, Has.Count.EqualTo(3));
            Assert.That(map["22_0505_003_SSL2_ND_Inj_27_01"], Is.EqualTo(@"D:\Data\22_0505_003_SSL2_ND_Inj_27_01.wiff"));
        }

        [Test]
        public void FactoryBuiltFileReportsSpectronautAndWritingIsNotSupported()
        {
            Assert.That(new SpectronautTricReportFile().Software, Is.EqualTo(Software.Spectronaut));
            Assert.Throws<NotSupportedException>(() => new SpectronautTricReportFile(FixturePath).WriteResults("x_TRICReport.tsv"));
            Assert.Throws<FileNotFoundException>(() => new SpectronautTricReportFile("missing_TRICReport.tsv").LoadResults());
        }
    }
}
