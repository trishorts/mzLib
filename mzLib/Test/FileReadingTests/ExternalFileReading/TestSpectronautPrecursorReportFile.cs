using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.ExternalFileReading
{
    /// <summary>
    /// Fixtures, all from public MassIVE deposits:
    /// <list type="bullet">
    /// <item>SpectronautPrecursor_pep_Report.tsv: 60 rows of MSV000098115's <c>_Report_pep.tsv</c>, one run. It is the
    /// thin precursor export deposits hold: EG.PrecursorId, EG.Qvalue (not filtered, up to 0.89) and quantities, but no
    /// retention time, m/z, charge column, PG.ProteinGroups or PG.Genes. Modifications <c>[Oxidation (M)]</c>,
    /// <c>[Acetyl (Protein N-term)]</c> and the siteless <c>[Biotin]</c>.</item>
    /// <item>Spectronaut10_fragment_Report.csv: every fragment row of 5 precursors from a Spectronaut 10 comma-separated
    /// export in MSV000087085. Mass-delta modifications in EG.ModifiedSequence; light and heavy acetyl labels in FG.Id.</item>
    /// <item>Spectronaut10_quoted_Report.csv: the same kind of export re-saved by R, every string quoted, plus a
    /// submitter's <c>sample_id</c> column.</item>
    /// </list>
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestSpectronautPrecursorReportFile
    {
        private static string Fixture(string name) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes", name);

        [Test]
        public void ThePepExportReadsWithChargeFromTheId()
        {
            var file = new SpectronautPrecursorReportFile(Fixture("SpectronautPrecursor_pep_Report.tsv"));
            Assert.That(file.Count(), Is.EqualTo(60));
            Assert.That(file.Software, Is.EqualTo(Software.Spectronaut));

            SpectronautPrecursor first = file.First();
            Assert.That(first.RunName, Is.EqualTo("240320_AKS_SD374566_YN1"));
            Assert.That(first.Condition, Is.EqualTo("YN"));
            Assert.That(first.PrecursorId, Is.EqualTo("_QLLGEHLTAILQK_.3"));
            Assert.That(first.Charge, Is.Null); // no FG.Charge column
            Assert.That(first.PrecursorCharge, Is.EqualTo(3));
            Assert.That(first.ChargeState, Is.EqualTo(3));
            Assert.That(first.ModifiedSequence, Is.EqualTo("_QLLGEHLTAILQK_"));
            Assert.That(first.FullSequence, Is.EqualTo("_QLLGEHLTAILQK_"));
            Assert.That(first.BaseSequence, Is.EqualTo("QLLGEHLTAILQK"));
            Assert.That(first.QValue, Is.EqualTo(0.12917055793935892));
            Assert.That(first.ProteinGroupQValue, Is.EqualTo(0));
            Assert.That(first.ProteinAccessions, Is.EqualTo("A2A432"));
            Assert.That(first.ProteinGroups, Is.Null);
            Assert.That(first.Quantity, Is.EqualTo(113.43197631835938)); // EG.TotalQuantity (Settings): no FG.Quantity here
            Assert.That(first.ProteinGroupInfos.Single().proteinAccessions, Is.EqualTo("A2A432"));
            Assert.That(double.IsNaN(first.RetentionTime));
            Assert.That(double.IsNaN(first.MonoisotopicMass));
        }

        [Test]
        public void ModificationsAreKeptAsWritten()
        {
            var file = new SpectronautPrecursorReportFile(Fixture("SpectronautPrecursor_pep_Report.tsv"));
            var biotin = file.Single(p => p.PrecursorId == "_ILNSFSVM[Oxidation (M)]PSPK[Biotin]_.2");
            Assert.That(biotin.FullSequence, Is.EqualTo("_ILNSFSVM[Oxidation (M)]PSPK[Biotin]_"));
            Assert.That(biotin.BaseSequence, Is.EqualTo("ILNSFSVMPSPK"));

            var acetyl = file.Single(p => p.PrecursorId == "_[Acetyl (Protein N-term)]GDDSEWLK_.2");
            Assert.That(acetyl.BaseSequence, Is.EqualTo("GDDSEWLK"));

            var singlyCharged = file.Single(p => p.PrecursorId == "_ILNSFSVM[Oxidation (M)]PSPK_.1");
            Assert.That(singlyCharged.PrecursorCharge, Is.EqualTo(1));

            // Not filtered at 1%: the reader keeps every q-value it is given
            Assert.That(file.Max(p => p.QValue), Is.EqualTo(0.8943779126152751));
            Assert.That(file.Count(p => p.ProteinGroupQuantityText == "NaN"), Is.EqualTo(4));
        }

        [Test]
        public void FragmentRowsCollapseToPrecursors()
        {
            var file = new SpectronautPrecursorReportFile(Fixture("Spectronaut10_fragment_Report.csv"));
            Assert.That(file.Count(), Is.EqualTo(5)); // 30 fragment rows of 5 precursors

            SpectronautPrecursor first = file.First();
            Assert.That(first.RunName, Is.EqualTo("20170514_Aging_DIA_D6_AcK_19_1"));
            Assert.That(first.FragmentGroupId, Is.EqualTo("_K[Acetyl (K)]K[Acetyl (K)]LLAK[Acetyl D3 (K)]DLR_.2"));
            Assert.That(first.ModifiedSequenceText, Is.EqualTo("_K[+42]K[+42]LLAK[+42]DLR_"));
            Assert.That(first.StrippedSequence, Is.EqualTo("KKLLAKDLR")); // EG.StrippedSequence in Spectronaut 10
            Assert.That(first.Charge, Is.EqualTo(2));
            Assert.That(first.PrecursorMz, Is.EqualTo(607.388916015625));
            Assert.That(first.QValue, Is.EqualTo(0.00285647621940142));
            Assert.That(first.Decoy, Is.False);
            Assert.That(first.ProteinGroups, Is.EqualTo("A0JNU3"));
            Assert.That(first.MonoisotopicMass, Is.EqualTo(607.388916015625 * 2 - 2 * 1.007276466812).Within(1e-6));

            var last = file.Last();
            Assert.That(last.FragmentGroupId, Does.Contain("M[+O]"));
            Assert.That(last.ModifiedSequenceText, Is.EqualTo("_K[+42]K[+42]DLK[+42]M[+16]SVLISK[+42]LNDTE_"));
            Assert.That(last.PrecursorCharge, Is.EqualTo(3));
        }

        /// <summary>
        /// EG.ModifiedSequence, written as mass deltas, is the same for light and heavy acetyl label channels. The four
        /// precursors of KKLLAKDLR differ in m/z, so the sequence FlashLFQ keys on has to come from FG.Id.
        /// </summary>
        [Test]
        public void LabelChannelsKeepDistinctSequences()
        {
            var file = new SpectronautPrecursorReportFile(Fixture("Spectronaut10_fragment_Report.csv"));
            var channels = file.Where(p => p.BaseSequence == "KKLLAKDLR").ToList();
            Assert.That(channels, Has.Count.EqualTo(4));
            Assert.That(channels.Select(p => p.ModifiedSequenceText).Distinct().Count(), Is.EqualTo(1));
            Assert.That(channels.Select(p => p.FullSequence).Distinct().Count(), Is.EqualTo(4));
            Assert.That(channels.Select(p => p.PrecursorMz).Distinct().Count(), Is.EqualTo(3)); // two channels share an m/z
        }

        [Test]
        public void TheRQuotedCsvParsesAndIgnoresTheExtraColumn()
        {
            var file = new SpectronautPrecursorReportFile(Fixture("Spectronaut10_quoted_Report.csv"));
            Assert.That(file.Count(), Is.EqualTo(2));
            Assert.That(file.First().FragmentGroupId, Is.EqualTo("_K[Acetyl (K)]K[Acetyl (K)]LLAK[Acetyl D3 (K)]DLR_.2"));
            Assert.That(file.First().Decoy, Is.False);
        }

        [Test]
        public void FlashLfqIsToldWhatIsMissing()
        {
            var pep = new SpectronautPrecursorReportFile(Fixture("SpectronautPrecursor_pep_Report.tsv"));
            var e = Assert.Throws<MzLibException>(() => pep.GetQuantifiableResults());
            Assert.That(e!.Message, Does.Contain("retention time"));
            Assert.That(e.Message, Does.Contain("FG.PrecMz"));

            // Spectronaut 10 exported m/z but no retention time
            var fragments = new SpectronautPrecursorReportFile(Fixture("Spectronaut10_fragment_Report.csv"));
            e = Assert.Throws<MzLibException>(() => fragments.GetQuantifiableResults());
            Assert.That(e!.Message, Does.Contain("retention time"));
            Assert.That(e.Message, Does.Not.Contain("FG.PrecMz"));
        }

        [Test]
        public void RunsMapToSpectraFiles()
        {
            var file = new SpectronautPrecursorReportFile(Fixture("SpectronautPrecursor_pep_Report.tsv"));
            var map = file.FileNameToFilePath(new List<string> { @"D:\data\240320_AKS_SD374566_YN1.raw", @"D:\data\other.raw" });
            Assert.That(map, Has.Count.EqualTo(1));
            Assert.That(map["240320_AKS_SD374566_YN1"], Is.EqualTo(@"D:\data\240320_AKS_SD374566_YN1.raw"));
        }

        [Test]
        public void FactoryBuiltFileReportsSpectronautAndWritingIsNotSupported()
        {
            Assert.That(new SpectronautPrecursorReportFile().Software, Is.EqualTo(Software.Spectronaut));
            var file = FileReader.ReadQuantifiableResultFile(Fixture("SpectronautPrecursor_pep_Report.tsv"));
            Assert.That(file, Is.TypeOf<SpectronautPrecursorReportFile>());
            Assert.Throws<NotSupportedException>(() => ((SpectronautPrecursorReportFile)file).WriteResults("x_Report.tsv"));
        }
    }
}
