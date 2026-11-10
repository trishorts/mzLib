using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.ExternalFileReading
{
    /// <summary>
    /// What the Spectronaut readers share: telling the report grains apart by header in .tsv, .xls and .csv; Spectronaut's
    /// value tokens; and its PG.* lists, which are not always one entry per member. The fixtures are cut from public
    /// MassIVE deposits. Their provenance is in FileReadingTests\ExternalFileTypes\SPECTRONAUT_PROVENANCE.md.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestSpectronautReports
    {
        private static string Fixture(string name) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes", name);

        [Test]
        [TestCase("Spectronaut19_PG_Report.tsv", SupportedFileType.SpectronautProteinGroupReport)]
        [TestCase("Spectronaut_PG_Report.xls", SupportedFileType.SpectronautProteinGroupReport)]
        [TestCase("SpectronautPrecursor_pep_Report.tsv", SupportedFileType.SpectronautPrecursorReport)]
        [TestCase("Spectronaut10_fragment_Report.csv", SupportedFileType.SpectronautPrecursorReport)]
        [TestCase("Spectronaut10_quoted_Report.csv", SupportedFileType.SpectronautPrecursorReport)]
        [TestCase("Spectronaut_TRICReport.tsv", SupportedFileType.SpectronautTricReport)]
        [TestCase("Spectronaut_PG_Pivot_Report.tsv", SupportedFileType.SpectronautPivotReport)]
        [TestCase("Spectronaut_custom_Pivot_Report.csv", SupportedFileType.SpectronautPivotReport)]
        [TestCase("Spectronaut_PTMSiteReport.tsv", SupportedFileType.SpectronautPtmSiteReport)]
        [TestCase("Spectronaut19_ConditionSetup.tsv", SupportedFileType.SpectronautConditionSetup)]
        public void EachGrainIsDetectedInEachExtensionSpectronautWrites(string fixture, SupportedFileType expected)
        {
            Assert.That(Fixture(fixture).ParseFileType(), Is.EqualTo(expected));
            Assert.That(FileReader.ReadResultFile(Fixture(fixture)).FileType, Is.EqualTo(expected));
        }

        /// <summary>The grain comes from the header, so a report keeps its type when renamed (the TRIC and PTM site names
        /// are also suffixes, which is why those two are renamed to a plain report name here).</summary>
        [Test]
        [TestCase("Spectronaut_TRICReport.tsv", SupportedFileType.SpectronautTricReport)]
        [TestCase("Spectronaut_PTMSiteReport.tsv", SupportedFileType.SpectronautPtmSiteReport)]
        [TestCase("SpectronautPrecursor_pep_Report.tsv", SupportedFileType.SpectronautPrecursorReport)]
        [TestCase("Spectronaut_PG_Pivot_Report.tsv", SupportedFileType.SpectronautPivotReport)]
        public void RenamedReportsAreDetectedByHeader(string fixture, SupportedFileType expected)
        {
            string renamed = Path.Combine(TestContext.CurrentContext.TestDirectory, $"renamed_{Guid.NewGuid():N}.tsv");
            File.Copy(Fixture(fixture), renamed);
            try
            {
                Assert.That(renamed.ParseFileType(), Is.EqualTo(expected));
            }
            finally
            {
                File.Delete(renamed);
            }
        }

        [Test]
        public void AnXlsThatIsNotSpectronautIsRefused()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"not_spectronaut_{Guid.NewGuid():N}.xls");
            File.WriteAllText(path, "Name\tValue\r\nA\t1\r\n");
            try
            {
                var e = Assert.Throws<MzLibException>(() => path.ParseFileType());
                Assert.That(e!.Message, Is.EqualTo("Xls file type not supported"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void AMissingXlsOrCsvSaysNotSupported()
        {
            Assert.That(Assert.Throws<MzLibException>(() => "missing_taco.xls".ParseFileType())!.Message, Is.EqualTo("Xls file type not supported"));
            Assert.That(Assert.Throws<MzLibException>(() => "missing_taco.csv".ParseFileType())!.Message, Is.EqualTo("Csv file type not supported"));
        }

        [Test]
        public void DiaNnAndSpectronautReportsAreNotTakenForEachOther()
        {
            string diaNn = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes\DiaNn_LongFormat_report.tsv");
            Assert.That(diaNn.ParseFileType(), Is.EqualTo(SupportedFileType.DiaNnReport));
            Assert.That(Fixture("SpectronautPrecursor_pep_Report.tsv").ParseFileType(), Is.Not.EqualTo(SupportedFileType.DiaNnReport));
        }

        [Test]
        public void AReaderRefusesAReportOfAnotherGrain()
        {
            var file = new SpectronautProteinGroupReportFile(Fixture("SpectronautPrecursor_pep_Report.tsv"));
            var e = Assert.Throws<MzLibException>(() => file.LoadResults());
            Assert.That(e!.Message, Does.Contain("is not a Spectronaut protein-group report"));
            Assert.That(e.Message, Does.Contain("precursor report"));
        }

        [Test]
        public void ValueTokens()
        {
            Assert.That(SpectronautReports.ParseDouble("71528.109375"), Is.EqualTo(71528.109375));
            Assert.That(SpectronautReports.ParseDouble("5.6152671614283445E-08"), Is.EqualTo(5.6152671614283445E-08));
            Assert.That(SpectronautReports.ParseDouble("NaN"), Is.Null);
            Assert.That(SpectronautReports.ParseDouble("Filtered"), Is.Null);
            Assert.That(SpectronautReports.ParseDouble(""), Is.Null);
            Assert.Throws<FormatException>(() => SpectronautReports.ParseDouble("abc"));

            Assert.That(SpectronautReports.StateOf("1.5"), Is.EqualTo(SpectronautValueState.Value));
            Assert.That(SpectronautReports.StateOf("NaN"), Is.EqualTo(SpectronautValueState.NotMeasured));
            Assert.That(SpectronautReports.StateOf(""), Is.EqualTo(SpectronautValueState.NotMeasured));
            Assert.That(SpectronautReports.StateOf("Filtered"), Is.EqualTo(SpectronautValueState.Filtered));

            Assert.That(SpectronautReports.ParsePercent("22.8%"), Is.EqualTo(0.228).Within(1e-12));
            Assert.That(SpectronautReports.ParsePercent("0%"), Is.EqualTo(0));
            Assert.That(SpectronautReports.Stated("Not Defined"), Is.Null);
            Assert.That(SpectronautReports.Stated(" Old "), Is.EqualTo("Old"));
        }

        [Test]
        [SetCulture("de-DE")]
        public void NumbersAreReadInvariantlyUnderACommaDecimalCulture()
        {
            Assert.That(SpectronautReports.ParseDouble("71528.109375"), Is.EqualTo(71528.109375));
            Assert.That(SpectronautReports.ParsePercent("22.8%"), Is.EqualTo(0.228).Within(1e-12));
            var file = new SpectronautProteinGroupReportFile(Fixture("Spectronaut19_PG_Report.tsv"));
            Assert.That(file.First().Quantity, Is.EqualTo(71528.109375));
        }

        [Test]
        public void MembersDropDuplicateIdsAndSplitGroups()
        {
            Assert.That(SpectronautReports.Members("A2A5Y0;Q61765;Q9D646;Q9D646"), Is.EqualTo(new[] { "A2A5Y0", "Q61765", "Q9D646" }));
            Assert.That(SpectronautReports.Members("P1;P2|P3"), Is.EqualTo(new[] { "P1", "P2", "P3" }));
            Assert.That(SpectronautReports.Members(""), Is.Empty);
        }

        [Test]
        public void PerMemberAttributesOnlyWhatItCan()
        {
            Assert.That(SpectronautReports.PerMember(";Krt31;", 3), Is.EqualTo(new string?[] { null, "Krt31", null }));
            Assert.That(SpectronautReports.PerMember("Akap7", 2), Is.EqualTo(new string?[] { "Akap7", "Akap7" }));
            Assert.That(SpectronautReports.PerMember("a;b", 3), Is.Null);
            Assert.That(SpectronautReports.PerMember("", 3), Is.Null);
        }

        [Test]
        public void PrecursorIdsSplitIntoSequenceAndCharge()
        {
            Assert.That(SpectronautReports.SplitPrecursorId("_GPLATGGIK_.1", out int? charge), Is.EqualTo("_GPLATGGIK_"));
            Assert.That(charge, Is.EqualTo(1));
            Assert.That(SpectronautReports.SplitPrecursorId("_M[+15.995]PEPTIDE_.3", out charge), Is.EqualTo("_M[+15.995]PEPTIDE_"));
            Assert.That(charge, Is.EqualTo(3));
            Assert.That(SpectronautReports.SplitPrecursorId("_PEPTIDE_", out charge), Is.EqualTo("_PEPTIDE_"));
            Assert.That(charge, Is.Null);
            Assert.That(SpectronautReports.StripModifications("_[Acetyl (Protein N-term)]M[Oxidation (M)]PEPK[Biotin]_"), Is.EqualTo("MPEPK"));
        }
    }
}
