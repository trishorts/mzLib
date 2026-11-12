using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.ExternalFileReading
{
    /// <summary>
    /// Spectronaut_PTMSiteReport.tsv: 23 rows of MSV000090777's <c>_-PTMSiteReport.tsv</c>, one run, its condition never set
    /// (<c>Not Defined</c>). Carbamidomethyl (C) sites, plus one Oxidation (M) site.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestSpectronautPtmSiteReportFile
    {
        private static string FixturePath =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes", "Spectronaut_PTMSiteReport.tsv");

        [Test]
        public void ReadsEverySite()
        {
            var file = new SpectronautPtmSiteReportFile(FixturePath);
            Assert.That(file.Count(), Is.EqualTo(23));
            Assert.That(file.FileType, Is.EqualTo(SupportedFileType.SpectronautPtmSiteReport));
            Assert.That(file.All(site => site.ConditionText == "Not Defined" && site.Condition is null));

            SpectronautPtmSite first = file.First();
            Assert.That(first.RunName, Is.EqualTo("22_0505_003_SSL2_ND_Inj_27_01"));
            Assert.That(first.ProteinGroups, Is.EqualTo("A0A067XG53;O70589;O70589-3;O70589-4;O70589-5"));
            Assert.That(first.ProteinId, Is.EqualTo("A0A067XG53"));
            Assert.That(first.Group, Is.EqualTo("_DVKPHC[Carbamidomethyl (C)]VLLASK_.3"));
            Assert.That(first.ModificationTitle, Is.EqualTo("Carbamidomethyl (C)"));
            Assert.That(first.SiteLocation, Is.EqualTo(120));
            Assert.That(first.SiteAminoAcid, Is.EqualTo("C"));
            Assert.That(first.SiteProbability, Is.EqualTo(1));
            Assert.That(first.Multiplicity, Is.EqualTo(1));
            Assert.That(first.Quantity, Is.EqualTo(2.9687089920043945));

            var oxidation = file.Single(site => site.ModificationTitle == "Oxidation (M)");
            Assert.That(oxidation.ProteinId, Is.EqualTo("A0A087WSN6"));
            Assert.That(oxidation.SiteLocation, Is.EqualTo(1282));
        }

        [Test]
        public void FactoryBuiltFileReportsSpectronautAndWritingIsNotSupported()
        {
            Assert.That(new SpectronautPtmSiteReportFile().Software, Is.EqualTo(Software.Spectronaut));
            Assert.Throws<NotSupportedException>(() => new SpectronautPtmSiteReportFile(FixturePath).WriteResults("x_-PTMSiteReport.tsv"));
            Assert.Throws<FileNotFoundException>(() => new SpectronautPtmSiteReportFile("missing_-PTMSiteReport.tsv").LoadResults());
        }
    }
}
