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
    /// Spectronaut19_PG_Report.tsv holds every row of 8 protein groups from the Spectronaut 19.2 protein-group report of
    /// MSV000096226 (mouse brain, 10 runs). The groups were chosen for the cases the reader has to get right:
    /// <list type="bullet">
    /// <item>a plain group;</item>
    /// <item>NaN in every run (groups are NaN in all runs or none);</item>
    /// <item>blank genes and blank existence;</item>
    /// <item>gene and existence lists collapsed to one value;</item>
    /// <item>a duplicate id;</item>
    /// <item>eight members with no genes.</item>
    /// </list>
    /// Spectronaut_PG_Report.xls is the 10-column variant from MSV000094212, tab text under an .xls name.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestSpectronautProteinGroupReportFile
    {
        private static string Fixture(string name) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes", name);

        private static SpectronautProteinGroupReportFile Report => new(Fixture("Spectronaut19_PG_Report.tsv"));

        [Test]
        public void LoadsEveryRowWithFileTypeAndSoftware()
        {
            var file = Report;
            Assert.That(file.Count(), Is.EqualTo(80));
            Assert.That(file.FileType, Is.EqualTo(SupportedFileType.SpectronautProteinGroupReport));
            Assert.That(file.Software, Is.EqualTo(Software.Spectronaut));
            Assert.That(file.Select(r => r.RunName).Distinct().Count(), Is.EqualTo(10));
        }

        [Test]
        public void FactoryBuiltFileReportsSpectronaut()
        {
            Assert.That(new SpectronautProteinGroupReportFile().Software, Is.EqualTo(Software.Spectronaut));
            var read = FileReader.ReadResultFile(Fixture("Spectronaut19_PG_Report.tsv"));
            Assert.That(read, Is.TypeOf<SpectronautProteinGroupReportFile>());
            Assert.That(((SpectronautProteinGroupReportFile)read).Software, Is.EqualTo(Software.Spectronaut));
        }

        [Test]
        public void FirstRowFieldByField()
        {
            SpectronautProteinGroup first = Report.First();
            Assert.That(first.ConditionText, Is.EqualTo("Old"));
            Assert.That(first.Condition, Is.EqualTo("Old"));
            Assert.That(first.RunName, Is.EqualTo("240418_AM_O1"));
            Assert.That(first.Replicate, Is.EqualTo(1));
            Assert.That(first.Genes, Is.EqualTo("Cul4b"));
            Assert.That(first.ProteinAccessions, Is.EqualTo("A2A432"));
            Assert.That(first.ProteinDescriptions, Is.EqualTo("Cullin-4B"));
            Assert.That(first.ProteinGroups, Is.EqualTo("A2A432"));
            Assert.That(first.ProteinNames, Is.EqualTo("CUL4B_MOUSE"));
            Assert.That(first.ProteinExistence, Is.EqualTo("1"));
            Assert.That(first.PrecursorsUsedForQuantification, Is.EqualTo(3));
            Assert.That(first.StrippedSequencesUsedForQuantification, Is.EqualTo(3));
            Assert.That(first.Quantity, Is.EqualTo(71528.109375));
            Assert.That(first.QuantityState, Is.EqualTo(SpectronautValueState.Value));
            Assert.That(first.Members, Is.EqualTo(new[] { "A2A432" }));
        }

        [Test]
        public void LastRowFieldByField()
        {
            SpectronautProteinGroup last = Report.Last();
            Assert.That(last.Condition, Is.EqualTo("RF"));
            Assert.That(last.RunName, Is.EqualTo("240418_AM_ODR_RF_4"));
            Assert.That(last.Replicate, Is.EqualTo(6));
            Assert.That(last.Genes, Is.Null.Or.Empty);
            Assert.That(last.Members, Has.Count.EqualTo(8));
            Assert.That(last.PerMember(last.ProteinExistence), Is.EqualTo(new[] { "1", "1", "1", "1", "1", "2", "2", "2" }));
            Assert.That(last.PerMember(last.ProteinNames)![7], Is.EqualTo("KV6AA_MOUSE"));
            Assert.That(last.PerMember(last.Genes), Is.Null);
            Assert.That(last.PrecursorsUsedForQuantification, Is.EqualTo(0));
            Assert.That(last.Quantity, Is.Null);
            Assert.That(last.QuantityState, Is.EqualTo(SpectronautValueState.NotMeasured));
        }

        [Test]
        public void NaNIsNotMeasuredAndAllOrNothingPerGroup()
        {
            var file = Report;
            Assert.That(file.Count(r => r.QuantityState == SpectronautValueState.NotMeasured), Is.EqualTo(50));
            var unquantified = file.Where(r => r.ProteinGroups == "A2A5Z6").ToList();
            Assert.That(unquantified, Has.Count.EqualTo(10));
            Assert.That(unquantified.All(r => r.Quantity is null && r.Genes == "Smurf2"));
        }

        [Test]
        public void ListsAreAttributedOnlyWhenTheyCanBe()
        {
            var file = Report;

            // A duplicate id in PG.ProteinGroups: three distinct members, and the deduplicated lists line up with them
            var duplicate = file.First(r => r.ProteinGroups == "A2A5Y0;Q61765;Q9D646;Q9D646");
            Assert.That(duplicate.Members, Is.EqualTo(new[] { "A2A5Y0", "Q61765", "Q9D646" }));
            Assert.That(duplicate.ProteinAccessions, Is.EqualTo("A2A5Y0;Q61765;Q9D646"));
            Assert.That(duplicate.PerMember(duplicate.Genes), Is.EqualTo(new string?[] { null, "Krt31", null }));
            Assert.That(duplicate.PerMember(duplicate.ProteinNames), Is.EqualTo(new[] { "A2A5Y0", "K1H1_MOUSE", "Q9D646" }));

            // One gene for two members: collapsed, so it applies to both
            var collapsedGene = file.First(r => r.ProteinGroups == "O55074;Q7TN79");
            Assert.That(collapsedGene.PerMember(collapsedGene.Genes), Is.EqualTo(new[] { "Akap7", "Akap7" }));

            // Two genes for two members, existence collapsed
            var twoGenes = file.First(r => r.ProteinGroups == "O08580;Q61539");
            Assert.That(twoGenes.PerMember(twoGenes.Genes), Is.EqualTo(new[] { "Esrra", "Esrrb" }));
            Assert.That(twoGenes.PerMember(twoGenes.ProteinExistence), Is.EqualTo(new[] { "1", "1" }));

            // Blank genes and blank existence read as absent
            var trypsin = file.First(r => r.ProteinGroups == "P00761");
            Assert.That(trypsin.PerMember(trypsin.Genes), Is.Null);
            Assert.That(trypsin.PerMember(trypsin.ProteinExistence), Is.Null);
            Assert.That(trypsin.Quantity, Is.EqualTo(7113235.5));
        }

        [Test]
        public void TheXlsVariantReadsAsTabText()
        {
            var file = new SpectronautProteinGroupReportFile(Fixture("Spectronaut_PG_Report.xls"));
            Assert.That(file.Count(), Is.EqualTo(18));
            Assert.That(file.First().ProteinGroups, Is.EqualTo("A1A5B6"));
            Assert.That(file.First().Quantity, Is.EqualTo(5512.521484375));
            Assert.That(file.Count(r => r.QuantityState == SpectronautValueState.NotMeasured), Is.EqualTo(3));
            Assert.That(file.All(r => r.PrecursorsUsedForQuantification is null)); // that variant does not export the column
        }

        [Test]
        public void AMissingFileThrowsFileNotFound()
        {
            var file = new SpectronautProteinGroupReportFile(Fixture("no_such_Report.tsv"));
            var e = Assert.Throws<FileNotFoundException>(() => file.LoadResults());
            Assert.That(e!.FileName, Does.EndWith("no_such_Report.tsv"));
        }

        [Test]
        public void ACorruptValueThrowsMzLibExceptionNamingTheFile()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"corrupt_{Guid.NewGuid():N}_Report.tsv");
            File.WriteAllText(path, "R.FileName\tPG.ProteinGroups\tPG.Quantity\r\nrun1\tP1\tnot-a-number\r\n");
            try
            {
                var e = Assert.Throws<MzLibException>(() => new SpectronautProteinGroupReportFile(path).LoadResults());
                Assert.That(e!.Message, Does.Contain(path));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void AnEmptyFileThrowsMzLibException()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"empty_{Guid.NewGuid():N}_Report.tsv");
            File.WriteAllText(path, "");
            try
            {
                Assert.That(Assert.Throws<MzLibException>(() => new SpectronautProteinGroupReportFile(path).LoadResults())!.Message,
                    Does.Contain("is empty"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void WritingIsNotSupported()
        {
            Assert.Throws<NotSupportedException>(() => Report.WriteResults("anything_Report.tsv"));
        }
    }
}
