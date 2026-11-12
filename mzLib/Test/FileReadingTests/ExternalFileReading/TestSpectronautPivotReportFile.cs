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
    /// Run Pivot (wide) reports from MSV000090777:
    /// <list type="bullet">
    /// <item>Spectronaut_PG_Pivot_Report.tsv: 19 rows of <c>_Report_Protein_Quant_Pivot.xls</c> (tab text, renamed). One
    /// row label, PG.ProteinGroups, and PG.Quantity for 8 runs. The rows include cells written <c>Filtered</c> and
    /// <c>NaN</c>.</item>
    /// <item>Spectronaut_custom_Pivot_Report.csv: 10 rows of <c>_Report_Birgit.csv</c>, a custom schema written comma-separated
    /// with 8 row labels and three cell values per run.</item>
    /// </list>
    /// Run names carry the <c>.wiff</c> extension, so the run is split from the column at the last <c>.PG.</c>.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestSpectronautPivotReportFile
    {
        private static string Fixture(string name) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes", name);

        [Test]
        public void ReadsTheLayoutAndRows()
        {
            var file = new SpectronautPivotReportFile(Fixture("Spectronaut_PG_Pivot_Report.tsv"));
            Assert.That(file.Count(), Is.EqualTo(19));
            var layout = file.Layout!;
            Assert.That(layout.RowLabelColumns, Is.EqualTo(new[] { "PG.ProteinGroups" }));
            Assert.That(layout.CellColumns, Is.EqualTo(new[] { "PG.Quantity" }));
            Assert.That(layout.Runs.Select(r => r.Index), Is.EqualTo(Enumerable.Range(1, 8)));
            Assert.That(layout.Runs[0].FileName, Is.EqualTo("22_0505_002_SSL2_ND_Non_23_01.wiff"));
            Assert.That(layout.Runs[0].Stem, Is.EqualTo("22_0505_002_SSL2_ND_Non_23_01"));

            SpectronautPivotRow first = file.First();
            Assert.That(first.Label("PG.ProteinGroups"), Is.EqualTo("A0A023T672"));
            Assert.That(first.Value(1, "PG.Quantity"), Is.EqualTo(6.992893218994141));
            Assert.That(first.Value(2, "PG.Quantity"), Is.EqualTo(19.109498977661133));
            Assert.That(first.State(1, "PG.Quantity"), Is.EqualTo(SpectronautValueState.Value));
            Assert.That(first.Label("PG.Genes"), Is.Null);
            Assert.That(first.Text(9, "PG.Quantity"), Is.Null);
        }

        [Test]
        public void FilteredAndNaNAreDistinctStates()
        {
            var file = new SpectronautPivotReportFile(Fixture("Spectronaut_PG_Pivot_Report.tsv"));
            var cells = file.ToLong().ToList();
            Assert.That(cells, Has.Count.EqualTo(19 * 8));
            Assert.That(cells.Count(c => c.State == SpectronautValueState.Filtered), Is.EqualTo(12));
            Assert.That(cells.Count(c => c.State == SpectronautValueState.NotMeasured), Is.EqualTo(7));
            Assert.That(cells.Where(c => c.State != SpectronautValueState.Value).All(c => c.Value is null));

            var filtered = file.First(r => r.Label("PG.ProteinGroups") == "A0A075B5N7");
            Assert.That(filtered.State(3, "PG.Quantity"), Is.EqualTo(SpectronautValueState.Filtered));
        }

        [Test]
        public void TheCustomCommaSeparatedPivotHasThreeCellValuesPerRun()
        {
            var file = new SpectronautPivotReportFile(Fixture("Spectronaut_custom_Pivot_Report.csv"));
            Assert.That(file.Count(), Is.EqualTo(10));
            var layout = file.Layout!;
            Assert.That(layout.RowLabelColumns, Has.Count.EqualTo(8));
            Assert.That(layout.CellColumns, Is.EqualTo(new[]
                { "PG.NrOfPrecursorsIdentified", "PG.NrOfPrecursorsUsedForQuantification", "PG.Quantity" }));
            Assert.That(layout.Runs, Has.Count.EqualTo(8));

            SpectronautPivotRow first = file.First();
            Assert.That(first.Label("PG.Genes"), Is.EqualTo("RBM8"));
            Assert.That(first.Label("PG.CellularComponent"), Is.EqualTo(""));
            Assert.That(first.Value(1, "PG.NrOfPrecursorsIdentified"), Is.EqualTo(3));
            Assert.That(first.Value(1, "PG.Quantity"), Is.EqualTo(6.992893218994141));
            Assert.That(file.ToLong().Count(), Is.EqualTo(10 * 8 * 3));
        }

        [Test]
        public void ARowOfTheWrongWidthIsRefused()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"ragged_{Guid.NewGuid():N}_Report.tsv");
            File.WriteAllText(path, "PG.ProteinGroups\t[1] a.raw.PG.Quantity\t[2] b.raw.PG.Quantity\r\nP1\t1\r\n");
            try
            {
                Assert.That(Assert.Throws<MzLibException>(() => new SpectronautPivotReportFile(path).LoadResults())!.Message,
                    Does.Contain("fields"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void ARunIndexWithTwoNamesIsRefused()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"twonames_{Guid.NewGuid():N}_Report.tsv");
            File.WriteAllText(path, "PG.ProteinGroups\t[1] a.raw.PG.Quantity\t[1] b.raw.PG.Qvalue\r\nP1\t1\t0\r\n");
            try
            {
                Assert.That(Assert.Throws<MzLibException>(() => new SpectronautPivotReportFile(path).LoadResults())!.Message,
                    Does.Contain("names run [1]"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void FactoryBuiltFileReportsSpectronautAndWritingIsNotSupported()
        {
            Assert.That(new SpectronautPivotReportFile().Software, Is.EqualTo(Software.Spectronaut));
            Assert.That(new SpectronautPivotReportFile().ToLong(), Is.Empty);
            Assert.Throws<NotSupportedException>(() => new SpectronautPivotReportFile().WriteResults("x_Report.tsv"));
            Assert.Throws<FileNotFoundException>(() => new SpectronautPivotReportFile("missing_Report.tsv").LoadResults());
        }
    }
}
