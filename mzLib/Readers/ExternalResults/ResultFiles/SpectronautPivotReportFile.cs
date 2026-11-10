using System.Globalization;
using System.Text.RegularExpressions;
using CsvHelper;
using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// A Spectronaut Run Pivot report: wide format, one row per row label (usually a protein group) and one column per
    /// run and cell value, headed <c>[&lt;index&gt;] &lt;run file name&gt;.&lt;column&gt;</c>, for example
    /// <c>[1] 22_0505_002_SSL2_ND_Non_23_01.wiff.PG.Quantity</c>.
    /// <para>
    /// Detected by its header: no R.FileName, and at least one column of that form. It may be <c>.tsv</c>, <c>.xls</c>
    /// (tab text) or <c>.csv</c>; the delimiter is the one its first line uses. The run name is split from the column at
    /// the last <c>.PG.</c>, <c>.PEP.</c>, <c>.EG.</c>, <c>.FG.</c> or <c>.F.</c>, because run file names carry dots
    /// of their own.
    /// </para>
    /// </summary>
    public class SpectronautPivotReportFile : ResultFile<SpectronautPivotRow>
    {
        public override SupportedFileType FileType => SupportedFileType.SpectronautPivotReport;
        public override Software Software { get; set; }

        /// <summary>The report's columns, read with its rows; null before loading.</summary>
        public SpectronautPivotLayout? Layout { get; private set; }

        public SpectronautPivotReportFile(string filePath) : base(filePath, Software.Spectronaut) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public SpectronautPivotReportFile() : base()
        {
            Software = Software.Spectronaut;
        }

        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="MzLibException">
        /// The file is empty or is not a Run Pivot report, one run index carries two file names, or a row has a different
        /// number of fields from the header; the message names the file.
        /// </exception>
        public override void LoadResults()
        {
            const string what = "Spectronaut Run Pivot report";
            if (!File.Exists(FilePath))
                throw new FileNotFoundException($"{what} not found: '{FilePath}'", FilePath);

            try
            {
                string[] header = SpectronautReports.ReadHeader(FilePath, out char delimiter);
                if (header.Length == 0)
                    throw new MzLibException($"{what} '{FilePath}' is empty.");
                SpectronautGrain grain = SpectronautReports.Sniff(header);
                if (grain != SpectronautGrain.Pivot)
                    throw new MzLibException($"'{FilePath}' is not a {what}: its header is that of a {SpectronautReports.Describe(grain)}.");

                var layout = ReadLayout(header);
                var rows = new List<SpectronautPivotRow>();
                using var stream = new FileStream(FilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
                using var reader = new StreamReader(stream);
                using var parser = new CsvParser(reader, SpectronautReports.Configuration(delimiter));
                parser.Read(); // the header
                while (parser.Read())
                {
                    string[] fields = parser.Record!;
                    if (fields.Length == 1 && fields[0].Length == 0)
                        continue;
                    if (fields.Length != header.Length)
                        throw new MzLibException($"{what} '{FilePath}' row {parser.Row} has {fields.Length} fields, but its header has {header.Length}.");
                    rows.Add(new SpectronautPivotRow(layout, fields.Select(field => field.Trim()).ToArray()));
                }

                Layout = layout;
                Results = rows;
            }
            catch (Exception e) when (e is not MzLibException)
            {
                throw new MzLibException($"Could not read {what} '{FilePath}': {e.Message}", e);
            }

            SpectronautPivotLayout ReadLayout(string[] header)
            {
                var labels = new List<string>();
                var labelIndex = new Dictionary<string, int>(StringComparer.Ordinal);
                var runs = new SortedDictionary<int, string>();
                var cellColumns = new List<string>();
                var cellIndex = new Dictionary<(int, string), int>();
                for (int i = 0; i < header.Length; i++)
                {
                    Match match = SpectronautReports.PivotColumn.Match(header[i]);
                    if (!match.Success)
                    {
                        labels.Add(header[i]);
                        labelIndex.TryAdd(header[i], i);
                        continue;
                    }
                    int runIndex = int.Parse(match.Groups[1].Value, CultureInfo.InvariantCulture);
                    string runName = match.Groups[2].Value;
                    string column = match.Groups[3].Value;
                    if (runs.TryGetValue(runIndex, out string? known) && known != runName)
                        throw new MzLibException($"{what} '{FilePath}' names run [{runIndex}] both '{known}' and '{runName}'.");
                    runs[runIndex] = runName;
                    if (!cellColumns.Contains(column))
                        cellColumns.Add(column);
                    cellIndex[(runIndex, column)] = i;
                }
                return new SpectronautPivotLayout(labels, runs.Select(run => new SpectronautPivotRun(run.Key, run.Value)).ToList(),
                    cellColumns, labelIndex, cellIndex);
            }
        }

        /// <summary>
        /// The report in long form: one cell per row, run and cell-value column, row by row, then run by run in index
        /// order, then column by column in header order.
        /// </summary>
        public IEnumerable<SpectronautPivotCell> ToLong()
        {
            var rows = Results; // loads, and so sets Layout
            if (Layout is null)
                yield break;
            foreach (var row in rows)
                foreach (var run in Layout.Runs)
                    foreach (var column in Layout.CellColumns)
                        if (Layout.CellIndex.ContainsKey((run.Index, column)))
                            yield return new SpectronautPivotCell(row, run, column);
        }

        /// <exception cref="NotSupportedException">Always.</exception>
        public override void WriteResults(string outputPath) =>
            throw new NotSupportedException("Writing a Spectronaut report is not supported; SpectronautPivotReportFile only reads it.");
    }
}
