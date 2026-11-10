using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// A Spectronaut protein-group report: a Normal (long) report whose deepest level is the protein group, one row per
    /// protein group per run. This is what deposits most often hold. Submitters put it under MassIVE's <c>quant/</c>,
    /// named <c>&lt;date&gt;_&lt;time&gt;_&lt;analysis&gt;_Report.tsv</c>, or <c>.xls</c> (tab text) from older versions.
    /// <para>
    /// Detected by its header, whatever its name: R.FileName and PG.Quantity with PG.ProteinGroups or PG.ProteinAccessions,
    /// and no precursor (EG.*, FG.*), PTM site or TRIC columns. Read with the delimiter its first line uses.
    /// </para>
    /// <para>
    /// There is no q-value filter to apply here: Spectronaut wrote the rows that passed the cutoffs set in its analysis
    /// settings (the report's <c>.setup.txt</c> sidecar). The quantity is the group quantity those settings define, already
    /// rolled up and normalized.
    /// </para>
    /// </summary>
    public class SpectronautProteinGroupReportFile : ResultFile<SpectronautProteinGroup>
    {
        public override SupportedFileType FileType => SupportedFileType.SpectronautProteinGroupReport;
        public override Software Software { get; set; }

        public SpectronautProteinGroupReportFile(string filePath) : base(filePath, Software.Spectronaut) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public SpectronautProteinGroupReportFile() : base()
        {
            Software = Software.Spectronaut;
        }

        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="MzLibException">
        /// The file is empty, is another Spectronaut grain or another tool's file, lacks a core column, or holds a value
        /// that cannot be read; the message names the file.
        /// </exception>
        public override void LoadResults()
        {
            var pool = new SpectronautReports.StringPool();
            Results = SpectronautReports.Read<SpectronautProteinGroup>(FilePath, SpectronautGrain.ProteinGroup, MissingCore,
                records => records.Select(record =>
                {
                    record.RunName = pool.Get(record.RunName)!;
                    record.ConditionText = pool.Get(record.ConditionText);
                    _ = record.Quantity; // a quantity that is neither a number nor a Spectronaut token fails the load, not a later read
                    return record;
                }).ToList(), out _);
        }

        private static string? MissingCore(IReadOnlySet<string> header)
        {
            if (!header.Contains("PG.ProteinGroups") && !header.Contains("PG.ProteinAccessions"))
                return "PG.ProteinGroups or PG.ProteinAccessions";
            return null; // R.FileName and PG.Quantity are already required by the grain
        }

        /// <exception cref="NotSupportedException">Always.</exception>
        public override void WriteResults(string outputPath) =>
            throw new NotSupportedException("Writing a Spectronaut report is not supported; SpectronautProteinGroupReportFile only reads it.");
    }
}
