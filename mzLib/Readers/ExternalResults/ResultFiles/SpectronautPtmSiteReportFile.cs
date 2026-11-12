using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// A Spectronaut PTM site report: one row per PTM site per run, for analyses run with Spectronaut's PTM workflow.
    /// Detected by its header (R.FileName, PTM.SiteLocation and PTM.Quantity) whatever its name. Deposits name it
    /// <c>&lt;analysis&gt;_-PTMSiteReport.tsv</c>.
    /// </summary>
    public class SpectronautPtmSiteReportFile : ResultFile<SpectronautPtmSite>
    {
        public override SupportedFileType FileType => SupportedFileType.SpectronautPtmSiteReport;
        public override Software Software { get; set; }

        public SpectronautPtmSiteReportFile(string filePath) : base(filePath, Software.Spectronaut) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public SpectronautPtmSiteReportFile() : base()
        {
            Software = Software.Spectronaut;
        }

        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="MzLibException">
        /// The file is empty, is not a PTM site report, or holds a value that cannot be read; the message names the file.
        /// </exception>
        public override void LoadResults()
        {
            var pool = new SpectronautReports.StringPool();
            Results = SpectronautReports.Read<SpectronautPtmSite>(FilePath, SpectronautGrain.PtmSite, _ => null,
                records => records.Select(record =>
                {
                    record.RunName = pool.Get(record.RunName)!;
                    record.ConditionText = pool.Get(record.ConditionText);
                    record.ProteinGroups = pool.Get(record.ProteinGroups);
                    _ = record.Quantity; // fail the load, not a later read, on a quantity that is neither a number nor a token
                    return record;
                }).ToList(), out _);
        }

        /// <exception cref="NotSupportedException">Always.</exception>
        public override void WriteResults(string outputPath) =>
            throw new NotSupportedException("Writing a Spectronaut report is not supported; SpectronautPtmSiteReportFile only reads it.");
    }
}
