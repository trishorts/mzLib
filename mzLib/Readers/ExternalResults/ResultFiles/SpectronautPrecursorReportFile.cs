using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// A Spectronaut Normal report whose deepest level is the precursor or the fragment: one record per precursor per
    /// run. A fragment-level report (Spectronaut 10 deposits run to 2.4 GB of comma-separated <c>.csv</c>) is reduced to
    /// its precursors while streaming, keeping the first row of each run and precursor.
    /// <para>
    /// Detected by its header: R.FileName with any of EG.PrecursorId, EG.ModifiedSequence or FG.Id, in <c>.tsv</c>,
    /// <c>.xls</c> (tab text) or <c>.csv</c>. Deposits name it <c>_Report_peptides.tsv</c>, <c>_Report_pep.tsv</c> or
    /// <c>_Report.csv</c>.
    /// </para>
    /// <para>
    /// FlashLFQ needs a retention time and a precursor mass, which many exports lack.
    /// <see cref="GetQuantifiableResults"/> throws, naming the missing columns, rather than hand FlashLFQ a record with
    /// no retention time or mass. Spectronaut's TRIC export (<see cref="SpectronautTricReportFile"/>) carries both.
    /// </para>
    /// </summary>
    public class SpectronautPrecursorReportFile : ResultFile<SpectronautPrecursor>, IQuantifiableResultFile
    {
        private string[] _header = Array.Empty<string>();

        public override SupportedFileType FileType => SupportedFileType.SpectronautPrecursorReport;
        public override Software Software { get; set; }

        public SpectronautPrecursorReportFile(string filePath) : base(filePath, Software.Spectronaut) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public SpectronautPrecursorReportFile() : base()
        {
            Software = Software.Spectronaut;
        }

        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="MzLibException">
        /// The file is empty, is another Spectronaut grain or another tool's file, or holds a value that cannot be read;
        /// the message names the file.
        /// </exception>
        public override void LoadResults()
        {
            var pool = new SpectronautReports.StringPool();
            Results = SpectronautReports.Read<SpectronautPrecursor>(FilePath, SpectronautGrain.Precursor, _ => null,
                records =>
                {
                    // Keep the first row of each run and precursor. A precursor-level report has one row per key already;
                    // a fragment-level report has one per fragment ion, and only the precursor columns are kept.
                    var seen = new HashSet<(string, string)>();
                    var kept = new List<SpectronautPrecursor>();
                    foreach (var record in records)
                    {
                        if (!seen.Add((record.RunName, record.PrecursorKey)))
                            continue;
                        record.RunName = pool.Get(record.RunName)!;
                        record.ConditionText = pool.Get(record.ConditionText);
                        record.ProteinGroups = pool.Get(record.ProteinGroups);
                        record.ProteinAccessions = pool.Get(record.ProteinAccessions);
                        record.Genes = pool.Get(record.Genes);
                        record.ProteinNames = pool.Get(record.ProteinNames);
                        record.ProteinDescriptions = pool.Get(record.ProteinDescriptions);
                        _ = record.Quantity; // fail the load, not a later read, on a quantity that is neither a number nor a token
                        kept.Add(record);
                    }
                    return kept;
                }, out _header);
        }

        /// <summary>
        /// The records, for FlashLFQ.
        /// </summary>
        /// <exception cref="MzLibException">
        /// The report has no retention-time column (EG.ApexRT, EG.MeanApexRT or EG.RTEmpirical) or no precursor m/z
        /// (FG.PrecMz), so FlashLFQ could not place its peaks. The message names what is missing. Re-export with those
        /// columns, or use the analysis's TRIC export.
        /// </exception>
        public IEnumerable<IQuantifiableRecord> GetQuantifiableResults()
        {
            var results = Results; // loads, and so reads the header
            var missing = new List<string>();
            if (!new[] { "EG.ApexRT", "EG.MeanApexRT", "EG.RTEmpirical" }.Any(_header.Contains))
                missing.Add("a retention time (EG.ApexRT, EG.MeanApexRT or EG.RTEmpirical)");
            if (!_header.Contains("FG.PrecMz"))
                missing.Add("a precursor m/z (FG.PrecMz)");
            if (missing.Count > 0)
                throw new MzLibException($"Spectronaut precursor report '{FilePath}' cannot be quantified by FlashLFQ: it has no " +
                    string.Join(" and no ", missing) + " column. Re-export it with those columns, or use the analysis's TRIC export.");
            return results;
        }

        /// <summary>Each run name (R.FileName) mapped to its spectra file, by <see cref="RunFileMatcher.Match"/>'s rules.</summary>
        public Dictionary<string, string> FileNameToFilePath(List<string> fullFilePath) =>
            RunFileMatcher.Match(Results.Select(record => record.FileName), fullFilePath);

        /// <exception cref="NotSupportedException">Always.</exception>
        public override void WriteResults(string outputPath) =>
            throw new NotSupportedException("Writing a Spectronaut report is not supported; SpectronautPrecursorReportFile only reads it.");
    }
}
