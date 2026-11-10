using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// Spectronaut's TRIC export (<c>&lt;analysis&gt;_TRICReport.tsv</c>): the OpenSWATH feature table TRIC aligns, written
    /// by Spectronaut. It carries up to four peak-group candidates per run and precursor, ranked; one record is kept per
    /// run and precursor, the rank-1 candidate.
    /// <para>
    /// This is the Spectronaut export that can feed FlashLFQ: retention time in minutes, precursor m/z, charge, a decoy
    /// flag and a q-value (<c>m_score</c>). Detected by its header (transition_group_id, PrecursorID, m_score and RT) as well as
    /// its name. The file is sorted by run and then m/z, so the rank-1 rows are chosen while streaming. A 595 MB export in
    /// MSV000090777 is the size to expect.
    /// </para>
    /// </summary>
    public class SpectronautTricReportFile : ResultFile<SpectronautTricPrecursor>, IQuantifiableResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.SpectronautTricReport;
        public override Software Software { get; set; }

        /// <summary>
        /// How many run-and-precursor keys had candidate rows but no rank-1 row. Those keys are left out. A TRIC export
        /// writes a rank 1 for every key, so this is 0 unless the file was cut.
        /// </summary>
        public int KeysWithoutRankOne { get; private set; }

        public SpectronautTricReportFile(string filePath) : base(filePath, Software.Spectronaut) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public SpectronautTricReportFile() : base()
        {
            Software = Software.Spectronaut;
        }

        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="MzLibException">
        /// The file is empty, is not a TRIC export, holds two rank-1 rows for one run and precursor, or holds a value that
        /// cannot be read; the message names the file.
        /// </exception>
        public override void LoadResults()
        {
            var pool = new SpectronautReports.StringPool();
            Results = SpectronautReports.Read<SpectronautTricPrecursor>(FilePath, SpectronautGrain.Tric,
                header => header.Contains("filename") ? null : "filename",
                records =>
                {
                    var rankOne = new Dictionary<(string, string), SpectronautTricPrecursor>();
                    var allKeys = new HashSet<(string, string)>();
                    foreach (var record in records)
                    {
                        var key = (record.RunId, record.PrecursorId);
                        allKeys.Add(key);
                        if (record.Rank != 1)
                            continue;
                        if (!rankOne.TryAdd(key, record))
                            throw new MzLibException($"Spectronaut TRIC report '{FilePath}' has two rank-1 peak groups for precursor " +
                                $"{record.PrecursorId} in run {record.RunId}.");
                        record.SpectraFileName = pool.Get(record.SpectraFileName)!;
                        record.RunId = pool.Get(record.RunId)!;
                        record.ConditionText = pool.Get(record.ConditionText);
                        record.ProteinName = pool.Get(record.ProteinName);
                    }
                    KeysWithoutRankOne = allKeys.Count - rankOne.Count;
                    return rankOne.Values.ToList();
                }, out _);
        }

        public IEnumerable<IQuantifiableRecord> GetQuantifiableResults() => Results;

        /// <summary>Each run (its spectra file name without extension) mapped to a spectra file, by <see cref="RunFileMatcher.Match"/>'s rules.</summary>
        public Dictionary<string, string> FileNameToFilePath(List<string> fullFilePath) =>
            RunFileMatcher.Match(Results.Select(record => record.FileName), fullFilePath);

        /// <exception cref="NotSupportedException">Always.</exception>
        public override void WriteResults(string outputPath) =>
            throw new NotSupportedException("Writing a Spectronaut TRIC report is not supported; SpectronautTricReportFile only reads it.");
    }
}
