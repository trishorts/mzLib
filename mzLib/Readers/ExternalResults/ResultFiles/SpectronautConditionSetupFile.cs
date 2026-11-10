using CsvHelper;
using MassSpectrometry;
using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// A Spectronaut condition setup (<c>&lt;analysis&gt;_ConditionSetup.tsv</c>): the experimental design an analysis was
    /// run with, one row per run. Deposits put it under MassIVE's <c>metadata/</c>. Recognized by its name, and read
    /// as tab text.
    /// <para>
    /// <see cref="ToExperimentalDesign"/> gives the design as mzLib quantification takes it. The 1-based Replicate and
    /// Fraction become mzLib's 0-based biological replicate and fraction, as the SDRF design does. The technical
    /// replicate is 0, because Spectronaut's setup has no such column.
    /// </para>
    /// </summary>
    public class SpectronautConditionSetupFile : ResultFile<SpectronautConditionSetupRow>
    {
        public override SupportedFileType FileType => SupportedFileType.SpectronautConditionSetup;
        public override Software Software { get; set; }

        public SpectronautConditionSetupFile(string filePath) : base(filePath, Software.Spectronaut) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public SpectronautConditionSetupFile() : base()
        {
            Software = Software.Spectronaut;
        }

        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="MzLibException">
        /// The file is empty, has no Run Label or Condition column, or holds a value that cannot be read; the message
        /// names the file.
        /// </exception>
        public override void LoadResults()
        {
            const string what = "Spectronaut condition setup";
            if (!File.Exists(FilePath))
                throw new FileNotFoundException($"{what} not found: '{FilePath}'", FilePath);

            try
            {
                string[] header = SpectronautReports.ReadHeader(FilePath, out char delimiter);
                if (header.Length == 0)
                    throw new MzLibException($"{what} '{FilePath}' is empty.");
                foreach (string required in new[] { "Run Label", "Condition" })
                    if (!header.Contains(required))
                        throw new MzLibException($"{what} '{FilePath}' has no {required} column.");

                using var stream = new FileStream(FilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
                using var reader = new StreamReader(stream);
                using var csv = new CsvReader(reader, SpectronautReports.Configuration(delimiter));
                Results = csv.GetRecords<SpectronautConditionSetupRow>().ToList();
            }
            catch (Exception e) when (e is not MzLibException)
            {
                throw new MzLibException($"Could not read {what} '{FilePath}': {e.Message}", e);
            }
        }

        /// <summary>
        /// The design as mzLib quantification takes it: one label-free <see cref="SpectraFileInfo"/> per run.
        /// </summary>
        /// <param name="spectraFilePaths">
        /// The spectra files to key the design on. Each run is matched by the exact stem of its Run Label, ignoring case,
        /// and then by its File Name. When null, each run is keyed on its Run Label as written.
        /// </param>
        /// <exception cref="MzLibException">
        /// A run has no condition (<c>Not Defined</c>) or no replicate, or matches no file or several, when
        /// <paramref name="spectraFilePaths"/> is given. The message names the run. A design is refused rather than
        /// built with a guessed entry.
        /// </exception>
        public SampleExperimentalDesign ToExperimentalDesign(IEnumerable<string>? spectraFilePaths = null)
        {
            List<string>? paths = spectraFilePaths?.ToList();
            var files = new List<SpectraFileInfo>();
            foreach (var row in Results)
            {
                string condition = row.Condition
                    ?? throw new MzLibException($"Run '{row.RunLabel}' in '{FilePath}' has no condition ('{row.ConditionText}').");
                if (row.Replicate is not int replicate || replicate < 1)
                    throw new MzLibException($"Run '{row.RunLabel}' in '{FilePath}' has no replicate number of 1 or more.");
                string path = paths is null ? row.RunLabel : PathFor(row, paths);
                int fraction = row.Fraction is int f && f >= 1 ? f - 1 : 0;
                files.Add(new SpectraFileInfo(path, condition, replicate - 1, 0, fraction));
            }
            return SampleExperimentalDesign.LabelFree(files);
        }

        private string PathFor(SpectronautConditionSetupRow row, List<string> paths)
        {
            foreach (string? stem in new[] { row.RunStem, row.FileName })
            {
                if (string.IsNullOrWhiteSpace(stem))
                    continue;
                var matches = paths.Where(p => RunFileMatcher.StripPathAndExtension(p).Equals(stem, StringComparison.OrdinalIgnoreCase))
                    .Distinct(StringComparer.OrdinalIgnoreCase).ToList();
                if (matches.Count == 1)
                    return matches[0];
                if (matches.Count > 1)
                    throw new MzLibException($"Run '{row.RunLabel}' in '{FilePath}' matches {matches.Count} spectra files: {string.Join(", ", matches)}.");
            }
            throw new MzLibException($"Run '{row.RunLabel}' in '{FilePath}' matches none of the given spectra files.");
        }

        /// <exception cref="NotSupportedException">Always.</exception>
        public override void WriteResults(string outputPath) =>
            throw new NotSupportedException("Writing a Spectronaut condition setup is not supported; SpectronautConditionSetupFile only reads it.");
    }
}
