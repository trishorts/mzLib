using System.Text;
using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// An SDRF-Proteomics (Sample and Data Relationship Format) document -- the HUPO-PSI standard
    /// for describing the sample-to-data-file relationship of a proteomics experiment. One row per
    /// data file; columns carry sample characteristics, acquisition and search metadata, and the
    /// experimental factor values. Specification v1.1.0 (2026-01).
    ///
    /// This reader is deliberately hand-rolled rather than CsvHelper-driven, which departs from
    /// every other reader in this project. The house pattern maps FIXED columns onto a record type
    /// with [Name(...)] attributes, and SDRF cannot be expressed that way: its column set is
    /// open-ended, its names are data, and names repeat within a single file. The nearest existing
    /// precedent in mzLib for dynamic headers is SpectrumMatchFromTsv's hand-rolled parsedHeader.
    ///
    /// Three further properties of real files were measured across the 1,236-file curated corpus at
    /// bigbio/sdrf-annotated-datasets and each one is load-bearing here:
    ///
    ///   1. THERE IS NO QUOTING OR ESCAPING. SDRF is strictly tab-delimited; a cell cannot contain
    ///      a tab or a newline, so no escape mechanism is needed or defined. Two corpus files
    ///      (PXD026824, PXD028735) contain literal double-quote characters INSIDE cell values, e.g.
    ///      "NT=Phospho;AC=UNIMOD:21;TA=S,T,Y;MT=Variable" with the quotes part of the data. A CSV
    ///      parser with default quoting would strip them and the file would not round-trip. Splitting
    ///      on '\t' with no escape handling is both simpler and more faithful.
    ///
    ///   2. LINE ENDINGS ARE CRLF, WITH A TRAILING NEWLINE. All 1,236 corpus files, without
    ///      exception, and none carry a byte-order mark. Writing LF would change every byte after
    ///      the first line, so <see cref="WriteResults"/> pins CRLF rather than inheriting the
    ///      platform default.
    ///
    ///   3. ROWS MAY BE SHORT. See <see cref="SdrfRow.Cells"/>.
    ///
    /// Together these make a byte-identical round trip achievable for every file in the corpus,
    /// which is the regression test this type exists to pass. Reporting a file's defects is the
    /// separate concern of the validator: a reader that repaired what it read could not round-trip,
    /// and silently rewriting a curated annotation is worse than reporting it.
    /// </summary>
    public class SdrfDocument : ResultFile<SdrfRow>, IResultFile
    {
        /// <summary>
        /// SDRF has no escape mechanism, so a cell may never contain a tab or a newline.
        /// Enforced on write, where a violation would silently corrupt the column structure.
        /// </summary>
        private const string CellSeparator = "\t";

        private SdrfHeader? _header;

        public override SupportedFileType FileType => SupportedFileType.Sdrf;

        public override Software Software { get; set; }

        public SdrfDocument(string filePath) : base(filePath, Software.Unspecified) { }

        /// <summary>
        /// Constructor used to initialize from the factory method. Required: FileReader.ReadResultFile
        /// creates result files through Activator.CreateInstance.
        /// </summary>
        public SdrfDocument() : base() { }

        /// <summary>
        /// Builds a document in memory, for writing. FilePath is left empty, which is what stops the
        /// base class lazy-loading over the top of these rows.
        /// </summary>
        public SdrfDocument(SdrfHeader header, IEnumerable<SdrfRow> rows) : base()
        {
            _header = header ?? throw new ArgumentNullException(nameof(header));
            Results = rows?.ToList() ?? throw new ArgumentNullException(nameof(rows));
        }

        /// <summary>
        /// The ordered column names. Triggers a load for the same reason Results does, so callers
        /// can inspect the columns without first touching a row.
        /// </summary>
        public SdrfHeader Header
        {
            get
            {
                if (_header is null && File.Exists(FilePath))
                    LoadResults();
                return _header ?? new SdrfHeader(Array.Empty<string>());
            }
        }

        public override void LoadResults()
        {
            // Read the whole file as text rather than streaming lines: SDRF documents are small
            // (the largest in the corpus is 5,798 rows) and reading the raw text is what lets the
            // line-ending and trailing-newline handling below stay explicit instead of being
            // silently normalized by StreamReader.ReadLine.
            string text = File.ReadAllText(FilePath, Encoding.UTF8);

            var lines = SplitLines(text);
            if (lines.Count == 0)
                throw new MzLibException($"SDRF file is empty: '{FilePath}'");

            _header = new SdrfHeader(lines[0].Split(CellSeparator));

            var rows = new List<SdrfRow>(lines.Count - 1);
            for (int i = 1; i < lines.Count; i++)
                rows.Add(new SdrfRow(_header, lines[i].Split(CellSeparator)));

            Results = rows;
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            var header = Header;

            // UTF8Encoding(false) -- no byte-order mark. Encoding.UTF8 emits one, and not one of the
            // 1,236 corpus files has it; adding a BOM would break a byte-identical round trip on the
            // very first byte.
            using var writer = new StreamWriter(File.Create(outputPath), new UTF8Encoding(false));
            writer.NewLine = "\r\n";

            writer.WriteLine(string.Join(CellSeparator, header));
            foreach (var row in Results)
            {
                foreach (var cell in row.Cells)
                {
                    // A tab or newline inside a cell has no representation in this format -- there is
                    // no quoting to fall back on -- so it would silently shift every following column.
                    // Fail loudly instead of writing a file that reads back wrong.
                    if (cell.Contains('\t') || cell.Contains('\n') || cell.Contains('\r'))
                        throw new MzLibException(
                            "SDRF cells cannot contain a tab or newline; the format defines no escape " +
                            $"mechanism. Offending value in '{outputPath}': '{cell}'.");
                }
                writer.WriteLine(string.Join(CellSeparator, row.Cells));
            }
        }

        /// <summary>
        /// Splits on CRLF or LF, dropping the final empty fragment produced by a trailing newline.
        ///
        /// Blank lines are NOT skipped in the middle of a document -- a blank line is a row of one
        /// empty cell, and dropping it would lose a row on write. Only the single trailing newline
        /// that every well-formed file ends with is absorbed.
        /// </summary>
        private static List<string> SplitLines(string text)
        {
            var lines = new List<string>();
            if (text.Length == 0)
                return lines;

            int start = 0;
            for (int i = 0; i < text.Length; i++)
            {
                if (text[i] != '\n')
                    continue;

                int end = i;
                if (end > start && text[end - 1] == '\r')
                    end--;
                lines.Add(text.Substring(start, end - start));
                start = i + 1;
            }

            // A file not ending in a newline still has a final line.
            if (start < text.Length)
                lines.Add(text.Substring(start));

            return lines;
        }
    }
}
