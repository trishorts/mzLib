using System.Globalization;
using System.Text.RegularExpressions;
using CsvHelper;
using CsvHelper.Configuration;
using CsvHelper.TypeConversion;
using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// Why a Spectronaut quantity cell holds no number. Spectronaut writes two different tokens. <c>NaN</c> means
    /// nothing was measured: a protein group with no quantified precursor, or a run where it was not seen.
    /// <c>Filtered</c> means a value was removed by the q-value filter. The Run Pivot report writes both.
    /// </summary>
    public enum SpectronautValueState
    {
        /// <summary>The cell holds a number.</summary>
        Value,

        /// <summary>The cell is blank, <c>NaN</c> or <c>NA</c>: nothing was measured.</summary>
        NotMeasured,

        /// <summary>The cell is <c>Filtered</c>: a value existed but did not pass the q-value filter.</summary>
        Filtered,
    }

    /// <summary>The export grains of a Spectronaut report, told apart by the report's header.</summary>
    internal enum SpectronautGrain
    {
        None,
        Tric,
        PtmSite,
        Precursor,
        ProteinGroup,
        Pivot,
    }

    /// <summary>
    /// What every Spectronaut reader shares: header sniffing, delimiter choice, Spectronaut's missing-value tokens,
    /// <c>;</c>-list handling and the load-failure contract.
    /// <para>
    /// Spectronaut report schemas can be configured column by column, so one deposit's report can differ from the
    /// next. Five real MassIVE deposits differed in grain (protein group, precursor, fragment, PTM site, wide pivot,
    /// TRIC), in extension (<c>.tsv</c>, <c>.xls</c> holding tab text, comma <c>.csv</c> from Spectronaut 10) and in
    /// which columns were exported. Readers therefore detect by header, choose the delimiter by content, require only
    /// a small core of columns and read the rest as optional.
    /// </para>
    /// </summary>
    internal static class SpectronautReports
    {
        /// <summary>A Run Pivot column: <c>[1] run.wiff.PG.Quantity</c>, a 1-based run index, the run's file name with
        /// its extension, then the hierarchy column. The run name is greedy so the split falls on the LAST
        /// <c>.PG.</c>/<c>.EG.</c>/... marker, because run names carry dots of their own.</summary>
        internal static readonly Regex PivotColumn = new(@"^\[(\d+)\] (.+)\.((?:PG|PEP|EG|FG|F)\.[^.].*)$", RegexOptions.Compiled);

        /// <summary>The charge suffix of a precursor id (<c>_PEPTIDE_.2</c>).</summary>
        private static readonly Regex ChargeSuffix = new(@"^(.*)\.(\d+)$", RegexOptions.Compiled);

        private static readonly Regex Bracketed = new(@"\[[^\]]*\]", RegexOptions.Compiled);

        /// <summary>Tokens Spectronaut writes for "no value". <c>Not Defined</c> is its mark for an unset condition.</summary>
        private static readonly HashSet<string> NotStatedTokens = new(StringComparer.OrdinalIgnoreCase)
            { "", "NaN", "NA", "N/A", "Filtered", "Not Defined" };

        /// <summary>The configuration every Spectronaut reader uses: invariant culture, the file's own delimiter, quotes honoured
        /// (Spectronaut 10 exports re-saved by R quote every string).</summary>
        internal static CsvConfiguration Configuration(char delimiter) => new(CultureInfo.InvariantCulture)
        {
            Delimiter = delimiter.ToString(),
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = TrimOptions.Trim,
            BadDataFound = null,
        };

        /// <summary>
        /// The first line's column names and the delimiter: a tab when the line holds one, otherwise a comma. Quotes are
        /// removed. Empty when the file is empty.
        /// </summary>
        internal static string[] ReadHeader(string filePath, out char delimiter)
        {
            using var stream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
            using var reader = new StreamReader(stream);
            string line = reader.ReadLine() ?? "";
            return SplitHeader(line, out delimiter);
        }

        internal static string[] SplitHeader(string line, out char delimiter)
        {
            delimiter = line.Contains('\t') ? '\t' : ',';
            if (line.Length == 0) return Array.Empty<string>();
            return line.Split(delimiter).Select(name => name.Trim().Trim('"')).ToArray();
        }

        /// <summary>
        /// The grain a header belongs to, or <see cref="SpectronautGrain.None"/>. The most specific grain is tested first,
        /// because the grains share columns: precursor and PTM site reports also carry the PG.* columns of a protein-group
        /// report.
        /// </summary>
        internal static SpectronautGrain Sniff(IReadOnlyCollection<string> header)
        {
            var names = header as ISet<string> ?? new HashSet<string>(header, StringComparer.Ordinal);
            if (names.Contains("transition_group_id") && names.Contains("PrecursorID") && names.Contains("m_score") && names.Contains("RT"))
                return SpectronautGrain.Tric;
            if (names.Contains("R.FileName") && names.Contains("PTM.SiteLocation") && names.Contains("PTM.Quantity"))
                return SpectronautGrain.PtmSite;
            if (names.Contains("R.FileName") && (names.Contains("EG.PrecursorId") || names.Contains("EG.ModifiedSequence") || names.Contains("FG.Id")))
                return SpectronautGrain.Precursor;
            if (names.Contains("R.FileName") && names.Contains("PG.Quantity") && (names.Contains("PG.ProteinGroups") || names.Contains("PG.ProteinAccessions")))
                return SpectronautGrain.ProteinGroup;
            if (!names.Contains("R.FileName") && names.Any(name => PivotColumn.IsMatch(name)))
                return SpectronautGrain.Pivot;
            return SpectronautGrain.None;
        }

        /// <summary>The grain of the file at <paramref name="filePath"/>; None for an empty file or a header no grain claims.</summary>
        internal static SpectronautGrain SniffFile(string filePath) => Sniff(ReadHeader(filePath, out _));

        internal static string Describe(SpectronautGrain grain) => grain switch
        {
            SpectronautGrain.Tric => "TRIC report",
            SpectronautGrain.PtmSite => "PTM site report",
            SpectronautGrain.Precursor => "precursor report",
            SpectronautGrain.ProteinGroup => "protein-group report",
            SpectronautGrain.Pivot => "Run Pivot report",
            _ => "file no Spectronaut reader recognizes",
        };

        /// <summary>True for a blank cell or one of Spectronaut's "no value" tokens.</summary>
        internal static bool IsNotStated(string? text) => text is null || NotStatedTokens.Contains(text.Trim());

        /// <summary>The cell's text as a value, or null for a "no value" token.</summary>
        internal static string? Stated(string? text) => IsNotStated(text) ? null : text!.Trim();

        /// <summary>Why a quantity cell holds no number, or <see cref="SpectronautValueState.Value"/> when it holds one.</summary>
        internal static SpectronautValueState StateOf(string? text)
        {
            if (text is not null && text.Trim().Equals("Filtered", StringComparison.OrdinalIgnoreCase))
                return SpectronautValueState.Filtered;
            return IsNotStated(text) ? SpectronautValueState.NotMeasured : SpectronautValueState.Value;
        }

        /// <summary>A number in invariant culture, or null for a "no value" token.</summary>
        /// <exception cref="FormatException">The text is neither a number nor a "no value" token.</exception>
        internal static double? ParseDouble(string? text)
        {
            if (IsNotStated(text)) return null;
            if (double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out double value))
                return value;
            throw new FormatException($"'{text}' is not a number");
        }

        /// <summary>A percentage written with its sign (<c>22.8%</c>) as a fraction (0.228); null for a "no value" token.</summary>
        /// <exception cref="FormatException">The text is not a percentage.</exception>
        internal static double? ParsePercent(string? text)
        {
            if (IsNotStated(text)) return null;
            string trimmed = text!.Trim();
            if (trimmed.EndsWith('%')) trimmed = trimmed[..^1];
            if (double.TryParse(trimmed, NumberStyles.Float, CultureInfo.InvariantCulture, out double value))
                return value / 100;
            throw new FormatException($"'{text}' is not a percentage");
        }

        /// <summary>
        /// The distinct members of a protein-group cell, in order. Spectronaut separates proteins within a group with
        /// <c>;</c> and groups with <c>|</c>. A protein can be listed twice (<c>Q8VED5;Q8VED5</c> in MSV000096226), so
        /// duplicates are dropped, as Spectronaut itself drops them from PG.ProteinAccessions.
        /// </summary>
        internal static IReadOnlyList<string> Members(string? proteinGroups) =>
            string.IsNullOrWhiteSpace(proteinGroups)
                ? Array.Empty<string>()
                : proteinGroups.Split(new[] { ';', '|' }, StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
                    .Distinct(StringComparer.Ordinal).ToArray();

        /// <summary>
        /// One value per member from a PG.* list column, or null when the list cannot be attributed. Spectronaut writes
        /// these lists "in the same order as the protein ids in PG.ProteinGroups", but collapses a list to one value when
        /// every member shares it, and leaves empty entries (<c>;Krt31;</c>). So a list as long as the member list gives one
        /// value per member, a single value applies to all, and any other length cannot be attributed. In MSV000096226
        /// 32 of the 44 multi-protein groups had a list of another length. An empty entry is null.
        /// </summary>
        internal static IReadOnlyList<string?>? PerMember(string? list, int memberCount)
        {
            if (string.IsNullOrEmpty(list) || memberCount == 0) return null;
            string[] parts = list.Split(';');
            if (parts.Length == memberCount)
                return parts.Select(part => string.IsNullOrWhiteSpace(part) ? null : part.Trim()).ToArray();
            if (parts.Length == 1)
                return Enumerable.Repeat<string?>(string.IsNullOrWhiteSpace(parts[0]) ? null : parts[0].Trim(), memberCount).ToArray();
            return null;
        }

        /// <summary>
        /// Splits a precursor id (<c>_PEPT[Phospho (STY)]IDE_.2</c>) into its modified sequence and charge. An id with no
        /// charge suffix is returned whole, with a null charge.
        /// </summary>
        internal static string SplitPrecursorId(string precursorId, out int? charge)
        {
            Match match = ChargeSuffix.Match(precursorId);
            if (match.Success && int.TryParse(match.Groups[2].Value, NumberStyles.None, CultureInfo.InvariantCulture, out int z))
            {
                charge = z;
                return match.Groups[1].Value;
            }
            charge = null;
            return precursorId;
        }

        /// <summary>The bare residues of a modified sequence: brackets and the <c>_</c> terminus anchors removed.</summary>
        internal static string StripModifications(string modifiedSequence) =>
            Bracketed.Replace(modifiedSequence, "").Replace("_", "");

        /// <summary>
        /// Reads a Spectronaut report into records with the newest reader failure contract. A missing file throws
        /// <see cref="FileNotFoundException"/>. A file of another grain, a missing core column or any parse failure throws
        /// <see cref="MzLibException"/> naming the file.
        /// </summary>
        /// <param name="filePath">The report.</param>
        /// <param name="expected">The grain the calling reader reads.</param>
        /// <param name="missingCore">Given the header, names the missing core columns, or returns null when none are missing.</param>
        /// <param name="collect">Turns the record stream into the results (de-duplicating, selecting, pooling).</param>
        /// <param name="header">The file's column names.</param>
        internal static List<T> Read<T>(string filePath, SpectronautGrain expected, Func<IReadOnlySet<string>, string?> missingCore,
            Func<IEnumerable<T>, List<T>> collect, out string[] header)
        {
            string what = "Spectronaut " + Describe(expected);
            if (!File.Exists(filePath))
                throw new FileNotFoundException($"{what} not found: '{filePath}'", filePath);

            try
            {
                header = ReadHeader(filePath, out char delimiter);
                if (header.Length == 0)
                    throw new MzLibException($"{what} '{filePath}' is empty.");
                var names = new HashSet<string>(header, StringComparer.Ordinal);
                SpectronautGrain grain = Sniff(names);
                if (grain != expected)
                    throw new MzLibException($"'{filePath}' is not a Spectronaut {Describe(expected)}: its header is that of a {Describe(grain)}.");
                string? missing = missingCore(names);
                if (missing is not null)
                    throw new MzLibException($"{what} '{filePath}' has no {missing} column.");

                using var stream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
                using var reader = new StreamReader(stream);
                using var csv = new CsvReader(reader, Configuration(delimiter));
                return collect(csv.GetRecords<T>());
            }
            catch (Exception e) when (e is not MzLibException)
            {
                throw new MzLibException($"Could not read {what} '{filePath}': {e.Message}", e);
            }
        }

        /// <summary>A string pool for the columns a report repeats on every row (run names, protein groups).</summary>
        internal sealed class StringPool
        {
            private readonly Dictionary<string, string> _pool = new(StringComparer.Ordinal);

            public string? Get(string? value)
            {
                if (value is null) return null;
                if (_pool.TryGetValue(value, out var shared)) return shared;
                _pool[value] = value;
                return value;
            }
        }
    }

    /// <summary>A double cell that may hold one of Spectronaut's "no value" tokens (<c>NaN</c>, <c>Filtered</c>, ...): null then.</summary>
    internal class SpectronautNullableDoubleConverter : DefaultTypeConverter
    {
        public override object? ConvertFromString(string? text, IReaderRow row, MemberMapData memberMapData)
        {
            if (SpectronautReports.IsNotStated(text)) return null;
            return double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out double value)
                ? value
                : base.ConvertFromString(text, row, memberMapData);
        }

        public override string? ConvertToString(object? value, IWriterRow row, MemberMapData memberMapData) =>
            value is double d ? d.ToString(CultureInfo.InvariantCulture) : "";
    }

    /// <summary>An integer cell that may hold a "no value" token: null then.</summary>
    internal class SpectronautNullableIntConverter : DefaultTypeConverter
    {
        public override object? ConvertFromString(string? text, IReaderRow row, MemberMapData memberMapData)
        {
            if (SpectronautReports.IsNotStated(text)) return null;
            return int.TryParse(text, NumberStyles.Integer, CultureInfo.InvariantCulture, out int value)
                ? value
                : base.ConvertFromString(text, row, memberMapData);
        }

        public override string? ConvertToString(object? value, IWriterRow row, MemberMapData memberMapData) =>
            value is int i ? i.ToString(CultureInfo.InvariantCulture) : "";
    }

    /// <summary>A boolean cell written <c>True</c>/<c>False</c> (reports), <c>TRUE</c>/<c>FALSE</c> (libraries) or <c>1</c>/<c>0</c>
    /// (TRIC); null for a "no value" token.</summary>
    internal class SpectronautNullableBoolConverter : DefaultTypeConverter
    {
        public override object? ConvertFromString(string? text, IReaderRow row, MemberMapData memberMapData)
        {
            if (SpectronautReports.IsNotStated(text)) return null;
            string trimmed = text!.Trim();
            if (bool.TryParse(trimmed, out bool flag)) return flag;
            if (trimmed == "1") return true;
            if (trimmed == "0") return false;
            return base.ConvertFromString(text, row, memberMapData);
        }

        public override string? ConvertToString(object? value, IWriterRow row, MemberMapData memberMapData) =>
            value is bool b ? (b ? "True" : "False") : "";
    }
}
