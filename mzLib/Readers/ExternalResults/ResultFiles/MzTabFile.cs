using System.Globalization;
using System.IO.Compression;
using System.Text;
using System.Text.RegularExpressions;
using MzLibUtil;
using Omics;
using Omics.SequenceConversion;
using Readers.ExternalResults.IndividualResultRecords;

namespace Readers.ExternalResults.ResultFiles;

/// <summary>
/// An mzTab 1.0 file (.mztab), or a gzip-compressed one (.mztab.gz) read without decompressing it to disk: the
/// metadata and the PSM section, from any producer (MaxQuant, Proteome Discoverer, Mascot, OpenMS and quantms,
/// PRIDE's and MassIVE's converters).
/// </summary>
/// <remarks>
/// <para>
/// mzTab-M (metabolomics, 2.0.0-M) is a different format and is refused, as is a pre-1.0 file with no mzTab-version.
/// Casanovo's de novo dialect has its own reader, <see cref="CasanovoMzTabFile"/>, which file-type detection picks.
/// </para>
/// <para>
/// One PSM is the rows that share PSM_ID, spectra_ref, sequence and modifications. <see cref="ResultFile{TResult}.Results"/>
/// merges such rows wherever they are, because producers scatter them: PRIDE's converter writes protein by protein, so
/// a PSM returns under each of its proteins, and Proteome Discoverer repeats exact copies up to thousands of rows
/// apart. <see cref="EnumeratePsms"/> streams instead, for files too large to hold, and merges only adjacent rows.
/// </para>
/// <para>
/// A PSM whose modifications do not all resolve is kept, with <see cref="MzTabPsm.ModificationsResolved"/> false and
/// an empty FullSequence. Only a row that is not a PSM at all goes to <see cref="SkippedPsms"/>.
/// </para>
/// <para>
/// CAREFUL with downloads: .NET's GZipStream reports corrupt data but not a truncated stream, and mzTab has no end
/// marker. A file cut off part-way, compressed or not, reads as a shorter file without error. Check that a downloaded
/// file is complete (its size or checksum) before reading it.
/// </para>
/// </remarks>
public class MzTabFile : ResultFile<MzTabPsm>, IResultFile
{
    private static readonly Regex SpectraRefPattern = new(@"^ms_run\[(\d+)\]:(.*)$", RegexOptions.Compiled);
    private static readonly Regex ScoreColumn = new(@"^search_engine_score\[(\d+)\]$", RegexOptions.IgnoreCase | RegexOptions.Compiled);
    private const string ScanColumnSuffix = "_cv_MS:1000776_scan_number_only_nativeID_format";

    public override SupportedFileType FileType =>
        FilePath.EndsWith(".gz", StringComparison.OrdinalIgnoreCase) ? SupportedFileType.MzTabGz : SupportedFileType.MzTab;

    public override Software Software { get; set; }

    private List<MzTabSkippedPsm> _skippedPsms = [];
    private string? _loadedFrom;
    private MzTabMetadata? _metadata;
    private string? _metadataFrom;

    public MzTabFile(string filePath) : base(filePath, Software.MzTab) { }

    /// <summary>Constructor used to initialize from the factory method.</summary>
    public MzTabFile() : base()
    {
        Software = Software.MzTab;
    }

    /// <summary>
    /// The metadata section. Reading it reads only the MTD lines at the top of the file, not the PSMs.
    /// </summary>
    /// <exception cref="FileNotFoundException">The file does not exist.</exception>
    /// <exception cref="MzLibException">The file could not be read.</exception>
    public MzTabMetadata Metadata
    {
        get
        {
            if (_metadata == null || _metadataFrom != FilePath)
            {
                using var reader = Open(FilePath);
                var entries = new List<(string, string)>();
                string? line;
                while ((line = ReadLine(reader, FilePath)) != null)
                {
                    if (!MzTabMetadata.TryReadMetadataLine(line, entries, out bool endOfBlock) && endOfBlock)
                    {
                        break;
                    }
                }

                _metadata = new MzTabMetadata(entries);
                _metadataFrom = FilePath;
            }

            return _metadata;
        }
    }

    /// <summary>
    /// The PSM rows that are not PSMs at all, and why: a row with no sequence. Reading this forces the load.
    /// </summary>
    public IReadOnlyList<MzTabSkippedPsm> SkippedPsms
    {
        get
        {
            _ = Results;
            return _skippedPsms;
        }
    }

    /// <summary>
    /// Reads every PSM, merging the rows of one PSM across the whole file. A file is read once per
    /// <see cref="ResultFile{TResult}.FilePath"/>: ResultFile reloads whenever Results is empty, and a file with no
    /// PSMs would otherwise be parsed again on every access.
    /// </summary>
    /// <exception cref="FileNotFoundException">The file does not exist.</exception>
    /// <exception cref="MzLibException">The file is mzTab-M or pre-1.0 mzTab, is malformed, or could not be read.</exception>
    public override void LoadResults()
    {
        if (_loadedFrom == FilePath)
        {
            return;
        }

        var skipped = new List<MzTabSkippedPsm>();
        var byKey = new Dictionary<(string, string, string, string), MzTabPsm>();
        var results = new List<MzTabPsm>();
        foreach (var psm in Enumerate(FilePath, skipped.Add))
        {
            if (byKey.TryGetValue(psm.Key, out var earlier))
            {
                earlier.AddProteins(psm.Proteins);
            }
            else
            {
                byKey.Add(psm.Key, psm);
                results.Add(psm);
            }
        }

        _skippedPsms = skipped;
        Results = results;
        _loadedFrom = FilePath;
    }

    /// <summary>
    /// Streams the PSMs without holding the file: the path for files too large to load (quantms writes mzTab files of
    /// over 100 GB). It neither reads nor fills <see cref="ResultFile{TResult}.Results"/>.
    /// </summary>
    /// <remarks>
    /// Only adjacent rows of one PSM are merged. That is exact for producers that write a PSM's rows together
    /// (MaxQuant, Mascot, MassIVE, OpenMS, quantms), but a PSM from PRIDE's converter or Proteome Discoverer can come
    /// out more than once, each time with the proteins of that block. Use Results when that matters.
    /// </remarks>
    /// <param name="onSkipped">Called for each row that is not a PSM at all.</param>
    /// <exception cref="FileNotFoundException">The file does not exist.</exception>
    /// <exception cref="MzLibException">The file is mzTab-M or pre-1.0 mzTab, is malformed, or could not be read.
    /// It is thrown when enumeration reaches the problem, so the PSMs before it have already been returned.</exception>
    public IEnumerable<MzTabPsm> EnumeratePsms(Action<MzTabSkippedPsm>? onSkipped = null) => Enumerate(FilePath, onSkipped);

    /// <summary>
    /// mzTab cannot be regenerated from these records, so writing is not supported.
    /// </summary>
    /// <exception cref="NotSupportedException">Always.</exception>
    public override void WriteResults(string outputPath) =>
        throw new NotSupportedException("Writing mzTab is not supported; MzTabFile only reads it.");

    /// <summary>
    /// The file type of an .mztab or .mztab.gz file that is not Casanovo's, from its mzTab-version.
    /// </summary>
    /// <exception cref="MzLibException">The file is mzTab-M, has no mzTab-version, or could not be read.</exception>
    internal static SupportedFileType DetectFileType(string path)
    {
        MzTabMetadata metadata;
        using (var reader = Open(path))
        {
            try
            {
                metadata = MzTabMetadata.Read(reader);
            }
            catch (Exception e) when (IsReadFailure(e))
            {
                throw new MzLibException($"Could not read mzTab file '{path}': {e.Message}", e);
            }
        }

        if (metadata.IsMzTabM)
        {
            throw new MzLibException($"MzTab file type not supported: mzTab-M ({metadata.Version}) is a different format");
        }

        if (metadata.Version == null)
        {
            throw new MzLibException("MzTab file type not supported: the metadata has no mzTab-version, as in files older than mzTab 1.0");
        }

        return path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase) ? SupportedFileType.MzTabGz : SupportedFileType.MzTab;
    }

    private static IEnumerable<MzTabPsm> Enumerate(string path, Action<MzTabSkippedPsm>? onSkipped)
    {
        using var reader = Open(path);
        var entries = new List<(string, string)>();
        MzTabMetadata? metadata = null;
        PsmColumns? columns = null;
        var lookup = new UnimodModificationLookup();
        var placedBySequence = new Dictionary<(string, string), MzTabModifications.Placed>();
        var group = new List<string[]>();
        (string, string, string, string) groupKey = default;
        int lineNumber = 0;

        string? line;
        while ((line = ReadLine(reader, path)) != null)
        {
            lineNumber++;
            if (metadata == null)
            {
                if (MzTabMetadata.TryReadMetadataLine(line, entries, out bool endOfBlock) || !endOfBlock)
                {
                    continue;
                }

                metadata = Validated(new MzTabMetadata(entries), path);
            }

            if (string.IsNullOrWhiteSpace(line))
            {
                continue;
            }

            var fields = line.Split('\t');
            string prefix = fields[0].Trim();
            if (prefix == "PSM")
            {
                if (columns == null)
                {
                    throw new MzLibException($"mzTab file '{path}', line {lineNumber}: a PSM row comes before the PSH header");
                }

                var key = columns.KeyOf(fields);
                if (columns.Value(fields, "sequence") == null)
                {
                    onSkipped?.Invoke(new MzTabSkippedPsm(key.Item1, columns.Value(fields, "spectra_ref"), null,
                        columns.Value(fields, "modifications"), "the row has no sequence"));
                    continue;
                }

                if (group.Count > 0 && key == groupKey)
                {
                    group.Add(fields);
                    continue;
                }

                if (group.Count > 0)
                {
                    yield return Build(group, groupKey, columns, metadata, lookup, placedBySequence);
                }

                group = [fields];
                groupKey = key;
                continue;
            }

            if (prefix == "COM")
            {
                continue;
            }

            // Any other line ends a run of PSM rows
            if (group.Count > 0)
            {
                yield return Build(group, groupKey, columns!, metadata, lookup, placedBySequence);
                group = [];
            }

            if (prefix == "PSH")
            {
                columns = new PsmColumns(fields);
            }
        }

        metadata ??= Validated(new MzTabMetadata(entries), path);
        if (group.Count > 0)
        {
            yield return Build(group, groupKey, columns!, metadata, lookup, placedBySequence);
        }
    }

    private static MzTabMetadata Validated(MzTabMetadata metadata, string path)
    {
        if (metadata.IsMzTabM)
        {
            throw new MzLibException($"'{path}' is mzTab-M ({metadata.Version}), a different format; MzTabFile reads mzTab 1.0");
        }

        if (metadata.Version == null)
        {
            throw new MzLibException($"'{path}' has no mzTab-version in its metadata. MzTabFile reads mzTab 1.0, not files older than it.");
        }

        return metadata;
    }

    private static MzTabPsm Build(List<string[]> rows, (string, string, string, string) key, PsmColumns columns,
        MzTabMetadata metadata, UnimodModificationLookup lookup, Dictionary<(string, string), MzTabModifications.Placed> placedBySequence)
    {
        var first = rows[0];
        string sequence = columns.Value(first, "sequence")!;
        string? modifications = columns.Value(first, "modifications");

        if (!placedBySequence.TryGetValue((sequence, modifications ?? ""), out var placed))
        {
            placed = MzTabModifications.Place(modifications, sequence, metadata, lookup);
            placedBySequence.Add((sequence, modifications ?? ""), placed);
        }

        string spectraRef = columns.Value(first, "spectra_ref") ?? "";
        string firstRef = spectraRef.Split('|')[0].Trim();
        var reference = SpectraRefPattern.Match(firstRef);
        int msRunIndex = reference.Success ? int.Parse(reference.Groups[1].Value, CultureInfo.InvariantCulture) : -1;
        string nativeId = reference.Success ? reference.Groups[2].Value : firstRef;
        string location = metadata.MsRuns.TryGetValue(msRunIndex, out var run) ? run.Location ?? "" : "";

        int scanNumber = columns.ScanColumn != null && Integer(columns.Value(first, columns.ScanColumn)) is int scan
            ? scan
            : SpectrumNativeId.OneBasedScanNumberOf(nativeId);

        var mods = new Dictionary<int, Omics.Modifications.Modification>(placed.Mods);
        var psm = new MzTabPsm
        {
            Key = key,
            PsmId = columns.Value(first, "PSM_ID") ?? "",
            BaseSequence = sequence,
            Modifications = modifications,
            AllModsOneIsNterminus = mods,
            UnresolvedModifications = placed.Unresolved,
            FullSequence = placed.Unresolved.Count == 0 ? IBioPolymerWithSetMods.DetermineFullSequence(sequence, mods) : "",
            Unique = columns.Value(first, "unique") switch { "1" => true, "0" => false, _ => null },
            Database = columns.Value(first, "database"),
            DatabaseVersion = columns.Value(first, "database_version"),
            SearchEngines = MzTabCvParam.ParseList(columns.Value(first, "search_engine")),
            SearchEngineScores = columns.Scores.ToDictionary(s => s.Key, s => Number(columns.Value(first, s.Value))),
            Reliability = Integer(columns.Value(first, "reliability")),
            RetentionTimeText = columns.Value(first, "retention_time"),
            RetentionTime = FirstRetentionTime(columns.Value(first, "retention_time")),
            Charge = Integer(columns.Value(first, "charge")),
            ExperimentalMassToCharge = Number(columns.Value(first, "exp_mass_to_charge")),
            CalculatedMassToCharge = Number(columns.Value(first, "calc_mass_to_charge")),
            Uri = columns.Value(first, "uri"),
            SpectraRef = spectraRef,
            MsRunIndex = msRunIndex,
            NativeId = nativeId,
            MsRunLocation = location,
            // Locations are written on the submitter's machine, so either separator ends a directory on any OS
            FileNameWithoutExtension = location.Length == 0 ? "" :
                PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(location[(location.LastIndexOfAny(['/', '\\']) + 1)..]),
            OneBasedScanNumber = scanNumber,
            OptionalColumns = columns.Optional.ToDictionary(c => c, c => columns.Value(first, c)),
            DecoyColumn = columns.DecoyColumn == null ? null : columns.Value(first, columns.DecoyColumn) switch
            {
                "1" => true,
                "0" => false,
                { } text when text.Equals("true", StringComparison.OrdinalIgnoreCase) => true,
                { } text when text.Equals("false", StringComparison.OrdinalIgnoreCase) => false,
                _ => null,
            },
        };
        psm.AddProteins(rows.SelectMany(row => ProteinsOf(row, columns)));
        return psm;
    }

    /// <summary>
    /// The proteins of one row. OpenMS and quantms join several in one cell with commas, and pre/post/start/end
    /// alongside; a per-protein value is taken only when its list is as long as the accessions, or has one value.
    /// </summary>
    private static IEnumerable<MzTabPsmProtein> ProteinsOf(string[] row, PsmColumns columns)
    {
        var accessions = SplitList(columns.Value(row, "accession"));
        var pre = SplitList(columns.Value(row, "pre"));
        var post = SplitList(columns.Value(row, "post"));
        var start = SplitList(columns.Value(row, "start"));
        var end = SplitList(columns.Value(row, "end"));
        for (int i = 0; i < accessions.Count; i++)
        {
            yield return new MzTabPsmProtein(accessions[i], Pick(pre, i, accessions.Count), Pick(post, i, accessions.Count),
                Integer(Pick(start, i, accessions.Count)), Integer(Pick(end, i, accessions.Count)));
        }
    }

    private static List<string> SplitList(string? cell) =>
        cell == null ? [] : cell.Split([',', ';']).Select(s => s.Trim()).Where(s => s.Length > 0).ToList();

    private static string? Pick(List<string> values, int i, int count)
    {
        string? value = values.Count == count ? values[i] : values.Count == 1 ? values[0] : null;
        return value == null || value.Equals("null", StringComparison.OrdinalIgnoreCase) ? null : value;
    }

    /// <summary>The first of several '|'-separated values, or the start of a Mascot range ("2345.46-2351.07").</summary>
    private static double? FirstRetentionTime(string? text)
    {
        if (text == null)
        {
            return null;
        }

        string first = text.Split('|')[0].Trim();
        int dash = first.IndexOf('-', 1);
        if (dash > 0 && first[dash - 1] is not ('e' or 'E'))
        {
            first = first[..dash];
        }

        return Number(first);
    }

    private static double? Number(string? token)
    {
        if (token == null)
        {
            return null;
        }

        if (token.Equals("nan", StringComparison.OrdinalIgnoreCase))
        {
            return double.NaN;
        }

        return double.TryParse(token, NumberStyles.Float, CultureInfo.InvariantCulture, out double value) ? value : null;
    }

    private static int? Integer(string? token)
    {
        if (token == null)
        {
            return null;
        }

        if (int.TryParse(token, NumberStyles.Integer, CultureInfo.InvariantCulture, out int value))
        {
            return value;
        }

        // "2.0", which some converters write for a charge
        return double.TryParse(token, NumberStyles.Float, CultureInfo.InvariantCulture, out double number)
               && number == Math.Floor(number) && Math.Abs(number) <= int.MaxValue
            ? (int)number
            : null;
    }

    private static StreamReader Open(string path)
    {
        try
        {
            Stream stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read);
            if (path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase))
            {
                stream = new GZipStream(stream, CompressionMode.Decompress);
            }

            return new StreamReader(stream, Encoding.UTF8, detectEncodingFromByteOrderMarks: true);
        }
        catch (Exception e) when (IsReadFailure(e))
        {
            throw new MzLibException($"Could not read mzTab file '{path}': {e.Message}", e);
        }
    }

    private static string? ReadLine(StreamReader reader, string path)
    {
        try
        {
            return reader.ReadLine();
        }
        catch (Exception e) when (IsReadFailure(e))
        {
            throw new MzLibException($"Could not read mzTab file '{path}': {e.Message}", e);
        }
    }

    /// <summary>
    /// Anything that goes wrong reading an existing file: a corrupt gzip, an I/O fault mid-read, or access denied. A
    /// missing file or directory is not one of these; it propagates as it does from the other result-file readers.
    /// </summary>
    private static bool IsReadFailure(Exception e) =>
        e is InvalidDataException or UnauthorizedAccessException
        || e is IOException and not FileNotFoundException and not DirectoryNotFoundException;

    /// <summary>
    /// The PSH header: each column by its trimmed name (MaxQuant and HUPO examples leave trailing spaces).
    /// </summary>
    private sealed class PsmColumns
    {
        private readonly Dictionary<string, int> _index = new(StringComparer.OrdinalIgnoreCase);

        public List<string> Optional { get; } = [];
        public Dictionary<int, string> Scores { get; } = [];
        public string? DecoyColumn { get; }
        public string? ScanColumn { get; }

        public PsmColumns(string[] header)
        {
            for (int i = 1; i < header.Length; i++)
            {
                string name = header[i].Trim();
                if (name.Length == 0 || !_index.TryAdd(name, i))
                {
                    continue;
                }

                var score = ScoreColumn.Match(name);
                if (score.Success)
                {
                    Scores[int.Parse(score.Groups[1].Value, CultureInfo.InvariantCulture)] = name;
                }

                if (name.StartsWith("opt_", StringComparison.OrdinalIgnoreCase))
                {
                    Optional.Add(name);
                    // opt_global_cv_MS:1002217_decoy_peptide, opt_cv_..., opt_ms_run[1]_cv_..., and the "decoy__peptide" typo
                    if (DecoyColumn == null && (name.EndsWith("decoy_peptide", StringComparison.OrdinalIgnoreCase)
                                                || name.EndsWith("decoy__peptide", StringComparison.OrdinalIgnoreCase)))
                    {
                        DecoyColumn = name;
                    }

                    if (ScanColumn == null && name.EndsWith(ScanColumnSuffix, StringComparison.OrdinalIgnoreCase))
                    {
                        ScanColumn = name;
                    }
                }
            }
        }

        /// <summary>The cell, trimmed; null when the column is missing, the row is short, or the cell is empty or null.</summary>
        public string? Value(string[] row, string column)
        {
            if (!_index.TryGetValue(column, out int i) || i >= row.Length)
            {
                return null;
            }

            string cell = row[i].Trim();
            return cell.Length == 0 || cell.Equals("null", StringComparison.OrdinalIgnoreCase) ? null : cell;
        }

        public (string, string, string, string) KeyOf(string[] row) =>
            (Value(row, "PSM_ID") ?? "", Value(row, "spectra_ref") ?? "", Value(row, "sequence") ?? "", Value(row, "modifications") ?? "");
    }
}
