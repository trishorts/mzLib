using System.Globalization;
using System.Text.RegularExpressions;
using MzLibUtil;

namespace Readers.ExternalResults.ResultFiles;

/// <summary>
/// One ms_run of an mzTab file: the spectra file a PSM's spectra_ref points into.
/// </summary>
/// <param name="Index">The n of ms_run[n].</param>
/// <param name="Location">The file as the producer wrote it, often a path on the submitter's machine
/// ("file:///home/user/BSA1_F1.mzML", "file://c:/users/user/....apl"). Null when the file gives none.</param>
/// <param name="Format">ms_run[n]-format, such as MS:1000584 mzML format.</param>
/// <param name="IdFormat">ms_run[n]-id_format, the nativeID format the spectra_refs use. Producers do not always
/// keep to it: MaxQuant 2.0.3 declares the Thermo format and writes "scan=".</param>
public sealed record MzTabMsRun(int Index, string? Location, CvParam? Format, CvParam? IdFormat);

/// <summary>
/// A fixed_mod[n] or variable_mod[n] declaration from the metadata.
/// </summary>
/// <param name="Index">The n.</param>
/// <param name="Param">The modification, as written. "No fixed modifications searched" is the MS:1002453 param.
/// MaxQuant writes its CHEMMOD in the name: "[, , CHEMMOD:57.0214637236,]".</param>
/// <param name="Site">The -site value: a residue letter, "N-term" or "C-term". Null when not given.</param>
/// <param name="Position">The -position value: "Anywhere", "Protein N-term", "Any N-term", .... Null when not given.</param>
public sealed record MzTabModificationDeclaration(int Index, CvParam? Param, string? Site, string? Position);

/// <summary>
/// The metadata (MTD) section of an mzTab 1.0 file. Reading it does not read the PSM table.
/// </summary>
public sealed class MzTabMetadata
{
    private static readonly Regex MsRunKey = new(@"^ms_run\[(\d+)\]-(location|format|id_format)$", RegexOptions.IgnoreCase | RegexOptions.Compiled);
    private static readonly Regex IndexedKey = new(@"^(software|psm_search_engine_score|fixed_mod|variable_mod)\[(\d+)\](-site|-position)?$", RegexOptions.IgnoreCase | RegexOptions.Compiled);

    private readonly Dictionary<string, string> _values = new(StringComparer.OrdinalIgnoreCase);

    /// <summary>
    /// Every metadata key with its value, as written. Keys compare without case, because producers differ
    /// (Proteome Discoverer writes "mzTab-Version"). A key written twice keeps its first value: Mascot repeats some.
    /// </summary>
    public IReadOnlyDictionary<string, string> Values => _values;

    /// <summary>mzTab-version: "1.0.0", "1.0" (Mascot, PRIDE), "1.0 rc5" (early PRIDE), "2.0.0-M" (mzTab-M).</summary>
    public string? Version => Get("mzTab-version");

    /// <summary>mzTab-mode: "Summary" or "Complete". Also read from Proteome Discoverer's misspelt "mzTag-mode".</summary>
    public string? Mode => Get("mzTab-mode") ?? Get("mzTag-mode");

    /// <summary>mzTab-type: "Identification" or "Quantification".</summary>
    public string? Type => Get("mzTab-type");

    /// <summary>mzTab-ID, when given.</summary>
    public string? Id => Get("mzTab-ID");

    /// <summary>The free-text description, when given.</summary>
    public string? Description => Get("description");

    /// <summary>True for mzTab-M (metabolomics, version 2.0.0-M), a different format that MzTabFile does not read.</summary>
    public bool IsMzTabM => Version?.Trim().EndsWith("-M", StringComparison.OrdinalIgnoreCase) == true;

    /// <summary>The ms_run entries, by index.</summary>
    public IReadOnlyDictionary<int, MzTabMsRun> MsRuns { get; }

    /// <summary>The software[n] params, in index order.</summary>
    public IReadOnlyList<CvParam> Software { get; }

    /// <summary>psm_search_engine_score[n], by n: what each PSM search_engine_score[n] column holds.</summary>
    public IReadOnlyDictionary<int, CvParam> PsmSearchEngineScores { get; }

    /// <summary>The fixed_mod[n] declarations, in index order.</summary>
    public IReadOnlyList<MzTabModificationDeclaration> FixedModifications { get; }

    /// <summary>The variable_mod[n] declarations, in index order.</summary>
    public IReadOnlyList<MzTabModificationDeclaration> VariableModifications { get; }

    /// <summary>
    /// True when the producer leaves its fixed modifications out of the PSM rows, against the spec. MaxQuant does, in
    /// every version seen (1.6.3.1, 1.6.3.4, 2.0.3.0): none of 5,581 Cys PSMs in PRIDE PXD011194 carries the declared
    /// fixed Carbamidomethyl. OpenMS, quantms and MassIVE's converter carry theirs.
    /// </summary>
    internal bool OmitsFixedModificationsFromPsms =>
        Software.Any(s => s.Accession == "MS:1001583" || s.Name.Contains("MaxQuant", StringComparison.OrdinalIgnoreCase));

    internal MzTabMetadata(IEnumerable<(string Key, string Value)> entries)
    {
        foreach (var (key, value) in entries)
        {
            _values.TryAdd(key, value);
        }

        var msRuns = new SortedDictionary<int, (string? Location, CvParam? Format, CvParam? IdFormat)>();
        var software = new SortedDictionary<int, CvParam>();
        var scores = new SortedDictionary<int, CvParam>();
        var fixedMods = new SortedDictionary<int, (CvParam? Param, string? Site, string? Position)>();
        var variableMods = new SortedDictionary<int, (CvParam? Param, string? Site, string? Position)>();

        foreach (var (key, value) in _values)
        {
            var run = MsRunKey.Match(key);
            if (run.Success)
            {
                int n = int.Parse(run.Groups[1].Value, CultureInfo.InvariantCulture);
                msRuns.TryGetValue(n, out var entry);
                switch (run.Groups[2].Value.ToLowerInvariant())
                {
                    case "location": entry.Location = value; break;
                    case "format": entry.Format = ParamOrNull(value); break;
                    default: entry.IdFormat = ParamOrNull(value); break;
                }
                msRuns[n] = entry;
                continue;
            }

            var indexed = IndexedKey.Match(key);
            if (!indexed.Success)
            {
                continue;
            }

            int index = int.Parse(indexed.Groups[2].Value, CultureInfo.InvariantCulture);
            string suffix = indexed.Groups[3].Value.ToLowerInvariant();
            switch (indexed.Groups[1].Value.ToLowerInvariant())
            {
                case "software" when suffix.Length == 0 && MzTabCvParam.TryParse(value, out var tool):
                    software[index] = tool;
                    break;
                case "psm_search_engine_score" when suffix.Length == 0 && MzTabCvParam.TryParse(value, out var score):
                    scores[index] = score;
                    break;
                case "fixed_mod":
                    fixedMods[index] = WithPart(fixedMods.GetValueOrDefault(index), suffix, value);
                    break;
                case "variable_mod":
                    variableMods[index] = WithPart(variableMods.GetValueOrDefault(index), suffix, value);
                    break;
            }
        }

        MsRuns = msRuns.ToDictionary(p => p.Key, p => new MzTabMsRun(p.Key, p.Value.Location, p.Value.Format, p.Value.IdFormat));
        Software = software.Values.ToList();
        PsmSearchEngineScores = new Dictionary<int, CvParam>(scores);
        FixedModifications = fixedMods.Select(p => new MzTabModificationDeclaration(p.Key, p.Value.Param, p.Value.Site, p.Value.Position)).ToList();
        VariableModifications = variableMods.Select(p => new MzTabModificationDeclaration(p.Key, p.Value.Param, p.Value.Site, p.Value.Position)).ToList();
    }

    /// <summary>
    /// Reads the metadata block: the MTD lines before the first line of another section. COM and blank lines are
    /// skipped. Nothing after the block is read.
    /// </summary>
    internal static MzTabMetadata Read(TextReader reader)
    {
        var entries = new List<(string, string)>();
        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            if (!TryReadMetadataLine(line, entries, out bool endOfBlock) && endOfBlock)
            {
                break;
            }
        }

        return new MzTabMetadata(entries);
    }

    /// <summary>
    /// Adds <paramref name="line"/> to <paramref name="entries"/> when it is an MTD line. A COM or blank line is
    /// skipped; any other line ends the metadata block.
    /// </summary>
    internal static bool TryReadMetadataLine(string line, List<(string, string)> entries, out bool endOfBlock)
    {
        endOfBlock = false;
        if (string.IsNullOrWhiteSpace(line))
        {
            return false;
        }

        var fields = line.Split('\t');
        switch (fields[0].Trim())
        {
            case "MTD":
                if (fields.Length > 1 && fields[1].Trim().Length > 0)
                {
                    entries.Add((fields[1].Trim(), fields.Length > 2 ? fields[2].Trim() : ""));
                }
                return true;
            case "COM":
                return false;
            default:
                endOfBlock = true;
                return false;
        }
    }

    private string? Get(string key) => _values.TryGetValue(key, out var value) && value.Length > 0 ? value : null;

    private static CvParam? ParamOrNull(string value) => MzTabCvParam.TryParse(value, out var param) ? param : null;

    private static (CvParam? Param, string? Site, string? Position) WithPart(
        (CvParam? Param, string? Site, string? Position) declaration, string suffix, string value)
    {
        string? text = value.Length > 0 ? value : null;
        return suffix switch
        {
            "-site" => (declaration.Param, text, declaration.Position),
            "-position" => (declaration.Param, declaration.Site, text),
            _ => (ParamOrNull(value), declaration.Site, declaration.Position),
        };
    }
}
