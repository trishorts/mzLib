using MzLibUtil;
using Omics.Modifications;

namespace Readers.ExternalResults.IndividualResultRecords;

/// <summary>
/// One peptide-spectrum match read from the PSM section of an mzTab 1.0 file by
/// <see cref="ResultFiles.MzTabFile"/>.
/// </summary>
/// <remarks>
/// <para>
/// A PSM is the rows that share PSM_ID, spectra_ref, sequence and modifications. Producers write one row per protein
/// and do not always keep those rows together, so <see cref="ResultFiles.MzTabFile"/>'s Results merges them across the
/// whole file. See <see cref="ResultFiles.MzTabFile.EnumeratePsms"/> for the streaming path, which does not.
/// </para>
/// <para>
/// Every cell is kept as written next to its typed reading, because producers fill the columns differently: MaxQuant
/// writes retention times in minutes and Mascot writes a range; Proteome Discoverer writes null.
/// </para>
/// </remarks>
public class MzTabPsm : ISpectralMatch
{
    private readonly List<MzTabPsmProtein> _proteins = [];

    /// <summary>PSM_ID as written. Not always an integer (Proteome Discoverer), and not unique on its own: MassIVE
    /// reuses one for different spectra, and PRIDE's converter for different candidate peptides of one spectrum.</summary>
    public string PsmId { get; internal set; } = string.Empty;

    /// <summary>The peptide sequence, as written.</summary>
    public string BaseSequence { get; internal set; } = string.Empty;

    /// <summary>
    /// The peptidoform in mzLib's notation, built from <see cref="AllModsOneIsNterminus"/>. Empty when
    /// <see cref="ModificationsResolved"/> is false, so a partial peptidoform never passes for a complete one.
    /// </summary>
    public string FullSequence { get; internal set; } = string.Empty;

    /// <summary>
    /// The modifications that resolved, keyed the mzLib way: 1 is the N-terminus, residue r (zero-based) is r + 2, the
    /// C-terminus is length + 2. For MaxQuant, which leaves fixed modifications out of the rows, the declared fixed
    /// modifications are added at every unmodified site.
    /// </summary>
    public Dictionary<int, Modification> AllModsOneIsNterminus { get; internal set; } = [];

    /// <summary>False when any entry of the modifications cell, or any fixed modification added for the producer,
    /// could not be placed on the peptide. <see cref="UnresolvedModifications"/> says which and why.</summary>
    public bool ModificationsResolved => UnresolvedModifications.Count == 0;

    /// <summary>The modifications that did not resolve, with reasons. Empty when all resolved.</summary>
    public IReadOnlyList<MzTabUnresolvedModification> UnresolvedModifications { get; internal set; } = [];

    /// <summary>The modifications cell as written, or null.</summary>
    public string? Modifications { get; internal set; }

    /// <summary>The proteins the PSM maps to, in the order the file first gives them.</summary>
    public IReadOnlyList<MzTabPsmProtein> Proteins => _proteins;

    /// <summary>The accessions of <see cref="Proteins"/>, in order.</summary>
    public IReadOnlyList<string> Accessions => _proteins.Select(p => p.Accession).ToList();

    /// <summary>
    /// The accessions joined with '|', as MzIdentMLRecord and psmtsv join them. Accessions themselves contain '|'
    /// ("sp|P21889|SYD_ECOLI"), so use <see cref="Accessions"/> to tell them apart.
    /// </summary>
    public string Accession { get; private set; } = string.Empty;

    /// <summary>
    /// True when the opt_ decoy column says 1. Without that column, true when every accession carries a decoy prefix
    /// (REV__ from MaxQuant, XXX_ from MS-GF+, DECOY_, dec_).
    /// </summary>
    public bool IsDecoy { get; private set; }

    /// <summary>unique: true for 1, false for 0, null when not given.</summary>
    public bool? Unique { get; internal set; }

    public string? Database { get; internal set; }

    public string? DatabaseVersion { get; internal set; }

    /// <summary>The search_engine params, in order.</summary>
    public IReadOnlyList<CvParam> SearchEngines { get; internal set; } = [];

    /// <summary>
    /// search_engine_score[n] by n; null where the cell is null. The metadata's psm_search_engine_score[n] names
    /// each one (<see cref="ResultFiles.MzTabMetadata.PsmSearchEngineScores"/>).
    /// </summary>
    public IReadOnlyDictionary<int, double?> SearchEngineScores { get; internal set; } = new Dictionary<int, double?>();

    /// <summary>reliability (1 high, 2 medium, 3 poor), when given.</summary>
    public int? Reliability { get; internal set; }

    /// <summary>The retention_time cell as written: one value, '|'-separated values, or a Mascot range "a-b".</summary>
    public string? RetentionTimeText { get; internal set; }

    /// <summary>
    /// The first value of <see cref="RetentionTimeText"/> (the start of a range), in the producer's unit. The spec
    /// says seconds; MaxQuant writes minutes.
    /// </summary>
    public double? RetentionTime { get; internal set; }

    public int? Charge { get; internal set; }

    public double? ExperimentalMassToCharge { get; internal set; }

    public double? CalculatedMassToCharge { get; internal set; }

    public string? Uri { get; internal set; }

    /// <summary>spectra_ref as written, possibly several '|'-separated references. The first one is used.</summary>
    public string SpectraRef { get; internal set; } = string.Empty;

    /// <summary>The n of the first reference's ms_run[n]; -1 when it names none.</summary>
    public int MsRunIndex { get; internal set; } = -1;

    /// <summary>The nativeID after "ms_run[n]:", such as "scan=5", "index=3508" or "spectrum=2458".</summary>
    public string NativeId { get; internal set; } = string.Empty;

    /// <summary>ms_run[n]-location for <see cref="MsRunIndex"/>, as written; empty when the metadata gives none.</summary>
    public string MsRunLocation { get; internal set; } = string.Empty;

    /// <summary>The file name in <see cref="MsRunLocation"/>, without directory or extension.</summary>
    public string FileNameWithoutExtension { get; internal set; } = string.Empty;

    /// <summary>
    /// The producer's scan number. A MaxQuant opt_..._cv_MS:1000776_scan_number_only_nativeID_format column wins,
    /// because MaxQuant's spectra_ref indexes an .apl peak list while that column has the raw scan. Otherwise it is
    /// read from <see cref="NativeId"/>: "scan=N" and "spectrum=N" give N, "index=N" gives N + 1, anything else -1.
    /// CAREFUL when comparing with mzLib: its Mzml reader numbers spectra by their position in the file, so for a
    /// cut-down or MS2-only mzML this is not the scan number mzLib gives the same spectrum. Look the spectrum up by
    /// <see cref="NativeId"/> (MsDataScan.NativeId) instead.
    /// </summary>
    public int OneBasedScanNumber { get; internal set; } = -1;

    /// <summary>Every opt_ column of the PSM's first row, by header name, as written (null for null).</summary>
    public IReadOnlyDictionary<string, string?> OptionalColumns { get; internal set; } = new Dictionary<string, string?>();

    /// <summary>The record key: PSM_ID, spectra_ref, sequence and modifications as written.</summary>
    internal (string PsmId, string SpectraRef, string Sequence, string Modifications) Key { get; set; }

    /// <summary>The decoy column's verdict for the PSM, or null when the file has no decoy column or leaves it null.</summary>
    internal bool? DecoyColumn { get; set; }

    private static readonly string[] DecoyPrefixes = ["REV__", "XXX_", "DECOY_", "dec_"];

    /// <summary>
    /// Adds the proteins that are not already listed, in order. An exact copy of a row (Proteome Discoverer repeats
    /// whole blocks) adds nothing.
    /// </summary>
    internal void AddProteins(IEnumerable<MzTabPsmProtein> proteins)
    {
        foreach (var protein in proteins)
        {
            if (!_proteins.Contains(protein))
            {
                _proteins.Add(protein);
            }
        }

        Accession = string.Join('|', _proteins.Select(p => p.Accession));
        IsDecoy = DecoyColumn ?? (_proteins.Count > 0 && _proteins.All(p =>
            DecoyPrefixes.Any(prefix => p.Accession.StartsWith(prefix, StringComparison.OrdinalIgnoreCase))));
    }
}
