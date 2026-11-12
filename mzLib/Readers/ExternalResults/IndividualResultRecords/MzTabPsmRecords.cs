namespace Readers.ExternalResults.IndividualResultRecords;

/// <summary>
/// One protein an mzTab PSM maps to, with the residues around the peptide in it. A producer may write each protein
/// on its own row (MaxQuant, PRIDE's converter) or join them in one cell with commas (OpenMS, quantms); either way
/// each protein is one of these.
/// </summary>
/// <param name="Accession">The accession as written, such as "sp|P21889|SYD_ECOLI" or "REV__Q06822".</param>
/// <param name="Pre">The residue before the peptide, "-" at the protein's N-terminus. Null when not given.</param>
/// <param name="Post">The residue after the peptide, "-" at the protein's C-terminus. Null when not given.</param>
/// <param name="Start">The one-based position of the peptide's first residue in the protein, when given.</param>
/// <param name="End">The one-based position of its last residue, when given.</param>
public sealed record MzTabPsmProtein(string Accession, string? Pre, string? Post, int? Start, int? End);

/// <summary>
/// An entry of a PSM's modifications cell that could not be placed on the peptide as an mzLib modification, and why.
/// A fixed modification that MzTabFile adds itself, for a producer that leaves them out, is reported as
/// "fixed_mod[n] at location L".
/// </summary>
/// <param name="Entry">The entry as written, such as "3-MOD:00719" or "1-CHEMMOD:+25".</param>
/// <param name="Reason">Why it is unresolved.</param>
public sealed record MzTabUnresolvedModification(string Entry, string Reason);

/// <summary>
/// A PSM row that <see cref="ResultFiles.MzTabFile"/> could not read as a PSM at all, such as one with no sequence.
/// A PSM whose modifications do not resolve is not skipped: it is kept, with
/// <see cref="MzTabPsm.ModificationsResolved"/> false.
/// </summary>
public sealed record MzTabSkippedPsm(string? PsmId, string? SpectraRef, string? Sequence, string? Modifications, string Reason);
