using System.Globalization;
using Chemistry;
using MzLibUtil;
using Omics.Modifications;
using Omics.SequenceConversion;
using Readers.ExternalResults.IndividualResultRecords;

namespace Readers.ExternalResults.ResultFiles;

/// <summary>
/// Places the entries of an mzTab PSM's modifications cell on its peptide as mzLib modifications.
/// </summary>
/// <remarks>
/// The grammar is "{position}{[param]}-{identifier}", comma-separated (", " for MaxQuant). Each entry is placed or
/// reported, never dropped silently:
/// <list type="table">
/// <item><term>null, NULL, 0, empty</term><description>no modifications</description></item>
/// <item><term>3-UNIMOD:35</term><description>resolved by Unimod id</description></item>
/// <item><term>3-CHEMMOD:15.9949, 3-CHEMMOD:+15.99</term><description>resolved by mass, when the mass has at least two
/// decimals and names exactly one Unimod id at that site</description></item>
/// <item><term>3-CHEMMOD:NH4</term><description>resolved by the formula's mass</description></item>
/// <item><term>9[MS, MS:1002380, ...]-UNIMOD:21</term><description>the score in the position is dropped</description></item>
/// <item><term>3-[MS, MS:1001524, fragment neutral loss, 63.99]</term><description>ignored: an annotation of the
/// spectrum, not a modification of the peptide</description></item>
/// <item><term>MOD:, SUBST:, 3|4-, no position, anything else</term><description>unresolved, with a reason</description></item>
/// </list>
/// </remarks>
internal static class MzTabModifications
{
    private const string NeutralLossAccession = "MS:1001524";
    private static readonly string[] IdentifierPrefixes = ["UNIMOD:", "MOD:", "CHEMMOD:", "SUBST:", "["];

    /// <summary>The placed modifications and the unresolved entries. Shared by every PSM with the same sequence and
    /// cell, so the reader copies <see cref="Mods"/> per record.</summary>
    internal sealed record Placed(Dictionary<int, Modification> Mods, IReadOnlyList<MzTabUnresolvedModification> Unresolved);

    /// <summary>
    /// The modifications of <paramref name="cell"/> placed on <paramref name="baseSequence"/>, and the entries that
    /// could not be. When <paramref name="metadata"/> says the producer leaves fixed modifications out of its rows, the
    /// declared fixed modifications are added at every site the row leaves unmodified.
    /// </summary>
    internal static Placed Place(string? cell, string baseSequence, MzTabMetadata metadata, UnimodModificationLookup lookup)
    {
        var mods = new Dictionary<int, Modification>();
        var unresolved = new List<MzTabUnresolvedModification>();

        if (!IsNone(cell))
        {
            foreach (var raw in MzTabCvParam.SplitTopLevel(cell!, ','))
            {
                string entry = raw.Trim();
                if (entry.Length == 0)
                {
                    continue;
                }

                if (PlaceEntry(entry, baseSequence, mods, lookup) is { } reason)
                {
                    unresolved.Add(new MzTabUnresolvedModification(entry, reason));
                }
            }
        }

        if (metadata.OmitsFixedModificationsFromPsms)
        {
            AddFixedModifications(metadata, baseSequence, mods, unresolved, lookup);
        }

        return new Placed(mods, unresolved.ToArray());
    }

    private static bool IsNone(string? cell)
    {
        string trimmed = cell?.Trim() ?? "";
        return trimmed.Length == 0 || trimmed.Equals("null", StringComparison.OrdinalIgnoreCase) || trimmed == "0";
    }

    /// <summary>Null when the entry was placed, or ignored as a neutral loss; otherwise why it was not placed.</summary>
    private static string? PlaceEntry(string entry, string baseSequence, Dictionary<int, Modification> mods,
        UnimodModificationLookup lookup)
    {
        int split = IdentifierStart(entry);
        string positionText = split < 0 ? "" : entry[..split];
        string identifier = split < 0 ? entry : entry[(split + 1)..];

        if (identifier.StartsWith('['))
        {
            if (MzTabCvParam.TryParse(identifier, out var param) && param.Accession == NeutralLossAccession)
            {
                return null;
            }

            return "a parameter stands where a modification accession should be";
        }

        // A score in the position ("9[MS, MS:1002380, false localization rate, 0.105]") is not part of it
        string position = StripParams(positionText).Trim();
        if (position.Length == 0)
        {
            return "the modification has no position";
        }

        if (position.Contains('|'))
        {
            return $"the position {position} is ambiguous";
        }

        if (!int.TryParse(position, NumberStyles.None, CultureInfo.InvariantCulture, out int location))
        {
            return $"the position {position} is not a number";
        }

        // Proteome Discoverer 2.5 writes N-terminal modifications at location 1, the first residue, never at 0: all
        // 3,944 rows of MSV000096151 have "1-UNIMOD:737" for the N-terminal TMT. A peptide starting with K then carries
        // two entries at 1, its side-chain label and its N-terminal one. The second is the N-terminus.
        bool firstResidueTaken = location == 1 && mods.ContainsKey(2);
        string? reason = Resolve(identifier, baseSequence, location, mods, lookup);
        if (reason != null && firstResidueTaken && Resolve(identifier, baseSequence, 0, mods, lookup) == null)
        {
            return null;
        }

        // A modification only allowed at the N-terminus belongs at the N-terminus, where mzLib keys it
        // (DigestionProduct puts "N-terminal." and "Peptide N-terminal." modifications at 1, not on the residue).
        if (reason == null && location == 1 && !mods.ContainsKey(1) && mods.TryGetValue(2, out var placed)
            && placed.LocationRestriction is "N-terminal." or "Peptide N-terminal.")
        {
            mods.Remove(2);
            mods[1] = placed;
        }

        return reason;
    }

    /// <summary>Null when <paramref name="identifier"/> was placed at <paramref name="location"/>; otherwise why not.</summary>
    private static string? Resolve(string identifier, string baseSequence, int location, Dictionary<int, Modification> mods,
        UnimodModificationLookup lookup)
    {
        if (identifier.StartsWith("UNIMOD:", StringComparison.OrdinalIgnoreCase))
        {
            if (!int.TryParse(identifier["UNIMOD:".Length..], NumberStyles.None, CultureInfo.InvariantCulture, out int unimodId))
            {
                return $"{identifier} is not a Unimod accession";
            }

            return PsiModificationLocation.TryAddUnimod(mods, baseSequence, location, identifier, unimodId, null, lookup, out string? reason)
                ? null
                : reason;
        }

        if (identifier.StartsWith("CHEMMOD:", StringComparison.OrdinalIgnoreCase))
        {
            return PlaceChemMod(identifier, identifier["CHEMMOD:".Length..].Trim(), baseSequence, location, mods, lookup);
        }

        if (identifier.StartsWith("MOD:", StringComparison.OrdinalIgnoreCase))
        {
            return $"{identifier} is a PSI-MOD accession, which this reader does not resolve yet";
        }

        if (identifier.StartsWith("SUBST:", StringComparison.OrdinalIgnoreCase))
        {
            return $"{identifier} is a substitution, not a modification";
        }

        return $"{identifier} is not a modification accession this reader knows";
    }

    private static string? PlaceChemMod(string identifier, string value, string baseSequence, int location,
        Dictionary<int, Modification> mods, UnimodModificationLookup lookup)
    {
        if (double.TryParse(value, NumberStyles.Float, CultureInfo.InvariantCulture, out double mass))
        {
            // MassIVE's converter writes nominal masses ("CHEMMOD:+1"). One Da of uncertainty spans many Unimod
            // entries, and +1 sits within 0.01 Da of Label:13C(1) (1.00335), so a whole number names nothing.
            int point = value.IndexOf('.');
            if (point < 0 || value.Length - point - 1 < 2)
            {
                return $"{identifier} at location {location} gives the mass with fewer than two decimals, too coarse to name one modification";
            }

            return PsiModificationLocation.TryAddByMass(mods, baseSequence, location, identifier, mass, lookup, out string? reason)
                ? null
                : reason;
        }

        // A formula, signed or not: "+NH4", "-H2O"
        int sign = value.StartsWith('-') ? -1 : 1;
        string formulaText = value.TrimStart('+', '-');
        ChemicalFormula formula;
        try
        {
            formula = ChemicalFormula.ParseFormula(formulaText);
        }
        catch (Exception e) when (e is MzLibException or FormatException or ArgumentException or KeyNotFoundException)
        {
            return $"{identifier} at location {location} is neither a mass nor a chemical formula";
        }

        return PsiModificationLocation.TryAddByMass(mods, baseSequence, location, identifier, sign * formula.MonoisotopicMass, lookup, out string? formulaReason)
            ? null
            : formulaReason;
    }

    private static void AddFixedModifications(MzTabMetadata metadata, string baseSequence,
        Dictionary<int, Modification> mods, List<MzTabUnresolvedModification> unresolved, UnimodModificationLookup lookup)
    {
        foreach (var declaration in metadata.FixedModifications)
        {
            // Only residue sites. No MaxQuant file seen declares a terminal fixed modification.
            if (declaration.Param == null || declaration.Site is not { Length: 1 } site || !char.IsLetter(site[0]))
            {
                continue;
            }

            if (declaration.Position != null && !declaration.Position.Equals("Anywhere", StringComparison.OrdinalIgnoreCase))
            {
                continue;
            }

            // MaxQuant writes its CHEMMOD in the param's name: "[, , CHEMMOD:57.0214637236,]"
            string identifier = new[] { declaration.Param.Accession, declaration.Param.Name }
                .FirstOrDefault(s => s.StartsWith("UNIMOD:", StringComparison.OrdinalIgnoreCase)
                                     || s.StartsWith("CHEMMOD:", StringComparison.OrdinalIgnoreCase)) ?? "";
            if (identifier.Length == 0)
            {
                continue;
            }

            for (int i = 0; i < baseSequence.Length; i++)
            {
                int location = i + 1;
                if (char.ToUpperInvariant(baseSequence[i]) != char.ToUpperInvariant(site[0]) || mods.ContainsKey(location + 1))
                {
                    continue;
                }

                if (PlaceEntry($"{location}-{identifier}", baseSequence, mods, lookup) is { } reason)
                {
                    unresolved.Add(new MzTabUnresolvedModification($"fixed_mod[{declaration.Index}] at location {location}", reason));
                }
            }
        }
    }

    /// <summary>
    /// The index of the '-' that separates the position from the identifier: the first one outside brackets that is
    /// followed by a known identifier prefix. 0 or more is that index; -1 means the entry has no position, as in
    /// "UNIMOD:35" or "CHEMMOD:+159.93".
    /// </summary>
    private static int IdentifierStart(string entry)
    {
        int depth = 0;
        for (int i = 0; i < entry.Length; i++)
        {
            char c = entry[i];
            if (c == '[')
            {
                depth++;
            }
            else if (c == ']' && depth > 0)
            {
                depth--;
            }
            else if (c == '-' && depth == 0 && IdentifierPrefixes.Any(p => entry.Length - (i + 1) >= p.Length
                         && string.Compare(entry, i + 1, p, 0, p.Length, StringComparison.OrdinalIgnoreCase) == 0))
            {
                return i;
            }
        }

        // An unknown identifier after a numeric position ("4-FOO:1") still splits there, so the reason names the
        // identifier rather than a missing position
        int dash = entry.IndexOf('-');
        return dash > 0 && entry[..dash].All(char.IsAsciiDigit) ? dash : -1;
    }

    private static string StripParams(string position)
    {
        var kept = new System.Text.StringBuilder();
        int depth = 0;
        foreach (char c in position)
        {
            if (c == '[')
            {
                depth++;
            }
            else if (c == ']' && depth > 0)
            {
                depth--;
            }
            else if (depth == 0)
            {
                kept.Append(c);
            }
        }

        return kept.ToString();
    }
}
