using System.Globalization;
using Omics.Modifications;
using Omics.SequenceConversion;

namespace Readers
{
    /// <summary>
    /// Places a modification given at a PSI location into an AllModsOneIsNterminus dictionary. mzIdentML and mzTab
    /// number locations the same way: 0 is the N-terminus, 1 to the peptide length a residue, length + 1 the
    /// C-terminus.
    /// </summary>
    internal static class PsiModificationLocation
    {
        /// <summary>
        /// Null when <paramref name="location"/> is on <paramref name="baseSequence"/>, from the N-terminus (0) to the
        /// C-terminus (length + 1); otherwise why it is not.
        /// </summary>
        internal static string? ReasonNotOnPeptide(string baseSequence, int location) =>
            location < 0 || location > baseSequence.Length + 1
                ? $"modification location {location} is not on the peptide"
                : null;

        /// <summary>
        /// Resolves the Unimod modification <paramref name="unimodAccession"/> at <paramref name="location"/> and adds
        /// it to <paramref name="mods"/>.
        /// </summary>
        /// <param name="mods">Modifications keyed the AllModsOneIsNterminus way. Unchanged when this returns false.</param>
        /// <param name="baseSequence">The unmodified peptide, not empty. A terminal modification's residue is taken
        /// from it.</param>
        /// <param name="location">The PSI location.</param>
        /// <param name="unimodAccession">The accession as written, such as "UNIMOD:35", used in the reason.</param>
        /// <param name="unimodId">The accession's number.</param>
        /// <param name="massDelta">The mass shift the file gives, if any.</param>
        /// <param name="lookup">The Unimod catalog to resolve against.</param>
        /// <param name="reason">Null when the modification was added; otherwise why it was not.</param>
        internal static bool TryAddUnimod(Dictionary<int, Modification> mods, string baseSequence, int location,
            string unimodAccession, int unimodId, double? massDelta, UnimodModificationLookup lookup, out string? reason)
        {
            reason = ReasonNotOnPeptide(baseSequence, location);
            if (reason != null)
            {
                return false;
            }

            // AllModsOneIsNterminus: 1 is the N-terminus, residue r (zero-based) is r + 2, the C-terminus length + 2.
            // Both numberings count the N-terminus as a position, so every key is the location plus one.
            var canonical = location switch
            {
                0 => CanonicalModification.AtNTerminus(unimodAccession, baseSequence[0], massDelta, unimodId: unimodId),
                _ when location == baseSequence.Length + 1 => CanonicalModification.AtCTerminus(unimodAccession, baseSequence[^1], massDelta, unimodId: unimodId),
                _ => CanonicalModification.AtResidue(location - 1, baseSequence[location - 1], unimodAccession, massDelta, unimodId: unimodId),
            };
            int key = location + 1;

            var resolved = lookup.TryResolve(canonical)?.MzLibModification;
            if (resolved == null)
            {
                reason = $"{unimodAccession} at location {location} does not resolve to a Unimod modification";
                return false;
            }

            if (!mods.TryAdd(key, resolved))
            {
                reason = $"more than one modification at location {location}";
                return false;
            }

            return true;
        }

        /// <summary>
        /// The tolerance <see cref="UnimodModificationLookup"/> matches masses with by default.
        /// </summary>
        internal const double MassTolerance = 0.01;

        /// <summary>
        /// Resolves a modification given only by its mass shift (mzTab's "CHEMMOD:15.9949") at <paramref name="location"/>
        /// and adds it to <paramref name="mods"/>.
        /// </summary>
        /// <remarks>
        /// The Unimod entries within <see cref="MassTolerance"/> of the mass that can sit at the location are found
        /// first: on the residue, or on X, and, at a terminus, restricted to that terminus when any such entry exists.
        /// They must name exactly one Unimod id, which is then resolved as <see cref="TryAddUnimod"/> resolves it.
        /// The lookup's own mass match is not used, because it stops at the first residue match: it would make an
        /// N-terminal acetyl on a peptide starting with S the side-chain "Acetyl on S".
        /// </remarks>
        /// <param name="mods">Modifications keyed the AllModsOneIsNterminus way. Unchanged when this returns false.</param>
        /// <param name="baseSequence">The unmodified peptide, not empty.</param>
        /// <param name="location">The PSI location.</param>
        /// <param name="representation">The modification as written, such as "CHEMMOD:15.9949", used in the reason.</param>
        /// <param name="mass">The monoisotopic mass shift.</param>
        /// <param name="lookup">The Unimod catalog to resolve against.</param>
        /// <param name="reason">Null when the modification was added; otherwise why it was not.</param>
        internal static bool TryAddByMass(Dictionary<int, Modification> mods, string baseSequence, int location,
            string representation, double mass, UnimodModificationLookup lookup, out string? reason)
        {
            reason = ReasonNotOnPeptide(baseSequence, location);
            if (reason != null)
            {
                return false;
            }

            bool nTerminus = location == 0;
            bool cTerminus = location == baseSequence.Length + 1;
            char residue = nTerminus ? baseSequence[0] : cTerminus ? baseSequence[^1] : baseSequence[location - 1];

            var near = Mods.UnimodModifications
                .Where(m => m.MonoisotopicMass is double delta && Math.Abs(delta - mass) <= MassTolerance)
                .Where(m => m.Target?.Motif == residue.ToString() || m.Target?.Motif == "X")
                .ToList();

            List<Modification> placed;
            if (nTerminus || cTerminus)
            {
                var terminal = near.Where(m => RestrictedTo(m, nTerminus ? "N-terminal" : "C-terminal")).ToList();
                placed = terminal.Count > 0 ? terminal : near.Where(AllowedAnywhere).ToList();
            }
            else
            {
                placed = near.Where(AllowedAnywhere).ToList();
                if (placed.Count == 0)
                {
                    // a residue-numbered terminal modification, such as pyro-Glu written at location 1
                    placed = near.Where(m => (location == 1 && RestrictedTo(m, "N-terminal"))
                                             || (location == baseSequence.Length && RestrictedTo(m, "C-terminal"))).ToList();
                }
            }

            var ids = placed.Select(UnimodIdOf).OfType<int>().Distinct().OrderBy(id => id).ToList();
            if (ids.Count == 0)
            {
                reason = $"{representation} at location {location} matches no Unimod modification within {MassTolerance} Da";
                return false;
            }

            if (ids.Count > 1)
            {
                reason = $"{representation} at location {location} matches more than one Unimod modification ("
                         + string.Join(", ", ids.Select(id => $"UNIMOD:{id}")) + ")";
                return false;
            }

            return TryAddUnimod(mods, baseSequence, location, $"UNIMOD:{ids[0]}", ids[0], mass, lookup, out reason);
        }

        private static bool AllowedAnywhere(Modification m) =>
            string.IsNullOrWhiteSpace(m.LocationRestriction)
            || m.LocationRestriction.Equals("Anywhere.", StringComparison.OrdinalIgnoreCase)
            || m.LocationRestriction.Equals("Unassigned.", StringComparison.OrdinalIgnoreCase);

        private static bool RestrictedTo(Modification m, string terminus) =>
            m.LocationRestriction?.Contains(terminus, StringComparison.OrdinalIgnoreCase) == true;

        /// <summary>
        /// The Unimod id an entry names for itself: the first value of its "Unimod" database reference, else a
        /// "UNIMOD:n" accession. Omics reads it the same way, internally (CanonicalModification.GetUnimodId).
        /// </summary>
        private static int? UnimodIdOf(Modification m)
        {
            var key = m.DatabaseReference?.Keys.FirstOrDefault(k => k.Equals("UNIMOD", StringComparison.OrdinalIgnoreCase));
            string? id = key != null && m.DatabaseReference![key].Count > 0 ? m.DatabaseReference[key][0] : m.Accession;
            if (string.IsNullOrWhiteSpace(id))
            {
                return null;
            }

            id = id.Replace("UNIMOD:", "", StringComparison.OrdinalIgnoreCase).Replace(":", "").Trim();
            return int.TryParse(id, NumberStyles.None, CultureInfo.InvariantCulture, out int parsed) ? parsed : null;
        }
    }
}
