using System.Collections.Generic;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// Pins the PSI-MS accessions in proteases.tsv against the embedded psi-ms.obo.
    ///
    /// Nothing checked these before, and three of them were wrong: subtilisin claimed MS:1001312,
    /// which is TrypChymo; singleN and singleC claimed MS:1001957 and MS:1001958, which are not
    /// cleavage agents at all but the literal regular-expression strings PSI-MS hangs off them
    /// (is_a MS:1001180, "Cleavage agent regular expression"). Each was individually plausible and
    /// silently reached every mzIdentML MetaMorpheus wrote.
    ///
    /// The check that catches all of them is comparing the NAME column against the vocabulary.
    /// Resolving the accession is not enough and was never going to be: all three wrong accessions
    /// ARE real PSI-MS terms, so an existence check passes on every one of them. What no one does
    /// is write a wrong accession and then write that accession's true name beside it, so the pair
    /// disagreeing is the signal.
    ///
    /// The trypsin pair is pinned by motif rather than by name because it moved twice: the file
    /// once used "|P" for the proline-RESTRICTED enzyme, the inverse of the /P convention every
    /// other search engine uses, and the accessions did not follow when that was repaired. Read
    /// the term's has_regexp against the Motif column, never the spelling of the row name.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestProteaseAccessions
    {
        private static IEnumerable<Protease> Accessioned => ProteaseDictionary.Dictionary.Values
            .Where(p => !string.IsNullOrWhiteSpace(p.PsiMsAccessionNumber));

        [Test]
        public static void EveryAccessionResolvesInThePinnedVocabulary()
        {
            var unresolved = Accessioned
                .Where(p => !ControlledVocabulary.PsiMs.TryGetByAccession(p.PsiMsAccessionNumber, out _))
                .Select(p => $"{p.Name} -> {p.PsiMsAccessionNumber}")
                .ToList();

            Assert.That(unresolved, Is.Empty,
                "proteases.tsv names accessions that are not in psi-ms.obo " +
                $"{ControlledVocabulary.PsiMs.Version}: {string.Join(", ", unresolved)}");
        }

        [Test]
        public static void EveryAccessionCarriesThatTermsOwnName()
        {
            var mismatched = new List<string>();

            foreach (var protease in Accessioned)
            {
                Assert.That(ControlledVocabulary.PsiMs.TryGetByAccession(
                    protease.PsiMsAccessionNumber, out CvParam term), Is.True,
                    $"{protease.Name}: {protease.PsiMsAccessionNumber} does not resolve");

                if (!string.Equals(protease.PsiMsName, term.Name, System.StringComparison.Ordinal))
                    mismatched.Add(
                        $"{protease.Name}: file says \"{protease.PsiMsName}\", " +
                        $"{protease.PsiMsAccessionNumber} is \"{term.Name}\"");
            }

            Assert.That(mismatched, Is.Empty,
                "the accession and the name disagree, so one of them is wrong: " +
                string.Join(" | ", mismatched));
        }

        /// <summary>
        /// Pinned against the motif each term actually describes, not against the row name. The
        /// bare name is the proline-restricted enzyme and /P is the variant that ignores the
        /// restriction, so the two trypsin rows take DIFFERENT accessions -- they carried the same
        /// one, MS:1001313, until this test.
        /// </summary>
        [TestCase("trypsin", "MS:1001251", "Trypsin", Description = "K[P]|,R[P]| blocks P, so Trypsin (?<=[KR])(?!P)")]
        [TestCase("trypsin/P", "MS:1001313", "Trypsin/P", Description = "K|,R| cleaves before P, so Trypsin/P (?<=[KR])")]
        [TestCase("Lys-C", "MS:1001309", "Lys-C", Description = "K[P]| blocks P, so Lys-C (?<=K)(?!P), not Lys-C/P MS:1001310")]
        [TestCase("chymotrypsin", "MS:1001306", "Chymotrypsin")]
        [TestCase("Lys-N", "MS:1003093", "Lys-N", Description = "|K is (?=K) exactly")]
        [TestCase("Glu-C", "MS:1001315", "V8-E", Description = "E| alone")]
        [TestCase("Glu-C (with asp)", "MS:1001314", "V8-DE", Description = "E|,D| -- the D is the whole difference")]
        [TestCase("peptidomics", "MS:1001955", "no cleavage")]
        [TestCase("top-down", "MS:1001955", "no cleavage")]
        [TestCase("non-specific", "MS:1001956", "unspecific cleavage")]
        public static void TheCleavageRuleAndTheTermAgree(string protease, string accession, string name)
        {
            var entry = ProteaseDictionary.Dictionary[protease];

            Assert.That(entry.PsiMsAccessionNumber, Is.EqualTo(accession));
            Assert.That(entry.PsiMsName, Is.EqualTo(name));
        }

        /// <summary>
        /// PSI-MS has no term for these agents, verified against the embedded vocabulary and
        /// against OLS4 on 2026-08-20. Blank is the correct entry: SDRF and mzIdentML both permit
        /// the name alone, and an accession that resolves to a different agent is worse than none
        /// because a consumer cannot tell it is wrong.
        ///
        /// A negative pinned deliberately. If PSI-MS later adds one of these, this fails and the
        /// mapping gets made on purpose rather than by whoever notices first.
        /// </summary>
        [TestCase("elastase", Description = "MS:1001915 is leukocyte elastase, (?<=[ALIV])(?!P); this motif is 14 residues")]
        [TestCase("subtilisin", Description = "MS:1001312 is TrypChymo, not Subtilisin")]
        [TestCase("tryptophan oxidation", Description = "MS:1001918 2-iodobenzoate cleaves after W, but is a specific reagent")]
        [TestCase("collagenase")]
        [TestCase("StcE-trypsin")]
        [TestCase("ProAlanase")]
        [TestCase("singleN", Description = "MS:1001957 is a regular-expression term, not a cleavage agent")]
        [TestCase("singleC", Description = "MS:1001958 is a regular-expression term, not a cleavage agent")]
        public static void AgentsWithNoPsiMsTermCarryNoAccession(string protease)
        {
            Assert.That(ProteaseDictionary.Dictionary[protease].PsiMsAccessionNumber,
                Is.Empty.Or.Null,
                $"{protease} acquired an accession; if PSI-MS added a term, verify its has_regexp " +
                "against the Motif column before accepting it");
        }
    }
}
