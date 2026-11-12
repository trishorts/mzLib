using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;
using Readers;

namespace Test.FileReadingTests.ReadersInfrastructure
{
    /// <summary>
    /// <see cref="PsiModificationLocation"/> called directly, with the numbering mzIdentML and mzTab share: 0 is the
    /// N-terminus, 1 to length a residue, length + 1 the C-terminus. MzIdentMLResultFile resolves every modification
    /// through it, so TestMzIdentMLResultFile covers it in context; these cases pin the shared entry point itself.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestPsiModificationLocation
    {
        private static readonly UnimodModificationLookup Lookup = new();

        [TestCase(0, null)]
        [TestCase(1, null)]
        [TestCase(7, null)]
        [TestCase(8, null)]
        [TestCase(-1, "modification location -1 is not on the peptide")]
        [TestCase(9, "modification location 9 is not on the peptide")]
        public void ALocationIsOnThePeptideFromTheNTerminusToTheCTerminus(int location, string expectedReason)
        {
            Assert.That(PsiModificationLocation.ReasonNotOnPeptide("PEPTIDE", location), Is.EqualTo(expectedReason));
        }

        /// <summary>
        /// Every key is the location plus one. The residue a terminal modification sits on comes from the sequence.
        /// </summary>
        [Test]
        public void EachLocationIsKeyedOnePastIt()
        {
            var mods = new Dictionary<int, Modification>();

            Assert.Multiple(() =>
            {
                Assert.That(PsiModificationLocation.TryAddUnimod(mods, "PEPTIDE", 0, "UNIMOD:1", 1, null, Lookup, out var nTerm), Is.True, nTerm);
                Assert.That(PsiModificationLocation.TryAddUnimod(mods, "PEPTIDE", 4, "UNIMOD:21", 21, 79.966331, Lookup, out var residue), Is.True, residue);
                Assert.That(PsiModificationLocation.TryAddUnimod(mods, "PEPTIDE", 8, "UNIMOD:2", 2, null, Lookup, out var cTerm), Is.True, cTerm);

                Assert.That(mods.Keys, Is.EquivalentTo(new[] { 1, 5, 9 }));
                Assert.That(mods[1].IdWithMotif, Is.EqualTo("Acetyl on X"));
                Assert.That(mods[5].IdWithMotif, Is.EqualTo("Phospho on T"));
                Assert.That(mods[9].IdWithMotif, Is.EqualTo("Amidated on X"));
                Assert.That(mods[9].LocationRestriction, Does.Contain("C-terminal"));
            });
        }

        /// <summary>
        /// Location length + 1 on a peptide ending in K is the C-terminal methyl, not the lysine one.
        /// </summary>
        [Test]
        public void TheLastLocationIsTheCTerminusEvenOnALysine()
        {
            var mods = new Dictionary<int, Modification>();

            Assert.That(PsiModificationLocation.TryAddUnimod(mods, "PEPTIDEK", 9, "UNIMOD:34", 34, null, Lookup, out _), Is.True);
            Assert.That(mods[10].IdWithMotif, Is.EqualTo("Methyl on X"));
            Assert.That(mods[10].LocationRestriction, Does.Contain("C-terminal"));
        }

        [TestCase(-1, "UNIMOD:35", 35, "modification location -1 is not on the peptide")]
        [TestCase(9, "UNIMOD:35", 35, "modification location 9 is not on the peptide")]
        [TestCase(2, "UNIMOD:999999", 999999, "UNIMOD:999999 at location 2 does not resolve to a Unimod modification")]
        public void AModificationThatCannotBePlacedIsNotAdded(int location, string accession, int unimodId, string expectedReason)
        {
            var mods = new Dictionary<int, Modification>();

            bool added = PsiModificationLocation.TryAddUnimod(mods, "PEPTIDE", location, accession, unimodId, null, Lookup, out var reason);

            Assert.That(added, Is.False);
            Assert.That(reason, Is.EqualTo(expectedReason));
            Assert.That(mods, Is.Empty);
        }

        [TestCase("PEPMTIDE", 4, 15.9949146221, 5, "Oxidation on M")]
        [TestCase("SAMPLER", 0, 42.0105646863, 1, "Acetyl on X")]           // not the side-chain Acetyl on S
        [TestCase("KAMPLER", 0, 42.0105646863, 1, "Acetyl on X")]           // nor on K
        [TestCase("PEPCTIDE", 4, 57.0214637236, 5, "Carbamidomethyl on C")]
        [TestCase("PEPTIDEK", 8, 42.0105646863, 9, "Acetyl on K")]          // a residue location keeps the side chain
        public void AMassResolvesToTheOneUnimodModificationThatFitsTheSite(string sequence, int location, double mass, int key, string expected)
        {
            var mods = new Dictionary<int, Modification>();

            bool added = PsiModificationLocation.TryAddByMass(mods, sequence, location, $"CHEMMOD:{mass}", mass, Lookup, out var reason);

            Assert.That(added, Is.True, reason);
            Assert.That(mods.Keys, Is.EqualTo(new[] { key }));
            Assert.That(mods[key].IdWithMotif, Is.EqualTo(expected));
            if (location == 0)
            {
                Assert.That(mods[key].LocationRestriction, Does.Contain("N-terminal"));
            }
        }

        [TestCase("PEPTIDE", 2, 1234.5678, "CHEMMOD:x at location 2 matches no Unimod modification within 0.01 Da")]
        [TestCase("PEPTIDE", 9, 15.9949, "modification location 9 is not on the peptide")]
        public void AMassThatCannotBePlacedIsNotAdded(string sequence, int location, double mass, string expectedReason)
        {
            var mods = new Dictionary<int, Modification>();

            bool added = PsiModificationLocation.TryAddByMass(mods, sequence, location, "CHEMMOD:x", mass, Lookup, out var reason);

            Assert.That(added, Is.False);
            Assert.That(reason, Is.EqualTo(expectedReason));
            Assert.That(mods, Is.Empty);
        }

        /// <summary>
        /// Dimethyl and Ethyl have one formula, C2H4, and both sit on K: the mass cannot tell them apart.
        /// </summary>
        [Test]
        public void AMassThatFitsTwoUnimodModificationsAtTheSiteIsNotGuessed()
        {
            var mods = new Dictionary<int, Modification>();

            bool added = PsiModificationLocation.TryAddByMass(mods, "PEPTIDEK", 8, "CHEMMOD:28.0313", 28.0313, Lookup, out var reason);

            Assert.That(added, Is.False);
            Assert.That(reason, Does.StartWith("CHEMMOD:28.0313 at location 8 matches more than one Unimod modification (UNIMOD:"));
            Assert.That(reason, Does.Contain("UNIMOD:36").And.Contain("UNIMOD:280"));
        }

        [Test]
        public void ASecondModificationAtOneLocationIsNotAddedAndTheFirstIsKept()
        {
            var mods = new Dictionary<int, Modification>();
            PsiModificationLocation.TryAddUnimod(mods, "PEPTIDE", 4, "UNIMOD:21", 21, 79.966331, Lookup, out _);

            bool added = PsiModificationLocation.TryAddUnimod(mods, "PEPTIDE", 4, "UNIMOD:1", 1, null, Lookup, out var reason);

            Assert.That(added, Is.False);
            Assert.That(reason, Is.EqualTo("more than one modification at location 4"));
            Assert.That(mods[5].IdWithMotif, Is.EqualTo("Phospho on T"));
        }
    }
}
