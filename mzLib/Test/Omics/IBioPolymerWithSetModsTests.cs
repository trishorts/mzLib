using NUnit.Framework;
using Omics;
using Omics.Modifications;
using System.Collections.Generic;

namespace Test.Omics;

[TestFixture]
public class IBioPolymerWithSetModsTests
{
    [Test]
    [TestCase("AC[Phospho]DEFGH", "ACDEFGH", '[', ']')]
    [TestCase("A[Oxidation]C[Phospho]DEFGH", "ACDEFGH", '[', ']')]
    [TestCase("A[Oxidation]C[Phospho]D[Carbamidomethyl]EFGH", "ACDEFGH", '[', ']')]
    [TestCase("A[Oxidation]C[Phospho]D[Carbamidomethyl]E[Acetyl]FGH", "ACDEFGH", '[', ']')]
    [TestCase("A[Oxidation]C[Phospho]D[Carbamidomethyl]E[Acetyl]F[Phospho]GH", "ACDEFGH", '[', ']')]
    [TestCase("A[Oxidation]C[Phospho]D[Carbamidomethyl]E[Acetyl]F[Phospho]G[Oxidation]H", "ACDEFGH", '[', ']')]
    [TestCase("AC{Phospho}DEFGH", "ACDEFGH", '{', '}')]
    [TestCase("A<Oxidation>C<Phospho>DEFGH", "ACDEFGH", '<', '>')]
    [TestCase("A(Oxidation)C(Phospho)DEFGH", "ACDEFGH", '(', ')')]
    public void TestGetBaseSequenceFromFullSequence(string fullSequence, string expectedBaseSequence, char startDelimiter, char endDelimiter)
    {
        string actualBaseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSequence, startDelimiter, endDelimiter);
        Assert.That(actualBaseSequence, Is.EqualTo(expectedBaseSequence));
    }

    [Test]
    [TestCase("ABC", new[] { 3 }, new[] { "Oxidation" }, "AB[Oxidation]C")]
    [TestCase("ALANVNIGSLICNVGAGGPAPAAGAAPAGGPAPSTAAAPAEEK", new[] { 13 }, new[] { "Carbamidomethyl" }, "ALANVNIGSLIC[Carbamidomethyl]NVGAGGPAPAAGAAPAGGPAPSTAAAPAEEK")]
    [TestCase("STSFRGGMGSGGLATGIAGGLAGMGGIQNEK", new[] { 6, 25 }, new[] { "Dimethylation", "Oxidation" }, "STSFR[Dimethylation]GGMGSGGLATGIAGGLAGM[Oxidation]GGIQNEK")]
    [TestCase("KDLYANTVLSGGTTMYPGIADR", new[] { 7, 16 }, new[] { "Deamidation", "Oxidation" }, "KDLYAN[Deamidation]TVLSGGTTM[Oxidation]YPGIADR")]
    [TestCase("MEFDLGAALEPTSQKPGVGAGHGGDPK", new[] { 1 }, new[] { "N-acetylmethionine" }, "[N-acetylmethionine]MEFDLGAALEPTSQKPGVGAGHGGDPK")]
    public static void ProFormaOutputWithModificationInBrackets(string sequence, int[] position, string[] modificationName, string expectedResult)
    {
        var mods = new Dictionary<int, Modification>();
        for (int i = 0; i < position.Length; i++)
        {
            mods.Add(position[i], new Modification(modificationName[i]));
        }
        Assert.That(IBioPolymerWithSetMods.DetermineProFormaCompliantSequence(sequence, mods), Is.EqualTo(expectedResult));



    }
}