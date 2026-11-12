using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;
using Readers.ExternalResults.IndividualResultRecords;
using Readers.ExternalResults.ResultFiles;

namespace Test.FileReadingTests.ReadersInfrastructure
{
    /// <summary>
    /// <see cref="MzTabModifications"/> on inline cells, one per form of the modifications column. The cells are
    /// SYNTHETIC, written to the grammar and to the habits of the producers named. TestMzTabFile covers the same forms
    /// in real rows where a licensed file has them.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestMzTabModifications
    {
        private static readonly UnimodModificationLookup Lookup = new();

        private static MzTabMetadata Metadata(params (string Key, string Value)[] entries) =>
            new(new[] { ("mzTab-version", "1.0.0") }.Concat(entries));

        private static readonly MzTabMetadata Plain = Metadata();

        private static MzTabModifications.Placed Place(string cell, string sequence, MzTabMetadata metadata = null) =>
            MzTabModifications.Place(cell, sequence, metadata ?? Plain, Lookup);

        [TestCase(null)]
        [TestCase("")]
        [TestCase("null")]
        [TestCase("NULL")]
        [TestCase("0")]          // Mascot
        public void NoModificationsAreWrittenSeveralWays(string cell)
        {
            var placed = Place(cell, "PEPMTIDEK");

            Assert.That(placed.Mods, Is.Empty);
            Assert.That(placed.Unresolved, Is.Empty);
        }

        [TestCase("4-UNIMOD:35")]
        [TestCase("4-CHEMMOD:15.9949146221")]   // MaxQuant: unsigned
        [TestCase("4-CHEMMOD:+15.99")]          // the spec: signed
        [TestCase("4-CHEMMOD:O")]               // a formula
        [TestCase("4[MS, MS:1002380, false localization rate, 0.105]-UNIMOD:35")] // a score in the position (quantms)
        public void OxidationOnTheMethionineResolvesFromEveryForm(string cell)
        {
            var placed = Place(cell, "PEPMTIDEK");

            Assert.That(placed.Unresolved, Is.Empty);
            Assert.That(placed.Mods.Keys, Is.EqualTo(new[] { 5 }));
            Assert.That(placed.Mods[5].IdWithMotif, Is.EqualTo("Oxidation on M"));
        }

        /// <summary>
        /// The lookup's own mass match would make this the side-chain "Acetyl on S"; the reader resolves by Unimod id.
        /// </summary>
        [Test]
        public void AMassAtTheNTerminusResolvesToTheNTerminalModification()
        {
            var placed = Place("0-CHEMMOD:42.0105646863", "SAMPLER");

            Assert.That(placed.Mods.Keys, Is.EqualTo(new[] { 1 }));
            Assert.That(placed.Mods[1].IdWithMotif, Is.EqualTo("Acetyl on X"));
            Assert.That(placed.Mods[1].LocationRestriction, Does.Contain("N-terminal"));
        }

        [Test]
        public void ANeutralLossIsIgnored()
        {
            var placed = Place("4-UNIMOD:35,4-[PSI-MS, MS:1001524, fragment neutral loss, 63.998285]", "PEPMTIDEK");

            Assert.That(placed.Unresolved, Is.Empty);
            Assert.That(placed.Mods.Keys, Is.EqualTo(new[] { 5 }));
        }

        [TestCase("4-MOD:00719", "MOD:00719 is a PSI-MOD accession, which this reader does not resolve yet")]
        [TestCase("4-SUBST:R", "SUBST:R is a substitution, not a modification")]
        [TestCase("3|4-UNIMOD:35", "the position 3|4 is ambiguous")]
        [TestCase("UNIMOD:35", "the modification has no position")]
        [TestCase("CHEMMOD:+15.99", "the modification has no position")]
        [TestCase("x-UNIMOD:35", "the position x is not a number")]
        [TestCase("4-FOO:1", "FOO:1 is not a modification accession this reader knows")]
        [TestCase("4-UNIMOD:x", "UNIMOD:x is not a Unimod accession")]
        [TestCase("4-UNIMOD:999999", "UNIMOD:999999 at location 4 does not resolve to a Unimod modification")]
        [TestCase("12-UNIMOD:35", "modification location 12 is not on the peptide")]
        [TestCase("4-CHEMMOD:+16", "CHEMMOD:+16 at location 4 gives the mass with fewer than two decimals, too coarse to name one modification")]
        [TestCase("4-CHEMMOD:+16.5", "CHEMMOD:+16.5 at location 4 gives the mass with fewer than two decimals, too coarse to name one modification")]
        [TestCase("4-CHEMMOD:+1234.5678", "CHEMMOD:+1234.5678 at location 4 matches no Unimod modification within 0.01 Da")]
        [TestCase("4-CHEMMOD:Qq", "CHEMMOD:Qq at location 4 is neither a mass nor a chemical formula")]
        [TestCase("4-[MS, MS:1001460, unknown modification, ]", "a parameter stands where a modification accession should be")]
        public void AnEntryThatCannotBePlacedIsReportedWithItsReason(string cell, string reason)
        {
            var placed = Place(cell, "PEPMTIDEK");

            Assert.That(placed.Mods, Is.Empty);
            Assert.That(placed.Unresolved.Single(), Is.EqualTo(new MzTabUnresolvedModification(cell, reason)));
        }

        [Test]
        public void TheEntriesThatResolveAreKeptBesideTheOnesThatDoNot()
        {
            var placed = Place("4-UNIMOD:35, 6-MOD:00412", "PEPMTIDEK");

            Assert.That(placed.Mods.Keys, Is.EqualTo(new[] { 5 }));
            Assert.That(placed.Unresolved.Select(u => u.Entry), Is.EqualTo(new[] { "6-MOD:00412" }));
        }

        [Test]
        public void TwoModificationsOnOneResidueAreReported()
        {
            var placed = Place("4-UNIMOD:35,4-UNIMOD:1", "PEPMTIDEK");

            Assert.That(placed.Mods.Keys, Is.EqualTo(new[] { 5 }));
            Assert.That(placed.Unresolved.Single().Reason, Is.EqualTo("more than one modification at location 4"));
        }

        /// <summary>
        /// Proteome Discoverer 2.5 writes the N-terminal TMT at location 1, and twice at 1 on a K-start peptide.
        /// </summary>
        [Test]
        public void AnNTerminalModificationWrittenAtTheFirstResidueGoesToTheNTerminus()
        {
            var aStart = Place("9-UNIMOD:737,1-UNIMOD:737", "AEPTIDEPK");
            Assert.That(aStart.Mods.Keys, Is.EquivalentTo(new[] { 1, 10 }));
            Assert.That(aStart.Mods[1].LocationRestriction, Is.EqualTo("Peptide N-terminal."));

            var kStart = Place("1-UNIMOD:737,1-UNIMOD:737", "KPEPTIDE");
            Assert.That(kStart.Unresolved, Is.Empty);
            Assert.That(kStart.Mods.Keys, Is.EquivalentTo(new[] { 1, 2 }));
            Assert.That(kStart.Mods[2].IdWithMotif, Is.EqualTo("TMT11 on K"));

            // a modification allowed anywhere stays on the first residue
            Assert.That(Place("1-UNIMOD:35", "MPEPTIDE").Mods.Keys, Is.EqualTo(new[] { 2 }));
        }

        private static readonly (string, string)[] MaxQuantFixedCarbamidomethyl =
        [
            ("software[1]", "[MS, MS:1001583, MaxQuant,1.6.3.1]"),
            ("fixed_mod[1]", "[, , CHEMMOD:57.0214637236,]"),
            ("fixed_mod[1]-site", "C"),
            ("fixed_mod[1]-position", "Anywhere"),
        ];

        [Test]
        public void MaxQuantsFixedModificationIsAddedAtEachSiteTheRowLeavesUnmodified()
        {
            var metadata = Metadata(MaxQuantFixedCarbamidomethyl);

            // C is residue 4 and residue 9: keys 5 and 10
            var none = Place("null", "PEPCTIDECK", metadata);
            Assert.That(none.Mods.Keys, Is.EquivalentTo(new[] { 5, 10 }));
            Assert.That(none.Mods.Values.Select(m => m.IdWithMotif), Is.All.EqualTo("Carbamidomethyl on C"));

            // a Cys that already carries a modification keeps it
            var carried = Place("4-UNIMOD:39", "PEPCTIDECK", metadata);
            Assert.That(carried.Mods[5].IdWithMotif, Is.EqualTo("Methylthio on C"));
            Assert.That(carried.Mods[10].IdWithMotif, Is.EqualTo("Carbamidomethyl on C"));
        }

        [Test]
        public void TheFixedModificationIsNotAddedForAnotherProducer()
        {
            var metadata = Metadata(MaxQuantFixedCarbamidomethyl.Skip(1).Append(("software[1]", "[MS, MS:1000752, TOPP software, 3.7.0]")).ToArray());

            Assert.That(Place("null", "PEPCTIDECK", metadata).Mods, Is.Empty);
        }

        [Test]
        public void AFixedModificationThatDoesNotResolveFlagsThePsm()
        {
            var metadata = Metadata(MaxQuantFixedCarbamidomethyl.Select(e => e.Item1 == "fixed_mod[1]" ? (e.Item1, "[, , CHEMMOD:+1,]") : e).ToArray());

            var placed = Place("null", "PEPCTIDE", metadata);

            Assert.That(placed.Unresolved.Single().Entry, Is.EqualTo("fixed_mod[1] at location 4"));
        }
    }

    /// <summary>
    /// <see cref="MzTabCvParam"/>: the parameter shapes real files write. SYNTHETIC strings in those shapes.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestMzTabCvParam
    {
        [TestCase("[MS, MS:1001207, Mascot, ]", "MS", "MS:1001207", "Mascot", "")]
        [TestCase("[MS, MS:1001583, MaxQuant,1.6.3.1]", "MS", "MS:1001583", "MaxQuant", "1.6.3.1")]
        [TestCase("[, , CHEMMOD:57.0214637236,]", "", "", "CHEMMOD:57.0214637236", "")]                // MaxQuant user param
        [TestCase("[PSI-MS, MS:1000650, Proteome Discoverer]", "PSI-MS", "MS:1000650", "Proteome Discoverer", "")] // 3 fields
        [TestCase("[MS, MS:1001207, Mascot, \"2.3, build 1\"]", "MS", "MS:1001207", "Mascot", "2.3, build 1")]   // quoted comma
        [TestCase("[MS, MS:1001456, analysis software, Mascot v2.3, build 1]", "MS", "MS:1001456", "analysis software", "Mascot v2.3, build 1")]
        [TestCase("  [MS, MS:1, x, y]  ", "MS", "MS:1", "x", "y")]
        public void AParameterIsRead(string text, string label, string accession, string name, string value)
        {
            Assert.That(MzTabCvParam.TryParse(text, out var param), Is.True);
            Assert.That(param, Is.EqualTo(new MzLibUtil.CvParam(label, accession, name, value)));
        }

        [TestCase(null)]
        [TestCase("")]
        [TestCase("MS, MS:1, x, y")]
        [TestCase("[MS, MS:1]")]
        [TestCase("[MS, MS:1, x, y")]
        public void SomethingElseIsNotAParameter(string text)
        {
            Assert.That(MzTabCvParam.TryParse(text, out _), Is.False);
        }

        [Test]
        public void AListOfParametersIsSplitOnPipesOutsideBrackets()
        {
            var list = MzTabCvParam.ParseList("[MS, MS:1001207, Mascot, ]|[MS, MS:1001583, MaxQuant, 1|2]");

            Assert.That(list.Select(p => p.Name), Is.EqualTo(new[] { "Mascot", "MaxQuant" }));
            Assert.That(list[1].Value, Is.EqualTo("1|2"));
        }
    }

    /// <summary>
    /// <see cref="MzTabMetadata"/> on SYNTHETIC entries in producers' shapes; TestMzTabFile reads real metadata.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestMzTabMetadataEntries
    {
        [Test]
        public void KeysCompareWithoutCaseAndTheFirstValueWins()
        {
            var metadata = new MzTabMetadata([("mzTab-Version", "1.0"), ("mzTag-mode", "Complete"), ("description", "first"), ("description", "second")]);

            Assert.That(metadata.Version, Is.EqualTo("1.0"));       // Proteome Discoverer's capital V
            Assert.That(metadata.Mode, Is.EqualTo("Complete"));     // and its "mzTag-mode"
            Assert.That(metadata.Description, Is.EqualTo("first"));
            Assert.That(metadata.Values["DESCRIPTION"], Is.EqualTo("first"));
        }

        [Test]
        public void MsRunsAreGatheredByIndex()
        {
            var metadata = new MzTabMetadata([
                ("ms_run[2]-location", "file:///b.mzML"),
                ("ms_run[1]-location", "file:///a.mzML"),
                ("ms_run[1]-format", "[MS, MS:1000584, mzML file, ]"),
                ("ms_run[1]-id_format", "[MS, MS:1000768, Thermo nativeID format, ]")]);

            Assert.That(metadata.MsRuns.Keys, Is.EquivalentTo(new[] { 1, 2 }));
            Assert.That(metadata.MsRuns[1], Is.EqualTo(new MzTabMsRun(1, "file:///a.mzML",
                new MzLibUtil.CvParam("MS", "MS:1000584", "mzML file", ""),
                new MzLibUtil.CvParam("MS", "MS:1000768", "Thermo nativeID format", ""))));
            Assert.That(metadata.MsRuns[2].Format, Is.Null);
        }

        [TestCase("2.0.0-M", true)]
        [TestCase("1.0.0", false)]
        [TestCase("1.0 rc5", false)]
        public void MzTabMIsRecognisedByItsVersion(string version, bool isMzTabM)
        {
            Assert.That(new MzTabMetadata([("mzTab-version", version)]).IsMzTabM, Is.EqualTo(isMzTabM));
        }

        [TestCase("[MS, MS:1001583, MaxQuant, 2.0.3.0]", true)]
        [TestCase("[, , MaxQuant, 1.6.3.4]", true)]
        [TestCase("[MS, MS:1000752, TOPP software, 3.7.0]", false)]
        public void MaxQuantIsTheProducerThatLeavesFixedModificationsOut(string software, bool omits)
        {
            Assert.That(new MzTabMetadata([("software[1]", software)]).OmitsFixedModificationsFromPsms, Is.EqualTo(omits));
        }
    }
}
