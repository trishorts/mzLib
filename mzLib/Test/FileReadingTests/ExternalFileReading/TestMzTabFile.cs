using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using Readers.ExternalResults.IndividualResultRecords;
using Readers.ExternalResults.ResultFiles;

namespace Test.FileReadingTests.ExternalFileReading
{
    /// <summary>
    /// MzTabFile over one real file per producer, in FileReadingTests\ExternalFileTypes\MzTab (sources, licences and
    /// selection rules in MZTAB_PROVENANCE.md there). Expected values were read out of each file's rows, not off the
    /// reader:
    /// <list type="bullet">
    /// <item>OpenMS 3.7 ProteomicsLFQ_1_out (whole; BSD-3): spectrum= refs, UNIMOD, comma-joined accessions;</item>
    /// <item>MaxQuant 1.6.3.1, PRIDE PXD011194 (CC0): CHEMMOD, one row per protein, the opt scan column, REV__ decoys,
    /// fixed modifications left out of the rows;</item>
    /// <item>Proteome Discoverer 2.5, MassIVE MSV000096151 (CC0): exact copies of a PSM far apart, N-terminal TMT at
    /// location 1;</item>
    /// <item>MassIVE's converter, MSV000078886 (CC0): nominal CHEMMOD masses, XXX_ decoys with no decoy column, a PSM_ID
    /// reused for two spectra;</item>
    /// <item>Mascot 2.3 through jmzTab, small.mztab (Apache-2.0): neutral-loss params, bare-number spectra_refs;</item>
    /// <item>PRIDE's XML converter through jmzTab, PRIDE_Example (Apache-2.0): PSI-MOD, a PSM split across proteins, two
    /// candidate peptides under one PSM_ID.</item>
    /// </list>
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestMzTabFile
    {
        private const string OpenMs = "OpenMS3.7_ProteomicsLFQ_1_out.mzTab";
        private const string MaxQuant = "MaxQuant1.6.3.1_PXD011194_excerpt.mzTab";
        private const string ProteomeDiscoverer = "ProteomeDiscoverer2.5_MSV000096151_excerpt.mzTab";
        private const string MassIve = "MassIVE_MSV000078886_excerpt.mzTab";
        private const string Mascot = "Mascot2.3_jmzTab_small_excerpt.mzTab";
        private const string PrideConverter = "PrideConverter_jmzTab_PRIDE_Example.mzTab";
        private const string MzTabM = "pyteomics_mzTabM_test_mztab2.mzTab";

        private static string Fixture(string name) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes", "MzTab", name);

        private static MzTabFile Read(string name) => new(Fixture(name));

        private static MzTabPsm Psm(MzTabFile file, string psmId, string sequence) =>
            file.Results.Single(p => p.PsmId == psmId && p.BaseSequence == sequence);

        private string _outputDirectory;

        [OneTimeSetUp]
        public void OneTimeSetUp()
        {
            _outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "MzTabFileTests");
            Directory.CreateDirectory(_outputDirectory);
        }

        [OneTimeTearDown]
        public void OneTimeTearDown()
        {
            if (Directory.Exists(_outputDirectory))
            {
                Directory.Delete(_outputDirectory, true);
            }
        }

        private string Gzipped(string name)
        {
            string path = Path.Combine(_outputDirectory, name + ".gz");
            using var source = File.OpenRead(Fixture(name));
            using var target = File.Create(path);
            using var gzip = new GZipStream(target, CompressionLevel.Optimal);
            source.CopyTo(gzip);
            return path;
        }

        /// <summary>
        /// Records are the distinct (PSM_ID, spectra_ref, sequence, modifications) keys; streamed records are the runs of
        /// adjacent rows with one key. Both counted from the rows. No fixture row lacks a sequence.
        /// </summary>
        [TestCase(OpenMs, 115, 115)]
        [TestCase(MaxQuant, 49, 49)]
        [TestCase(ProteomeDiscoverer, 17, 19)]
        [TestCase(MassIve, 41, 41)]
        [TestCase(Mascot, 36, 36)]
        [TestCase(PrideConverter, 212, 220)]
        public void ResultsMergeAPsmsRowsAcrossTheFileAndStreamingMergesOnlyAdjacentRows(string name, int records, int streamed)
        {
            var file = Read(name);

            Assert.That(file.Results, Has.Count.EqualTo(records));
            Assert.That(file.SkippedPsms, Is.Empty);
            Assert.That(file.EnumeratePsms().Count(), Is.EqualTo(streamed));
            Assert.That(file.Software, Is.EqualTo(Software.MzTab));
            Assert.That(file.FileType, Is.EqualTo(SupportedFileType.MzTab));
        }

        [Test]
        public void TheFirstOpenMsPsmIsReadFieldByField()
        {
            var psm = Read(OpenMs).Results[0];

            Assert.Multiple(() =>
            {
                Assert.That(psm.PsmId, Is.EqualTo("0"));
                Assert.That(psm.BaseSequence, Is.EqualTo("DGDIEAEISR"));
                Assert.That(psm.FullSequence, Is.EqualTo("DGDIEAEISR"));
                Assert.That(psm.Modifications, Is.Null);
                Assert.That(psm.ModificationsResolved, Is.True);
                Assert.That(psm.AllModsOneIsNterminus, Is.Empty);
                Assert.That(psm.Accessions, Is.EqualTo(new[] { "tr|A9EY18|A9EY18_SORC5" }));
                Assert.That(psm.Accession, Is.EqualTo("tr|A9EY18|A9EY18_SORC5"));
                Assert.That(psm.Proteins.Single(), Is.EqualTo(new MzTabPsmProtein("tr|A9EY18|A9EY18_SORC5", "K", "E", 90, 99)));
                Assert.That(psm.Unique, Is.True);
                Assert.That(psm.Database, Is.EqualTo("18Protein_SoCe_Tr_detergents_trace_target_decoy"));
                Assert.That(psm.DatabaseVersion, Is.Null);
                Assert.That(psm.SearchEngines.Single(), Is.EqualTo(new CvParam("", "", "OMSSA", "2.1.4")));
                Assert.That(psm.SearchEngineScores[1], Is.EqualTo(8.620689655172414e-03));
                Assert.That(psm.RetentionTimeText, Is.EqualTo("1517.88525390625"));
                Assert.That(psm.RetentionTime, Is.EqualTo(1517.88525390625));
                Assert.That(psm.Charge, Is.EqualTo(3));
                Assert.That(psm.ExperimentalMassToCharge, Is.EqualTo(368.832366943358977));
                Assert.That(psm.CalculatedMassToCharge, Is.EqualTo(368.843772909670918));
                Assert.That(psm.SpectraRef, Is.EqualTo("ms_run[2]:spectrum=2311"));
                Assert.That(psm.MsRunIndex, Is.EqualTo(2));
                Assert.That(psm.NativeId, Is.EqualTo("spectrum=2311"));
                Assert.That(psm.MsRunLocation, Is.EqualTo("file:///home/user/OpenMS/share/OpenMS/examples/FRACTIONS/BSA2_F1.mzML"));
                Assert.That(psm.FileNameWithoutExtension, Is.EqualTo("BSA2_F1"));
                Assert.That(psm.OneBasedScanNumber, Is.EqualTo(2311));
                Assert.That(psm.IsDecoy, Is.False);
                Assert.That(psm.OptionalColumns["opt_global_cv_MS:1002217_decoy_peptide"], Is.EqualTo("0"));
                Assert.That(psm.OptionalColumns["opt_global_cv_MS:1000889_peptidoform_sequence"], Is.EqualTo("DGDIEAEISR"));
            });
        }

        [Test]
        public void TheLastOpenMsPsmIsReadFieldByField()
        {
            var psm = Read(OpenMs).Results[^1];

            Assert.Multiple(() =>
            {
                Assert.That(psm.PsmId, Is.EqualTo("114"));
                Assert.That(psm.BaseSequence, Is.EqualTo("LVVSTQTALA"));
                Assert.That(psm.Proteins.Single(), Is.EqualTo(new MzTabPsmProtein("P02769|ALBU_BOVIN", "K", "-", 598, 607)));
                Assert.That(psm.RetentionTime, Is.EqualTo(2341.02001953125));
                Assert.That(psm.Charge, Is.EqualTo(2));
                Assert.That(psm.ExperimentalMassToCharge, Is.EqualTo(501.794799804688012));
                Assert.That(psm.FileNameWithoutExtension, Is.EqualTo("BSA2_F2"));
                Assert.That(psm.OneBasedScanNumber, Is.EqualTo(3304));
            });
        }

        /// <summary>
        /// OpenMS joins a peptide's proteins in one cell, with pre/post/start/end joined alongside (rows 13 and 75).
        /// </summary>
        [Test]
        public void CommaJoinedProteinsAreSplitWithTheirFlankingResidues()
        {
            var file = Read(OpenMs);

            Assert.That(Psm(file, "12", "LSSPATLNSR").Proteins, Is.EqualTo(new[]
            {
                new MzTabPsmProtein("P06871|TRY1_CANFA", "K", "V", 113, 122),
                new MzTabPsmProtein("P00761|TRYP_PIG", "K", "V", 98, 107),
            }));

            var seven = Psm(file, "74", "LAADDFR");
            Assert.That(seven.Accessions, Has.Count.EqualTo(7));
            Assert.That(seven.Proteins[1], Is.EqualTo(new MzTabPsmProtein("Q14532|K1H2_HUMAN", "K", "A", 176, 182)));
            Assert.That(seven.Accession, Is.EqualTo(string.Join('|', seven.Accessions)));
        }

        [Test]
        public void UnimodModificationsResolveToAFullSequence()
        {
            var file = Read(OpenMs);

            var one = Psm(file, "2", "SHCIAEVEK");
            Assert.That(one.AllModsOneIsNterminus.Keys, Is.EqualTo(new[] { 4 }));
            Assert.That(one.AllModsOneIsNterminus[4].IdWithMotif, Is.EqualTo("Carbamidomethyl on C"));
            Assert.That(one.FullSequence, Is.EqualTo("SHC[Unimod:Carbamidomethyl on C]IAEVEK"));

            var two = Psm(file, "4", "CCTESLVNR");
            Assert.That(two.FullSequence, Is.EqualTo("C[Unimod:Carbamidomethyl on C]C[Unimod:Carbamidomethyl on C]TESLVNR"));
        }

        /// <summary>
        /// MaxQuant writes one row per protein. PSM_ID 1 is ms_run[16]:index=3508, whose opt scan column says 995: the
        /// index is a position in an .apl peak list, the opt column the raw scan.
        /// </summary>
        [Test]
        public void MaxQuantsRowsPerProteinAreOnePsmAndItsScanComesFromTheOptColumn()
        {
            var psm = Psm(Read(MaxQuant), "1", "AAAALAGGK");

            Assert.Multiple(() =>
            {
                Assert.That(psm.Accessions, Is.EqualTo(new[] { "Q3E792", "P0C0T4" }));
                Assert.That(psm.Accession, Is.EqualTo("Q3E792|P0C0T4"));
                Assert.That(psm.Unique, Is.False);
                Assert.That(psm.SpectraRef, Is.EqualTo("ms_run[16]:index=3508"));
                Assert.That(psm.NativeId, Is.EqualTo("index=3508"));
                Assert.That(psm.OneBasedScanNumber, Is.EqualTo(995));
                Assert.That(psm.FileNameWithoutExtension, Is.EqualTo("allspectra.cid.itms.iso_0"));
                Assert.That(psm.SearchEngines.Single(), Is.EqualTo(new CvParam("MS", "MS:1002338", "Andromeda:score", "1.6.3.1")));
                Assert.That(psm.SearchEngineScores[1], Is.EqualTo(85.2649238978181));
                // MaxQuant writes minutes, against the spec's seconds; the value is kept as written
                Assert.That(psm.RetentionTime, Is.EqualTo(8.10443210601807));
            });
        }

        [Test]
        public void MaxQuantsChemModsResolveByMassAtTheirSite()
        {
            var file = Read(MaxQuant);

            Assert.Multiple(() =>
            {
                Assert.That(Psm(file, "374", "ADFLLRPIK").FullSequence, Is.EqualTo("[Unimod:Acetyl on X]ADFLLRPIK"));
                Assert.That(Psm(file, "374", "ADFLLRPIK").AllModsOneIsNterminus[1].LocationRestriction, Does.Contain("N-terminal"));
                // ", " separates MaxQuant's entries
                Assert.That(Psm(file, "1279", "AHLMNPMVPGLAQGGK").FullSequence,
                    Is.EqualTo("AHLM[Unimod:Oxidation on M]NPM[Unimod:Oxidation on M]VPGLAQGGK"));
                // an Oxidation on the first residue stays on it, and the N-terminal acetyl takes the N-terminus
                Assert.That(Psm(file, "26723", "MAAIRDYK").FullSequence,
                    Is.EqualTo("[Unimod:Acetyl on X]M[Unimod:Oxidation on M]AAIRDYK"));
            });
        }

        /// <summary>
        /// MaxQuant leaves its fixed modifications out of the rows: fixed_mod[1] is CHEMMOD:57.0214637236 on C, and
        /// AACLVQNGIATWFPLAVTK's modifications cell is null.
        /// </summary>
        [Test]
        public void MaxQuantsDeclaredFixedModificationIsAddedAtEveryUnmodifiedSite()
        {
            var file = Read(MaxQuant);

            var noMods = Psm(file, "58", "AACLVQNGIATWFPLAVTK");
            Assert.That(noMods.Modifications, Is.Null);
            Assert.That(noMods.ModificationsResolved, Is.True);
            Assert.That(noMods.FullSequence, Is.EqualTo("AAC[Unimod:Carbamidomethyl on C]LVQNGIATWFPLAVTK"));

            var withOxidation = Psm(file, "299", "ACDDKIMYVDYK");
            Assert.That(withOxidation.FullSequence, Is.EqualTo("AC[Unimod:Carbamidomethyl on C]DDKIM[Unimod:Oxidation on M]YVDYK"));
        }

        [Test]
        public void DecoysAreReadFromTheDecoyColumn()
        {
            var file = Read(MaxQuant);

            Assert.That(file.Results.Where(p => p.IsDecoy).Select(p => p.PsmId), Is.EquivalentTo(new[] { "5223", "12687" }));
            Assert.That(Psm(file, "12687", "GSVEIENEALIR").Proteins.Single(), Is.EqualTo(new MzTabPsmProtein("REV__P43561", null, null, null, null)));
        }

        /// <summary>
        /// MassIVE's converter writes no decoy column, and marks decoy proteins with MS-GF+'s XXX_ prefix.
        /// </summary>
        [Test]
        public void WithoutADecoyColumnADecoyPrefixOnEveryAccessionMakesADecoy()
        {
            var file = Read(MassIve);

            Assert.That(file.Results.Where(p => p.IsDecoy).Select(p => p.PsmId), Is.EquivalentTo(new[] { "48", "87" }));
        }

        /// <summary>
        /// MassIVE writes pyro-Glu, an N-terminal-only modification, at the first residue ("1-UNIMOD:28"); mzLib keys it at
        /// the N-terminus.
        /// </summary>
        [Test]
        public void AnNTerminalOnlyModificationAtTheFirstResidueIsKeyedAtTheNTerminus()
        {
            var psm = Psm(Read(MassIve), "31", "QTISRPK");

            Assert.That(psm.AllModsOneIsNterminus.Keys, Is.EqualTo(new[] { 1 }));
            Assert.That(psm.FullSequence, Is.EqualTo("[Unimod:Gln->pyro-Glu on Q]QTISRPK"));
        }

        /// <summary>
        /// PSM_ID 790 is written on two adjacent rows for two spectra and two peptides; PSM_ID 1082 for two spectra of one
        /// peptide. Each spectrum is its own PSM.
        /// </summary>
        [Test]
        public void APsmIdReusedForAnotherSpectrumIsAnotherPsm()
        {
            var file = Read(MassIve);

            Assert.That(file.Results.Where(p => p.PsmId == "790").Select(p => p.NativeId), Is.EqualTo(new[] { "index=39184", "index=39185" }));
            Assert.That(file.Results.Where(p => p.PsmId == "1082").Select(p => p.OneBasedScanNumber), Is.EqualTo(new[] { 46121, 46131 }));
        }

        /// <summary>
        /// HCYYQIRAK carries "2-UNIMOD:4,9-CHEMMOD:-11": a Carbamidomethyl that resolves and a nominal mass that names
        /// nothing. The PSM is kept, flagged, with the resolved modification and no FullSequence.
        /// </summary>
        [Test]
        public void APsmWithAnUnresolvedModificationIsKeptFlaggedWithNoFullSequence()
        {
            var file = Read(MassIve);
            var psm = Psm(file, "48", "HCYYQIRAK");

            Assert.Multiple(() =>
            {
                Assert.That(file.SkippedPsms, Is.Empty);
                Assert.That(psm.ModificationsResolved, Is.False);
                Assert.That(psm.FullSequence, Is.Empty);
                Assert.That(psm.AllModsOneIsNterminus.Keys, Is.EqualTo(new[] { 3 }));
                Assert.That(psm.AllModsOneIsNterminus[3].IdWithMotif, Is.EqualTo("Carbamidomethyl on C"));
                Assert.That(psm.UnresolvedModifications.Single(), Is.EqualTo(new MzTabUnresolvedModification("9-CHEMMOD:-11",
                    "CHEMMOD:-11 at location 9 gives the mass with fewer than two decimals, too coarse to name one modification")));
                Assert.That(psm.FileNameWithoutExtension, Is.EqualTo("specs_ms"));
                Assert.That(psm.OneBasedScanNumber, Is.EqualTo(4197));
            });
        }

        /// <summary>
        /// Proteome Discoverer 2.5 writes the N-terminal TMT at location 1 on every row (here "11-UNIMOD:737,1-UNIMOD:737"),
        /// and twice at 1 on a peptide that starts with K.
        /// </summary>
        [Test]
        public void ProteomeDiscoverersNTerminalLabelAtLocationOneGoesToTheNTerminus()
        {
            var file = Read(ProteomeDiscoverer);
            var psm = Psm(file, "1", "AFITNIPFDVK");

            Assert.Multiple(() =>
            {
                Assert.That(psm.AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { 1, 12 }));
                Assert.That(psm.AllModsOneIsNterminus[1].IdWithMotif, Is.EqualTo("TMT11 on X"));
                Assert.That(psm.AllModsOneIsNterminus[1].LocationRestriction, Is.EqualTo("Peptide N-terminal."));
                Assert.That(psm.AllModsOneIsNterminus[12].IdWithMotif, Is.EqualTo("TMT11 on K"));
                // the K-start peptide, "1-UNIMOD:737,10-UNIMOD:737,1-UNIMOD:737": the lysine and the N-terminus
                Assert.That(Psm(file, "359", "KDAEAWFNEK").AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { 1, 2, 11 }));
                Assert.That(psm.NativeId, Is.EqualTo("scan=23239"));
                Assert.That(psm.OneBasedScanNumber, Is.EqualTo(23239));
                Assert.That(psm.RetentionTime, Is.Null);
                Assert.That(psm.Unique, Is.Null);
                Assert.That(file.Results.All(p => p.ModificationsResolved), Is.True);
            });
        }

        /// <summary>
        /// PSM_ID 4 is rows 4-6 (three proteins), and the same three rows again at 20-22 and 27-29.
        /// </summary>
        [Test]
        public void ExactCopiesOfAPsmFarApartAddNothingToResultsButStreamAgain()
        {
            var file = Read(ProteomeDiscoverer);

            var merged = file.Results.Where(p => p.PsmId == "4").ToList();
            Assert.That(merged, Has.Count.EqualTo(1));
            Assert.That(merged[0].Accessions, Is.EqualTo(new[] { "P02533", "Q04695", "P08779" }));

            var streamed = file.EnumeratePsms().Where(p => p.PsmId == "4").ToList();
            Assert.That(streamed, Has.Count.EqualTo(3));
            Assert.That(streamed.Select(p => p.Accessions.Count), Is.All.EqualTo(3));
        }

        /// <summary>
        /// PRIDE's converter writes protein by protein: PSM_ID 759 is row 27 under Q5FWB7_MOUSE and row 126 under
        /// A6ZI47_MOUSE. PSM_ID 466 holds two candidate peptides of one spectrum, I and L at one position.
        /// </summary>
        [Test]
        public void PrideConverterPsmsMergeAcrossProteinsAndKeepEachCandidatePeptide()
        {
            var file = Read(PrideConverter);

            Assert.That(Psm(file, "759", "QLLLTADDRVNPCIGGVILFHETLYQK").Accessions, Is.EqualTo(new[] { "Q5FWB7_MOUSE", "A6ZI47_MOUSE" }));
            Assert.That(file.EnumeratePsms().Count(p => p.PsmId == "759"), Is.EqualTo(2));
            Assert.That(file.Results.Where(p => p.PsmId == "466").Select(p => p.BaseSequence),
                Is.EquivalentTo(new[] { "DLYANTVLSGGTTMYPGIADR", "DLYANTVLSGGTTMYPGLADR" }));
        }

        [Test]
        public void PsiModAccessionsAreUnresolvedForNow()
        {
            var psm = Read(PrideConverter).Results[0];

            Assert.That(psm.Modifications, Is.EqualTo("11-MOD:00397"));
            Assert.That(psm.ModificationsResolved, Is.False);
            Assert.That(psm.FullSequence, Is.Empty);
            Assert.That(psm.UnresolvedModifications.Single().Reason, Is.EqualTo("MOD:00397 is a PSI-MOD accession, which this reader does not resolve yet"));
            Assert.That(psm.OneBasedScanNumber, Is.EqualTo(391));
            Assert.That(psm.FileNameWithoutExtension, Is.EqualTo("PRIDE_Example"));
        }

        /// <summary>
        /// Mascot through jmzTab writes "5-UNIMOD:35,5-[PSI-MS, MS:1001524, fragment neutral loss, 0]", and a bare number
        /// after ms_run[1]:, whose meaning the file does not say.
        /// </summary>
        [Test]
        public void ANeutralLossIsIgnoredAndABareNumberGivesNoScan()
        {
            var file = Read(Mascot);
            var psm = Psm(file, "2", "STAAMSTYTGIFTDQVLSVLK");

            Assert.That(psm.ModificationsResolved, Is.True);
            Assert.That(psm.FullSequence, Is.EqualTo("STAAM[Unimod:Oxidation on M]STYTGIFTDQVLSVLK"));
            Assert.That(psm.NativeId, Is.EqualTo("244"));
            Assert.That(psm.OneBasedScanNumber, Is.EqualTo(-1));
            Assert.That(psm.FileNameWithoutExtension, Is.EqualTo("small"));
        }

        /// <summary>
        /// What the reader is for: comparing a deposited search with one run on the same spectra. The 24 PSMs of
        /// ms_run[1] each name a spectrum of BSA1_F1.mzML (the OpenMS input, whole) by nativeID, and it is theirs: MS2,
        /// precursor m/z within 0.01, retention time within 0.5 s.
        /// CAREFUL: the files' ids are not positions. The mzML's first spectrum is spectrum=1011, and mzLib numbers
        /// spectra by position, so spectrum=2458 is mzLib's scan 303. The lookup has to go through NativeId.
        /// </summary>
        [Test]
        public void EveryPsmOfARunFindsItsSpectrumInThatRunsMzmlByNativeId()
        {
            var run = Read(OpenMs).Results.Where(p => p.FileNameWithoutExtension == "BSA1_F1").ToList();
            var mzml = MsDataFileReader.GetDataFile(Fixture("BSA1_F1.mzML"));
            mzml.LoadAllStaticData();
            var byNativeId = mzml.GetAllScansList().ToDictionary(s => s.NativeId);

            Assert.That(run, Has.Count.EqualTo(24));
            Assert.Multiple(() =>
            {
                foreach (var psm in run)
                {
                    Assert.That(byNativeId.TryGetValue(psm.NativeId, out var scan), Is.True, psm.NativeId);
                    Assert.That(scan!.MsnOrder, Is.EqualTo(2), psm.NativeId);
                    Assert.That(scan.SelectedIonMZ!.Value, Is.EqualTo(psm.ExperimentalMassToCharge!.Value).Within(0.01), psm.NativeId);
                    Assert.That(scan.RetentionTime * 60, Is.EqualTo(psm.RetentionTime!.Value).Within(0.5), psm.NativeId);
                }

                var example = run.Single(p => p.PsmId == "2");
                Assert.That(example.NativeId, Is.EqualTo("spectrum=2458"));
                Assert.That(example.OneBasedScanNumber, Is.EqualTo(2458));
                Assert.That(byNativeId["spectrum=2458"].OneBasedScanNumber, Is.EqualTo(303));
            });
        }

        [Test]
        public void TheMetadataIsReadWithoutThePsms()
        {
            var metadata = Read(MaxQuant).Metadata;

            Assert.Multiple(() =>
            {
                Assert.That(metadata.Version, Is.EqualTo("1.0.0"));
                Assert.That(metadata.Mode, Is.EqualTo("Complete"));
                Assert.That(metadata.Type, Is.EqualTo("Identification"));
                Assert.That(metadata.IsMzTabM, Is.False);
                Assert.That(metadata.Software.Single(), Is.EqualTo(new CvParam("MS", "MS:1001583", "MaxQuant", "1.6.3.1")));
                Assert.That(metadata.PsmSearchEngineScores[1].Name, Is.EqualTo("Andromeda:score"));
                Assert.That(metadata.FixedModifications.Single(), Is.EqualTo(new MzTabModificationDeclaration(1,
                    new CvParam("", "", "CHEMMOD:57.0214637236", ""), "C", "Anywhere")));
                Assert.That(metadata.VariableModifications.Select(m => m.Param!.Name),
                    Is.EqualTo(new[] { "CHEMMOD:15.9949146221", "CHEMMOD:42.0105646863" }));
                Assert.That(metadata.MsRuns[16].Location, Does.EndWith("/allspectra.cid.itms.iso_0.apl"));
                Assert.That(metadata.MsRuns[16].Format, Is.Not.Null);
            });
        }

        [Test]
        public void TheGzippedFileGivesTheSamePsms()
        {
            string gz = Gzipped(OpenMs);

            Assert.That(gz.ParseFileType(), Is.EqualTo(SupportedFileType.MzTabGz));
            var plain = Read(OpenMs).Results;
            var compressed = new MzTabFile(gz);
            Assert.That(compressed.FileType, Is.EqualTo(SupportedFileType.MzTabGz));
            Assert.That(compressed.Results.Select(p => (p.PsmId, p.FullSequence, p.Accession, p.OneBasedScanNumber)),
                Is.EqualTo(plain.Select(p => (p.PsmId, p.FullSequence, p.Accession, p.OneBasedScanNumber))));
            Assert.That(FileReader.ReadResultFile(gz), Is.TypeOf<MzTabFile>());
        }

        /// <summary>
        /// Streaming returns PSMs as it reaches them: a gzip whose compressed data is corrupt 90% of the way in gives its
        /// first PSMs, and fails only when enumeration reaches the damage.
        /// </summary>
        [Test]
        public void StreamingReturnsThePsmsBeforeAProblemLaterInTheFile()
        {
            byte[] bytes = File.ReadAllBytes(Gzipped(OpenMs));
            int at = bytes.Length * 9 / 10;
            for (int i = at; i < at + 64; i++)
            {
                bytes[i] = 0xFF;
            }

            string path = Path.Combine(_outputDirectory, "corrupt_late.mztab.gz");
            File.WriteAllBytes(path, bytes);
            var file = new MzTabFile(path);

            Assert.That(file.EnumeratePsms().First().PsmId, Is.EqualTo("0"));
            Assert.Throws<MzLibException>(() => file.EnumeratePsms().ToList());
        }

        /// <summary>
        /// CAREFUL, and why the class says so: .NET's GZipStream does not report a truncated gzip. A copy cut off at 90%,
        /// as an interrupted download leaves it, reads as a shorter file without any error.
        /// </summary>
        [Test]
        public void ATruncatedGzipReadsAsAShorterFileWithoutError()
        {
            byte[] whole = File.ReadAllBytes(Gzipped(OpenMs));
            string path = Path.Combine(_outputDirectory, "truncated.mztab.gz");
            File.WriteAllBytes(path, whole[..(whole.Length * 9 / 10)]);

            Assert.That(new MzTabFile(path).Results.Count, Is.LessThan(115));
        }

        /// <summary>
        /// The file is deleted after enumerating, so Results can only be non-empty if enumerating had filled it.
        /// </summary>
        [Test]
        public void EnumeratingDoesNotFillResults()
        {
            string path = Path.Combine(_outputDirectory, "enumerate_only.mztab");
            File.Copy(Fixture(OpenMs), path, true);
            var file = new MzTabFile(path);

            Assert.That(file.EnumeratePsms().Count(), Is.EqualTo(115));
            File.Delete(path);
            Assert.That(file.Results, Is.Empty);
        }

        /// <summary>
        /// A PSM row before the PSH header cannot be read: there are no columns to read it by. SYNTHETIC.
        /// </summary>
        [Test]
        public void APsmRowBeforeTheHeaderIsAnMzLibException()
        {
            string path = Path.Combine(_outputDirectory, "row_before_header.mztab");
            var lines = File.ReadAllLines(Fixture(OpenMs)).ToList();
            int psh = lines.FindIndex(l => l.StartsWith("PSH\t"));
            lines.Insert(psh, lines[psh + 1]);
            File.WriteAllLines(path, lines);

            var e = Assert.Throws<MzLibException>(() => new MzTabFile(path).EnumeratePsms().First());
            Assert.That(e!.Message, Does.Contain($"line {psh + 1}: a PSM row comes before the PSH header"));
        }

        /// <summary>
        /// A file with no PSMs leaves Results empty, which would make ResultFile reload it on every access.
        /// </summary>
        [Test]
        public void AFileWithNoPsmsIsReadOnce()
        {
            string path = Path.Combine(_outputDirectory, "metadata_only.mztab");
            File.WriteAllLines(path, File.ReadAllLines(Fixture(OpenMs)).TakeWhile(l => !l.StartsWith("PSH\t")));
            var file = new MzTabFile(path);

            Assert.That(file.Results, Is.Empty);
            File.Delete(path);
            Assert.That(file.Results, Is.Empty, "a second access must not reread the file");
            Assert.DoesNotThrow(() => file.LoadResults());
        }

        [Test]
        public void ARowWithNoSequenceIsSkippedWithItsReason()
        {
            string path = Path.Combine(_outputDirectory, "no_sequence.mztab");
            var lines = File.ReadAllLines(Fixture(OpenMs)).ToList();
            int first = lines.FindIndex(l => l.StartsWith("PSM\t"));
            lines[first] = "PSM\tnull" + lines[first][lines[first].IndexOf('\t', 4)..]; // SYNTHETIC: sequence nulled
            File.WriteAllLines(path, lines);

            var file = new MzTabFile(path);

            Assert.That(file.Results, Has.Count.EqualTo(114));
            Assert.That(file.SkippedPsms.Single(), Is.EqualTo(new MzTabSkippedPsm("0", "ms_run[2]:spectrum=2311", null, null, "the row has no sequence")));
        }

        [Test]
        public void MzTabMIsRefusedByDetectionAndByTheReader()
        {
            var e = Assert.Throws<MzLibException>(() => Fixture(MzTabM).ParseFileType());
            Assert.That(e!.Message, Is.EqualTo("MzTab file type not supported: mzTab-M (2.0.0-M) is a different format"));

            var file = new MzTabFile(Fixture(MzTabM));
            Assert.Throws<MzLibException>(() => _ = file.Results);
            Assert.That(file.Metadata.IsMzTabM, Is.True);
        }

        [Test]
        public void AFileWithNoMzTabVersionIsRefused()
        {
            string path = Path.Combine(_outputDirectory, "pre_1_0.mztab");
            File.WriteAllLines(path, File.ReadAllLines(Fixture(OpenMs)).Where(l => !l.StartsWith("MTD\tmzTab-version"))); // SYNTHETIC

            var e = Assert.Throws<MzLibException>(() => path.ParseFileType());
            Assert.That(e!.Message, Is.EqualTo("MzTab file type not supported: the metadata has no mzTab-version, as in files older than mzTab 1.0"));
            Assert.Throws<MzLibException>(() => _ = new MzTabFile(path).Results);
        }

        [Test]
        public void ACorruptGzipIsAnMzLibException()
        {
            string path = Path.Combine(_outputDirectory, "corrupt.mztab.gz");
            File.WriteAllBytes(path, [0x1f, 0x8b, 0x08, 0x00, 0xde, 0xad, 0xbe, 0xef]);

            Assert.Throws<MzLibException>(() => path.ParseFileType());
            Assert.Throws<MzLibException>(() => _ = new MzTabFile(path).Results);
        }

        [Test]
        public void AMissingFileIsFileNotFound()
        {
            var file = new MzTabFile(Path.Combine(_outputDirectory, "missing.mztab"));

            Assert.That(file.Results, Is.Empty);
            Assert.Throws<FileNotFoundException>(() => file.EnumeratePsms().First());
            Assert.Throws<FileNotFoundException>(() => _ = file.Metadata);
        }

        [Test]
        public void WritingIsNotSupported()
        {
            Assert.Throws<NotSupportedException>(() => Read(OpenMs).WriteResults(Path.Combine(_outputDirectory, "out.mztab")));
        }

        [Test]
        public void TheFactoryBuildsAnMzTabFile()
        {
            var file = FileReader.ReadResultFile(Fixture(MaxQuant));

            Assert.That(file, Is.TypeOf<MzTabFile>());
            Assert.That(file.Software, Is.EqualTo(Software.MzTab));
            Assert.That(((MzTabFile)file).Results, Has.Count.EqualTo(49));
        }
    }
}
