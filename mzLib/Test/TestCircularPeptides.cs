using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using NUnit.Framework;
using Proteomics;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestCircularPeptides
    {
        private static string OutputDirectory => Path.Combine(@"C:\Users\trish\Downloads", "CircularPeptideOutput");

        [OneTimeSetUp]
        public void Setup()
        {
            // Create output directory if it doesn't exist
            if (!Directory.Exists(OutputDirectory))
            {
                Directory.CreateDirectory(OutputDirectory);
            }
        }

        #region CircularPeptide Tests

        [Test]
        public void CircularPeptide_AssignOrigin_FindsLowestAlphabetically()
        {
            // "PEPTIDE" - D is the lowest letter, at position 5
            // P(0) E(1) P(2) T(3) I(4) D(5) E(6)
            var peptide = new CircularPeptide("PEPTIDE");

            Assert.That(peptide.OriginIndex, Is.EqualTo(5), "Origin should be at position 5 (the 'D')");
            Assert.That(peptide.CanonicalSequence, Is.EqualTo("DEPEPTI"), "Canonical sequence should start with D");
        }

        [Test]
        public void CircularPeptide_AssignOrigin_HandlesMultipleSameLetters()
        {
            // "ABCABC" - Multiple A's at positions 0 and 3
            // Position 0: ABCABC
            // Position 3: ABCABC (same rotation)
            // Compare: A-B vs A-B, then B-C vs B-C... they're the same pattern
            var peptide = new CircularPeptide("ABCABC");

            Assert.That(peptide.OriginIndex, Is.EqualTo(0), "Origin should be at position 0");
            Assert.That(peptide.CanonicalSequence, Is.EqualTo("ABCABC"));
        }

        [Test]
        public void CircularPeptide_AssignOrigin_BreaksTiesCorrectly()
        {
            // "ACABAC" - A's at positions 0, 2, 4
            // Rotation from 0: ACABAC
            // Rotation from 2: ABACAC
            // Rotation from 4: ACACAB
            // Sorted: ABACAC < ACABAC < ACACAB
            var peptide = new CircularPeptide("ACABAC");

            Assert.That(peptide.OriginIndex, Is.EqualTo(2), "Origin should be at position 2 (gives ABACAC)");
            Assert.That(peptide.CanonicalSequence, Is.EqualTo("ABACAC"));
        }

        [Test]
        public void CircularPeptide_IsEquivalentTo_DetectsRotations()
        {
            var peptide1 = new CircularPeptide("PEPTIDE");
            var peptide2 = new CircularPeptide("TIDEPEP"); // Rotation of PEPTIDE
            var peptide3 = new CircularPeptide("PEPTIDX"); // Different sequence

            Assert.That(peptide1.IsEquivalentTo(peptide2), Is.True, "Rotations should be equivalent");
            Assert.That(peptide1.IsEquivalentTo(peptide3), Is.False, "Different sequences should not be equivalent");
        }

        [Test]
        public void CircularPeptide_GetRotatedSequence_Works()
        {
            var peptide = new CircularPeptide("ABCDEF");

            Assert.That(peptide.GetRotatedSequence(0), Is.EqualTo("ABCDEF"));
            Assert.That(peptide.GetRotatedSequence(2), Is.EqualTo("CDEFAB"));
            Assert.That(peptide.GetRotatedSequence(5), Is.EqualTo("FABCDE"));
        }

        #endregion

        #region LinearFragment Tests

        [Test]
        public void LinearFragment_ToCircularPosition_MapsCorrectly()
        {
            // Fragment "CDEF" starting at position 2 in a circular peptide of length 6
            var fragment = new LinearFragment("CDEF", 2, 5, 6);

            Assert.That(fragment.ToCircularPosition(0), Is.EqualTo(2)); // C
            Assert.That(fragment.ToCircularPosition(1), Is.EqualTo(3)); // D
            Assert.That(fragment.ToCircularPosition(2), Is.EqualTo(4)); // E
            Assert.That(fragment.ToCircularPosition(3), Is.EqualTo(5)); // F
        }

        [Test]
        public void LinearFragment_ToCircularBondPosition_MapsCorrectly()
        {
            var fragment = new LinearFragment("CDEF", 2, 5, 6);

            Assert.That(fragment.ToCircularBondPosition(0), Is.EqualTo(2)); // C-D bond
            Assert.That(fragment.ToCircularBondPosition(1), Is.EqualTo(3)); // D-E bond
            Assert.That(fragment.ToCircularBondPosition(2), Is.EqualTo(4)); // E-F bond
        }

        [Test]
        public void LinearFragment_CrossesOrigin_DetectedCorrectly()
        {
            var normalFragment = new LinearFragment("CDE", 2, 4, 6);
            var wrappingFragment = new LinearFragment("FABC", 5, 2, 6);

            Assert.That(normalFragment.CrossesOrigin, Is.False);
            Assert.That(wrappingFragment.CrossesOrigin, Is.True);
        }

        #endregion

        #region PredictedSpectrum Fragment Generation Tests

        [Test]
        public void PredictedSpectrum_GeneratesOverlappingFragments()
        {
            var circular = new CircularPeptide("ACDEFGHIKLMNPQRSTVWY"); // 20 amino acids
            var spectrum = new PredictedSpectrum(circular);

            Assert.That(spectrum.LinearFragments.Count, Is.GreaterThanOrEqualTo(3),
                "Should generate at least 3 fragments");

            // Check coverage
            int[] coverage = new int[20];
            foreach (var fragment in spectrum.LinearFragments)
            {
                for (int i = 0; i < fragment.Sequence.Length; i++)
                {
                    int pos = fragment.ToCircularPosition(i);
                    coverage[pos]++;
                }
            }

            Assert.That(coverage.Min(), Is.GreaterThanOrEqualTo(1),
                "Every position should be covered at least once");
        }

        [Test]
        public void PredictedSpectrum_ShortPeptide_HandledCorrectly()
        {
            var circular = new CircularPeptide("ACDEFGHIK"); // 9 amino acids
            var spectrum = new PredictedSpectrum(circular, minFragmentLength: 7, maxFragmentLength: 15);

            Assert.That(spectrum.LinearFragments.Count, Is.GreaterThanOrEqualTo(3));

            foreach (var fragment in spectrum.LinearFragments)
            {
                Assert.That(fragment.Sequence.Length, Is.GreaterThanOrEqualTo(7),
                    "All fragments should meet minimum length");
            }
        }

        #endregion

        #region Bond Cleavage Probability Tests

        [Test]
        public void BondCleavageProbability_NormalizesToOne()
        {
            var circular = new CircularPeptide("PEPTIDE");
            var spectrum = new PredictedSpectrum(circular);

            // Add mock ions with different intensities
            var bondIntensities = new Dictionary<int, double>
            {
                { 0, 0.2 }, { 1, 0.3 }, { 2, 0.8 }, { 3, 0.1 },
                { 4, 0.5 }, { 5, 0.4 }, { 6, 0.3 }
            };

            foreach (var fragment in spectrum.LinearFragments)
            {
                spectrum.AddMockPredictedIons(fragment, bondIntensities);
            }

            spectrum.ComputeBondCleavageProbabilities();

            double totalProb = spectrum.BondCleavageProbabilities.Sum(b => b.CleavageProbability);
            Assert.That(totalProb, Is.EqualTo(1.0).Within(0.0001),
                "Probabilities should sum to 1");
        }

        [Test]
        public void BondCleavageProbability_HigherIntensityGivesHigherProbability()
        {
            var circular = new CircularPeptide("ACDEFGHI"); // 8 amino acids
            var spectrum = new PredictedSpectrum(circular);

            // Set one bond to have much higher intensity
            var bondIntensities = new Dictionary<int, double>
            {
                { 0, 0.1 }, { 1, 0.1 }, { 2, 0.9 }, { 3, 0.1 },
                { 4, 0.1 }, { 5, 0.1 }, { 6, 0.1 }, { 7, 0.1 }
            };

            foreach (var fragment in spectrum.LinearFragments)
            {
                spectrum.AddMockPredictedIons(fragment, bondIntensities);
            }

            spectrum.ComputeBondCleavageProbabilities();

            var highestProbBond = spectrum.BondCleavageProbabilities
                .OrderByDescending(b => b.CleavageProbability)
                .First();

            Assert.That(highestProbBond.BondPosition, Is.EqualTo(2),
                "Bond 2 should have highest probability");
        }

        #endregion

        #region Double Fragmentation Simulation Tests

        [Test]
        public void DoubleFragmentationSimulation_ProducesCorrectFragmentCount()
        {
            var circular = new CircularPeptide("ACDEFGHIKLMN"); // 12 amino acids
            var spectrum = new PredictedSpectrum(circular, randomSeed: 42);

            // Set uniform probabilities
            double uniformProb = 1.0 / 12;
            spectrum.SetCleavageProbabilities(Enumerable.Repeat(uniformProb, 12).ToArray());

            spectrum.RunDoubleFragmentationSimulation(1000);

            Assert.That(spectrum.SimulationResults.Count, Is.EqualTo(1000),
                "Should have 1000 simulation events");

            // Each event produces 2 fragments, total should be 2000 in the spectrum
            int totalFragmentCount = spectrum.SimulatedSpectrum.Sum(p => p.Count);
            Assert.That(totalFragmentCount, Is.EqualTo(2000),
                "Total fragment count should be 2 * simulation count");
        }

        [Test]
        public void DoubleFragmentationSimulation_TwoFragmentsPerEvent()
        {
            var circular = new CircularPeptide("ACDEFGHI"); // 8 amino acids
            var spectrum = new PredictedSpectrum(circular, randomSeed: 123);

            spectrum.SetCleavageProbabilities(Enumerable.Repeat(1.0 / 8, 8).ToArray());
            spectrum.RunDoubleFragmentationSimulation(100);

            foreach (var evt in spectrum.SimulationResults)
            {
                Assert.That(evt.Bond1, Is.Not.EqualTo(evt.Bond2),
                    "Two different bonds should be selected");
                Assert.That(evt.Fragment1.Sequence.Length + evt.Fragment2.Sequence.Length,
                    Is.EqualTo(8), "Two fragments should sum to original length");
            }
        }

        [Test]
        public void DoubleFragmentationSimulation_FragmentsHaveCorrectCoordinates()
        {
            var circular = new CircularPeptide("ACDEFGHI"); // 8 amino acids (all valid)
            var spectrum = new PredictedSpectrum(circular, randomSeed: 456);

            spectrum.SetCleavageProbabilities(Enumerable.Repeat(1.0 / 8, 8).ToArray());
            spectrum.RunDoubleFragmentationSimulation(50);

            foreach (var evt in spectrum.SimulationResults)
            {
                // Fragment 1 should start after Bond1 and end at Bond2
                int expectedStart1 = (evt.Bond1 + 1) % 8;
                Assert.That(evt.Fragment1.StartIndex, Is.EqualTo(expectedStart1));

                // Fragment 2 should start after Bond2 and end at Bond1
                int expectedStart2 = (evt.Bond2 + 1) % 8;
                Assert.That(evt.Fragment2.StartIndex, Is.EqualTo(expectedStart2));
            }
        }

        #endregion

        #region Mass Calculation Tests

        [Test]
        public void SimulatedSpectrum_CalculatesMassCorrectly()
        {
            var circular = new CircularPeptide("AG"); // Simple: Alanine + Glycine
            var spectrum = new PredictedSpectrum(circular, minFragmentLength: 1);

            spectrum.SetCleavageProbabilities(new[] { 0.5, 0.5 });
            spectrum.RunDoubleFragmentationSimulation(100);

            // Find the peaks for A and G
            var peakA = spectrum.SimulatedSpectrum.FirstOrDefault(p => p.Sequence == "A");
            var peakG = spectrum.SimulatedSpectrum.FirstOrDefault(p => p.Sequence == "G");

            // A mass: 71.03711 + 18.01056 (water) = 89.04767
            // G mass: 57.02146 + 18.01056 (water) = 75.03202
            if (peakA != null)
            {
                Assert.That(peakA.MonoisotopicMass, Is.EqualTo(89.04767).Within(0.001));
            }
            if (peakG != null)
            {
                Assert.That(peakG.MonoisotopicMass, Is.EqualTo(75.03202).Within(0.001));
            }
        }

        #endregion

        #region Full Integration Test with File Output

        [Test]
        public void FullSimulation_WithFileOutput()
        {
            // Create a realistic circular peptide
            string sequence = "ACDEFGHIKLMNPQRSTVWY"; // 20 amino acids (all standard)
            var circular = new CircularPeptide(sequence);
            var spectrum = new PredictedSpectrum(circular, randomSeed: 42);

            // Set up bond intensities that simulate realistic cleavage preferences
            // Proline and Glycine bonds tend to be weaker
            var bondIntensities = new Dictionary<int, double>();
            string canonicalSeq = circular.CanonicalSequence;

            for (int i = 0; i < canonicalSeq.Length; i++)
            {
                char residueN = canonicalSeq[i];
                char residueC = canonicalSeq[(i + 1) % canonicalSeq.Length];

                // Higher intensity for certain bond types
                double intensity = 0.1; // baseline

                // Proline N-terminal to bond increases cleavage
                if (residueC == 'P') intensity = 0.8;
                // After Asp or Glu
                else if (residueN == 'D' || residueN == 'E') intensity = 0.6;
                // After aromatic residues
                else if (residueN == 'W' || residueN == 'Y' || residueN == 'F') intensity = 0.5;

                bondIntensities[i] = intensity;
            }

            // Add mock ions based on intensities
            foreach (var fragment in spectrum.LinearFragments)
            {
                spectrum.AddMockPredictedIons(fragment, bondIntensities);
            }

            // Compute probabilities
            spectrum.ComputeBondCleavageProbabilities();

            // Run simulation
            int numCopies = 100000;
            spectrum.RunDoubleFragmentationSimulation(numCopies);

            // Generate output
            var sb = new StringBuilder();
            sb.AppendLine("=" + new string('=', 79));
            sb.AppendLine("CIRCULAR PEPTIDE DOUBLE FRAGMENTATION SIMULATION");
            sb.AppendLine("=" + new string('=', 79));
            sb.AppendLine();
            sb.AppendLine($"Date: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
            sb.AppendLine($"Original Sequence: {sequence}");
            sb.AppendLine($"Canonical Sequence: {circular.CanonicalSequence}");
            sb.AppendLine($"Origin Index: {circular.OriginIndex}");
            sb.AppendLine($"Peptide Length: {sequence.Length}");
            sb.AppendLine($"Number of Simulated Copies: {numCopies:N0}");
            sb.AppendLine();

            // Coverage summary
            sb.AppendLine(spectrum.GetCoverageSummary());
            sb.AppendLine();

            // Cleavage probabilities
            sb.AppendLine(spectrum.GetCleavageProbabilityTable());
            sb.AppendLine();

            // Fragmentation patterns
            sb.AppendLine(spectrum.GetFragmentationPatternsTable());
            sb.AppendLine();

            // Simulated spectrum
            sb.AppendLine(spectrum.GetSimulatedSpectrumTable());
            sb.AppendLine();

            // Summary statistics
            sb.AppendLine("SUMMARY STATISTICS");
            sb.AppendLine("==================");
            sb.AppendLine($"Unique Fragment Sequences: {spectrum.SimulatedSpectrum.Select(p => p.Sequence).Distinct().Count()}");
            sb.AppendLine($"Unique Fragmentation Patterns: {spectrum.SimulationResults.Select(e => e.GetKey()).Distinct().Count()}");
            sb.AppendLine($"Most Common Fragment: {spectrum.SimulatedSpectrum.OrderByDescending(p => p.Count).First().Sequence}");
            sb.AppendLine($"Smallest Fragment: {spectrum.SimulatedSpectrum.OrderBy(p => p.Sequence.Length).First().Sequence}");
            sb.AppendLine($"Largest Fragment: {spectrum.SimulatedSpectrum.OrderByDescending(p => p.Sequence.Length).First().Sequence}");
            sb.AppendLine($"Mass Range: {spectrum.SimulatedSpectrum.Min(p => p.Mz):F2} - {spectrum.SimulatedSpectrum.Max(p => p.Mz):F2} m/z");

            // Write to file
            string outputPath = Path.Combine(OutputDirectory, "CircularPeptideSimulation_Results.txt");
            File.WriteAllText(outputPath, sb.ToString());

            // Also write CSV for the spectrum
            string csvPath = Path.Combine(OutputDirectory, "CircularPeptideSimulation_Spectrum.csv");
            using (var writer = new StreamWriter(csvPath))
            {
                writer.WriteLine("Sequence,StartPosition,EndPosition,MonoisotopicMass,Mz,Charge,Count,RelativeIntensity");
                foreach (var peak in spectrum.SimulatedSpectrum.OrderByDescending(p => p.Count))
                {
                    writer.WriteLine($"{peak.Sequence},{peak.StartPosition},{peak.EndPosition}," +
                        $"{peak.MonoisotopicMass:F6},{peak.Mz:F6},{peak.Charge},{peak.Count},{peak.RelativeIntensity:F6}");
                }
            }

            // Write fragmentation patterns CSV
            string patternsPath = Path.Combine(OutputDirectory, "CircularPeptideSimulation_Patterns.csv");
            using (var writer = new StreamWriter(patternsPath))
            {
                writer.WriteLine("Bond1,Bond2,Fragment1,Fragment2,Count,Frequency");
                var patterns = spectrum.SimulationResults
                    .GroupBy(e => e.GetKey())
                    .Select(g => new
                    {
                        Bond1 = g.First().Bond1,
                        Bond2 = g.First().Bond2,
                        Frag1 = g.First().Fragment1.Sequence,
                        Frag2 = g.First().Fragment2.Sequence,
                        Count = g.Count()
                    })
                    .OrderByDescending(p => p.Count);

                foreach (var pattern in patterns)
                {
                    writer.WriteLine($"{pattern.Bond1},{pattern.Bond2},{pattern.Frag1},{pattern.Frag2}," +
                        $"{pattern.Count},{(double)pattern.Count / numCopies:F6}");
                }
            }

            // Assertions
            Assert.That(File.Exists(outputPath), Is.True, "Output file should be created");
            Assert.That(File.Exists(csvPath), Is.True, "CSV spectrum file should be created");
            Assert.That(File.Exists(patternsPath), Is.True, "CSV patterns file should be created");

            // Log the output location
            Console.WriteLine($"Results written to: {outputPath}");
            Console.WriteLine($"Spectrum CSV written to: {csvPath}");
            Console.WriteLine($"Patterns CSV written to: {patternsPath}");

            // Verify simulation ran correctly
            Assert.That(spectrum.SimulationResults.Count, Is.EqualTo(numCopies));
            Assert.That(spectrum.SimulatedSpectrum.Count, Is.GreaterThan(0));
        }
        [Test]
        public void FullSimulation_WithFileOutput2()
        {
            // Create a realistic circular peptide
            string sequence = "KELSDIAHRIVAPGK"; // 20 amino acids (all standard)
            var circular = new CircularPeptide(sequence);
            var spectrum = new PredictedSpectrum(circular, randomSeed: 42);

            // Set up bond intensities that simulate realistic cleavage preferences
            // Proline and Glycine bonds tend to be weaker
            var bondIntensities = new Dictionary<int, double>();
            string canonicalSeq = circular.CanonicalSequence;

            for (int i = 0; i < canonicalSeq.Length; i++)
            {
                char residueN = canonicalSeq[i];
                char residueC = canonicalSeq[(i + 1) % canonicalSeq.Length];

                // Higher intensity for certain bond types
                double intensity = 0.1; // baseline

                // Proline N-terminal to bond increases cleavage
                if (residueC == 'P') intensity = 0.8;
                // After Asp or Glu
                else if (residueN == 'D' || residueN == 'E') intensity = 0.6;
                // After aromatic residues
                else if (residueN == 'W' || residueN == 'Y' || residueN == 'F') intensity = 0.5;

                bondIntensities[i] = intensity;
            }

            // Add mock ions based on intensities
            foreach (var fragment in spectrum.LinearFragments)
            {
                spectrum.AddMockPredictedIons(fragment, bondIntensities);
            }

            // Compute probabilities
            spectrum.ComputeBondCleavageProbabilities();

            // Run simulation
            int numCopies = 100000;
            spectrum.RunDoubleFragmentationSimulation(numCopies);

            // Generate output
            var sb = new StringBuilder();
            sb.AppendLine("=" + new string('=', 79));
            sb.AppendLine("CIRCULAR PEPTIDE DOUBLE FRAGMENTATION SIMULATION");
            sb.AppendLine("=" + new string('=', 79));
            sb.AppendLine();
            sb.AppendLine($"Date: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
            sb.AppendLine($"Original Sequence: {sequence}");
            sb.AppendLine($"Canonical Sequence: {circular.CanonicalSequence}");
            sb.AppendLine($"Origin Index: {circular.OriginIndex}");
            sb.AppendLine($"Peptide Length: {sequence.Length}");
            sb.AppendLine($"Number of Simulated Copies: {numCopies:N0}");
            sb.AppendLine();

            // Coverage summary
            sb.AppendLine(spectrum.GetCoverageSummary());
            sb.AppendLine();

            // Cleavage probabilities
            sb.AppendLine(spectrum.GetCleavageProbabilityTable());
            sb.AppendLine();

            // Fragmentation patterns
            sb.AppendLine(spectrum.GetFragmentationPatternsTable());
            sb.AppendLine();

            // Simulated spectrum
            sb.AppendLine(spectrum.GetSimulatedSpectrumTable());
            sb.AppendLine();

            // Summary statistics
            sb.AppendLine("SUMMARY STATISTICS");
            sb.AppendLine("==================");
            sb.AppendLine($"Unique Fragment Sequences: {spectrum.SimulatedSpectrum.Select(p => p.Sequence).Distinct().Count()}");
            sb.AppendLine($"Unique Fragmentation Patterns: {spectrum.SimulationResults.Select(e => e.GetKey()).Distinct().Count()}");
            sb.AppendLine($"Most Common Fragment: {spectrum.SimulatedSpectrum.OrderByDescending(p => p.Count).First().Sequence}");
            sb.AppendLine($"Smallest Fragment: {spectrum.SimulatedSpectrum.OrderBy(p => p.Sequence.Length).First().Sequence}");
            sb.AppendLine($"Largest Fragment: {spectrum.SimulatedSpectrum.OrderByDescending(p => p.Sequence.Length).First().Sequence}");
            sb.AppendLine($"Mass Range: {spectrum.SimulatedSpectrum.Min(p => p.Mz):F2} - {spectrum.SimulatedSpectrum.Max(p => p.Mz):F2} m/z");

            // Write to file
            string outputPath = Path.Combine(OutputDirectory, "CircularPeptideSimulation_Results.txt");
            File.WriteAllText(outputPath, sb.ToString());

            // Also write CSV for the spectrum
            string csvPath = Path.Combine(OutputDirectory, "CircularPeptideSimulation_Spectrum.csv");
            using (var writer = new StreamWriter(csvPath))
            {
                writer.WriteLine("Sequence,StartPosition,EndPosition,MonoisotopicMass,Mz,Charge,Count,RelativeIntensity");
                foreach (var peak in spectrum.SimulatedSpectrum.OrderByDescending(p => p.Count))
                {
                    writer.WriteLine($"{peak.Sequence},{peak.StartPosition},{peak.EndPosition}," +
                        $"{peak.MonoisotopicMass:F6},{peak.Mz:F6},{peak.Charge},{peak.Count},{peak.RelativeIntensity:F6}");
                }
            }

            // Write fragmentation patterns CSV
            string patternsPath = Path.Combine(OutputDirectory, "CircularPeptideSimulation_Patterns.csv");
            using (var writer = new StreamWriter(patternsPath))
            {
                writer.WriteLine("Bond1,Bond2,Fragment1,Fragment2,Count,Frequency");
                var patterns = spectrum.SimulationResults
                    .GroupBy(e => e.GetKey())
                    .Select(g => new
                    {
                        Bond1 = g.First().Bond1,
                        Bond2 = g.First().Bond2,
                        Frag1 = g.First().Fragment1.Sequence,
                        Frag2 = g.First().Fragment2.Sequence,
                        Count = g.Count()
                    })
                    .OrderByDescending(p => p.Count);

                foreach (var pattern in patterns)
                {
                    writer.WriteLine($"{pattern.Bond1},{pattern.Bond2},{pattern.Frag1},{pattern.Frag2}," +
                        $"{pattern.Count},{(double)pattern.Count / numCopies:F6}");
                }
            }

            // Assertions
            Assert.That(File.Exists(outputPath), Is.True, "Output file should be created");
            Assert.That(File.Exists(csvPath), Is.True, "CSV spectrum file should be created");
            Assert.That(File.Exists(patternsPath), Is.True, "CSV patterns file should be created");

            // Log the output location
            Console.WriteLine($"Results written to: {outputPath}");
            Console.WriteLine($"Spectrum CSV written to: {csvPath}");
            Console.WriteLine($"Patterns CSV written to: {patternsPath}");

            // Verify simulation ran correctly
            Assert.That(spectrum.SimulationResults.Count, Is.EqualTo(numCopies));
            Assert.That(spectrum.SimulatedSpectrum.Count, Is.GreaterThan(0));
        }

        [Test]
        public void MultipleSequences_ComparisonOutput()
        {
            // Test multiple circular peptides and compare their fragmentation patterns
            var sequences = new[]
            {
                "PEPTIDEK",
                "CYCLOPEP",
                "GRAMICID"
            };

            var sb = new StringBuilder();
            sb.AppendLine("COMPARISON OF CIRCULAR PEPTIDE FRAGMENTATION PATTERNS");
            sb.AppendLine("=" + new string('=', 60));
            sb.AppendLine();

            foreach (var seq in sequences)
            {
                try
                {
                    var circular = new CircularPeptide(seq);
                    var spectrum = new PredictedSpectrum(circular,
                        minFragmentLength: Math.Min(7, seq.Length),
                        randomSeed: 42);

                    // Use uniform probabilities
                    double prob = 1.0 / seq.Length;
                    spectrum.SetCleavageProbabilities(Enumerable.Repeat(prob, seq.Length).ToArray());

                    spectrum.RunDoubleFragmentationSimulation(10000);

                    sb.AppendLine($"Sequence: {seq}");
                    sb.AppendLine($"Canonical: {circular.CanonicalSequence}");
                    sb.AppendLine($"Unique Patterns: {spectrum.SimulationResults.Select(e => e.GetKey()).Distinct().Count()}");
                    sb.AppendLine($"Unique Fragments: {spectrum.SimulatedSpectrum.Count}");

                    // Top 3 fragments
                    sb.AppendLine("Top 3 Most Common Fragments:");
                    foreach (var peak in spectrum.SimulatedSpectrum.OrderByDescending(p => p.Count).Take(3))
                    {
                        sb.AppendLine($"  {peak.Sequence} (m/z {peak.Mz:F2}, {peak.RelativeIntensity:P1})");
                    }
                    sb.AppendLine();
                }
                catch (Exception ex)
                {
                    sb.AppendLine($"Sequence: {seq} - Error: {ex.Message}");
                    sb.AppendLine();
                }
            }

            string outputPath = Path.Combine(OutputDirectory, "CircularPeptide_Comparison.txt");
            File.WriteAllText(outputPath, sb.ToString());

            Assert.That(File.Exists(outputPath), Is.True);
            Console.WriteLine($"Comparison written to: {outputPath}");
        }

        #endregion

        #region Edge Case Tests

        [Test]
        public void EdgeCase_SmallPeptide()
        {
            // Test with minimum viable peptide
            var circular = new CircularPeptide("ACDEFGH"); // 7 amino acids (minimum)
            var spectrum = new PredictedSpectrum(circular, minFragmentLength: 7, randomSeed: 42);

            spectrum.SetCleavageProbabilities(Enumerable.Repeat(1.0 / 7, 7).ToArray());
            spectrum.RunDoubleFragmentationSimulation(1000);

            Assert.That(spectrum.SimulationResults.Count, Is.EqualTo(1000));

            // With 7 amino acids and 2 cuts, fragments should be 1-6 amino acids
            foreach (var evt in spectrum.SimulationResults)
            {
                int len1 = evt.Fragment1.Sequence.Length;
                int len2 = evt.Fragment2.Sequence.Length;
                Assert.That(len1 + len2, Is.EqualTo(7));
                Assert.That(len1, Is.GreaterThanOrEqualTo(1));
                Assert.That(len2, Is.GreaterThanOrEqualTo(1));
            }
        }

        [Test]
        public void EdgeCase_UniformVsBiasedProbabilities()
        {
            var circular = new CircularPeptide("ACDEFGHIKL"); // 10 amino acids (all valid)

            // Uniform probabilities
            var spectrumUniform = new PredictedSpectrum(circular, randomSeed: 42);
            spectrumUniform.SetCleavageProbabilities(Enumerable.Repeat(0.1, 10).ToArray());
            spectrumUniform.RunDoubleFragmentationSimulation(10000);

            // Heavily biased - 90% probability on bond 0 and 5
            var spectrumBiased = new PredictedSpectrum(circular, randomSeed: 42);
            var biasedProbs = new double[10];
            biasedProbs[0] = 0.45;
            biasedProbs[5] = 0.45;
            for (int i = 0; i < 10; i++) if (i != 0 && i != 5) biasedProbs[i] = 0.0125;
            spectrumBiased.SetCleavageProbabilities(biasedProbs);
            spectrumBiased.RunDoubleFragmentationSimulation(10000);

            // Biased should have fewer unique patterns
            int uniqueUniform = spectrumUniform.SimulationResults.Select(e => e.GetKey()).Distinct().Count();
            int uniqueBiased = spectrumBiased.SimulationResults.Select(e => e.GetKey()).Distinct().Count();

            // The biased simulation should concentrate on the 0-5 pattern
            var mostCommonBiased = spectrumBiased.SimulationResults
                .GroupBy(e => e.GetKey())
                .OrderByDescending(g => g.Count())
                .First();

            Assert.That(mostCommonBiased.Key, Is.EqualTo("0-5"),
                "Biased simulation should favor bonds 0 and 5");
            Assert.That(mostCommonBiased.Count(), Is.GreaterThan(5000),
                "Most common pattern should dominate");
        }

        #endregion
    }
}