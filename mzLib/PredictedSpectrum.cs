using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace Proteomics
{
    /// <summary>
    /// Represents a linear fragment extracted from a circular peptide,
    /// with tracking of its origin coordinates in the parent circular peptide.
    /// </summary>
    public class LinearFragment
    {
        /// <summary>
        /// The amino acid sequence of this linear fragment.
        /// </summary>
        public string Sequence { get; }

        /// <summary>
        /// The starting index in the parent circular peptide (0-based, relative to canonical origin).
        /// </summary>
        public int StartIndex { get; }

        /// <summary>
        /// The ending index (inclusive) in the parent circular peptide.
        /// May wrap around (i.e., EndIndex < StartIndex for fragments crossing the origin).
        /// </summary>
        public int EndIndex { get; }

        /// <summary>
        /// The length of the parent circular peptide.
        /// </summary>
        public int ParentLength { get; }

        /// <summary>
        /// Indicates whether this fragment crosses the origin (wraps around).
        /// </summary>
        public bool CrossesOrigin => EndIndex < StartIndex;

        /// <summary>
        /// Predicted fragment ions from Prosit HCD model.
        /// Populated after calling PredictFragmentation().
        /// </summary>
        public List<PredictedFragmentIon> PredictedIons { get; set; }

        public LinearFragment(string sequence, int startIndex, int endIndex, int parentLength)
        {
            Sequence = sequence ?? throw new ArgumentNullException(nameof(sequence));
            StartIndex = startIndex;
            EndIndex = endIndex;
            ParentLength = parentLength;
            PredictedIons = new List<PredictedFragmentIon>();
        }

        /// <summary>
        /// Converts a position in this linear fragment to the corresponding position
        /// in the parent circular peptide.
        /// </summary>
        public int ToCircularPosition(int linearPosition)
        {
            if (linearPosition < 0 || linearPosition >= Sequence.Length)
            {
                throw new ArgumentOutOfRangeException(nameof(linearPosition));
            }
            return (StartIndex + linearPosition) % ParentLength;
        }

        /// <summary>
        /// Converts a bond position in this linear fragment to the corresponding bond position
        /// in the parent circular peptide.
        /// </summary>
        public int ToCircularBondPosition(int linearBondPosition)
        {
            if (linearBondPosition < 0 || linearBondPosition >= Sequence.Length - 1)
            {
                throw new ArgumentOutOfRangeException(nameof(linearBondPosition));
            }
            return (StartIndex + linearBondPosition) % ParentLength;
        }

        public override string ToString()
        {
            string wrapIndicator = CrossesOrigin ? " [wraps]" : "";
            return $"{Sequence} (positions {StartIndex}-{EndIndex}{wrapIndicator})";
        }
    }

    /// <summary>
    /// Represents a predicted fragment ion from the Prosit HCD model.
    /// </summary>
    public class PredictedFragmentIon
    {
        public string IonType { get; set; }
        public int IonNumber { get; set; }
        public int Charge { get; set; }
        public double Mz { get; set; }
        public double Intensity { get; set; }
        public int CircularBondPosition { get; set; }

        public override string ToString()
        {
            return $"{IonType}{IonNumber}+{Charge} @ {Mz:F4} ({Intensity:P1}) [bond {CircularBondPosition}]";
        }
    }

    /// <summary>
    /// Represents the cleavage probability for a single peptide bond in the circular peptide.
    /// </summary>
    public class BondCleavageProbability
    {
        public int BondPosition { get; }
        public char ResidueN { get; }
        public char ResidueC { get; }
        public List<double> BIonIntensities { get; } = new List<double>();
        public List<double> YIonIntensities { get; } = new List<double>();

        /// <summary>
        /// The normalized probability of cleavage at this bond (0-1, sums to 1 across all bonds).
        /// </summary>
        public double CleavageProbability { get; set; }

        public double AverageIntensity
        {
            get
            {
                var all = BIonIntensities.Concat(YIonIntensities).ToList();
                return all.Count > 0 ? all.Average() : 0.0;
            }
        }

        public double MaxIntensity
        {
            get
            {
                var all = BIonIntensities.Concat(YIonIntensities).ToList();
                return all.Count > 0 ? all.Max() : 0.0;
            }
        }

        public int ObservationCount => BIonIntensities.Count + YIonIntensities.Count;
        public int FragmentCount { get; set; }

        public BondCleavageProbability(int bondPosition, char residueN, char residueC)
        {
            BondPosition = bondPosition;
            ResidueN = residueN;
            ResidueC = residueC;
        }

        public override string ToString()
        {
            return $"Bond {BondPosition} ({ResidueN}-{ResidueC}): prob={CleavageProbability:F4}, avg={AverageIntensity:F4}";
        }
    }

    /// <summary>
    /// Represents a double fragmentation event producing two linear fragments.
    /// </summary>
    public class DoubleFragmentationEvent
    {
        /// <summary>
        /// First bond position where cleavage occurred.
        /// </summary>
        public int Bond1 { get; }

        /// <summary>
        /// Second bond position where cleavage occurred.
        /// </summary>
        public int Bond2 { get; }

        /// <summary>
        /// First resulting linear fragment.
        /// </summary>
        public LinearFragment Fragment1 { get; }

        /// <summary>
        /// Second resulting linear fragment.
        /// </summary>
        public LinearFragment Fragment2 { get; }

        public DoubleFragmentationEvent(int bond1, int bond2, LinearFragment fragment1, LinearFragment fragment2)
        {
            Bond1 = Math.Min(bond1, bond2);
            Bond2 = Math.Max(bond1, bond2);
            Fragment1 = fragment1;
            Fragment2 = fragment2;
        }

        /// <summary>
        /// Creates a unique key for this fragmentation pattern.
        /// </summary>
        public string GetKey() => $"{Bond1}-{Bond2}";

        public override string ToString()
        {
            return $"Bonds [{Bond1}, {Bond2}]: {Fragment1.Sequence} + {Fragment2.Sequence}";
        }
    }

    /// <summary>
    /// Represents a unique fragment in the simulated spectrum.
    /// </summary>
    public class SimulatedFragmentPeak
    {
        public string Sequence { get; set; }
        public int StartPosition { get; set; }
        public int EndPosition { get; set; }
        public double MonoisotopicMass { get; set; }
        public double Mz { get; set; }
        public int Charge { get; set; }
        public int Count { get; set; }
        public double RelativeIntensity { get; set; }

        public override string ToString()
        {
            return $"{Sequence} [{StartPosition}-{EndPosition}] m/z={Mz:F4} Count={Count} ({RelativeIntensity:P2})";
        }
    }

    /// <summary>
    /// Handles fragmentation simulation and spectrum prediction for circular peptides.
    /// </summary>
    public class PredictedSpectrum
    {
        // Amino acid monoisotopic masses (average masses could also be used)
        private static readonly Dictionary<char, double> AminoAcidMasses = new Dictionary<char, double>
        {
            {'A', 71.03711}, {'R', 156.10111}, {'N', 114.04293}, {'D', 115.02694},
            {'C', 103.00919}, {'E', 129.04259}, {'Q', 128.05858}, {'G', 57.02146},
            {'H', 137.05891}, {'I', 113.08406}, {'L', 113.08406}, {'K', 128.09496},
            {'M', 131.04049}, {'F', 147.06841}, {'P', 97.05276}, {'S', 87.03203},
            {'T', 101.04768}, {'W', 186.07931}, {'Y', 163.06333}, {'V', 99.06841}
        };

        private const double WaterMass = 18.01056;
        private const double ProtonMass = 1.00727;

        public CircularPeptide CircularPeptide { get; }
        public List<LinearFragment> LinearFragments { get; private set; }
        public List<BondCleavageProbability> BondCleavageProbabilities { get; private set; }
        public List<DoubleFragmentationEvent> SimulationResults { get; private set; }
        public List<SimulatedFragmentPeak> SimulatedSpectrum { get; private set; }

        public int MinFragmentLength { get; }
        public int MaxFragmentLength { get; }
        public int MinFragmentCount { get; }
        public float CollisionEnergy { get; set; } = 30f;
        public int PrecursorCharge { get; set; } = 2;

        private Random _random;

        public PredictedSpectrum(
            CircularPeptide circularPeptide,
            int minFragmentLength = 7,
            int maxFragmentLength = 15,
            int minFragmentCount = 3,
            int? randomSeed = null)
        {
            CircularPeptide = circularPeptide ?? throw new ArgumentNullException(nameof(circularPeptide));
            MinFragmentLength = minFragmentLength;
            MaxFragmentLength = maxFragmentLength;
            MinFragmentCount = minFragmentCount;
            LinearFragments = new List<LinearFragment>();
            BondCleavageProbabilities = new List<BondCleavageProbability>();
            SimulationResults = new List<DoubleFragmentationEvent>();
            SimulatedSpectrum = new List<SimulatedFragmentPeak>();
            _random = randomSeed.HasValue ? new Random(randomSeed.Value) : new Random();

            ValidateParameters();
            GenerateLinearFragments();
            InitializeBondCleavageProbabilities();
        }

        private void InitializeBondCleavageProbabilities()
        {
            BondCleavageProbabilities.Clear();
            string sequence = CircularPeptide.CanonicalSequence;
            int length = sequence.Length;

            for (int i = 0; i < length; i++)
            {
                char residueN = sequence[i];
                char residueC = sequence[(i + 1) % length];
                BondCleavageProbabilities.Add(new BondCleavageProbability(i, residueN, residueC));
            }
        }

        /// <summary>
        /// Normalizes bond cleavage intensities to probabilities that sum to 1.
        /// </summary>
        public void NormalizeCleavageProbabilities()
        {
            double totalIntensity = BondCleavageProbabilities.Sum(b => b.AverageIntensity);

            if (totalIntensity <= 0)
            {
                // If no intensities, use uniform distribution
                double uniformProb = 1.0 / BondCleavageProbabilities.Count;
                foreach (var bond in BondCleavageProbabilities)
                {
                    bond.CleavageProbability = uniformProb;
                }
            }
            else
            {
                foreach (var bond in BondCleavageProbabilities)
                {
                    bond.CleavageProbability = bond.AverageIntensity / totalIntensity;
                }
            }

            // Verify normalization
            double sum = BondCleavageProbabilities.Sum(b => b.CleavageProbability);
            if (Math.Abs(sum - 1.0) > 0.0001)
            {
                // Adjust for floating point errors
                var last = BondCleavageProbabilities.Last();
                last.CleavageProbability += (1.0 - sum);
            }
        }

        /// <summary>
        /// Sets custom cleavage probabilities (must sum to 1).
        /// </summary>
        public void SetCleavageProbabilities(double[] probabilities)
        {
            if (probabilities.Length != BondCleavageProbabilities.Count)
            {
                throw new ArgumentException($"Expected {BondCleavageProbabilities.Count} probabilities, got {probabilities.Length}");
            }

            double sum = probabilities.Sum();
            if (Math.Abs(sum - 1.0) > 0.0001)
            {
                throw new ArgumentException($"Probabilities must sum to 1, got {sum}");
            }

            for (int i = 0; i < probabilities.Length; i++)
            {
                BondCleavageProbabilities[i].CleavageProbability = probabilities[i];
            }
        }

        /// <summary>
        /// Samples a bond position based on the cleavage probabilities.
        /// </summary>
        private int SampleBondPosition(int? excludeBond = null)
        {
            double totalProb = BondCleavageProbabilities
                .Where((b, i) => i != excludeBond)
                .Sum(b => b.CleavageProbability);

            double r = _random.NextDouble() * totalProb;
            double cumulative = 0;

            for (int i = 0; i < BondCleavageProbabilities.Count; i++)
            {
                if (i == excludeBond) continue;

                cumulative += BondCleavageProbabilities[i].CleavageProbability;
                if (r <= cumulative)
                {
                    return i;
                }
            }

            // Fallback (shouldn't happen)
            return BondCleavageProbabilities.Count - 1;
        }

        /// <summary>
        /// Creates two linear fragments from a circular peptide broken at two bond positions.
        /// </summary>
        private (LinearFragment, LinearFragment) CreateFragmentsFromBreaks(int bond1, int bond2)
        {
            string sequence = CircularPeptide.CanonicalSequence;
            int length = sequence.Length;

            // Ensure bond1 < bond2
            if (bond1 > bond2)
            {
                (bond1, bond2) = (bond2, bond1);
            }

            // Fragment 1: from position (bond1 + 1) to position bond2 (inclusive)
            // This is the arc from after bond1 to bond2
            int start1 = (bond1 + 1) % length;
            int end1 = bond2;
            int length1 = bond2 - bond1;
            string seq1 = ExtractCircularSubstring(sequence, start1, length1);

            // Fragment 2: from position (bond2 + 1) to position bond1 (inclusive), wrapping around
            // This is the arc from after bond2 back to bond1
            int start2 = (bond2 + 1) % length;
            int end2 = bond1;
            int length2 = length - length1;
            string seq2 = ExtractCircularSubstring(sequence, start2, length2);

            var frag1 = new LinearFragment(seq1, start1, end1, length);
            var frag2 = new LinearFragment(seq2, start2, end2, length);

            return (frag1, frag2);
        }

        /// <summary>
        /// Runs the double fragmentation simulation.
        /// </summary>
        /// <param name="numberOfCopies">Number of circular peptide copies to simulate (default 100,000).</param>
        public void RunDoubleFragmentationSimulation(int numberOfCopies = 100000)
        {
            SimulationResults.Clear();

            // Ensure probabilities are normalized
            if (BondCleavageProbabilities.All(b => b.CleavageProbability == 0))
            {
                NormalizeCleavageProbabilities();
            }

            for (int i = 0; i < numberOfCopies; i++)
            {
                // Sample first bond
                int bond1 = SampleBondPosition();

                // Sample second bond (different from first)
                int bond2 = SampleBondPosition(excludeBond: bond1);

                // Create the two fragments
                var (frag1, frag2) = CreateFragmentsFromBreaks(bond1, bond2);

                // Record the event
                var evt = new DoubleFragmentationEvent(bond1, bond2, frag1, frag2);
                SimulationResults.Add(evt);
            }

            // Build the simulated spectrum
            BuildSimulatedSpectrum();
        }

        /// <summary>
        /// Builds the simulated spectrum by counting fragment occurrences and calculating m/z.
        /// </summary>
        private void BuildSimulatedSpectrum()
        {
            SimulatedSpectrum.Clear();

            // Group fragments by sequence (accounting for both fragments from each event)
            var fragmentCounts = new Dictionary<string, (string Sequence, int Start, int End, int Count)>();

            foreach (var evt in SimulationResults)
            {
                // Process Fragment 1
                string key1 = $"{evt.Fragment1.Sequence}_{evt.Fragment1.StartIndex}_{evt.Fragment1.EndIndex}";
                if (fragmentCounts.ContainsKey(key1))
                {
                    var entry = fragmentCounts[key1];
                    fragmentCounts[key1] = (entry.Sequence, entry.Start, entry.End, entry.Count + 1);
                }
                else
                {
                    fragmentCounts[key1] = (evt.Fragment1.Sequence, evt.Fragment1.StartIndex, evt.Fragment1.EndIndex, 1);
                }

                // Process Fragment 2
                string key2 = $"{evt.Fragment2.Sequence}_{evt.Fragment2.StartIndex}_{evt.Fragment2.EndIndex}";
                if (fragmentCounts.ContainsKey(key2))
                {
                    var entry = fragmentCounts[key2];
                    fragmentCounts[key2] = (entry.Sequence, entry.Start, entry.End, entry.Count + 1);
                }
                else
                {
                    fragmentCounts[key2] = (evt.Fragment2.Sequence, evt.Fragment2.StartIndex, evt.Fragment2.EndIndex, 1);
                }
            }

            // Calculate m/z and relative intensities
            int totalFragments = SimulationResults.Count * 2;

            foreach (var kvp in fragmentCounts)
            {
                var (sequence, start, end, count) = kvp.Value;
                double mass = CalculatePeptideMass(sequence);
                int charge = 1; // Default to singly charged

                var peak = new SimulatedFragmentPeak
                {
                    Sequence = sequence,
                    StartPosition = start,
                    EndPosition = end,
                    MonoisotopicMass = mass,
                    Mz = (mass + charge * ProtonMass) / charge,
                    Charge = charge,
                    Count = count,
                    RelativeIntensity = (double)count / totalFragments
                };

                SimulatedSpectrum.Add(peak);
            }

            // Sort by m/z
            SimulatedSpectrum = SimulatedSpectrum.OrderBy(p => p.Mz).ToList();
        }

        /// <summary>
        /// Calculates the monoisotopic mass of a peptide sequence.
        /// </summary>
        private double CalculatePeptideMass(string sequence)
        {
            double mass = WaterMass; // Add water for the termini
            foreach (char aa in sequence)
            {
                if (AminoAcidMasses.TryGetValue(aa, out double aaMass))
                {
                    mass += aaMass;
                }
                else
                {
                    throw new ArgumentException($"Unknown amino acid: {aa}");
                }
            }
            return mass;
        }

        /// <summary>
        /// Gets a summary table of unique double fragmentation patterns.
        /// </summary>
        public string GetFragmentationPatternsTable()
        {
            var sb = new System.Text.StringBuilder();
            sb.AppendLine("Double Fragmentation Pattern Summary");
            sb.AppendLine("=====================================");
            sb.AppendLine();

            // Group by bond pair
            var patterns = SimulationResults
                .GroupBy(e => e.GetKey())
                .Select(g => new
                {
                    Key = g.Key,
                    Bond1 = g.First().Bond1,
                    Bond2 = g.First().Bond2,
                    Fragment1 = g.First().Fragment1.Sequence,
                    Fragment2 = g.First().Fragment2.Sequence,
                    Count = g.Count(),
                    Frequency = (double)g.Count() / SimulationResults.Count
                })
                .OrderByDescending(p => p.Count)
                .ToList();

            sb.AppendLine($"Total simulated events: {SimulationResults.Count:N0}");
            sb.AppendLine($"Unique fragmentation patterns: {patterns.Count}");
            sb.AppendLine();
            sb.AppendLine("Bond1 | Bond2 | Fragment 1          | Fragment 2          | Count    | Frequency");
            sb.AppendLine(new string('-', 95));

            foreach (var pattern in patterns)
            {
                string frag1Display = pattern.Fragment1.Length > 18
                    ? pattern.Fragment1.Substring(0, 15) + "..."
                    : pattern.Fragment1.PadRight(18);
                string frag2Display = pattern.Fragment2.Length > 18
                    ? pattern.Fragment2.Substring(0, 15) + "..."
                    : pattern.Fragment2.PadRight(18);

                sb.AppendLine($"  {pattern.Bond1,3} |   {pattern.Bond2,3} | {frag1Display} | {frag2Display} | {pattern.Count,8:N0} | {pattern.Frequency,8:P2}");
            }

            return sb.ToString();
        }

        /// <summary>
        /// Gets the simulated spectrum as a formatted table.
        /// </summary>
        public string GetSimulatedSpectrumTable()
        {
            var sb = new System.Text.StringBuilder();
            sb.AppendLine("Simulated Fragmentation Spectrum");
            sb.AppendLine("=================================");
            sb.AppendLine();
            sb.AppendLine($"Circular Peptide: {CircularPeptide.CanonicalSequence}");
            sb.AppendLine($"Total Events: {SimulationResults.Count:N0}");
            sb.AppendLine($"Unique Fragments: {SimulatedSpectrum.Count}");
            sb.AppendLine();
            sb.AppendLine("  m/z       | Mass      | Charge | Count    | Rel.Int. | Pos.     | Sequence");
            sb.AppendLine(new string('-', 110));

            foreach (var peak in SimulatedSpectrum.OrderByDescending(p => p.Count))
            {
                string seqDisplay = peak.Sequence.Length > 25
                    ? peak.Sequence.Substring(0, 22) + "..."
                    : peak.Sequence;

                sb.AppendLine($" {peak.Mz,10:F4} | {peak.MonoisotopicMass,9:F4} |   {peak.Charge,2}   | {peak.Count,8:N0} | {peak.RelativeIntensity,7:P2} | {peak.StartPosition,2}-{peak.EndPosition,2}   | {seqDisplay}");
            }

            return sb.ToString();
        }

        /// <summary>
        /// Gets the cleavage probability distribution table.
        /// </summary>
        public string GetCleavageProbabilityTable()
        {
            var sb = new System.Text.StringBuilder();
            sb.AppendLine("Bond Cleavage Probabilities");
            sb.AppendLine("===========================");
            sb.AppendLine();
            sb.AppendLine("Position | Bond | Probability | Avg Intensity");
            sb.AppendLine(new string('-', 50));

            foreach (var bond in BondCleavageProbabilities)
            {
                sb.AppendLine($"   {bond.BondPosition,3}   | {bond.ResidueN}-{bond.ResidueC} |   {bond.CleavageProbability,8:F4}  |   {bond.AverageIntensity,8:F4}");
            }

            sb.AppendLine(new string('-', 50));
            sb.AppendLine($"   Total Probability: {BondCleavageProbabilities.Sum(b => b.CleavageProbability):F4}");

            return sb.ToString();
        }

        #region Existing methods (unchanged)

        private void ValidateParameters()
        {
            int length = CircularPeptide.Sequence.Length;
            if (length < MinFragmentLength)
                throw new ArgumentException($"Circular peptide length ({length}) is less than minimum fragment length ({MinFragmentLength}).");
            if (MinFragmentLength > MaxFragmentLength)
                throw new ArgumentException($"Minimum fragment length ({MinFragmentLength}) cannot exceed maximum ({MaxFragmentLength}).");
            if (MinFragmentCount < 3)
                throw new ArgumentException($"Minimum fragment count must be at least 3 for proper coverage.");
        }

        public void GenerateLinearFragments()
        {
            LinearFragments.Clear();
            int circularLength = CircularPeptide.Sequence.Length;
            string canonicalSequence = CircularPeptide.CanonicalSequence;

            if (circularLength <= MaxFragmentLength)
            {
                GenerateFragmentsForShortPeptide(canonicalSequence, circularLength);
                return;
            }

            int fragmentLength = CalculateOptimalFragmentLength(circularLength);
            int fragmentCount = CalculateFragmentCount(circularLength, fragmentLength);
            double stepSize = (double)circularLength / fragmentCount;

            for (int i = 0; i < fragmentCount; i++)
            {
                int startIndex = (int)Math.Round(i * stepSize) % circularLength;
                int endIndex = (startIndex + fragmentLength - 1) % circularLength;
                string fragmentSequence = ExtractCircularSubstring(canonicalSequence, startIndex, fragmentLength);
                LinearFragments.Add(new LinearFragment(fragmentSequence, startIndex, endIndex, circularLength));
            }

            EnsureCompleteCoverage(canonicalSequence);
        }

        private void GenerateFragmentsForShortPeptide(string sequence, int length)
        {
            int fragmentLength = Math.Max(Math.Min(length, MaxFragmentLength), MinFragmentLength);
            if (length <= MinFragmentLength) fragmentLength = length;

            int fragmentCount = Math.Max(MinFragmentCount, length);
            double stepSize = (double)length / fragmentCount;

            for (int i = 0; i < fragmentCount && i < length; i++)
            {
                int startIndex = (int)Math.Round(i * stepSize) % length;
                int actualLength = Math.Min(fragmentLength, length);
                int endIndex = (startIndex + actualLength - 1) % length;
                string fragmentSequence = ExtractCircularSubstring(sequence, startIndex, actualLength);

                if (fragmentSequence.Length >= MinFragmentLength)
                    LinearFragments.Add(new LinearFragment(fragmentSequence, startIndex, endIndex, length));
            }

            if (LinearFragments.Count < MinFragmentCount && length >= MinFragmentLength)
            {
                for (int i = LinearFragments.Count; i < MinFragmentCount; i++)
                {
                    int startIndex = (i * length / MinFragmentCount) % length;
                    int endIndex = (startIndex + fragmentLength - 1) % length;
                    string fragmentSequence = ExtractCircularSubstring(sequence, startIndex, fragmentLength);
                    LinearFragments.Add(new LinearFragment(fragmentSequence, startIndex, endIndex, length));
                }
            }
        }

        private int CalculateOptimalFragmentLength(int circularLength) =>
            Math.Max(Math.Min(circularLength / 2, MaxFragmentLength), MinFragmentLength);

        private int CalculateFragmentCount(int circularLength, int fragmentLength) =>
            Math.Max((int)Math.Ceiling(circularLength * 1.5 / fragmentLength), MinFragmentCount);

        private string ExtractCircularSubstring(string sequence, int startIndex, int length)
        {
            int seqLength = sequence.Length;
            startIndex = startIndex % seqLength;
            if (startIndex + length <= seqLength)
                return sequence.Substring(startIndex, length);
            return sequence.Substring(startIndex) + sequence.Substring(0, length - (seqLength - startIndex));
        }

        private void EnsureCompleteCoverage(string sequence)
        {
            int length = sequence.Length;
            bool[] covered = new bool[length];

            foreach (var fragment in LinearFragments)
                for (int i = 0; i < fragment.Sequence.Length; i++)
                    covered[fragment.ToCircularPosition(i)] = true;

            for (int i = 0; i < length; i++)
            {
                if (!covered[i])
                {
                    int fragmentLength = Math.Max(Math.Min(MaxFragmentLength, length), MinFragmentLength);
                    int endIndex = (i + fragmentLength - 1) % length;
                    string fragmentSequence = ExtractCircularSubstring(sequence, i, fragmentLength);
                    LinearFragments.Add(new LinearFragment(fragmentSequence, i, endIndex, length));
                    for (int j = 0; j < fragmentLength; j++)
                        covered[(i + j) % length] = true;
                }
            }
        }

        public void ComputeBondCleavageProbabilities()
        {
            foreach (var bond in BondCleavageProbabilities)
            {
                bond.BIonIntensities.Clear();
                bond.YIonIntensities.Clear();
                bond.FragmentCount = 0;
            }

            int circularLength = CircularPeptide.CanonicalSequence.Length;
            int[] bondCoverage = new int[circularLength];

            foreach (var fragment in LinearFragments)
                for (int localBondPos = 0; localBondPos < fragment.Sequence.Length - 1; localBondPos++)
                    bondCoverage[fragment.ToCircularBondPosition(localBondPos)]++;

            for (int i = 0; i < circularLength; i++)
                BondCleavageProbabilities[i].FragmentCount = bondCoverage[i];

            foreach (var fragment in LinearFragments)
                ProcessFragmentIons(fragment);

            NormalizeCleavageProbabilities();
        }

        private void ProcessFragmentIons(LinearFragment fragment)
        {
            int fragmentLength = fragment.Sequence.Length;

            foreach (var ion in fragment.PredictedIons)
            {
                int localBondPosition = -1;

                if (ion.IonType.Equals("b", StringComparison.OrdinalIgnoreCase))
                    localBondPosition = ion.IonNumber - 1;
                else if (ion.IonType.Equals("y", StringComparison.OrdinalIgnoreCase))
                    localBondPosition = fragmentLength - ion.IonNumber - 1;

                if (localBondPosition >= 0 && localBondPosition < fragmentLength - 1)
                {
                    int circularBondPosition = fragment.ToCircularBondPosition(localBondPosition);
                    ion.CircularBondPosition = circularBondPosition;
                    var bond = BondCleavageProbabilities[circularBondPosition];

                    if (ion.IonType.Equals("b", StringComparison.OrdinalIgnoreCase))
                        bond.BIonIntensities.Add(ion.Intensity);
                    else if (ion.IonType.Equals("y", StringComparison.OrdinalIgnoreCase))
                        bond.YIonIntensities.Add(ion.Intensity);
                }
            }
        }

        public void AddMockPredictedIons(LinearFragment fragment, Dictionary<int, double> bondIntensities)
        {
            int fragmentLength = fragment.Sequence.Length;
            for (int localBondPos = 0; localBondPos < fragmentLength - 1; localBondPos++)
            {
                int circularBondPos = fragment.ToCircularBondPosition(localBondPos);
                double intensity = bondIntensities.TryGetValue(circularBondPos, out var val) ? val : 0.1;

                fragment.PredictedIons.Add(new PredictedFragmentIon
                {
                    IonType = "b",
                    IonNumber = localBondPos + 1,
                    Charge = 1,
                    Intensity = intensity,
                    CircularBondPosition = circularBondPos
                });

                fragment.PredictedIons.Add(new PredictedFragmentIon
                {
                    IonType = "y",
                    IonNumber = fragmentLength - localBondPos - 1,
                    Charge = 1,
                    Intensity = intensity,
                    CircularBondPosition = circularBondPos
                });
            }
        }

        public string GetCoverageSummary()
        {
            if (LinearFragments.Count == 0) return "No fragments generated.";
            int[] coverage = new int[CircularPeptide.Sequence.Length];
            foreach (var fragment in LinearFragments)
                for (int i = 0; i < fragment.Sequence.Length; i++)
                    coverage[fragment.ToCircularPosition(i)]++;
            return $"Fragments: {LinearFragments.Count}, Coverage: min={coverage.Min()}x, max={coverage.Max()}x, avg={coverage.Average():F1}x";
        }

        public override string ToString() => $"PredictedSpectrum for {CircularPeptide.CanonicalSequence} - {GetCoverageSummary()}";

        #endregion
    }
}