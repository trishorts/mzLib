using System;

namespace Proteomics
{
    /// <summary>
    /// Provides utilities for working with circular peptides.
    /// Circular peptides have no defined start/end, so we need a consistent
    /// coordinate system to compare and identify them.
    /// </summary>
    public class CircularPeptide
    {
        /// <summary>
        /// Gets or sets the original sequence of the circular peptide.
        /// </summary>
        public string Sequence { get; set; }

        /// <summary>
        /// Gets or sets the index in the original sequence where the canonical origin is located.
        /// </summary>
        public int OriginIndex { get; private set; }

        /// <summary>
        /// Gets the canonical (normalized) sequence starting from the assigned origin.
        /// </summary>
        public string CanonicalSequence => GetRotatedSequence(OriginIndex);

        public CircularPeptide(string sequence)
        {
            Sequence = sequence ?? throw new ArgumentNullException(nameof(sequence));
            OriginIndex = AssignOrigin(sequence);
        }

        /// <summary>
        /// Assigns the origin (position 0) of a circular peptide by finding the 
        /// lexicographically smallest rotation of the sequence.
        /// 
        /// Algorithm:
        /// 1. Find all positions of the alphabetically lowest amino acid
        /// 2. If there's only one, that's the origin
        /// 3. If there are multiple, compare subsequent amino acids to break ties
        /// 4. Continue until a unique minimum rotation is found
        /// </summary>
        /// <param name="sequence">The amino acid sequence of the circular peptide.</param>
        /// <returns>The 0-based index in the original sequence where the origin should be.</returns>
        public static int AssignOrigin(string sequence)
        {
            if (string.IsNullOrEmpty(sequence))
            {
                return 0;
            }

            int length = sequence.Length;

            if (length == 1)
            {
                return 0;
            }

            // Find the lexicographically smallest rotation using Booth's algorithm concept
            // but implemented in a straightforward way for clarity

            // Step 1: Find the minimum character in the sequence
            char minChar = sequence[0];
            for (int i = 1; i < length; i++)
            {
                if (sequence[i] < minChar)
                {
                    minChar = sequence[i];
                }
            }

            // Step 2: Find all positions where the minimum character occurs
            var candidatePositions = new System.Collections.Generic.List<int>();
            for (int i = 0; i < length; i++)
            {
                if (sequence[i] == minChar)
                {
                    candidatePositions.Add(i);
                }
            }

            // Step 3: If only one candidate, we're done
            if (candidatePositions.Count == 1)
            {
                return candidatePositions[0];
            }

            // Step 4: Compare rotations starting from each candidate position
            // to find the lexicographically smallest one
            int bestPosition = candidatePositions[0];

            for (int i = 1; i < candidatePositions.Count; i++)
            {
                int currentPosition = candidatePositions[i];

                if (CompareRotations(sequence, currentPosition, bestPosition) < 0)
                {
                    bestPosition = currentPosition;
                }
            }

            return bestPosition;
        }

        /// <summary>
        /// Compares two rotations of a circular sequence lexicographically.
        /// </summary>
        /// <param name="sequence">The original sequence.</param>
        /// <param name="startA">Starting index of rotation A.</param>
        /// <param name="startB">Starting index of rotation B.</param>
        /// <returns>
        /// Negative if rotation A is smaller, positive if rotation B is smaller, 0 if equal.
        /// </returns>
        private static int CompareRotations(string sequence, int startA, int startB)
        {
            int length = sequence.Length;

            for (int offset = 0; offset < length; offset++)
            {
                char charA = sequence[(startA + offset) % length];
                char charB = sequence[(startB + offset) % length];

                if (charA < charB)
                {
                    return -1;
                }
                if (charA > charB)
                {
                    return 1;
                }
            }

            return 0; // Rotations are identical (sequence has a repeating pattern)
        }

        /// <summary>
        /// Gets the sequence rotated to start at the specified index.
        /// </summary>
        /// <param name="startIndex">The index to start from.</param>
        /// <returns>The rotated sequence.</returns>
        public string GetRotatedSequence(int startIndex)
        {
            if (string.IsNullOrEmpty(Sequence))
            {
                return Sequence;
            }

            int length = Sequence.Length;
            startIndex = ((startIndex % length) + length) % length; // Handle negative indices

            return Sequence.Substring(startIndex) + Sequence.Substring(0, startIndex);
        }

        /// <summary>
        /// Checks if two circular peptides are equivalent (same sequence when rotated).
        /// </summary>
        /// <param name="other">The other circular peptide to compare.</param>
        /// <returns>True if the peptides are equivalent rotations of each other.</returns>
        public bool IsEquivalentTo(CircularPeptide other)
        {
            if (other == null || Sequence.Length != other.Sequence.Length)
            {
                return false;
            }

            return CanonicalSequence == other.CanonicalSequence;
        }

        public override string ToString()
        {
            return $"CircularPeptide: {CanonicalSequence} (Origin at index {OriginIndex} of original: {Sequence})";
        }
    }
}