﻿using System;

namespace MassSpectrometry
{
    [Serializable]
    public class IndexedMassSpectralPeak : IIndexedPeak
    {
        public int ZeroBasedScanIndex { get; init; }
        public double Mz { get; init; }
        public double M => Mz;
        public double RetentionTime { get; init; }
        public double Intensity { get; init; }
        public IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedScanIndex, double retentionTime)
        {
            this.Mz = mz;
            this.ZeroBasedScanIndex = zeroBasedScanIndex;
            this.RetentionTime = retentionTime;
            this.Intensity = intensity;
        }
        public override bool Equals(object obj)
        {
            return Equals((IndexedMassSpectralPeak)obj);
        }

        public bool Equals(IIndexedPeak other)
        {
            return Equals((IndexedMassSpectralPeak)other);
        }

        public bool Equals(IndexedMassSpectralPeak other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return other.ZeroBasedScanIndex == this.ZeroBasedScanIndex
                   && Math.Abs(other.Mz - this.Mz) < 1e-9;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(Mz, ZeroBasedScanIndex);
        }

        public override string ToString()
        {
            return Mz.ToString("F3") + "; " + ZeroBasedScanIndex;
        }
    }
}