using MathNet.Numerics;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.Globalization;
using System.Linq;
using System.Runtime.ConstrainedExecution;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.MzSpectra
{
    public class CrossCorrelation
    {
        public CrossCorrelation(MzSpectrum experimentalSpectrum, MzSpectrum theoreticalSpectrum, SpectrumNormalizationScheme scheme, double toleranceInPpm, bool allPeaks, double filterOutBelowThisMz = 300)
        {
            ExperimentalYArray = Normalize(FilterOutIonsBelowThisMz(experimentalSpectrum.XArray, experimentalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            ExperimentalXArray = FilterOutIonsBelowThisMz(experimentalSpectrum.XArray, experimentalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            TheoreticalYArray = Normalize(FilterOutIonsBelowThisMz(theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            TheoreticalXArray = FilterOutIonsBelowThisMz(theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            LocalPpmTolerance = toleranceInPpm;
            normalizationScheme = scheme;
            _intensityPairs = IntensityPairs(allPeaks, toleranceInPpm);
        }

        public CrossCorrelation(MzSpectrum experimentalSpectrum, double[] theoreticalX, double[] theoreticalY, SpectrumNormalizationScheme scheme, double toleranceInPpm, bool allPeaks, double filterOutBelowThisMz = 300)
        {
            ExperimentalYArray = Normalize(FilterOutIonsBelowThisMz(experimentalSpectrum.XArray, experimentalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            ExperimentalXArray = FilterOutIonsBelowThisMz(experimentalSpectrum.XArray, experimentalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            TheoreticalYArray = Normalize(FilterOutIonsBelowThisMz(theoreticalX, theoreticalY, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            TheoreticalXArray = FilterOutIonsBelowThisMz(theoreticalX, theoreticalY, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            LocalPpmTolerance = toleranceInPpm;
            normalizationScheme = scheme;
            _intensityPairs = IntensityPairs(allPeaks, toleranceInPpm);
        }

        public CrossCorrelation(double[] P_XArray, double[] P_YArray, double[] Q_XArray, double[] Q_YArray, SpectrumNormalizationScheme scheme, double toleranceInPpm, bool allPeaks, double filterOutBelowThisMz = 300)
        {
            ExperimentalYArray = Normalize(FilterOutIonsBelowThisMz(P_XArray, P_YArray, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            ExperimentalXArray = FilterOutIonsBelowThisMz(P_XArray, P_YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            TheoreticalYArray = Normalize(FilterOutIonsBelowThisMz(Q_XArray, Q_YArray, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            TheoreticalXArray = FilterOutIonsBelowThisMz(Q_XArray, Q_YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            LocalPpmTolerance = toleranceInPpm;
            _intensityPairs = IntensityPairs(allPeaks, toleranceInPpm);
        }
        public double[] ExperimentalYArray { get; private set; }
        public double[] ExperimentalXArray { get; private set; }
        public double[] TheoreticalYArray { get; private set; }
        public double[] TheoreticalXArray { get; private set; }

        private SpectrumNormalizationScheme normalizationScheme;

        private double LocalPpmTolerance;

        private readonly List<(double, double)> _intensityPairs = new();

        public List<(double, double)> intensityPairs
        { get { return _intensityPairs; } }


        /// <summary>
        /// All peaks with mz less than the cutOff will be filtered out. default to zero to remove an mz values that are accidentally negative. this is an unexpected error. 
        private static List<(double, double)> FilterOutIonsBelowThisMz(double[] spectrumX, double[] spectrumY, double filterOutBelowThisMz = 0)
        {
            if (spectrumY.Length == 0)
            {
                throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Empty YArray in spectrum."));
            }
            if (spectrumY.Sum() == 0)
            {
                throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Spectrum has no intensity."));
            }

            List<(double, double)> spectrumWithMzCutoff = new List<(double, double)>();
            for (int i = 0; i < spectrumX.Length; i++)
            {
                //second conditional to avoid getting an accidental negative intensities
                if (spectrumX[i] >= filterOutBelowThisMz && spectrumY[i] >= 0)
                {
                    spectrumWithMzCutoff.Add((spectrumX[i], spectrumY[i]));
                }
            }
            return spectrumWithMzCutoff;
        }

        /// <summary>
        /// Every spectrum gets normalized when the SpectralSimilarity object gets created. This methods sends the spectra to the appropriate normalization.
        /// </summary>
        /// <param name="spectrum"></param>
        /// <param name="scheme"></param>
        /// <returns></returns>
        private double[] Normalize(double[] spectrum, SpectrumNormalizationScheme scheme)
        {
            if (spectrum.Length == 0)
            {
                return null;
            }

            return scheme switch
            {
                SpectrumNormalizationScheme.mostAbundantPeak => NormalizeMostAbundantPeak(spectrum),
                SpectrumNormalizationScheme.spectrumSum => NormalizeSpectrumSum(spectrum),
                SpectrumNormalizationScheme.squareRootSpectrumSum => NormalizeSquareRootSpectrumSum(spectrum),
                _ => spectrum,
            };
        }

        /// <summary>
        /// Intensity Pairs a computed immediately upon creation of the SpectralSimilarity object. That way they can be used in all the methods without being recomputed.
        /// We loop throught the secondaryXArray under the assumption that it is the shorter of the two arrays (i.e. typically the theoretical spectrum).
        /// Experimental spectrum defaults to 200 peaks and is therefore usually longer.
        /// We sort intensities in descending order so that when we make peak pairs, we're choosing pairs with the highest intensity so long as they are with mz range. 
        /// Sometimes you could have two peaks in mz range and I don't think you want to pair the lesser intensity peak first just because it is closer in mass.
        /// 
        /// Intensity Pair: (Experimental Intensity , Theoretical Intensity)
        /// 
        /// </summary>

        private List<(double, double)> IntensityPairs(bool allPeaks, double localPpmTolerance, double[] experimentalYArray = null, double[] theoreticalYArray = null)
        {
            if (experimentalYArray == null) experimentalYArray = ExperimentalYArray;
            if (theoreticalYArray == null) theoreticalYArray = TheoreticalYArray;

            if (experimentalYArray == null || theoreticalYArray == null)
            {
                //when all mz of theoretical peaks or experimental peaks are less than mz cut off , it is treated as no corresponding library spectrum is found and later the similarity score will be assigned as null.
                return new List<(double, double)> { (-1, -1) };
            }

            List<(double, double)> intensityPairs = new();
            List<(double, double)> experimental = new();
            List<(double, double)> theoretical = new();

            for (int i = 0; i < ExperimentalXArray.Length; i++)
            {
                experimental.Add((ExperimentalXArray[i], experimentalYArray[i]));
            }
            for (int i = 0; i < TheoreticalXArray.Length; i++)
            {
                theoretical.Add((TheoreticalXArray[i], theoreticalYArray[i]));
            }

            experimental = experimental.OrderByDescending(i => i.Item2).ToList();
            theoretical = theoretical.OrderByDescending(i => i.Item2).ToList();

            foreach ((double, double) xyPair in theoretical)
            {
                int index = 0;
                while (experimental.Count > 0 && index < experimental.Count)
                {
                    if (Within(experimental[index].Item1, xyPair.Item1, localPpmTolerance))
                    {
                        intensityPairs.Add((experimental[index].Item2, xyPair.Item2));
                        experimental.RemoveAt(index);
                        index = -1;
                        break;
                    }
                    index++;
                }
                if (experimental.Count == 0)
                {
                    index++;
                }
                if (index > 0)
                {
                    //didn't find a experimental mz in range
                    intensityPairs.Add((0, xyPair.Item2));
                }
            }

            //If we're keeping all experimental and theoretical peaks, then we add intensity pairs for all unpaired experimental peaks here.
            if (experimental.Count > 0 && allPeaks)
            {
                foreach ((double, double) xyPair in experimental)
                {
                    intensityPairs.Add((xyPair.Item2, 0));
                }
            }
            return intensityPairs;
        }

        public bool DoubleWithinToleranceBool(double[] array, double value, double tolerance)
        {
            if (array.Length == 1 && !Within(array[0], value, tolerance))
            {
                return false;
            }
            else
            {
                int mid = array.Length / 2;
                if (Within(array[mid], value, tolerance))
                {
                    return true;
                }
                else
                {
                    if (value > array[mid] && value <= array.Last())
                    {
                        return DoubleWithinToleranceBool(array.TakeLast(array.Length - mid).ToArray(), value, tolerance);
                    }
                    else
                    {
                        return DoubleWithinToleranceBool(array.Take(mid).ToArray(), value, tolerance);
                    }
                }
            }
        }

        public double? DoubleWithinToleranceValue(double[] array, double value, double tolerance)
        {
            if (array.Length == 1 && !Within(array[0], value, tolerance))
            {
                return null;
            }
            else
            {
                int mid = array.Length / 2;
                if (Within(array[mid], value, tolerance))
                {
                    return array[mid];
                }
                else
                {
                    if (value > array[mid] && value <= array.Last())
                    {
                        return DoubleWithinToleranceValue(array.TakeLast(array.Length - mid).ToArray(), value, tolerance);
                    }
                    else
                    {
                        return DoubleWithinToleranceValue(array.Take(mid).ToArray(), value, tolerance);
                    }
                }
            }
        }

        #region normalization

        public static double[] NormalizeSquareRootSpectrumSum(double[] spectrum)
        {
            double sqrtSum = spectrum.Select(y => Math.Sqrt(y)).Sum();
            double[] normalizedSpectrum = new double[spectrum.Length];

            for (int i = 0; i < spectrum.Length; i++)
            {
                normalizedSpectrum[i] = Math.Sqrt(spectrum[i]) / sqrtSum;
            }
            return normalizedSpectrum;
        }

        public static double[] NormalizeMostAbundantPeak(double[] spectrum)
        {
            double max = spectrum.Max();
            double[] normalizedSpectrum = new double[spectrum.Length];

            for (int i = 0; i < spectrum.Length; i++)
            {
                normalizedSpectrum[i] = spectrum[i] / max;
            }
            return normalizedSpectrum;
        }

        public static double[] NormalizeSpectrumSum(double[] spectrum)
        {
            double sum = spectrum.Sum();
            double[] normalizedSpectrum = new double[spectrum.Length];

            for (int i = 0; i < spectrum.Length; i++)
            {
                normalizedSpectrum[i] = spectrum[i] / sum;
            }
            return normalizedSpectrum;
        }

        #endregion normalization

        public double? Xcorr(IEnumerable<double> arrayOfValues1, IEnumerable<double> arrayOfValues2)
        {
            if (arrayOfValues1.Count() != arrayOfValues2.Count())
            {
                return null;
            }
            else
            {
                return Correlation.Pearson(arrayOfValues1, arrayOfValues2);
            }
        }


        private bool Within(double mz1, double mz2, double localPpmTolerance)
        {
            return ((Math.Abs(mz1 - mz2) / Math.Max(mz1, mz2) * 1000000.0) < localPpmTolerance);
        }
        public enum SpectrumNormalizationScheme
        {
            squareRootSpectrumSum,
            spectrumSum,
            mostAbundantPeak,
            unnormalized
        }
    }
}
