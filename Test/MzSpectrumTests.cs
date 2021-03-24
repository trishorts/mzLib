using MassSpectrometry;
using NUnit.Framework;
using System;
using System.IO;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class MzSpectrumTests
    {
        [Test]
        public static void MyTest()
        {
            string spectrumFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\PEPTIDE_MS1_mz_and_intensities.txt");
            MzSpectrum mzs = ReadMzAndIntensity(spectrumFilePath);
            var isotopeEnvelopes = mzs.Deconvolute(mzs.Range, 1, 1, 5, 100).ToList();
            Assert.AreEqual(1, 0);
        }

        public static MzSpectrum ReadMzAndIntensity(string file)
        {
            string[] mzIntensityInput = File.ReadAllLines(file);
            double[] mzValues = new double[mzIntensityInput.Length];
            double[] intensities = new double[mzIntensityInput.Length];
            for (int i = 0; i < mzIntensityInput.Length; i++)
            {
                string[] pair = mzIntensityInput[i].Split('\t');
                mzValues[i] = Convert.ToDouble(pair[0]);
                intensities[i] = Convert.ToDouble(pair[1]);
            }
            return new MzSpectrum(mzValues, intensities, true);
        }
    }
}