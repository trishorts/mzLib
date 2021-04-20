using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;
using Chemistry;
using Proteomics.AminoAcidPolymer;
using MassSpectrometry;
using MzLibUtil;

namespace Test
{
    [TestFixture]
    public static class AverageineTest
    {
        [Test]
        public static void MyTest()
        {
            List<string> peptides = new List<string> { "PEPTIDE" };

            List<Protein> prots = ProteinDbLoader.LoadProteinFasta(@"C:\Users\Michael Shortreed\Downloads\MetaMorpheusVignette\FASTA\uniprot-mouse-reviewed-2-5-2018.fasta", true, DecoyType.None, false, out var a,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex);


            //SequenceSpecificAverageine j = new SequenceSpecificAverageine(prots.Select(s=>s.BaseSequence).ToList());

            List<string> sequences = prots.Select(s => s.BaseSequence).ToList();

            foreach (string sequence in sequences)
            {
                Peptide p = new Peptide(sequence);
                double monoIsotopicMass = p.MonoisotopicMass;

                int chargeFor1000mz = (int)Math.Round(monoIsotopicMass / 999.0, 0);

                double massResolution = monoIsotopicMass / 1000000.0 * 20.0;

                //var isotopeDistribution = IsotopicDistribution.GetDistribution(p.GetChemicalFormula(), 0.125, 1e-08, massResolution);
                double fineResolution = 0.125;
                double minProbability = 1e-06;
                var isotopeDistribution = IsotopicDistribution.GetDistribution(p.GetChemicalFormula(), fineResolution, minProbability);
                List<string> myout = new List<string>();
                List<double> masses = isotopeDistribution.Masses.ToList();
                List<double> intensities = isotopeDistribution.Intensities.ToList();
                for (int i = 0; i < isotopeDistribution.Masses.Count(); i++)
                {
                    myout.Add(Convert.ToString(masses[i].ToMz(1)) + "\t" + intensities[i]);
                }
                File.WriteAllLines(@"E:\junk\isotopeDist.txt", myout);

                if (chargeFor1000mz < 31)
                {
                    List<IsotopicEnvelope> isotopicEnvelopes = RetrieveIsotopeEnvelopes(sequence, chargeFor1000mz, massResolution, fineResolution, minProbability).Where(m=>m.TotalIntensity > 0.01).ToList();
                    if(isotopicEnvelopes.Count > 0)
                    {
                        List<double> mm = isotopicEnvelopes.Select(m => m.MonoisotopicMass).ToList();
                    }
                }
                
            }

            Assert.AreEqual(1, 0);
        }

        public static List<IsotopicEnvelope> RetrieveIsotopeEnvelopes(string sequence, int charge, double massResolution, double fineResolution, double minProbability)
        {
            Peptide p = new Peptide(sequence);
            double mim = p.MonoisotopicMass;
            var isotopeDistribution = IsotopicDistribution.GetDistribution(p.GetChemicalFormula(), fineResolution, minProbability);

            List<double> myMasses = isotopeDistribution.Masses.Select(m=>m.ToMz(charge + 1)).ToList();
            myMasses.AddRange(isotopeDistribution.Masses.Select(m => m.ToMz(charge)).ToList());
            myMasses.AddRange(isotopeDistribution.Masses.Select(m => m.ToMz(charge - 1)).ToList());
            List<double> myIntensities = isotopeDistribution.Intensities.ToList();
            myIntensities.AddRange(isotopeDistribution.Intensities);
            myIntensities.AddRange(isotopeDistribution.Intensities);

            double[] mz = myMasses.ToArray();
            double[] intensities = myIntensities.ToArray();
            bool shouldCopy = false;
            MzSpectrum massSpecrum = new MzSpectrum(mz, intensities, shouldCopy);
            List<string> myNewOut = new List<string>();
            for (int i = 0; i < mz.Count(); i++)
            {
                myNewOut.Add(Convert.ToString(mz[i]) + "\t" + intensities[i]);
            }
            File.WriteAllLines(@"E:\junk\spectrum.txt", myNewOut);
            double minMZ = mz.Min();
            double maxMZ = mz.Max();
            MzRange scanWindowRange = new MzRange(minMZ, maxMZ);
            double? isolationMZ = null;
            double? isolationWidth = null;
            MsDataScan j = new MsDataScan(massSpecrum, 1, 1, true, Polarity.Positive, 1, scanWindowRange, "", MZAnalyzerType.Orbitrap, 1, null, null, null, null, null, null, isolationMZ, isolationWidth, null, null, null, null);
            int minAssumedChargeState = Math.Max(1,charge -5);
            int maxAssumedChargeState = charge + 5;
            double deconvolutionTolerancePpm = 20;
            double intensityRatioLimit = 3;
            return j.MassSpectrum.Deconvolute(new MzRange(0, double.PositiveInfinity), minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatioLimit).ToList();
        }
    }
}
