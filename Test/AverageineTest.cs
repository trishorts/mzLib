using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;

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


            Averagine j = new Averagine(prots.Select(s=>s.BaseSequence).ToList());

            Chemistry.ChemicalFormula cf = j.GlobalChemicalFormula;

            Assert.AreEqual(1, 0);
        }
    }
}
