using Chemistry;
using Proteomics.AminoAcidPolymer;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    public class Averagine
    {
        public ChemicalFormula GlobalChemicalFormula { get; private set; }
        public Dictionary<string, double> SequenceSpecificAveragine { get; private set; }

        public Averagine(List<string> sequences)
        {
            ComputeGlobalChemicalFormula(sequences);
        }

        private void ComputeGlobalChemicalFormula(List<string> sequences)
        {
            GlobalChemicalFormula = new ChemicalFormula();
            foreach (string sequence in sequences)
            {
                if(!sequence.Contains("X") && !sequence.Contains("B"))
                {
                    Peptide baseSequence = new Peptide(sequence);
                    GlobalChemicalFormula.Add(baseSequence.GetChemicalFormula());
                }
                
            }
        }
    }
}