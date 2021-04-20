using Chemistry;
using Proteomics.AminoAcidPolymer;
using System.Collections.Generic;

namespace Proteomics
{
    public class SequenceSpecificAverageine
    {
        private ChemicalFormula GlobalChemicalFormula { get; set; }
        private int TotalSequenceLength;
        public Dictionary<string, double> OneAA_Averagine { get; private set; }
        

        public SequenceSpecificAverageine(List<string> sequences)
        {
            ComputeGlobalChemicalFormula(sequences);
            ComputeSequenceSpecificAveragine();
        }

        private void ComputeSequenceSpecificAveragine()
        {
            OneAA_Averagine = new Dictionary<string, double>();
            foreach (Element element in GlobalChemicalFormula.Elements.Keys)
            {
                OneAA_Averagine.Add(element.ToString(), (double)GlobalChemicalFormula.Elements[element] / (double)TotalSequenceLength);
            }
        }

        private void ComputeGlobalChemicalFormula(List<string> sequences)
        {
            GlobalChemicalFormula = new ChemicalFormula();
            foreach (string sequence in sequences)
            {
                string sequenceWithNoUnusualAAs = sequence.Replace("X", "").Replace("B", "").Replace("U", "");
                TotalSequenceLength += sequenceWithNoUnusualAAs.Length;
                GlobalChemicalFormula.Add(new Peptide(sequenceWithNoUnusualAAs).GetChemicalFormula());
            }
        }
    }
}