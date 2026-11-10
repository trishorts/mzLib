using CsvHelper.Configuration.Attributes;

namespace Readers
{
    /// <summary>
    /// One PTM site in one run, from a Spectronaut PTM site report (<c>&lt;analysis&gt;_-PTMSiteReport.tsv</c>), which
    /// Spectronaut writes for analyses run with its PTM workflow. Protein-group columns repeat on every site.
    /// </summary>
    public class SpectronautPtmSite
    {
        [Name("R.FileName")] public string RunName { get; set; } = "";
        [Name("R.Condition")] [Optional] public string? ConditionText { get; set; }

        [Name("R.Replicate")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? Replicate { get; set; }

        [Name("PG.ProteinGroups")] [Optional] public string? ProteinGroups { get; set; }
        [Name("PG.ProteinNames")] [Optional] public string? ProteinNames { get; set; }
        [Name("PG.Organisms")] [Optional] public string? Organisms { get; set; }
        [Name("PG.Genes")] [Optional] public string? Genes { get; set; }

        [Name("PG.Qvalue")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? ProteinGroupQValue { get; set; }

        [Name("PG.QValue (Run-Wise)")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? ProteinGroupRunWiseQValue { get; set; }

        /// <summary>The protein the site is located on.</summary>
        [Name("PTM.ProteinId")] [Optional] public string? ProteinId { get; set; }

        [Name("PTM.CollapseKey")] [Optional] public string? CollapseKey { get; set; }

        /// <summary>The precursors collapsed into this site (<c>_DVKPHC[Carbamidomethyl (C)]VLLASK_.3</c>).</summary>
        [Name("PTM.Group")] [Optional] public string? Group { get; set; }

        /// <summary>Spectronaut's modification name (<c>Carbamidomethyl (C)</c>), as written.</summary>
        [Name("PTM.ModificationTitle")] [Optional] public string? ModificationTitle { get; set; }

        [Name("PTM.Multiplicity")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? Multiplicity { get; set; }

        [Name("PTM.NrOfCollapsedPeptides")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? CollapsedPeptides { get; set; }

        /// <summary>The 1-based position of the site on <see cref="ProteinId"/>.</summary>
        [Name("PTM.SiteLocation")] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? SiteLocation { get; set; }

        [Name("PTM.SiteAA")] [Optional] public string? SiteAminoAcid { get; set; }

        [Name("PTM.SiteProbability")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? SiteProbability { get; set; }

        [Name("PTM.FlankingRegion")] [Optional] public string? FlankingRegion { get; set; }

        [Name("PTM.Quantity")] public string QuantityText { get; set; } = "";

        /// <summary>The run's condition, or null for Spectronaut's <c>Not Defined</c>.</summary>
        [Ignore] public string? Condition => SpectronautReports.Stated(ConditionText);

        [Ignore] public double? Quantity => SpectronautReports.ParseDouble(QuantityText);

        [Ignore] public SpectronautValueState QuantityState => SpectronautReports.StateOf(QuantityText);
    }
}
