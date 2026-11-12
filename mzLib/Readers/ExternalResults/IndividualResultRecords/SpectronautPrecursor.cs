using Chemistry;
using CsvHelper.Configuration.Attributes;

namespace Readers
{
    /// <summary>
    /// One precursor in one run, from a Spectronaut Normal report whose deepest level is the precursor (elution
    /// group / fragment group) or the fragment. A fragment-level report has one row per fragment ion. Its reader keeps
    /// the first row of each run and precursor, and this record holds only the precursor-level columns.
    /// <para>
    /// Deposited precursor reports are often thin. MassIVE's <c>_Report_peptides.tsv</c> exports carry only
    /// EG.PrecursorId (<c>_PEPTIDE_.2</c>), EG.Qvalue and a quantity: no retention time, no m/z, no charge column. The
    /// charge is then read from the id's <c>.z</c> suffix, and the record cannot feed FlashLFQ.
    /// <see cref="SpectronautPrecursorReportFile.GetQuantifiableResults"/> says which columns are missing.
    /// </para>
    /// <para>
    /// Modified sequences are kept as Spectronaut wrote them, with no modification lookup. Spectronaut 10 wrote mass
    /// deltas (<c>[+42]</c>, <c>[+16]</c>) in EG.ModifiedSequence and names, or a formula (<c>[+O]</c>), in FG.Id. Newer
    /// versions write names (<c>[Oxidation (M)]</c>), sometimes without a site list (<c>[Biotin]</c>).
    /// </para>
    /// </summary>
    public class SpectronautPrecursor : IQuantifiableRecord
    {
        [Name("R.FileName")] public string RunName { get; set; } = "";
        [Name("R.Condition")] [Optional] public string? ConditionText { get; set; }

        [Name("R.Replicate")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? Replicate { get; set; }

        [Name("PG.ProteinGroups")] [Optional] public string? ProteinGroups { get; set; }
        [Name("PG.ProteinAccessions")] [Optional] public string? ProteinAccessions { get; set; }
        [Name("PG.Genes")] [Optional] public string? Genes { get; set; }
        [Name("PG.ProteinNames")] [Optional] public string? ProteinNames { get; set; }
        [Name("PG.ProteinDescriptions")] [Optional] public string? ProteinDescriptions { get; set; }

        [Name("PG.Qvalue")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? ProteinGroupQValue { get; set; }

        [Name("PG.Quantity")] [Optional] public string? ProteinGroupQuantityText { get; set; }

        /// <summary>Spectronaut 10 named this column EG.StrippedSequence.</summary>
        [Name("PEP.StrippedSequence", "EG.StrippedSequence")] [Optional] public string? StrippedSequence { get; set; }

        /// <summary>The modified sequence and charge: <c>_PEPTIDE_.2</c>.</summary>
        [Name("EG.PrecursorId")] [Optional] public string? PrecursorId { get; set; }

        [Name("EG.ModifiedSequence")] [Optional] public string? ModifiedSequenceText { get; set; }
        [Name("EG.ModifiedPeptide")] [Optional] public string? ModifiedPeptide { get; set; }

        /// <summary>The fragment group id. Label-free it is the precursor id; Spectronaut 10 wrote modification names here.</summary>
        [Name("FG.Id")] [Optional] public string? FragmentGroupId { get; set; }

        [Name("FG.Charge")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? Charge { get; set; }

        [Name("FG.PrecMz")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? PrecursorMz { get; set; }

        /// <summary>The precursor q-value. It is per run, and deposited reports are often not filtered on it.</summary>
        [Name("EG.Qvalue")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? QValue { get; set; }

        [Name("EG.PEP")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? PosteriorErrorProbability { get; set; }

        [Name("EG.Cscore")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? CScore { get; set; }

        [Name("EG.IsDecoy")] [Optional] [TypeConverter(typeof(SpectronautNullableBoolConverter))]
        public bool? Decoy { get; set; }

        /// <summary>Retention time in minutes: EG.ApexRT, or EG.MeanApexRT / EG.RTEmpirical where a version wrote those instead.</summary>
        [Name("EG.ApexRT", "EG.MeanApexRT", "EG.RTEmpirical")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? ApexRetentionTime { get; set; }

        [Name("FG.Quantity")] [Optional] public string? QuantityText { get; set; }
        [Name("EG.TotalQuantity (Settings)")] [Optional] public string? TotalQuantityText { get; set; }

        #region Derived

        /// <summary>The run's condition, or null for Spectronaut's <c>Not Defined</c>.</summary>
        [Ignore] public string? Condition => SpectronautReports.Stated(ConditionText);

        /// <summary>The precursor's identity within a run: EG.PrecursorId, else FG.Id, else the modified sequence and charge.</summary>
        [Ignore] public string PrecursorKey => PrecursorId ?? FragmentGroupId ?? $"{ModifiedSequenceText}.{Charge}";

        /// <summary>
        /// The modified sequence of this precursor: EG.PrecursorId, or else FG.Id, without its charge suffix. Falls back to
        /// EG.ModifiedSequence when neither id was exported. The ids come first because they distinguish label channels
        /// that EG.ModifiedSequence merges. In a Spectronaut 10 stoichiometry export (MSV000087085), four precursors of
        /// m/z 607.39 to 610.41 (light and heavy acetyl, FG.Id <c>_K[Acetyl (K)]K[Acetyl D3 (K)]...</c>) all have the
        /// EG.ModifiedSequence <c>_K[+42]K[+42]LLAK[+42]DLR_</c>.
        /// </summary>
        [Ignore]
        public string ModifiedSequence => (PrecursorId ?? FragmentGroupId) is string id
            ? SpectronautReports.SplitPrecursorId(id, out _)
            : ModifiedSequenceText ?? "";

        /// <summary>FG.Charge, or the charge suffix of EG.PrecursorId or FG.Id when FG.Charge was not exported.</summary>
        [Ignore]
        public int? PrecursorCharge
        {
            get
            {
                if (Charge.HasValue) return Charge;
                SpectronautReports.SplitPrecursorId(PrecursorId ?? FragmentGroupId ?? "", out int? charge);
                return charge;
            }
        }

        /// <summary>FG.Quantity, or EG.TotalQuantity (Settings) when FG.Quantity was not exported; null when neither holds a number.</summary>
        [Ignore] public double? Quantity => SpectronautReports.ParseDouble(QuantityText ?? TotalQuantityText);

        [Ignore] public SpectronautValueState QuantityState => SpectronautReports.StateOf(QuantityText ?? TotalQuantityText);

        #endregion

        #region IQuantifiableRecord

        [Ignore] public string FileName => RunName;

        [Ignore] public string FullSequence => ModifiedSequence;

        /// <summary>The stripped sequence column, or the modified sequence with its modifications and anchors removed.</summary>
        [Ignore] public string BaseSequence => StrippedSequence ?? SpectronautReports.StripModifications(ModifiedSequence);

        /// <summary>Minutes; NaN when the report has no retention-time column (the file refuses to hand such records to FlashLFQ).</summary>
        [Ignore] public double RetentionTime => ApexRetentionTime ?? double.NaN;

        [Ignore] public int ChargeState => PrecursorCharge ?? 0;

        [Ignore] public bool IsDecoy => Decoy ?? false;

        /// <summary>From FG.PrecMz and the charge; NaN when either is absent (the file refuses to hand such records to FlashLFQ).</summary>
        [Ignore]
        public double MonoisotopicMass => PrecursorMz is double mz && PrecursorCharge is int z && z != 0 ? mz.ToMass(z) : double.NaN;

        [Ignore] private List<(string, string, string)>? _proteinGroupInfos;

        /// <summary>
        /// One tuple per distinct member of the group. A gene is assigned only when the group has one member and one gene,
        /// as for DIA-NN: a gene list longer than one cannot be tied to a member beyond doubt. Organism is not reported per
        /// member.
        /// </summary>
        [Ignore]
        public List<(string proteinAccessions, string geneName, string organism)> ProteinGroupInfos =>
            _proteinGroupInfos ??= BuildProteinGroupInfos();

        private List<(string, string, string)> BuildProteinGroupInfos()
        {
            var members = SpectronautReports.Members(ProteinGroups ?? ProteinAccessions);
            var genes = SpectronautReports.Members(Genes);
            string gene = members.Count == 1 && genes.Count == 1 ? genes[0] : "";
            if (members.Count == 0)
                return new List<(string, string, string)> { ("", gene, "") };
            return members.Select(member => (member, gene, "")).ToList();
        }

        #endregion
    }
}
