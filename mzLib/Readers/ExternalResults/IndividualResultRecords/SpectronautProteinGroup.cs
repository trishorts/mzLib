using CsvHelper.Configuration.Attributes;

namespace Readers
{
    /// <summary>
    /// One row of a Spectronaut protein-group report (a Normal report with only R.* and PG.* columns): one protein
    /// group in one run. Spectronaut's columns are configurable, so every column but the core (R.FileName, PG.Quantity,
    /// and PG.ProteinGroups or PG.ProteinAccessions) is optional.
    /// <para>
    /// The PG.* list columns are kept as written. They are lists over the group's members, but Spectronaut collapses a
    /// list to one value when every member shares it, keeps duplicate ids in PG.ProteinGroups while removing them from
    /// PG.ProteinAccessions, and leaves empty entries. Use <see cref="PerMember"/> to read one value per member; it
    /// refuses a list it cannot attribute rather than pairing values with the wrong protein.
    /// </para>
    /// </summary>
    public class SpectronautProteinGroup
    {
        /// <summary>The run, without directory or extension.</summary>
        [Name("R.FileName")] public string RunName { get; set; } = "";

        /// <summary>The run's condition, as written. <see cref="Condition"/> gives null for Spectronaut's <c>Not Defined</c>.</summary>
        [Name("R.Condition")] [Optional] public string? ConditionText { get; set; }

        [Name("R.Replicate")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? Replicate { get; set; }

        [Name("R.Fraction")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? Fraction { get; set; }

        /// <summary>The group's protein ids: <c>;</c> within a group, <c>|</c> between groups.</summary>
        [Name("PG.ProteinGroups")] [Optional] public string? ProteinGroups { get; set; }

        [Name("PG.ProteinAccessions")] [Optional] public string? ProteinAccessions { get; set; }
        [Name("PG.Genes")] [Optional] public string? Genes { get; set; }
        [Name("PG.ProteinNames")] [Optional] public string? ProteinNames { get; set; }
        [Name("PG.ProteinDescriptions")] [Optional] public string? ProteinDescriptions { get; set; }
        [Name("PG.UniProtIds")] [Optional] public string? UniProtIds { get; set; }
        [Name("PG.Organisms")] [Optional] public string? Organisms { get; set; }
        [Name("PG.OrganismId")] [Optional] public string? OrganismId { get; set; }

        /// <summary>UniProt protein-existence levels, one per member when not collapsed.</summary>
        [Name("PG.Protein Existence")] [Optional] public string? ProteinExistence { get; set; }

        /// <summary>The experiment-wide protein-group q-value.</summary>
        [Name("PG.Qvalue")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? QValue { get; set; }

        [Name("PG.QValue (Run-Wise)")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? RunWiseQValue { get; set; }

        /// <summary>Sequence coverage as written, with its percent sign (<c>22.8%</c>).</summary>
        [Name("PG.Coverage")] [Optional] public string? CoverageText { get; set; }

        [Name("PG.IsSingleHit")] [Optional] [TypeConverter(typeof(SpectronautNullableBoolConverter))]
        public bool? IsSingleHit { get; set; }

        [Name("PG.NrOfPrecursorsUsedForQuantification")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? PrecursorsUsedForQuantification { get; set; }

        [Name("PG.NrOfStrippedSequencesUsedForQuantification")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? StrippedSequencesUsedForQuantification { get; set; }

        /// <summary>
        /// The quantity as written. Spectronaut's PG.Quantity is the group quantity its settings define: usually rolled up
        /// from peptides and normalized across runs, not a raw intensity.
        /// </summary>
        [Name("PG.Quantity")] public string QuantityText { get; set; } = "";

        /// <summary>The run's condition, or null when Spectronaut wrote none (<c>Not Defined</c>, <c>NA</c>, blank).</summary>
        [Ignore] public string? Condition => SpectronautReports.Stated(ConditionText);

        /// <summary>The quantity, or null when the cell is <c>NaN</c> or <c>Filtered</c>; <see cref="QuantityState"/> says which.</summary>
        [Ignore] public double? Quantity => SpectronautReports.ParseDouble(QuantityText);

        [Ignore] public SpectronautValueState QuantityState => SpectronautReports.StateOf(QuantityText);

        /// <summary>Sequence coverage as a fraction (22.8% gives 0.228), or null.</summary>
        [Ignore] public double? Coverage => SpectronautReports.ParsePercent(CoverageText);

        /// <summary>The group's distinct proteins, in order, from PG.ProteinGroups (or PG.ProteinAccessions when that is absent).</summary>
        [Ignore] public IReadOnlyList<string> Members => _members ??= SpectronautReports.Members(ProteinGroups ?? ProteinAccessions);
        [Ignore] private IReadOnlyList<string>? _members;

        /// <summary>
        /// One value per member of <see cref="Members"/> from one of this row's PG.* list columns (for example
        /// <see cref="Genes"/>), or null when the list cannot be attributed: absent, or a length other than the member
        /// count or one. An empty entry is null.
        /// </summary>
        public IReadOnlyList<string?>? PerMember(string? list) => SpectronautReports.PerMember(list, Members.Count);
    }
}
