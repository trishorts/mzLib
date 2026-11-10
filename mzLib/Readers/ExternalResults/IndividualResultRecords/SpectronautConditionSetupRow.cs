using CsvHelper.Configuration.Attributes;

namespace Readers
{
    /// <summary>
    /// One run of a Spectronaut condition setup (<c>&lt;analysis&gt;_ConditionSetup.tsv</c>): the run, its condition,
    /// replicate and fraction.
    /// <para>
    /// The Label and Color columns are not read. Label repeats the condition or cuts it short (MSV000090777 writes
    /// <c>Non-i...</c>), and Color is display only. <see cref="SampleEvidenceExtractor"/> reads the same table for SDRF
    /// evidence under the same rules.
    /// </para>
    /// <para>
    /// <see cref="IsReferenceCondition"/> is Spectronaut's "Reference" column. It marks the reference condition that
    /// ratios are taken against. It is NOT a pooled reference sample, which is what <c>SpectraFileInfo.IsReference</c>
    /// means, and is never mapped to it.
    /// </para>
    /// </summary>
    public class SpectronautConditionSetupRow
    {
        [Name("#")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? Number { get; set; }

        /// <summary>The ratio reference condition flag. Not a pooled-reference sample flag.</summary>
        [Name("Reference")] [Optional] [TypeConverter(typeof(SpectronautNullableBoolConverter))]
        public bool? IsReferenceCondition { get; set; }

        /// <summary>The run as Spectronaut lists it, with its extension (<c>240418_AM_O1.raw</c>).</summary>
        [Name("Run Label")] public string RunLabel { get; set; } = "";

        [Name("Condition")] public string ConditionText { get; set; } = "";

        [Name("Fraction")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? Fraction { get; set; }

        /// <summary>The biological replicate, 1-based as Spectronaut writes it.</summary>
        [Name("Replicate")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? Replicate { get; set; }

        [Name("Quantity Correction Factor")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? QuantityCorrectionFactor { get; set; }

        /// <summary>The run's file name without extension, as the report's R.FileName gives it.</summary>
        [Name("File Name")] [Optional] public string? FileName { get; set; }

        /// <summary>The condition, or null for Spectronaut's <c>Not Defined</c>, <c>NA</c> or a blank.</summary>
        [Ignore] public string? Condition => SpectronautReports.Stated(ConditionText);

        /// <summary>The run's stem: Run Label without directory or extension.</summary>
        [Ignore] public string RunStem => RunFileMatcher.StripPathAndExtension(RunLabel);
    }
}
