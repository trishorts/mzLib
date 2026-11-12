using Chemistry;
using CsvHelper.Configuration.Attributes;

namespace Readers
{
    /// <summary>
    /// One peak-group candidate of one precursor in one run, from Spectronaut's TRIC export: the OpenSWATH feature
    /// table that TRIC (<c>feature_alignment.py</c>) reads, written by Spectronaut. Spectronaut writes up to four
    /// candidates per run and precursor, ranked by <see cref="Rank"/>. Its reader keeps rank 1.
    /// <para>
    /// This is the one Spectronaut export seen in deposits that carries everything FlashLFQ needs: retention time
    /// (minutes), precursor m/z, charge, a decoy flag and a q-value (<c>m_score</c>).
    /// </para>
    /// </summary>
    public class SpectronautTricPrecursor : IQuantifiableRecord
    {
        [Name("transition_group_id")] public string TransitionGroupId { get; set; } = "";
        [Name("peptide_group_label")] [Optional] public string? PeptideGroupLabel { get; set; }

        /// <summary>The modified sequence and charge, as EG.PrecursorId writes it: <c>_YTGVPDR_.2</c>.</summary>
        [Name("PrecursorID")] public string PrecursorId { get; set; } = "";

        [Name("FullPeptideName")] [Optional] public string? FullPeptideName { get; set; }
        [Name("Sequence")] [Optional] public string? Sequence { get; set; }

        [Name("m/z")] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? PrecursorMz { get; set; }

        [Name("Charge")] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? Charge { get; set; }

        /// <summary>The group's protein ids, <c>;</c>-separated.</summary>
        [Name("ProteinName")] [Optional] public string? ProteinName { get; set; }

        [Name("ProteinGroupQvalue")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? ProteinGroupQValue { get; set; }

        [Name("ProteinGroupCScore")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? ProteinGroupCScore { get; set; }

        /// <summary>A 0-based run index.</summary>
        [Name("run_id")] public string RunId { get; set; } = "";

        /// <summary>The run's spectra file name, with its extension (<c>22_0505_002_SSL2_ND_Non_23_01.wiff</c>).</summary>
        [Name("filename")] public string SpectraFileName { get; set; } = "";

        [Name("Condition")] [Optional] public string? ConditionText { get; set; }

        [Name("Replicate")] [Optional] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? Replicate { get; set; }

        [Name("decoy")] [TypeConverter(typeof(SpectronautNullableBoolConverter))]
        public bool? Decoy { get; set; }

        [Name("usedForPeptideQuantity")] [Optional] [TypeConverter(typeof(SpectronautNullableBoolConverter))]
        public bool? UsedForPeptideQuantity { get; set; }

        [Name("usedForProteinQuantity")] [Optional] [TypeConverter(typeof(SpectronautNullableBoolConverter))]
        public bool? UsedForProteinQuantity { get; set; }

        /// <summary>The library (assay) retention time, in minutes.</summary>
        [Name("assay_rt")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? AssayRetentionTime { get; set; }

        [Name("id")] [Optional] public string? Id { get; set; }

        /// <summary>The candidate's rank among this precursor's peak groups in this run; 1 is the chosen peak.</summary>
        [Name("rank")] [TypeConverter(typeof(SpectronautNullableIntConverter))]
        public int? Rank { get; set; }

        [Name("isImputed")] [Optional] [TypeConverter(typeof(SpectronautNullableBoolConverter))]
        public bool? IsImputed { get; set; }

        [Name("delta_rt")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? DeltaRetentionTime { get; set; }

        /// <summary>The peak apex retention time, in minutes.</summary>
        [Name("RT")] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? ApexRetentionTime { get; set; }

        [Name("StartRT")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? StartRetentionTime { get; set; }

        [Name("EndRT")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? EndRetentionTime { get; set; }

        /// <summary>The discriminant score; higher is better.</summary>
        [Name("d_score")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? DScore { get; set; }

        /// <summary>The q-value of the peak group. Not filtered at 1%: values above 0.1 occur.</summary>
        [Name("m_score")] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? MScore { get; set; }

        [Name("Intensity")] [Optional] [TypeConverter(typeof(SpectronautNullableDoubleConverter))]
        public double? Intensity { get; set; }

        #region IQuantifiableRecord

        /// <summary>The spectra file name without directory or extension.</summary>
        [Ignore] public string FileName => _fileName ??= RunFileMatcher.StripPathAndExtension(SpectraFileName);
        [Ignore] private string? _fileName;

        [Ignore] public string FullSequence => FullPeptideName ?? SpectronautReports.SplitPrecursorId(PrecursorId, out _);

        [Ignore] public string BaseSequence => Sequence ?? SpectronautReports.StripModifications(FullSequence);

        /// <summary>The peak apex, in minutes.</summary>
        [Ignore] public double RetentionTime => ApexRetentionTime ?? double.NaN;

        [Ignore] public int ChargeState => Charge ?? 0;

        [Ignore] public bool IsDecoy => Decoy ?? false;

        [Ignore]
        public double MonoisotopicMass => PrecursorMz is double mz && Charge is int z && z != 0 ? mz.ToMass(z) : double.NaN;

        [Ignore] private List<(string, string, string)>? _proteinGroupInfos;

        /// <summary>One tuple per protein in ProteinName. The TRIC export reports no genes or organisms.</summary>
        [Ignore]
        public List<(string proteinAccessions, string geneName, string organism)> ProteinGroupInfos =>
            _proteinGroupInfos ??= SpectronautReports.Members(ProteinName) is { Count: > 0 } members
                ? members.Select(member => (member, "", "")).ToList()
                : new List<(string, string, string)> { ("", "", "") };

        #endregion
    }
}
