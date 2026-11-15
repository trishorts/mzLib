using System;
using System.Collections.Generic;
using System.Linq;
using StatisticalModels;

namespace Quantification.Differential;

/// <summary>
/// The share models (QuantProject <c>DEF-DIFF-SHARE-MODEL</c>, ruling S1): how much more, or less, a member's marker
/// peptides change than its reference peptides, in log2. The reference peptides are in every form of the gene, so it is
/// the change in log2 of the named form's share of the gene (for an <c>exon_peptide</c> row, of exon inclusion). Each
/// sample's marker values are compared with its own reference values, so the sample's shift, its normalization and each
/// peptide's own digestion efficiency cancel from it.
/// </summary>
public static class ShareModel
{
    /// <summary>A single's <c>model_used</c>.</summary>
    public const string ShareMixedModel = "share_mixed_model";
    /// <summary>A share set's <c>model_used</c>.</summary>
    public const string ShareSetMixedModel = "share_set_mixed_model";

    /// <summary>
    /// A member alone (a single):
    /// <c>log2 value ~ peptide + (the stratum's fixed design) + marker:(the stratum's fixed design) + random terms</c>, over
    /// its marker and reference peptides. The random terms are STAT1's protein model's for the stratum (GR-6): one intercept
    /// per sample, or per individual when individuals are given. The share effect is the contrast on the <c>marker:</c>
    /// terms.
    /// </summary>
    /// <param name="setup">The contrast's setup; must be run.</param>
    /// <param name="member">The member, <c>tested</c> (<see cref="ShareMemberStates"/>), with its peptides classed on <paramref name="normalized"/>.</param>
    /// <param name="normalized">The stratum's normalized observation table.</param>
    /// <param name="notRun">Why no estimate came out, when it is null.</param>
    /// <param name="robust">STAT1's robust setting for the stratum (GR-8), as <see cref="SetModel.FitSet"/> takes it.</param>
    public static SetEstimate? FitSingle(SetContrastSetup setup, ShareMemberResult member, ObservationTable normalized,
        out string? notRun, bool robust = true)
    {
        ArgumentNullException.ThrowIfNull(member);
        return Fit(setup, new[] { member }, normalized, set: false, out notRun, robust);
    }

    /// <summary>
    /// A share set, one model extending <c>DEF-DIFF-SET-MODEL</c>'s:
    /// <c>log2 value ~ peptide + member:(design) + member:marker:(design) + (1|sample) + (1|sample:member)</c>. β̄ is the
    /// equal-weight average of the members' <c>marker:</c> contrasts; its se comes from the fit and its df is Satterthwaite's.
    /// </summary>
    /// <param name="setup">The contrast's setup; must be run.</param>
    /// <param name="members">The tested members, each once, in order.</param>
    /// <param name="normalized">The stratum's normalized observation table.</param>
    /// <param name="notRun">Why no estimate came out, when it is null.</param>
    /// <param name="robust">STAT1's robust setting for the stratum (GR-8).</param>
    public static SetEstimate? FitSet(SetContrastSetup setup, IReadOnlyList<ShareMemberResult> members, ObservationTable normalized,
        out string? notRun, bool robust = true)
    {
        ArgumentNullException.ThrowIfNull(members);
        if (members.Select(m => m.Assignment.Member).Distinct(StringComparer.Ordinal).Count() != members.Count)
            throw new ArgumentException("A member enters a share set at most once.", nameof(members));
        return Fit(setup, members, normalized, set: true, out notRun, robust);
    }

    private static SetEstimate? Fit(SetContrastSetup setup, IReadOnlyList<ShareMemberResult> members, ObservationTable normalized,
        bool set, out string? notRun, bool robust)
    {
        ArgumentNullException.ThrowIfNull(setup);
        ArgumentNullException.ThrowIfNull(normalized);
        if (setup.NotRun is not null || setup.Design is null || setup.Contrast is null)
            throw new ArgumentException("The setup is not run, so there is no model to fit.", nameof(setup));
        if (members.Any(m => !m.InEstimate || m.Peptides is null))
            throw new ArgumentException("Only tested members, with their peptides classed, enter a share model.", nameof(members));

        var design = setup.Design;
        double[,] x = design.FittedMatrix;
        double[] w = design.FittedWeights(setup.Contrast);
        int p = design.FittedColumns.Count, m = members.Count;
        // marker:intercept is the sum of the marker peptides' own terms, so it is aliased and left out; a level contrast
        // never weighs the intercept, so the share effect does not need it.
        if (design.FittedColumns[0] != 0 || !string.Equals(design.ColumnNames[0], "intercept", StringComparison.Ordinal) || w[0] != 0)
            throw new ArgumentException("A share model needs a design whose first column is the intercept and a contrast that does not weigh it.", nameof(setup));
        int width = 2 * p - 1; // design, then marker:(design without the intercept)

        var sampleRow = design.Samples.Select((s, i) => (s.SampleId, i)).ToDictionary(t => t.SampleId, t => t.i, StringComparer.Ordinal);
        var tableColumn = normalized.Samples.Select((s, i) => (s, i)).ToDictionary(t => t.s, t => t.i, StringComparer.Ordinal);
        var rowOf = new Dictionary<ObservedPeptide, int>(ReferenceEqualityComparer.Instance);
        for (int k = 0; k < normalized.Peptides.Count; k++) rowOf[normalized.Peptides[k]] = k;

        var response = new List<double>();
        var commonRows = new List<double[]>();
        var nuisanceIndex = new List<int>(); // -1: the member's first peptide (no column)
        var sampleOf = new List<string>();
        var individualOf = new List<string?>();
        var memberOf = new List<string>();
        var peptideOf = new List<string>();
        int nuisanceColumns = 0, peptides = 0;
        for (int g = 0; g < m; g++)
        {
            var classes = members[g].Peptides!;
            var used = classes.Markers.Select(pp => (Peptide: pp, Marker: true)).Concat(classes.References.Select(pp => (Peptide: pp, Marker: false)))
                .Select(t => (t.Peptide, t.Marker, Row: rowOf.TryGetValue(t.Peptide, out int r) ? r
                    : throw new ArgumentException($"{t.Peptide.FullSequence} is not a peptide of the table.", nameof(members))))
                .Where(t => design.Samples.Any(s => tableColumn.TryGetValue(s.SampleId, out int c) && double.IsFinite(normalized.Log2Intensity(t.Row, c))))
                .OrderBy(t => t.Peptide.FullSequence, StringComparer.Ordinal)
                .ToList();
            for (int u = 0; u < used.Count; u++)
            {
                int column = u == 0 ? -1 : nuisanceColumns++;
                peptides++;
                foreach (var sample in design.Samples)
                {
                    if (!tableColumn.TryGetValue(sample.SampleId, out int c))
                        throw new ArgumentException($"The table has no sample '{sample.SampleId}'.", nameof(normalized));
                    double value = normalized.Log2Intensity(used[u].Row, c);
                    if (!double.IsFinite(value)) continue;
                    var row = new double[m * width];
                    int r = sampleRow[sample.SampleId];
                    for (int j = 0; j < p; j++) row[g * width + j] = x[r, j];
                    if (used[u].Marker)
                        for (int j = 1; j < p; j++) row[g * width + p + j - 1] = x[r, j];
                    response.Add(value);
                    commonRows.Add(row);
                    nuisanceIndex.Add(column);
                    sampleOf.Add(sample.SampleId);
                    individualOf.Add(sample.Individual);
                    memberOf.Add(members[g].Assignment.Member);
                    peptideOf.Add(used[u].Peptide.FullSequence);
                }
            }
        }

        int n = response.Count;
        var common = new double[n, m * width];
        var nuisance = nuisanceColumns > 0 ? new double[n, nuisanceColumns] : null;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m * width; j++) common[i, j] = commonRows[i][j];
            if (nuisanceIndex[i] >= 0) nuisance![i, nuisanceIndex[i]] = 1;
        }
        var fitted = design.FittedColumns.Select(j => design.ColumnNames[j]).ToList();
        var names = members.SelectMany(t => fitted.Select(c => $"member={t.Assignment.Member}:{c}")
            .Concat(fitted.Skip(1).Select(c => $"member={t.Assignment.Member}:marker:{c}"))).ToList();
        // A single's random terms are STAT1's protein model's: the sample, or the individual when given. A set's are the set
        // model's: (1|sample) shared by its members and (1|sample:member).
        var problem = set
            ? new NestedMixedModelProblem(response, common, nuisance, memberOf, sampleOf)
            : new NestedMixedModelProblem(response, common, nuisance, sampleOf, individualOf.Any(v => v is not null) ? individualOf : null);
        var fit = NestedMixedModel.Fit(new[] { problem }, names, maxThreads: 1, robust: robust);
        if (fit.Status[0] != FeatureFitStatus.Fitted)
        {
            notRun = fit.Status[0] == FeatureFitStatus.NotConverged
                ? SetNotRunReason.NotConverged
                : SetNotRunReason.NotConverged + ":" + SnakeCase(fit.Status[0].ToString());
            return null;
        }

        // The share effect: each member's contrast on its marker: terms, weighted 1/m.
        var weights = new double[m * width];
        for (int g = 0; g < m; g++)
            for (int j = 1; j < p; j++)
                weights[g * width + p + j - 1] = w[j] / m;
        double estimate = 0, quad = 0;
        for (int i = 0; i < m * width; i++)
        {
            estimate += weights[i] * fit.Coefficient(0, i);
            for (int j = 0; j < m * width; j++) quad += weights[i] * fit.UnscaledCovariance(0, i, j) * weights[j];
        }
        notRun = null;
        var structure = fit.Structure[0];
        var valueWeights = fit.Weights(0)!;
        return new SetEstimate
        {
            Log2Effect = estimate,
            StandardError = Math.Sqrt(quad * fit.ResidualVariance[0]),
            Df = fit.SatterthwaiteDf(0, weights),
            DfMethod = "satterthwaite",
            ModelUsed = set ? ShareSetMixedModel : ShareMixedModel,
            RandomTerms = RandomTermsOf[(set, structure!.Value)],
            Terms = m,
            Peptides = peptides,
            Observations = n,
            ResidualVariance = fit.ResidualVariance[0],
            SampleVariance = fit.IndividualVariance[0],
            SampleMemberVariance = set && structure == RandomStructure.IndividualAndSample ? fit.SampleVariance[0] : double.NaN,
            Robust = robust,
            RobustFellBack = fit.RobustFellBack[0],
            Weights = Enumerable.Range(0, n).Select(i => new SetValueWeight(memberOf[i], peptideOf[i], sampleOf[i], valueWeights[i])).ToList(),
        };
    }

    /// <summary>
    /// The random terms each structure means. A set's "sample" is the sample and its "individual" the member term
    /// (<see cref="SetModel.FitSet"/>'s layout); a single's are STAT1's.
    /// </summary>
    private static readonly IReadOnlyDictionary<(bool Set, RandomStructure Structure), string> RandomTermsOf =
        new Dictionary<(bool, RandomStructure), string>
        {
            [(true, RandomStructure.IndividualAndSample)] = "sample;sample:member",
            [(true, RandomStructure.Individual)] = "sample",
            [(true, RandomStructure.Sample)] = "sample",
            [(false, RandomStructure.IndividualAndSample)] = "individual;individual:sample",
            [(false, RandomStructure.Individual)] = "individual",
            [(false, RandomStructure.Sample)] = "sample",
        };

    private static string SnakeCase(string name) =>
        string.Concat(name.Select((ch, i) => char.IsUpper(ch) && i > 0 ? "_" + char.ToLowerInvariant(ch) : char.ToLowerInvariant(ch).ToString()));
}

/// <summary>Runs one share part on one dataset contrast (<c>DEF-DIFF-SHARE-MODEL</c>, <c>-PART</c>; <c>DEF-DIFF-SET-VERDICT</c>).</summary>
public static class ShareTest
{
    /// <summary>
    /// Tests one share part on one dataset contrast: each member's share state (<see cref="ShareMemberStates"/>), the share
    /// single or share set model, the prediction's test on the share effect, and <c>DEF-DIFF-SET-VERDICT</c>'s guards. The
    /// outcome's <see cref="SetContrastOutcome.ShareMembers"/> carries the members; pooling and the family are as for a
    /// protein part.
    /// </summary>
    /// <param name="part">A share part (<see cref="SetHypothesisPart.Measure"/> is <see cref="SetTestFreeze.ShareMeasure"/>).</param>
    /// <param name="contrast">One of its dataset contrasts.</param>
    /// <param name="setup">The contrast's setup (<see cref="SetContrastSetup.Build"/>).</param>
    /// <param name="resolved">Its members, resolved (<see cref="SetMemberResolver.Resolve"/>) on the same samples.</param>
    /// <param name="normalized">The stratum's normalized observation table.</param>
    /// <param name="sequences">The search's databases' sequences.</param>
    /// <param name="genes">The same databases' gene names.</param>
    /// <param name="robust">STAT1's robust setting for the stratum (GR-8).</param>
    public static SetContrastOutcome Run(SetHypothesisPart part, SetDatasetContrast contrast, SetContrastSetup setup,
        IReadOnlyList<SetMemberAssignment> resolved, ObservationTable normalized, ShareSequenceIndex sequences, SetGeneIndex genes,
        bool robust = true)
    {
        ArgumentNullException.ThrowIfNull(part);
        ArgumentNullException.ThrowIfNull(contrast);
        ArgumentNullException.ThrowIfNull(setup);
        if (!string.Equals(part.Measure, SetTestFreeze.ShareMeasure, StringComparison.Ordinal))
            throw new ArgumentException($"{part.HypothesisId} {SetTestFreeze.PredictionName(part.Prediction)} is a {part.Measure} part; SetTest.Run tests it.", nameof(part));
        var outcome = new SetContrastOutcome { Part = part, Contrast = contrast, Setup = setup };
        if (setup.NotRun is not null)
            return outcome with { Untested = SetVerdictName.NotRun, UntestedReason = setup.NotRun };

        var members = ShareMemberStates.Apply(contrast, part.Prediction, resolved, setup, normalized, sequences, genes);
        outcome = outcome with { ShareMembers = members };
        var tested = members.Where(x => x.InEstimate).ToList();

        double? margin = SetTest.OneMargin(contrast.Rows.Select(r => r.MarginLog2), "margin_log2", contrast);
        double? attenuation = SetTest.OneMargin(contrast.Rows.Select(r => r.AttenuationFraction), "attenuation_fraction", contrast);
        double? mrnaMean = null;
        if ((part.Prediction is SetPrediction.AttenuatedDown or SetPrediction.AttenuatedUp) && tested.Count > 0)
        {
            // m̄ over the tested members, each weighted 1/m as β̄ weighs them (v3.25's rule; a share member is its own unit).
            if (tested.Any(x => x.Row.MrnaLog2 is null))
                return outcome with { Untested = SetVerdictName.NotRun, UntestedReason = SetNotRunReason.MrnaNotUsablePrefix + "missing" };
            mrnaMean = tested.Average(x => x.Row.MrnaLog2!.Value);
            if (part.Prediction == SetPrediction.AttenuatedDown ? mrnaMean >= 0 : mrnaMean <= 0)
                return outcome with { Untested = SetVerdictName.NotRun, UntestedReason = SetNotRunReason.MrnaNotUsablePrefix + "direction" };
        }

        SetEstimate? estimate;
        string? notRun;
        if (!part.IsSet)
        {
            if (members.Count != 1)
                throw new ArgumentException(
                    $"{part.HypothesisId} is tested as a single but lists {members.Count} members in {contrast.Dataset}.", nameof(contrast));
            if (tested.Count == 0)
                return outcome with { Untested = SetVerdictName.NotTestable, UntestedReason = "no_tested_member" };
            estimate = ShareModel.FitSingle(setup, tested[0], normalized, out notRun, robust);
        }
        else
        {
            if (tested.Count < SetTest.MinimumTestedMembers)
                return outcome with { Untested = SetVerdictName.NotTestable, UntestedReason = "fewer_than_3_tested" };
            estimate = ShareModel.FitSet(setup, tested, normalized, out notRun, robust);
        }
        if (estimate is null)
            return outcome with { Untested = SetVerdictName.NotRun, UntestedReason = notRun };

        var test = SetVerdict.Test(part.Prediction, estimate, margin, attenuation, mrnaMean);
        var (blocked, agree, disagree) = SetTest.Guards(part.Prediction, members.Select(x => (x.Assignment.Member, x.Assignment.State)));
        return outcome with { Test = test, BlockedBy = blocked, OneSideAgree = agree, OneSideDisagree = disagree };
    }
}
