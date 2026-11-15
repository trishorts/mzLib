using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// SETTEST1 S4b (QuantProject <c>DEF-DIFF-SHARE-MODEL</c>, <c>-PART</c>): the share single and share set models against
/// lmerTest, plain and robust (<c>ReferenceData/make_settest_share_model_fixtures.R</c>, R never run here); a part split by
/// measure with its flags read over the whole original part (ruling S5); and <see cref="ShareTest.Run"/>'s outcomes.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SetTestShareModelTests
{
    private const double Fit = 1e-6;
    private const double Df = 1e-4;

    private static void Close(double actual, double expected, double tolerance, string what) =>
        Assert.That(Math.Abs(actual - expected), Is.LessThanOrEqualTo(tolerance * Math.Max(1, Math.Abs(expected))),
            $"{what}: {actual:R} vs {expected:R}");

    private static (string[] Header, List<string[]> Rows) Tsv(string file)
    {
        var lines = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "Quantification", "Differential",
            "ReferenceData", file)).Where(l => l.Length > 0).ToArray();
        return (lines[0].Split('\t'), lines.Skip(1).Select(l => l.Split('\t')).ToList());
    }

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private const string Header =
        "hypothesis_id\tkind\tprediction\tmember\tspecies\tdataset\tfactor\ttest_level\treference_level\tcontexts\tmrna_log2" +
        "\tindependence\tpost_hoc\tmeasure\tisoform_accession\tmarker_region";

    private static SetTestFreeze Freeze(params string[] rows) =>
        SetTestFreeze.Read(Encoding.UTF8.GetBytes(string.Join('\n', new[] { Header }.Concat(rows)) + "\n"), "t");

    private static string Row(string member, string kind = "single", string prediction = "down", string contexts = "",
        string measure = "exon_peptide", string mrna = "", string independence = "independent", string postHoc = "no",
        string id = "H1", string dataset = "D1") =>
        string.Join('\t', id, kind, prediction, member, "mouse", dataset, "age group", "old", "young", contexts, mrna, independence,
            postHoc, measure, measure == "protein" ? "" : "A_" + member, measure == "protein" ? "" : "A_" + member + ":1-5");

    private static readonly DesignSpec AgeSpec = new() { Factors = new[] { new FactorSpec("age group", "young") } };

    // ── the share models against lmerTest ───────────────────────────────────

    private sealed record Scenario(SetContrastSetup Setup, ObservationTable Table, List<ShareMemberResult> Members,
        Func<string, string, double> Result, List<string[]> RobustWeights);

    private static Scenario Load(string name, bool individuals = false, string? extraSample = null)
    {
        bool set = name.StartsWith("set", StringComparison.Ordinal);
        var samples = Tsv($"settest_share_model_{name}_samples.tsv").Rows
            .Select(c => new AnalysisSample(c[0])
            {
                Factors = new Dictionary<string, string> { ["age group"] = c[1], ["treatment"] = c[2] },
                Individual = individuals ? "mouse-" + c[0] : null,
            })
            .ToList();
        var designSamples = extraSample is null ? samples
            : samples.Append(new AnalysisSample(extraSample) { Factors = new Dictionary<string, string> { ["age group"] = "old", ["treatment"] = "none" } }).ToList();
        var values = Tsv($"settest_share_model_{name}_values.tsv").Rows.Select(c => (Member: c[0], Peptide: c[1], Kind: c[2], Sample: c[3], Y: D(c[4]))).ToList();
        var runs = samples.Select(s => new ObservationRun(s.SampleId + ".raw", s.SampleId)).ToList();
        var peptides = values.GroupBy(v => (v.Member, v.Peptide))
            .Select(g => new PeptideRunValues(g.Key.Peptide, g.Key.Peptide, new[] { g.Key.Member },
                g.ToDictionary(v => v.Sample + ".raw", v => new RunValue(Math.Pow(2, v.Y), ObservationState.Quantified)))).ToList();
        var members = values.Select(v => v.Member).Distinct().Order(StringComparer.Ordinal).ToList();
        var table = ObservationTable.Build(runs, runs.Select(r => r.FileName), peptides,
            members.Select(m => new ProteinGroupInfo(m, 0.001, false, false)), QuantBasis.MsmsOnly);
        string contexts = name.EndsWith("context", StringComparison.Ordinal) ? "treatment=a|b" : "";
        var part = Freeze(members.Select(m => Row(m, set ? "set" : "single", contexts: contexts)).ToArray()).Parts.Single();
        var contrast = part.Contrasts.Single();
        var kindOf = values.GroupBy(v => v.Peptide).ToDictionary(g => g.Key, g => g.First().Kind);
        var results = members.Select(m => new ShareMemberResult
        {
            Assignment = new SetMemberAssignment { Member = m, ProteinGroup = m, State = SetMemberStateFromRow.Tested },
            Row = contrast.Rows.Single(r => r.Member == m),
            Peptides = new SharePeptides
            {
                Markers = table.Peptides.Where(p => p.ProteinGroups[0] == m && kindOf[p.FullSequence] == "marker").ToList(),
                References = table.Peptides.Where(p => p.ProteinGroups[0] == m && kindOf[p.FullSequence] == "reference").ToList(),
            },
        }).ToList();
        var (rh, rr) = Tsv($"settest_share_model_{name}_result.tsv");
        var (bh, br) = Tsv($"settest_share_model_{name}_robust_result.tsv");
        double Result(string which, string column) => which == "robust"
            ? D(br.Single()[Array.IndexOf(bh, column)])
            : D(rr.Single()[Array.IndexOf(rh, column)]);
        return new Scenario(SetContrastSetup.Build(contrast, designSamples, AgeSpec), table, results, Result,
            Tsv($"settest_share_model_{name}_robust_weights.tsv").Rows);
    }

    private static SetEstimate FitScenario(Scenario s, bool robust)
    {
        var estimate = s.Members.Count == 1
            ? ShareModel.FitSingle(s.Setup, s.Members[0], s.Table, out string? notRun, robust)
            : ShareModel.FitSet(s.Setup, s.Members, s.Table, out notRun, robust);
        Assert.That(notRun, Is.Null);
        return estimate!;
    }

    [TestCase("single_balanced")]
    [TestCase("single_missing")]
    [TestCase("single_context")]
    [TestCase("set_balanced")]
    [TestCase("set_context")]
    public void TheShareModel_MatchesLmerTest(string name)
    {
        var s = Load(name);
        Assert.That(s.Setup.NotRun, Is.Null);

        var e = FitScenario(s, robust: false);

        Close(e.Log2Effect, s.Result("plain", "estimate"), Fit, name + " share effect");
        Close(e.StandardError, s.Result("plain", "se"), Fit, name + " se");
        Close(e.Df, s.Result("plain", "df"), Df, name + " Satterthwaite df");
        Close(e.ResidualVariance, s.Result("plain", "sigma2"), Fit, name + " sigma^2");
        Close(e.SampleVariance, s.Result("plain", "var_sample"), Fit, name + " var (1|sample)");
        bool set = s.Members.Count > 1;
        if (set) Close(e.SampleMemberVariance, s.Result("plain", "var_sample_member"), Fit, name + " var (1|sample:member)");
        else Assert.That(double.IsNaN(e.SampleMemberVariance));
        Assert.That(e.Terms, Is.EqualTo((int)s.Result("plain", "members")));
        Assert.That(e.Peptides, Is.EqualTo((int)s.Result("plain", "peptides")));
        Assert.That(e.Observations, Is.EqualTo((int)s.Result("plain", "observations")));
        Assert.That(e.ModelUsed, Is.EqualTo(set ? ShareModel.ShareSetMixedModel : ShareModel.ShareMixedModel));
        Assert.That(e.RandomTerms, Is.EqualTo(set ? "sample;sample:member" : "sample"));
        Assert.That(e.DfMethod, Is.EqualTo("satterthwaite"));
        Assert.That(e.Robust, Is.False);
    }

    [TestCase("single_balanced")]
    [TestCase("single_missing")]
    [TestCase("single_context")]
    [TestCase("set_balanced")]
    [TestCase("set_context")]
    public void TheRobustShareModel_MatchesLmerTestUnderMsqrob2sWeights(string name)
    {
        var s = Load(name);

        var e = FitScenario(s, robust: true);

        Close(e.Log2Effect, s.Result("robust", "estimate"), Fit, name + " robust share effect");
        Close(e.StandardError, s.Result("robust", "se"), Fit, name + " robust se");
        Close(e.Df, s.Result("robust", "df"), Df, name + " robust df");
        Close(e.ResidualVariance, s.Result("robust", "sigma2"), Fit, name + " robust sigma^2");
        Assert.That(e.Robust, Is.True);
        Assert.That(e.RobustFellBack, Is.False);
        Assert.That(e.Weights.Count(v => v.Weight < 1), Is.EqualTo((int)s.Result("robust", "down_weighted")));
        var mine = e.Weights.ToDictionary(v => (v.ProteinGroup, v.Peptide, v.Sample), v => v.Weight);
        Assert.That(mine, Has.Count.EqualTo(s.RobustWeights.Count));
        foreach (var w in s.RobustWeights)
            Close(mine[(w[0], w[1], w[2])], D(w[3]), Fit, $"{name} weight {w[1]} {w[2]}");
    }

    [Test]
    public void TheShareEffect_IsTheMarkersChange_LessTheReferencesChange()
    {
        // Two references and a marker over o1, o2 (old) and y1, y2 (young): every reference rises 1 with age, the marker
        // 0.4; noise-free values give the share effect -0.6 exactly, whatever the sample shifts.
        var shift = new Dictionary<string, double> { ["o1"] = 0.3, ["o2"] = -0.2, ["y1"] = 0.1, ["y2"] = -0.15 };
        var rows = new List<(string Peptide, string Sample, double Y)>();
        foreach (var (peptide, baseline, ageChange) in new[] { ("R1", 20.0, 1.0), ("R2", 22.0, 1.0), ("M1", 18.0, 0.4) })
            foreach (var sample in shift.Keys)
                rows.Add((peptide, sample, baseline + shift[sample] + (sample.StartsWith('o') ? ageChange : 0) + (sample == "o1" && peptide == "R2" ? 0.01 : 0)));
        var samples = shift.Keys.Select(s => new AnalysisSample(s) { Factors = new Dictionary<string, string> { ["age group"] = s.StartsWith('o') ? "old" : "young" } }).ToList();
        var runs = samples.Select(s => new ObservationRun(s.SampleId + ".raw", s.SampleId)).ToList();
        var table = ObservationTable.Build(runs, runs.Select(r => r.FileName),
            rows.GroupBy(r => r.Peptide).Select(g => new PeptideRunValues(g.Key, g.Key, new[] { "G" },
                g.ToDictionary(r => r.Sample + ".raw", r => new RunValue(Math.Pow(2, r.Y), ObservationState.Quantified)))).ToList(),
            new[] { new ProteinGroupInfo("G", 0.001, false, false) }, QuantBasis.MsmsOnly);
        var contrast = Freeze(Row("G")).Parts.Single().Contrasts.Single();
        var member = new ShareMemberResult
        {
            Assignment = new SetMemberAssignment { Member = "G", ProteinGroup = "G", State = SetMemberStateFromRow.Tested },
            Row = contrast.Rows.Single(),
            Peptides = new SharePeptides { Markers = table.Peptides.Where(p => p.FullSequence == "M1").ToList(), References = table.Peptides.Where(p => p.FullSequence != "M1").ToList() },
        };

        var e = ShareModel.FitSingle(SetContrastSetup.Build(contrast, samples, AgeSpec), member, table, out string? notRun, robust: false);

        Assert.That(notRun, Is.Null);
        Assert.That(e!.Log2Effect, Is.EqualTo(-0.6 - 0.0025).Within(1e-9),
            "R2's 0.01 in o1 raises the reference's change by 0.01/4, so the share effect falls by it");
    }

    /// <summary>With one sample per individual, <c>(1|individual)</c> is <c>(1|sample)</c>: the same fit, named as STAT1 names it.</summary>
    [Test]
    public void ASingleWithIndividuals_TakesStat1sIndividualTerm()
    {
        var s = Load("single_balanced", individuals: true);
        var e = FitScenario(s, robust: false);
        Close(e.Log2Effect, s.Result("plain", "estimate"), Fit, "share effect");
        Close(e.StandardError, s.Result("plain", "se"), Fit, "se");
        Assert.That(e.RandomTerms, Is.EqualTo("individual"));
    }

    [Test]
    public void APeptideWithNoValueInTheContrast_IsLeftOutOfTheModel()
    {
        var s = Load("single_balanced");
        var empty = new ObservedPeptide("EMPTYPEPTIDE", "EMPTYPEPTIDE", new[] { "G1" });
        var table = ObservationTable.Build(s.Table.Samples.Select(x => new ObservationRun(x + ".raw", x)).ToList(),
            s.Table.Samples.Select(x => x + ".raw"),
            s.Table.Peptides.Select((p, k) => new PeptideRunValues(p.FullSequence, p.BaseSequence, p.ProteinGroups,
                    s.Table.Samples.Select((x, c) => (x, v: s.Table.Log2Intensity(k, c))).Where(t => double.IsFinite(t.v))
                        .ToDictionary(t => t.x + ".raw", t => new RunValue(Math.Pow(2, t.v), ObservationState.Quantified))))
                .Append(new PeptideRunValues("EMPTYPEPTIDE", "EMPTYPEPTIDE", new[] { "G1" }, new Dictionary<string, RunValue>())).ToList(),
            new[] { new ProteinGroupInfo("G1", 0.001, false, false) }, QuantBasis.MsmsOnly);
        var member = s.Members[0] with
        {
            Peptides = new SharePeptides
            {
                Markers = table.Peptides.Where(p => p.FullSequence.Contains('M')).ToList(),
                References = table.Peptides.Where(p => !p.FullSequence.Contains('M') || p.FullSequence == "EMPTYPEPTIDE").ToList(),
            },
        };
        var e = ShareModel.FitSingle(s.Setup, member, table, out _, robust: false);
        Assert.That(e!.Peptides, Is.EqualTo((int)s.Result("plain", "peptides")));
        Close(e.Log2Effect, s.Result("plain", "estimate"), Fit, "the same fit");
    }

    /// <summary>Each sample holds one value, so <c>(1|sample)</c> cannot be told from the residual: not run, naming why.</summary>
    [Test]
    public void AShareModelThatCannotBeFormed_IsNotRunNamingWhy()
    {
        var t = Toy(new[] { "Ga" }, "single", "down", value: (m, p, s) =>
            (p == Marker(m) && s is "o1" or "o2" or "y1" or "y2") || (p.StartsWith("AAPEG") && s is "o3" or "o4" or "y3" or "y4")
            || (p.StartsWith("LLN") && s is "o5" or "y5") ? null : double.NaN);
        var o = RunToy(t);
        Assert.That(o.ShareMembers.Single().Assignment.State, Is.EqualTo(SetMemberStateFromRow.Tested));
        Assert.That(o.Untested, Is.EqualTo(SetVerdictName.NotRun));
        Assert.That(o.UntestedReason, Does.StartWith(SetNotRunReason.NotConverged + ":"));
        Assert.That(o.Test, Is.Null);
    }

    // ── refusals ─────────────────────────────────────────────────────────────

    [Test]
    public void AMemberNotTested_ASecondCopy_OrAPeptideOutsideTheTable_IsRefused()
    {
        var s = Load("set_balanced");
        var untested = s.Members[0] with { Assignment = s.Members[0].Assignment with { State = ShareMemberStates.TooFewMarkers } };
        Assert.That(() => ShareModel.FitSet(s.Setup, new[] { untested, s.Members[1], s.Members[2] }, s.Table, out _),
            Throws.ArgumentException.With.Message.Contains("Only tested members"));
        Assert.That(() => ShareModel.FitSet(s.Setup, new[] { s.Members[0], s.Members[0], s.Members[1] }, s.Table, out _),
            Throws.ArgumentException.With.Message.Contains("at most once"));
        var stranger = s.Members[0] with
        {
            Peptides = s.Members[0].Peptides! with { Markers = new[] { new ObservedPeptide("NOTINTABLE", "NOTINTABLE", new[] { "G1" }) } },
        };
        Assert.That(() => ShareModel.FitSingle(s.Setup, stranger, s.Table, out _),
            Throws.ArgumentException.With.Message.Contains("is not a peptide of the table"));
        Assert.That(() => ShareModel.FitSingle(new SetContrastSetup { NotRun = "spans_strata" }, s.Members[0], s.Table, out _),
            Throws.ArgumentException.With.Message.Contains("not run"));
        var weighsIntercept = s.Setup with { Contrast = new Contrast("c", "c", new[] { 1.0, 0.0 }, null, null) };
        Assert.That(() => ShareModel.FitSingle(weighsIntercept, s.Members[0], s.Table, out _),
            Throws.ArgumentException.With.Message.Contains("does not weigh it"));
        var extra = Load("single_balanced", extraSample: "s99");
        Assert.That(() => ShareModel.FitSingle(extra.Setup, extra.Members[0], extra.Table, out _),
            Throws.ArgumentException.With.Message.Contains("no sample 's99'"));
    }

    [Test]
    public void AnIncompleteMember_EntersTheShareModel()
    {
        var s = Load("set_balanced");
        var incomplete = s.Members[0] with { Assignment = s.Members[0].Assignment with { State = SetMemberStateFromRow.Incomplete } };
        var withIt = ShareModel.FitSet(s.Setup, new[] { incomplete, s.Members[1], s.Members[2] }, s.Table, out _, robust: false);
        Close(withIt!.Log2Effect, s.Result("plain", "estimate"), Fit, "same model, an incomplete member included (ruled 2026-11-15)");
    }

    // ── a part split by measure (ruling S5) ─────────────────────────────────

    /// <summary>
    /// H15's case: one part names an exon (EIF4H) and a protein (SFPQ). It is split into a protein part and a share part,
    /// each with its own members and contrasts; the flags are read over the whole original part, so SFPQ's post-hoc row
    /// makes EIF4H's D1 post hoc, and SFPQ's reproduction row makes EIF4H's D2 a reproduction.
    /// </summary>
    [Test]
    public void APartMixingMeasures_IsSplit_WithItsFlagsReadOverTheWholePart()
    {
        var freeze = Freeze(
            Row("EIF4H", "set", dataset: "D1"), Row("SFPQ", "set", measure: "protein", postHoc: "yes", dataset: "D1"),
            Row("EIF4H", "set", dataset: "D2"), Row("SFPQ", "set", measure: "protein", independence: "reproduces [9]", dataset: "D2"),
            Row("EIF4H", "set", dataset: "D3"), Row("SFPQ", "set", measure: "protein", dataset: "D3"));

        Assert.That(freeze.Parts.Select(p => (p.HypothesisId, p.Measure)), Is.EqualTo(new[] { ("H1", "share"), ("H1", "protein") }));
        var share = freeze.Parts.Single(p => p.Measure == SetTestFreeze.ShareMeasure);
        var protein = freeze.Parts.Single(p => p.Measure == SetTestFreeze.ProteinMeasure);
        Assert.That(share.Members, Is.EqualTo(new[] { "EIF4H" }));
        Assert.That(protein.Members, Is.EqualTo(new[] { "SFPQ" }));
        Assert.That(share.Contrasts.SelectMany(c => c.Rows).All(r => r.Member == "EIF4H"));
        Assert.That(share.Contrasts.Select(c => (c.Dataset, c.PostHoc, c.Independence)), Is.EqualTo(new[]
        {
            ("D1", true, SetIndependence.Independent), ("D2", false, SetIndependence.Reproduction), ("D3", false, SetIndependence.Independent),
        }));
        Assert.That(protein.Contrasts.Select(c => (c.PostHoc, c.Independence)),
            Is.EqualTo(share.Contrasts.Select(c => (c.PostHoc, c.Independence))));
        Assert.That(share.Kind, Is.EqualTo("set"));
    }

    [Test]
    public void APartOfOneMeasure_IsOnePart_AndAnUnknownMeasureIsRefused()
    {
        Assert.That(Freeze(Row("G1", measure: "protein"), Row("G2", measure: "protein")).Parts.Single().Measure, Is.EqualTo("protein"));
        Assert.That(Freeze(Row("G1"), Row("G2", measure: "isoform_peptide")).Parts.Single().Measure, Is.EqualTo("share"));
        Assert.That(() => Freeze(Row("G1", measure: "rna_peptide")), Throws.TypeOf<FormatException>().With.Message.Contains("Line 2: measure 'rna_peptide'"));
    }

    // ── ShareTest.Run ────────────────────────────────────────────────────────
    //
    // Each member has two forms carrying its own residue t (W, Y, H or F), so its peptides are its own:
    //   C1: MStRK AAPEGtK VDLtEWNR QQFtSGR LLNtAMK TTPtWR   (the exon QQFtSGR at 21-27)
    //   I1: the same without the exon.
    // Markers QQFtSGR; references AAPEGtK and LLNtAMK. Five samples a side, so 2 of 5 is incomplete but not too few.
    private static readonly string[] Old = { "o1", "o2", "o3", "o4", "o5" };
    private static readonly string[] Young = { "y1", "y2", "y3", "y4", "y5" };
    private static readonly string[] Members4 = { "Ga", "Gb", "Gc", "Gd" };

    private static string Tag(string member) => "WYHF"[Array.IndexOf(Members4, member)].ToString();
    private static string Form(string member, bool exon) => "MS" + Tag(member) + "RK" + "AAPEG" + Tag(member) + "K" + "VDL" + Tag(member) + "EWNR"
        + (exon ? "QQF" + Tag(member) + "SGR" : "") + "LLN" + Tag(member) + "AMK" + "TTP" + Tag(member) + "WR";
    private static string Marker(string member) => "QQF" + Tag(member) + "SGR";

    private static string ShareRow(string member, string kind, string prediction, string mrna) =>
        string.Join('\t', "H1", kind, prediction, member, "mouse", "D1", "age group", "old", "young", "", mrna, "independent", "no",
            "exon_peptide", member + "_C1", member + "_C1:21-27");

    /// <summary>
    /// A toy contrast: <paramref name="value"/> gives each (member, peptide, sample) value, NaN for none; by default the
    /// references rise 1 with age and the marker 0.2, so the share effect is -0.8.
    /// </summary>
    private static (SetHypothesisPart Part, SetContrastSetup Setup, ObservationTable Table, ShareSequenceIndex Sequences, SetGeneIndex Genes,
        IReadOnlyList<SetMemberAssignment> Resolved) Toy(string[] members, string kind, string prediction, string mrna = "",
            Func<string, string, string, double?>? value = null)
    {
        var part = SetTestFreeze.Read(Encoding.UTF8.GetBytes(Header + "\n" + string.Join("\n", members.Select(m => ShareRow(m, kind, prediction, mrna))) + "\n"), "t")
            .Parts.Single();
        var samples = Old.Select(s => new AnalysisSample(s) { Factors = new Dictionary<string, string> { ["age group"] = "old" } })
            .Concat(Young.Select(s => new AnalysisSample(s) { Factors = new Dictionary<string, string> { ["age group"] = "young" } })).ToList();
        var setup = SetContrastSetup.Build(part.Contrasts.Single(), samples, AgeSpec);
        var sequences = ShareSequenceIndex.FromSequences(members.SelectMany(m => new[] { (m + "_C1", Form(m, true)), (m + "_I1", Form(m, false)) }));
        var genes = SetGeneIndex.FromPrimaryNames(members.SelectMany(m => new[] { (m + "_C1", (IEnumerable<string>)new[] { m }), (m + "_I1", new[] { m }) }));
        var rng = new Random(7);
        var rows = new List<PeptideRunValues>();
        foreach (var m in members)
            foreach (var peptide in new[] { "AAPEG" + Tag(m) + "K", "LLN" + Tag(m) + "AMK", Marker(m) })
            {
                bool marker = peptide == Marker(m);
                var cells = new Dictionary<string, RunValue>();
                foreach (var s in Old.Concat(Young))
                {
                    double noise = rng.NextDouble() * 0.1;
                    double v = value?.Invoke(m, peptide, s) ?? 20 + (s[0] == 'o' ? (marker ? 0.2 : 1.0) : 0) + noise;
                    if (double.IsFinite(v)) cells[s + ".raw"] = new RunValue(Math.Pow(2, v), ObservationState.Quantified);
                }
                rows.Add(new PeptideRunValues(peptide, peptide, new[] { m + "_C1" }, cells));
            }
        var runs = Old.Concat(Young).Select(s => new ObservationRun(s + ".raw", s)).ToList();
        var table = ObservationTable.Build(runs, runs.Select(r => r.FileName), rows,
            members.Select(m => new ProteinGroupInfo(m + "_C1", 0.001, false, false)), QuantBasis.MsmsOnly);
        return (part, setup, table, sequences, genes, members.Select(m => new SetMemberAssignment { Member = m, ProteinGroup = m + "_C1" }).ToList());
    }

    private static SetContrastOutcome RunToy((SetHypothesisPart Part, SetContrastSetup Setup, ObservationTable Table, ShareSequenceIndex Sequences,
        SetGeneIndex Genes, IReadOnlyList<SetMemberAssignment> Resolved) t) =>
        ShareTest.Run(t.Part, t.Part.Contrasts.Single(), t.Setup, t.Resolved, t.Table, t.Sequences, t.Genes);

    [Test]
    public void AShareSingle_IsTestedOnItsShareEffect()
    {
        var o = RunToy(Toy(new[] { "Ga" }, "single", "down"));

        Assert.That(o.Untested, Is.Null);
        Assert.That(o.Members, Is.Empty);
        var member = o.ShareMembers.Single();
        Assert.That(member.Assignment.State, Is.EqualTo(SetMemberStateFromRow.Tested));
        Assert.That(member.Peptides!.MarkerBaseSequences, Is.EqualTo(new[] { "QQFWSGR" }));
        Assert.That(o.Test!.Estimate.ModelUsed, Is.EqualTo(ShareModel.ShareMixedModel));
        Assert.That(o.Test.Estimate.Log2Effect, Is.EqualTo(-0.8).Within(0.1), "the marker rises 0.2 where the reference rises 1");
        Assert.That(o.Test.Estimate.Robust, Is.True, "robust by default, as STAT1");
    }

    [Test]
    public void AShareSet_NeedsThreeTestedMembers_AndASingle_ItsOne()
    {
        var o = RunToy(Toy(new[] { "Ga", "Gb" }, "set", "down"));
        Assert.That((o.Untested, o.UntestedReason), Is.EqualTo((SetVerdictName.NotTestable, "fewer_than_3_tested")));
        Assert.That(o.ShareMembers, Has.Count.EqualTo(2), "each member's own result is kept");

        var tested = RunToy(Toy(new[] { "Ga", "Gb", "Gc" }, "set", "down"));
        Assert.That(tested.Test!.Estimate.ModelUsed, Is.EqualTo(ShareModel.ShareSetMixedModel));
        Assert.That(tested.Test.Estimate.Terms, Is.EqualTo(3));
        Assert.That(tested.Test.Estimate.Log2Effect, Is.EqualTo(-0.8).Within(0.1));

        var untested = RunToy(Toy(new[] { "Ga" }, "single", "down", value: (m, p, s) => p == Marker(m) ? double.NaN : 20));
        Assert.That((untested.Untested, untested.UntestedReason), Is.EqualTo((SetVerdictName.NotTestable, "no_tested_member")));
        Assert.That(untested.ShareMembers.Single().Assignment.State, Is.EqualTo(ShareMemberStates.TooFewMarkers));
    }

    /// <summary>Gd's marker is seen in old samples only: counted as agreeing with up, never folded into the effect.</summary>
    [Test]
    public void AOneSidedMember_IsCounted_NotFolded()
    {
        var o = RunToy(Toy(Members4, "set", "up",
            value: (m, p, s) => m == "Gd" && p == Marker(m) && s[0] == 'y' ? double.NaN : null));
        Assert.That(o.ShareMembers.Single(x => x.Assignment.Member == "Gd").Assignment.State, Is.EqualTo(SetMemberStateFromRow.SeenNumeratorOnly));
        Assert.That((o.OneSideAgree, o.OneSideDisagree), Is.EqualTo((1, 0)));
        Assert.That(o.Test!.Estimate.Terms, Is.EqualTo(3));
        Assert.That(o.BlockedBy, Is.Empty);
    }

    /// <summary>Gd's marker has values in 2 of 5 young samples: incomplete for unchanged, so in the effect and blocking (2026-11-15).</summary>
    [Test]
    public void AnIncompleteMember_StaysInTheEffect_AndBlocksSupported()
    {
        var o = RunToy(Toy(Members4, "set", "unchanged",
            value: (m, p, s) => m == "Gd" && p == Marker(m) && s is "y3" or "y4" or "y5" ? double.NaN : null));
        Assert.That(o.ShareMembers.Single(x => x.Assignment.Member == "Gd").Assignment.State, Is.EqualTo(SetMemberStateFromRow.Incomplete));
        Assert.That(o.Test!.Estimate.Terms, Is.EqualTo(4));
        Assert.That(o.BlockedBy, Is.EqualTo(new[] { "Gd" }));
    }

    [Test]
    public void AnAttenuatedPrediction_NeedsTheMembersMrna_OnItsSide()
    {
        var o = RunToy(Toy(new[] { "Ga" }, "single", "attenuated_down"));
        Assert.That((o.Untested, o.UntestedReason), Is.EqualTo((SetVerdictName.NotRun, "mrna_not_usable:missing")));

        o = RunToy(Toy(new[] { "Ga" }, "single", "attenuated_down", mrna: "1.5"));
        Assert.That((o.Untested, o.UntestedReason), Is.EqualTo((SetVerdictName.NotRun, "mrna_not_usable:direction")));
        o = RunToy(Toy(new[] { "Ga" }, "single", "attenuated_up", mrna: "-1.5"));
        Assert.That((o.Untested, o.UntestedReason), Is.EqualTo((SetVerdictName.NotRun, "mrna_not_usable:direction")));

        o = RunToy(Toy(new[] { "Ga", "Gb", "Gc" }, "set", "attenuated_down", mrna: "-2"));
        Assert.That(o.Test!.MrnaMean, Is.EqualTo(-2));
    }

    [Test]
    public void EachRunnerTestsOnlyItsOwnMeasure()
    {
        var t = Toy(new[] { "Ga" }, "single", "down");
        var protein = Freeze(Row("Ga", measure: "protein")).Parts.Single();
        Assert.That(() => ShareTest.Run(protein, protein.Contrasts.Single(), t.Setup, t.Resolved, t.Table, t.Sequences, t.Genes),
            Throws.ArgumentException.With.Message.Contains("is a protein part"));
        Assert.That(() => SetTest.Run(t.Part, t.Part.Contrasts.Single(), t.Setup, t.Resolved, t.Table, new Dictionary<string, DifferentialResult>()),
            Throws.ArgumentException.With.Message.Contains("is a share part"));
        var two = Toy(new[] { "Ga", "Gb" }, "single", "down");
        Assert.That(() => RunToy(two), Throws.ArgumentException.With.Message.Contains("is tested as a single but lists 2 members"));
        var notRun = ShareTest.Run(t.Part, t.Part.Contrasts.Single(), new SetContrastSetup { NotRun = "spans_strata" }, t.Resolved, t.Table, t.Sequences, t.Genes);
        Assert.That((notRun.Untested, notRun.UntestedReason), Is.EqualTo((SetVerdictName.NotRun, "spans_strata")));
    }
}
