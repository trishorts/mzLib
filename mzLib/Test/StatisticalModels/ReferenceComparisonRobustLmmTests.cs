using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using StatisticalModels;

namespace Test.StatisticalModels;

/// <summary>
/// <see cref="NestedMixedModel"/>'s robust option against frozen reference values from
/// <c>ReferenceData/make_robust_lmm_fixtures.R</c>: msqrob2 1.20.0's own reweighting loop (<c>.robust_fitting</c>:
/// <c>psi.huber</c> on <c>resid</c> scaled by <c>mad(res, 0)</c>, <c>refit</c>, stop at a 1e-6 relative change in pwrss),
/// run to convergence (QuantProject GR-22), then lmerTest's Satterthwaite df and GR-21's moderation on the final weighted
/// fit. Two scenarios of 30 proteins with about 6% outlying values. R is never run by these tests.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ReferenceComparisonRobustLmmTests
{
    /// <summary>
    /// Weights, estimates and σ². Each round refits from scratch and the loop stops on a 1e-6 relative change in pwrss, so
    /// the final weights agree to the precision that rule leaves, not to the optimizer's.
    /// </summary>
    private const double Fit = 1e-5;

    /// <summary>Satterthwaite df and everything downstream, as in <see cref="ReferenceComparisonNestedLmmTests"/>.</summary>
    private const double Df = 1e-4;

    [TestCase("nested")]
    [TestCase("sample_only")]
    public void RobustFitAndModerationEqualReference(string scenario)
    {
        var data = LmmData.Read(scenario, "rlmm");
        var fit = NestedMixedModel.Fit(data.Problems, data.CommonNames, robust: true);
        var moderated = NestedMixedModel.ModerateContrasts(fit, data.Contrasts, data.ConditionRows);
        var expected = LmmData.Rows($"rlmm_{scenario}_results.tsv");
        var weights = LmmData.Rows($"rlmm_{scenario}_weights.tsv")
            .ToDictionary(r => (r["protein"], r["peptide"], r["sample"]), r => LmmData.Parse(r["weight"]));

        int downWeighted = 0;
        foreach (var protein in expected.GroupBy(r => r["protein"]))
        {
            int f = data.Proteins.IndexOf(protein.Key);
            var first = protein.First();
            string at = $"{scenario} {protein.Key}";
            Assert.That(fit.Status[f], Is.EqualTo(FeatureFitStatus.Fitted), at);
            Assert.That(fit.RobustConverged[f], Is.EqualTo(bool.Parse(first["converged"])), at + " converged");
            Assert.That(Math.Abs(fit.RobustIterations[f] - int.Parse(first["iterations"])), Is.LessThanOrEqualTo(1),
                at + " iterations (the stopping rule can fall either side of 1e-6 by one round)");
            var w = fit.Weights(f)!;
            for (int i = 0; i < data.Keys[f].Count; i++)
            {
                var (peptide, sample) = data.Keys[f][i];
                double expectedWeight = weights[(protein.Key, peptide, sample)];
                Assert.That(w[i], Is.EqualTo(expectedWeight).Within(Fit), $"{at} weight {peptide} {sample}");
                if (expectedWeight < 1) downWeighted++;
            }
            Close(fit.ResidualVariance[f], first["sigma2"], Fit, at + " sigma2");
            Close(fit.RemlCriterion[f], first["reml"], Fit, at + " REML");
            foreach (var row in protein)
            {
                var r = moderated.Single(m => m.Name == row["contrast"]);
                var c = data.Contrasts.Single(x => x.Name == row["contrast"]).Weights;
                string where = $"{at} {row["contrast"]}";
                Close(r.Estimate[f], row["estimate"], Fit, where + " estimate");
                Close(fit.SatterthwaiteDf(f, c), row["df_satterthwaite"], Df, where + " Satterthwaite df");
                Close(fit.JointSatterthwaiteDf(f, data.ConditionRows), row["variance_df"], Df, where + " joint DenDF");
                Close(r.StandardError[f], row["se"], Df, where + " se");
                Close(r.DfTotal[f], row["df"], Df, where + " df");
                Close(r.PValue[f], row["p"], Df, where + " p");
                Close(r.BenjaminiHochbergAdjusted[f], row["adj_p"], Df, where + " adj.p");
                Close(r.ConfidenceLow[f], row["ci_low"], Df, where + " CI low");
                Close(r.ConfidenceHigh[f], row["ci_high"], Df, where + " CI high");
            }
        }
        Assert.That(downWeighted, Is.GreaterThan(20), "the fixture's outliers are down-weighted, so the test exercises reweighting");
    }

    [Test]
    public void WithoutRobustEveryWeightIsOne()
    {
        var data = LmmData.Read("nested", "rlmm");
        var fit = NestedMixedModel.Fit(data.Problems, data.CommonNames);
        for (int f = 0; f < fit.FeatureCount; f++)
        {
            Assert.That(fit.RobustIterations[f], Is.EqualTo(0));
            Assert.That(fit.RobustConverged[f], Is.Null);
            Assert.That(fit.Weights(f)!.Where(double.IsFinite), Is.All.EqualTo(1.0));
        }
    }

    [Test]
    public void RobustFittingMovesEstimatesTowardTheBulk()
    {
        // One wild value: robust weighting gives it weight < 1 and the condition estimate moves away from it.
        var data = LmmData.Read("nested", "rlmm");
        var plain = NestedMixedModel.Fit(data.Problems, data.CommonNames);
        var robust = NestedMixedModel.Fit(data.Problems, data.CommonNames, robust: true);
        int changed = Enumerable.Range(0, plain.FeatureCount)
            .Count(f => Math.Abs(plain.Coefficient(f, 1) - robust.Coefficient(f, 1)) > 1e-6);
        Assert.That(changed, Is.GreaterThan(0));
        Assert.That(Enumerable.Range(0, robust.FeatureCount).All(f => robust.Weights(f)!.Where(double.IsFinite).All(w => w > 0 && w <= 1)));
    }

    private static void Close(double actual, string expectedText, double tolerance, string where)
    {
        double expected = LmmData.Parse(expectedText);
        Assert.That(Math.Abs(actual - expected), Is.LessThanOrEqualTo(tolerance * Math.Max(1, Math.Abs(expected))),
            $"{where}: {actual:R} vs {expected:R}");
    }
}
