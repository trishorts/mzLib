using NUnit.Framework;
using StatisticalModels;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.StatisticalModels;

/// <summary>
/// A small feed-forward network (tanh hidden layers, sigmoid output, cross-entropy, Adam) and an averaged ensemble of
/// them. It gives the target-decoy rescorer a non-linear model, as DIA-NN's classifier does.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class MultilayerPerceptronTests
{
    private static double Gaussian(Random random) =>
        Math.Sqrt(-2 * Math.Log(1 - random.NextDouble())) * Math.Cos(2 * Math.PI * random.NextDouble());

    /// <summary>Positives inside a ring, negatives outside: no straight line separates them.</summary>
    private static (List<double[]> X, List<bool> Y) Ring(int n, int seed)
    {
        var random = new Random(seed);
        var x = new List<double[]>();
        var y = new List<bool>();
        for (int i = 0; i < n; i++)
        {
            double a = 3 * Gaussian(random), b = 3 * Gaussian(random);
            x.Add([a, b]);
            y.Add(Math.Sqrt(a * a + b * b) < 3);
        }
        return (x, y);
    }

    private static double Accuracy(Func<double[], double> predict, List<double[]> x, List<bool> y) =>
        x.Select((v, i) => (predict(v) > 0.5) == y[i]).Count(ok => ok) / (double)x.Count;

    [Test]
    public void TheNetworkLearnsWhatNoLineCan()
    {
        var (x, y) = Ring(3000, 1);
        var (tx, ty) = Ring(1000, 2);

        var net = MultilayerPerceptron.Train(x, y, hiddenLayers: [16, 8], epochs: 60, seed: 7);
        var line = LinearDiscriminant.Fit(x, y);

        Assert.That(Accuracy(net.Predict, tx, ty), Is.GreaterThan(0.93), "network on held-out data");
        Assert.That(Accuracy(v => line.Score(v) > 0 ? 1 : 0, tx, ty), Is.LessThan(0.75), "a line cannot do it");
    }

    [Test]
    public void TrainingIsDeterministicForASeed()
    {
        var (x, y) = Ring(500, 3);

        var a = MultilayerPerceptron.Train(x, y, [8, 4], epochs: 5, seed: 11);
        var b = MultilayerPerceptron.Train(x, y, [8, 4], epochs: 5, seed: 11);
        var c = MultilayerPerceptron.Train(x, y, [8, 4], epochs: 5, seed: 12);

        Assert.That(x.Select(a.Predict), Is.EqualTo(x.Select(b.Predict)));
        Assert.That(x.Select(a.Predict), Is.Not.EqualTo(x.Select(c.Predict)));
    }

    /// <summary>Outputs are probabilities, and the ensemble is the mean of its members.</summary>
    [Test]
    public void TheEnsembleAveragesItsMembers()
    {
        var (x, y) = Ring(500, 4);

        var ensemble = MultilayerPerceptron.TrainEnsemble(x, y, members: 3, [8, 4], epochs: 5, seed: 5);

        var v = x[17];
        Assert.That(ensemble.Members, Has.Count.EqualTo(3));
        Assert.That(ensemble.Predict(v), Is.EqualTo(ensemble.Members.Average(m => m.Predict(v))).Within(1e-12));
        Assert.That(x.Select(ensemble.Predict), Is.All.InRange(0.0, 1.0));
    }

    [Test]
    public void ArgumentsAreChecked()
    {
        var (x, y) = Ring(50, 5);
        Assert.Throws<ArgumentNullException>(() => MultilayerPerceptron.Train(null!, y, [4], 1, 1));
        Assert.Throws<ArgumentException>(() => MultilayerPerceptron.Train(x, y.Take(10).ToList(), [4], 1, 1));
        Assert.Throws<ArgumentException>(() => MultilayerPerceptron.Train(x, x.Select(_ => true).ToList(), [4], 1, 1), "both classes needed");
        Assert.Throws<ArgumentOutOfRangeException>(() => MultilayerPerceptron.Train(x, y, [0], 1, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => MultilayerPerceptron.Train(x, y, [4], 0, 1));
        Assert.Throws<ArgumentException>(() => MultilayerPerceptron.Train(x, y, [4], 1, 1).Predict([1.0]), "feature count");
    }
}
