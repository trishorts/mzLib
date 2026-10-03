using System;
using System.Collections.Generic;
using System.Linq;

namespace StatisticalModels
{
    /// <summary>
    /// A small feed-forward binary classifier: standardised inputs, tanh hidden layers and a sigmoid output, trained on
    /// cross-entropy with Adam in mini-batches. Training is deterministic for a seed. It gives the target-decoy rescorer a
    /// non-linear model, as DIA-NN's classifier ensemble does (Demichev et al. 2020, Nat. Methods 17:41).
    /// </summary>
    public sealed class MultilayerPerceptron
    {
        private readonly double[] _mean;
        private readonly double[] _sd;
        private readonly double[][,] _weights; // layer l: [inputs, outputs]
        private readonly double[][] _biases;

        private MultilayerPerceptron(double[] mean, double[] sd, double[][,] weights, double[][] biases)
        {
            _mean = mean;
            _sd = sd;
            _weights = weights;
            _biases = biases;
        }

        /// <summary>The probability that <paramref name="features"/> belong to the positive class.</summary>
        /// <exception cref="ArgumentException">The feature count differs from training.</exception>
        public double Predict(double[] features) => Sigmoid(PredictLogit(features));

        /// <summary>The log-odds that <paramref name="features"/> belong to the positive class: the output before the sigmoid,
        /// which keeps its resolution where the probability has rounded to 0 or 1.</summary>
        /// <exception cref="ArgumentException">The feature count differs from training.</exception>
        public double PredictLogit(double[] features)
        {
            ArgumentNullException.ThrowIfNull(features);
            if (features.Length != _mean.Length)
                throw new ArgumentException($"Expected {_mean.Length} features; got {features.Length}.", nameof(features));
            Span<double> first = stackalloc double[Width];
            Span<double> second = stackalloc double[Width];
            return Forward(features, first, second);
        }

        /// <summary>Features the network was trained on.</summary>
        internal int InputCount => _mean.Length;

        /// <summary>The widest layer, inputs included: the size of each buffer <see cref="Forward"/> needs.</summary>
        internal int Width => Math.Max(_mean.Length, _weights.Max(w => w.GetLength(1)));

        /// <summary>
        /// <see cref="PredictLogit(double[])"/> on two caller buffers of at least <see cref="Width"/>, allocating nothing:
        /// prediction is the hot loop of a DIA rescoring. Same arithmetic, in the same order.
        /// </summary>
        internal double Forward(double[] features, Span<double> first, Span<double> second)
        {
            int width = _mean.Length;
            for (int j = 0; j < width; j++)
                first[j] = (features[j] - _mean[j]) / _sd[j];
            bool inFirst = true;
            for (int l = 0; l < _weights.Length; l++)
            {
                var w = _weights[l];
                var bias = _biases[l];
                bool last = l == _weights.Length - 1;
                Span<double> input = (inFirst ? first : second)[..width];
                Span<double> output = inFirst ? second : first;
                int outputs = w.GetLength(1);
                for (int o = 0; o < outputs; o++)
                {
                    double z = bias[o];
                    for (int i = 0; i < width; i++)
                        z += input[i] * w[i, o];
                    output[o] = last ? z : Math.Tanh(z);
                }
                width = outputs;
                inFirst = !inFirst;
            }
            return (inFirst ? first : second)[0];
        }

        /// <param name="hiddenLayers">Units in each hidden layer, input side first.</param>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">Lengths disagree, rows are ragged or not finite, or a class is missing.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A layer has no units, or a count or rate is not positive.</exception>
        public static MultilayerPerceptron Train(IReadOnlyList<double[]> features, IReadOnlyList<bool> isPositive, IReadOnlyList<int> hiddenLayers,
            int epochs, int seed, int batchSize = 50, double learningRate = 0.003)
        {
            ArgumentNullException.ThrowIfNull(features);
            ArgumentNullException.ThrowIfNull(isPositive);
            ArgumentNullException.ThrowIfNull(hiddenLayers);
            if (features.Count != isPositive.Count)
                throw new ArgumentException($"There are {features.Count} rows but {isPositive.Count} labels.", nameof(isPositive));
            if (features.Count == 0 || isPositive.All(b => b) || !isPositive.Any(b => b))
                throw new ArgumentException("Both classes must be present.", nameof(isPositive));
            if (hiddenLayers.Any(units => units < 1))
                throw new ArgumentOutOfRangeException(nameof(hiddenLayers), "Every hidden layer needs at least one unit.");
            ArgumentOutOfRangeException.ThrowIfLessThan(epochs, 1);
            ArgumentOutOfRangeException.ThrowIfLessThan(batchSize, 1);
            if (!(learningRate > 0))
                throw new ArgumentOutOfRangeException(nameof(learningRate), learningRate, "The learning rate must be positive.");
            int p = features[0].Length;
            if (features.Any(row => row is null || row.Length != p || !row.All(double.IsFinite)))
                throw new ArgumentException("Every row needs the same number of finite features.", nameof(features));

            // Standardisation from the training rows; a constant feature is centred only
            var mean = new double[p];
            var sd = new double[p];
            for (int j = 0; j < p; j++)
            {
                mean[j] = features.Average(row => row[j]);
                double variance = features.Sum(row => (row[j] - mean[j]) * (row[j] - mean[j])) / Math.Max(1, features.Count - 1);
                sd[j] = variance > 1e-24 ? Math.Sqrt(variance) : 1;
            }

            var random = new Random(seed);
            int[] sizes = [p, .. hiddenLayers, 1];
            var weights = new double[sizes.Length - 1][,];
            var biases = new double[sizes.Length - 1][];
            for (int l = 0; l < weights.Length; l++)
            {
                double limit = Math.Sqrt(6.0 / (sizes[l] + sizes[l + 1])); // Xavier/Glorot uniform
                weights[l] = new double[sizes[l], sizes[l + 1]];
                for (int i = 0; i < sizes[l]; i++)
                    for (int o = 0; o < sizes[l + 1]; o++)
                        weights[l][i, o] = (2 * random.NextDouble() - 1) * limit;
                biases[l] = new double[sizes[l + 1]];
            }
            var net = new MultilayerPerceptron(mean, sd, weights, biases);
            net.Fit(features, isPositive, epochs, batchSize, learningRate, random);
            return net;
        }

        /// <summary>An ensemble of <paramref name="members"/> networks, each seeded differently; it predicts their mean.</summary>
        public static MultilayerPerceptronEnsemble TrainEnsemble(IReadOnlyList<double[]> features, IReadOnlyList<bool> isPositive, int members,
            IReadOnlyList<int> hiddenLayers, int epochs, int seed, int batchSize = 50, double learningRate = 0.003)
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(members, 1);
            // Members train in parallel; each has its own seed, so the result is the same as training them in turn
            var trained = new MultilayerPerceptron[members];
            System.Threading.Tasks.Parallel.For(0, members, m => trained[m] = Train(features, isPositive, hiddenLayers, epochs, seed + 7919 * m, batchSize, learningRate));
            return new MultilayerPerceptronEnsemble(trained);
        }

        private void Fit(IReadOnlyList<double[]> features, IReadOnlyList<bool> isPositive, int epochs, int batchSize, double learningRate, Random random)
        {
            int layers = _weights.Length;
            var mW = _weights.Select(w => new double[w.GetLength(0), w.GetLength(1)]).ToArray();
            var vW = _weights.Select(w => new double[w.GetLength(0), w.GetLength(1)]).ToArray();
            var mB = _biases.Select(b => new double[b.Length]).ToArray();
            var vB = _biases.Select(b => new double[b.Length]).ToArray();
            var gW = _weights.Select(w => new double[w.GetLength(0), w.GetLength(1)]).ToArray();
            var gB = _biases.Select(b => new double[b.Length]).ToArray();
            const double beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8;
            int step = 0;
            int[] order = Enumerable.Range(0, features.Count).ToArray();
            // Buffers reused for every row: activations per layer, and the backpropagated error at each layer's output
            var activations = new double[layers + 1][];
            activations[0] = new double[_mean.Length];
            for (int l = 0; l < layers; l++)
                activations[l + 1] = new double[_weights[l].GetLength(1)];
            var deltas = activations.Select(a => new double[a.Length]).ToArray();

            for (int epoch = 0; epoch < epochs; epoch++)
            {
                for (int i = order.Length - 1; i > 0; i--) // Fisher–Yates with the seeded generator
                {
                    int k = random.Next(i + 1);
                    (order[i], order[k]) = (order[k], order[i]);
                }
                for (int start = 0; start < order.Length; start += batchSize)
                {
                    int end = Math.Min(order.Length, start + batchSize);
                    foreach (var g in gW) Array.Clear(g);
                    foreach (var g in gB) Array.Clear(g);
                    for (int r = start; r < end; r++)
                    {
                        int row = order[r];
                        StandardiseInto(features[row], activations[0]);
                        for (int l = 0; l < layers; l++)
                            LayerInto(activations[l], l, last: l == layers - 1, activations[l + 1]);
                        // Cross-entropy on the sigmoid output: dLoss/dz = p - y
                        double[] delta = deltas[layers];
                        delta[0] = Sigmoid(activations[layers][0]) - (isPositive[row] ? 1 : 0);
                        for (int l = layers - 1; l >= 0; l--)
                        {
                            double[] input = activations[l];
                            for (int o = 0; o < delta.Length; o++)
                            {
                                gB[l][o] += delta[o];
                                for (int i = 0; i < input.Length; i++)
                                    gW[l][i, o] += input[i] * delta[o];
                            }
                            if (l == 0)
                                break;
                            double[] previous = deltas[l];
                            for (int i = 0; i < input.Length; i++)
                            {
                                double sum = 0;
                                for (int o = 0; o < delta.Length; o++)
                                    sum += _weights[l][i, o] * delta[o];
                                previous[i] = sum * (1 - input[i] * input[i]); // tanh'
                            }
                            delta = previous;
                        }
                    }

                    step++;
                    double n = end - start;
                    double correction1 = 1 - Math.Pow(beta1, step), correction2 = 1 - Math.Pow(beta2, step);
                    for (int l = 0; l < layers; l++)
                    {
                        for (int i = 0; i < _weights[l].GetLength(0); i++)
                            for (int o = 0; o < _weights[l].GetLength(1); o++)
                            {
                                double g = gW[l][i, o] / n;
                                mW[l][i, o] = beta1 * mW[l][i, o] + (1 - beta1) * g;
                                vW[l][i, o] = beta2 * vW[l][i, o] + (1 - beta2) * g * g;
                                _weights[l][i, o] -= learningRate * (mW[l][i, o] / correction1) / (Math.Sqrt(vW[l][i, o] / correction2) + epsilon);
                            }
                        for (int o = 0; o < _biases[l].Length; o++)
                        {
                            double g = gB[l][o] / n;
                            mB[l][o] = beta1 * mB[l][o] + (1 - beta1) * g;
                            vB[l][o] = beta2 * vB[l][o] + (1 - beta2) * g * g;
                            _biases[l][o] -= learningRate * (mB[l][o] / correction1) / (Math.Sqrt(vB[l][o] / correction2) + epsilon);
                        }
                    }
                }
            }
        }

        private void StandardiseInto(double[] features, double[] x)
        {
            for (int j = 0; j < x.Length; j++)
                x[j] = (features[j] - _mean[j]) / _sd[j];
        }

        private void LayerInto(double[] input, int l, bool last, double[] output)
        {
            var w = _weights[l];
            for (int o = 0; o < output.Length; o++)
            {
                double z = _biases[l][o];
                for (int i = 0; i < input.Length; i++)
                    z += input[i] * w[i, o];
                output[o] = last ? z : Math.Tanh(z);
            }
        }

        private static double Sigmoid(double z) => 1 / (1 + Math.Exp(-z));
    }

    /// <summary>An averaged ensemble of <see cref="MultilayerPerceptron"/>s.</summary>
    public sealed class MultilayerPerceptronEnsemble
    {
        internal MultilayerPerceptronEnsemble(IReadOnlyList<MultilayerPerceptron> members) => Members = members;

        public IReadOnlyList<MultilayerPerceptron> Members { get; }

        /// <summary>The members' mean probability.</summary>
        public double Predict(double[] features) => Members.Average(m => m.Predict(features));

        /// <summary>
        /// The log-odds of the members' mean probability, computed from their logits so that it never saturates: a confident
        /// member's probability rounds to exactly 1 in double precision, and every confident row would tie.
        /// </summary>
        public double PredictLogit(double[] features)
        {
            ArgumentNullException.ThrowIfNull(features);
            if (Members.Count == 0 || features.Length != Members[0].InputCount)
                return LogitOfMeanProbability(Members.Select(m => m.PredictLogit(features)).ToArray()); // the members report the error
            int width = Members.Max(m => m.Width);
            Span<double> first = stackalloc double[width];
            Span<double> second = stackalloc double[width];
            Span<double> logits = stackalloc double[Members.Count];
            for (int m = 0; m < Members.Count; m++)
                logits[m] = Members[m].Forward(features, first, second);
            return LogitOfMeanProbability(logits);
        }

        /// <summary>log(mean sigmoid(z)) - log(mean sigmoid(-z)), each a log-sum-exp of log-sigmoids.</summary>
        internal static double LogitOfMeanProbability(IReadOnlyList<double> logits) =>
            LogSumExp(logits.Select(z => -Softplus(-z))) - LogSumExp(logits.Select(z => -Softplus(z)));

        /// <summary><see cref="LogitOfMeanProbability(IReadOnlyList{double})"/> without allocating, summed in the same order.</summary>
        internal static double LogitOfMeanProbability(ReadOnlySpan<double> logits) => LogSumExpOfLogSigmoid(logits, -1) - LogSumExpOfLogSigmoid(logits, 1);

        /// <summary>log sum exp(-softplus(sign * z)) over the logits.</summary>
        private static double LogSumExpOfLogSigmoid(ReadOnlySpan<double> logits, double sign)
        {
            double max = double.NegativeInfinity;
            foreach (double z in logits)
                max = Math.Max(max, -Softplus(sign * z));
            double sum = 0;
            foreach (double z in logits)
                sum += Math.Exp(-Softplus(sign * z) - max);
            return max + Math.Log(sum);
        }

        private static double Softplus(double t) => Math.Max(t, 0) + Math.Log(1 + Math.Exp(-Math.Abs(t)));

        private static double LogSumExp(IEnumerable<double> values)
        {
            double[] v = values.ToArray();
            double max = v.Max();
            return max + Math.Log(v.Sum(x => Math.Exp(x - max)));
        }
    }
}
