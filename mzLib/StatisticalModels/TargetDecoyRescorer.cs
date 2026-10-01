using System;
using System.Collections.Generic;
using System.Linq;

namespace StatisticalModels
{
    /// <summary>How a target-decoy rescoring ended.</summary>
    public enum RescoreStatus
    {
        /// <summary>Every fold trained a discriminant.</summary>
        Rescored,

        /// <summary>There were no decoys, so nothing could be scored; every score is NaN.</summary>
        NoDecoys,

        /// <summary>There were no targets, so nothing could be scored; every score is NaN.</summary>
        NoTargets,

        /// <summary>
        /// At least one fold's training rows lacked positives or decoys to train on. That fold was scored by the best single
        /// feature instead.
        /// </summary>
        FoldStarved,
    }

    /// <summary>The model the rescorer fits within each fold.</summary>
    public enum RescoreModel
    {
        /// <summary>A ridge Fisher linear discriminant (<see cref="LinearDiscriminant"/>), refit each iteration.</summary>
        LinearDiscriminant,

        /// <summary>
        /// The linear iterations choose the training rows, then an averaged ensemble of small tanh networks
        /// (<see cref="MultilayerPerceptron"/>) is trained on them, all targets against all decoys (as DIA-NN does), and scores
        /// by the logit of its mean probability.
        /// </summary>
        NeuralNetworkEnsemble,
    }

    /// <summary>The combined score of each candidate, the fold that scored it, and how the rescoring ended.</summary>
    public sealed class RescoreResult
    {
        internal RescoreResult(double[] scores, int[] folds, RescoreStatus status)
        {
            Scores = scores;
            Folds = folds;
            Status = status;
        }

        /// <summary>
        /// One score per candidate, in input order; higher is better. Each fold's scores are normalized on its own training
        /// rows (0 at that fold's q-value cutoff, −1 at its median decoy), so scores from different folds are comparable.
        /// </summary>
        public double[] Scores { get; }

        /// <summary>The fold each candidate belonged to; that fold's model never trained on it.</summary>
        public int[] Folds { get; }

        public RescoreStatus Status { get; }
    }

    /// <summary>
    /// Semi-supervised target-decoy rescoring in the style of mProphet and Percolator: combine many features into one
    /// score for target-decoy q-values, trained only on the search's own targets and decoys.
    /// <para>
    /// Candidates are split into folds by group (for example one peptide's charge states), and a group never straddles
    /// folds. Fold assignment uses only the group keys, never the labels. Each fold is scored by a model trained on the
    /// other folds, so no candidate is scored by a model that saw it. Within the training rows, the first iteration ranks by
    /// the best single feature (either sign). Each later iteration fits a <see cref="LinearDiscriminant"/> of targets
    /// passing the q-value cutoff against all decoys, then re-ranks.
    /// </para>
    /// </summary>
    public static class TargetDecoyRescorer
    {
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">Lengths disagree, rows are ragged, or a value is not finite.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A count or the q-value cutoff is out of range.</exception>
        public static RescoreResult Score(IReadOnlyList<double[]> features, IReadOnlyList<bool> isDecoy, IReadOnlyList<string> groupKeys,
            int folds = 3, int iterations = 3, double positiveQValue = 0.01, IReadOnlyList<int>? candidateGroups = null,
            RescoreModel model = RescoreModel.LinearDiscriminant, int? maxNetworkTrainingRows = null)
        {
            ArgumentNullException.ThrowIfNull(features);
            ArgumentNullException.ThrowIfNull(isDecoy);
            ArgumentNullException.ThrowIfNull(groupKeys);
            if (isDecoy.Count != features.Count)
                throw new ArgumentException($"There are {features.Count} candidates but {isDecoy.Count} decoy labels.", nameof(isDecoy));
            if (groupKeys.Count != features.Count)
                throw new ArgumentException($"There are {features.Count} candidates but {groupKeys.Count} group keys.", nameof(groupKeys));
            if (candidateGroups is not null && candidateGroups.Count != features.Count)
                throw new ArgumentException($"There are {features.Count} candidates but {candidateGroups.Count} candidate groups.", nameof(candidateGroups));
            if (folds < 2)
                throw new ArgumentOutOfRangeException(nameof(folds), folds, "There must be at least 2 folds.");
            if (iterations < 1)
                throw new ArgumentOutOfRangeException(nameof(iterations), iterations, "There must be at least 1 iteration.");
            if (maxNetworkTrainingRows is < 2)
                throw new ArgumentOutOfRangeException(nameof(maxNetworkTrainingRows), maxNetworkTrainingRows, "The network needs at least 2 training rows.");
            if (!(positiveQValue > 0 && positiveQValue < 1))
                throw new ArgumentOutOfRangeException(nameof(positiveQValue), positiveQValue, "The q-value cutoff must be in (0, 1).");
            int n = features.Count;
            int p = n == 0 ? 0 : features[0].Length;
            foreach (var row in features)
                if (row is null || row.Length != p || !row.All(double.IsFinite))
                    throw new ArgumentException("Every candidate needs the same number of finite features.", nameof(features));

            int[] fold = AssignFolds(groupKeys, folds);
            var scores = Enumerable.Repeat(double.NaN, n).ToArray();
            if (!isDecoy.Any(d => d))
                return new RescoreResult(scores, fold, RescoreStatus.NoDecoys);
            if (isDecoy.All(d => d))
                return new RescoreResult(scores, fold, RescoreStatus.NoTargets);

            // Folds are independent (each writes only its own held-out rows, with its own seed), so they run in parallel
            var starved = new bool[folds];
            System.Threading.Tasks.Parallel.For(0, folds, f =>
            {
                int[] train = Enumerable.Range(0, n).Where(i => fold[i] != f).ToArray();
                int[] test = Enumerable.Range(0, n).Where(i => fold[i] == f).ToArray();
                if (test.Length == 0)
                    return;

                var seed = BestSingleFeature(features, isDecoy, train, positiveQValue, candidateGroups);
                Func<double[], double> scorer = x => seed.Sign * x[seed.Feature];
                bool trained = false;
                for (int iteration = 0; iteration < iterations; iteration++)
                {
                    // Only each candidate group's top row under the current model trains it (pyProphet)
                    int[] active = TopPerGroup(train, candidateGroups, i => scorer(features[i]));
                    double[] trainScores = active.Select(i => scorer(features[i])).ToArray();
                    double[] q = QValues(trainScores, active.Select(i => isDecoy[i]).ToArray());
                    double cutoff = TrainingCutoff(q, active.Select(i => isDecoy[i]).ToArray(), positiveQValue);
                    var rows = new List<double[]>();
                    var positive = new List<bool>();
                    for (int t = 0; t < active.Length; t++)
                    {
                        bool decoy = isDecoy[active[t]];
                        if (decoy || q[t] <= cutoff)
                        {
                            rows.Add(features[active[t]]);
                            positive.Add(!decoy);
                        }
                    }
                    if (positive.Count(b => b) < 2 || positive.All(b => b))
                        break;
                    var fit = LinearDiscriminant.Fit(rows, positive);
                    scorer = x => fit.Score(x);
                    trained = true;
                }
                if (!trained)
                    starved[f] = true;

                if (model == RescoreModel.NeuralNetworkEnsemble && trained)
                {
                    int[] rows = TopPerGroup(train, candidateGroups, i => scorer(features[i]));
                    if (maxNetworkTrainingRows is int cap && rows.Length > cap)
                    {
                        // A random subsample of the fold's own training rows, seeded by the fold. Not the top rows by the
                        // linear score: where the line misses the signal, its top rows are the wrong ones.
                        var sampler = new Random(31 + f);
                        rows = rows.OrderBy(_ => sampler.Next()).Take(cap).Order().ToArray();
                    }
                    var ensemble = MultilayerPerceptron.TrainEnsemble(rows.Select(i => features[i]).ToList(), rows.Select(i => !isDecoy[i]).ToList(),
                        NetworkMembers, NetworkLayers, NetworkEpochs, seed: 17 + f);
                    scorer = x => ensemble.PredictLogit(x);
                }

                // Normalize on the training rows so that folds are comparable when pooled
                int[] finalActive = TopPerGroup(train, candidateGroups, i => scorer(features[i]));
                double[] finalTrain = finalActive.Select(i => scorer(features[i])).ToArray();
                bool[] trainDecoy = finalActive.Select(i => isDecoy[i]).ToArray();
                double[] finalQ = QValues(finalTrain, trainDecoy);
                double[] decoyScores = finalTrain.Where((_, t) => trainDecoy[t]).Order().ToArray();
                double medianDecoy = decoyScores.Length == 0 ? 0 : decoyScores[decoyScores.Length / 2];
                double[] passing = finalTrain.Where((_, t) => !trainDecoy[t] && finalQ[t] <= positiveQValue).ToArray();
                double threshold = passing.Length > 0 ? passing.Min() : (decoyScores.Length > 0 ? decoyScores[^1] : medianDecoy + 1);
                double scale = threshold > medianDecoy ? threshold - medianDecoy : 1;
                foreach (int i in test)
                    scores[i] = (scorer(features[i]) - threshold) / scale;
            });
            return new RescoreResult(scores, fold, starved.Any(s => s) ? RescoreStatus.FoldStarved : RescoreStatus.Rescored);
        }

        /// <summary>
        /// Groups sorted by key (ordinal) and dealt round-robin, so a group stays whole and the assignment depends only on the
        /// keys, never on the labels or the row order.
        /// </summary>
        internal static int[] AssignFolds(IReadOnlyList<string> groupKeys, int folds)
        {
            var foldOfGroup = groupKeys.Distinct().Order(StringComparer.Ordinal)
                .Select((key, index) => (key, index))
                .ToDictionary(g => g.key, g => g.index % folds, StringComparer.Ordinal);
            return groupKeys.Select(key => foldOfGroup[key]).ToArray();
        }

        private const int NetworkMembers = 5;
        private const int NetworkEpochs = 10;
        private static readonly int[] NetworkLayers = [25, 20, 15, 10, 5]; // DIA-NN 2020's architecture

        /// <summary>Fewer training positives than this and the training cutoff is relaxed.</summary>
        internal const int MinimumPositives = 10;

        private static readonly double[] RelaxedCutoffs = [0.05, 0.10, 0.25, 0.50];

        /// <summary>
        /// The q-value cutoff for choosing training positives: the requested one, relaxed step by step only while fewer than
        /// <see cref="MinimumPositives"/> targets pass. Weak features rarely pass targets at 1% on the first iteration, so this
        /// seeds a first model, as Percolator and mokapot do. It affects only which training rows are positives; the reported
        /// q-values are the caller's.
        /// </summary>
        internal static double TrainingCutoff(double[] q, bool[] isDecoy, double requested)
        {
            int Passing(double cutoff) => Enumerable.Range(0, q.Length).Count(t => !isDecoy[t] && q[t] <= cutoff);
            if (Passing(requested) >= MinimumPositives)
                return requested;
            foreach (double relaxed in RelaxedCutoffs.Where(c => c > requested))
                if (Passing(relaxed) >= MinimumPositives)
                    return relaxed;
            return RelaxedCutoffs[^1];
        }

        /// <summary>
        /// The single feature, and sign, that passes the most training targets. It uses the requested cutoff, or the first
        /// relaxed one at which some feature passes <see cref="MinimumPositives"/>. Ties, including the case where nothing
        /// passes anywhere, go to the larger standardized target-minus-decoy mean difference. That difference always carries the
        /// right sign, so a direction is never chosen by the order the features happen to be in.
        /// </summary>
        internal static (int Feature, int Sign) BestSingleFeature(IReadOnlyList<double[]> features, IReadOnlyList<bool> isDecoy, int[] train, double cutoff,
            IReadOnlyList<int>? candidateGroups = null)
        {
            int p = features.Count == 0 ? 0 : features[0].Length;
            var candidates = new List<(int Feature, int Sign, double Separation, double[] Q, bool[] Decoy)>();
            for (int j = 0; j < p; j++)
            {
                foreach (int sign in new[] { 1, -1 })
                {
                    // With candidate groups, each feature and sign judges only the group rows it ranks top
                    int[] rows = TopPerGroup(train, candidateGroups, i => sign * features[i][j]);
                    bool[] rowDecoy = rows.Select(i => isDecoy[i]).ToArray();
                    double[] values = rows.Select(i => features[i][j]).ToArray();
                    double separation = StandardizedMeanDifference(values, rowDecoy);
                    candidates.Add((j, sign, sign * separation, QValues(values.Select(v => sign * v).ToArray(), rowDecoy), rowDecoy));
                }
            }

            int PassingOf((int Feature, int Sign, double Separation, double[] Q, bool[] Decoy) k, double c) =>
                Enumerable.Range(0, k.Q.Length).Count(t => !k.Decoy[t] && k.Q[t] <= c);
            double chosenCutoff = new[] { cutoff }.Concat(RelaxedCutoffs.Where(c => c > cutoff))
                .FirstOrDefault(c => candidates.Any(k => PassingOf(k, c) >= MinimumPositives), cutoff);
            var best = candidates
                .OrderByDescending(k => PassingOf(k, chosenCutoff))
                .ThenByDescending(k => k.Separation)
                .ThenBy(k => k.Feature).ThenByDescending(k => k.Sign)
                .First();
            return (best.Feature, best.Sign);
        }

        /// <summary>
        /// The rows that train: all of them, or, with candidate groups, each group's top-scoring row (ties go to the earlier row),
        /// in input order.
        /// </summary>
        private static int[] TopPerGroup(int[] rows, IReadOnlyList<int>? candidateGroups, Func<int, double> score) =>
            candidateGroups is null ? rows
                : rows.GroupBy(i => candidateGroups[i]).Select(g => g.OrderByDescending(score).ThenBy(i => i).First()).Order().ToArray();

        /// <summary>(mean of targets − mean of decoys) / pooled SD; 0 when either class is missing or the feature is constant.</summary>
        internal static double StandardizedMeanDifference(double[] values, bool[] isDecoy)
        {
            double[] targets = values.Where((_, t) => !isDecoy[t]).ToArray();
            double[] decoys = values.Where((_, t) => isDecoy[t]).ToArray();
            if (targets.Length == 0 || decoys.Length == 0)
                return 0;
            double mean = values.Average();
            double sd = Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / Math.Max(1, values.Length - 1));
            return sd > LinearDiscriminant.ConstantTolerance ? (targets.Average() - decoys.Average()) / sd : 0;
        }

        /// <summary>Target-decoy q-values in input order, from the shared <see cref="TargetDecoyQValues"/>.</summary>
        internal static double[] QValues(double[] scores, bool[] isDecoy) => TargetDecoyQValues.Compute(scores, isDecoy);
    }
}
