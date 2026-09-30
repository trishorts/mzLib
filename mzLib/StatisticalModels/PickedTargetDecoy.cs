using System;
using System.Collections.Generic;

namespace StatisticalModels
{
    public static class TargetDecoyQValues
    {
        public static double[] Compute(IReadOnlyList<double> scores, IReadOnlyList<bool> isDecoy) => throw new NotImplementedException();
    }

    public sealed class PickedResult
    {
        internal PickedResult(bool[] kept, double[] qValues)
        {
            Kept = kept;
            QValues = qValues;
        }

        public bool[] Kept { get; }
        public double[] QValues { get; }
    }

    public static class PickedTargetDecoy
    {
        public static PickedResult Compete(IReadOnlyList<string> pairKeys, IReadOnlyList<double> scores, IReadOnlyList<bool> isDecoy) =>
            throw new NotImplementedException();
    }
}
