using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMzSpectra
    {
        [Test]
        public void SpectrCroumCount()
        {
            CrossCorrelation c = new CrossCorrelation(new double[1], new double[1], new double[1], new double[1],
                CrossCorrelation.SpectrumNormalizationScheme.unnormalized, 5.0, true);


            double[] array = new double[3] { 1, 2, 3 };
            double value = 2.0;
            double tolerance = 0.5;

            Assert.IsTrue(c.DoubleWithinToleranceBool(array,value,tolerance));

        }
    }
}
