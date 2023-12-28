using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMzSpectra
    {
        [Test]
        public void BinarySearchArrayBool()
        {
            CrossCorrelation c = new CrossCorrelation(new double[1]{1}, new double[1] { 1 }, new double[1] { 1 }, new double[1] { 1 },
                CrossCorrelation.SpectrumNormalizationScheme.unnormalized, 5.0, true);


            double[] array = new double[3] { 1, 2, 3 };
            double value = 2.0;
            PpmTolerance ppmTolerance = new PpmTolerance(1);

            Assert.IsTrue(c.DoubleWithinToleranceBool(array,value,ppmTolerance));

            array = new double[5] { 1, 2, 3,4, 5 };
            Assert.IsTrue(c.DoubleWithinToleranceBool(array, value, ppmTolerance));

            value = 2.5;
            Assert.IsFalse(c.DoubleWithinToleranceBool(array, value, ppmTolerance));

            value = 5.5;
            Assert.IsFalse(c.DoubleWithinToleranceBool(array, value, ppmTolerance));

            value = 0.5;
            Assert.IsFalse(c.DoubleWithinToleranceBool(array, value, ppmTolerance));

            value = 1.0 - 1e-6;
            Assert.IsFalse(c.DoubleWithinToleranceBool(array, value, ppmTolerance));

            value = 1 + 1e-6;
            Assert.IsTrue(c.DoubleWithinToleranceBool(array, value, ppmTolerance));

            value = 5.0 - 1e-6;
            Assert.IsTrue(c.DoubleWithinToleranceBool(array, value, ppmTolerance));

            value = 5.0 + 1e-6;
            Assert.IsTrue(c.DoubleWithinToleranceBool(array, value, ppmTolerance));
        }
        [Test]
        public void BinarySearchArrayValue()
        {
            CrossCorrelation c = new CrossCorrelation(new double[1] { 1 }, new double[1] { 1 }, new double[1] { 1 }, new double[1] { 1 },
                CrossCorrelation.SpectrumNormalizationScheme.unnormalized, 5.0, true);


            double[] array = new double[3] { 1, 2, 3 };
            double value = 2.0;
            PpmTolerance ppmTolerance = new PpmTolerance(1);

            Assert.That(value, Is.EqualTo(c.DoubleWithinToleranceValue(array, value, ppmTolerance)).Within(1e-6));

            array = new double[5] { 1, 2, 3, 4, 5 };
            Assert.That(value, Is.EqualTo(c.DoubleWithinToleranceValue(array, value, ppmTolerance)).Within(1e-6));

            value = 2.5;
            Assert.IsNull(c.DoubleWithinToleranceValue(array, value, ppmTolerance));

            value = 5.5;
            Assert.IsNull(c.DoubleWithinToleranceValue(array, value, ppmTolerance));

            value = 0.5;
            Assert.IsNull(c.DoubleWithinToleranceValue(array, value, ppmTolerance));

            value = 1.0 - 1e-6;
            Assert.IsNull(c.DoubleWithinToleranceValue(array, value, ppmTolerance));

            value = 1 + 1e-6;
            Assert.That(value, Is.EqualTo(c.DoubleWithinToleranceValue(array, value, ppmTolerance)).Within(1e-6));

            value = 5.0 - 1e-6;
            Assert.That(value, Is.EqualTo(c.DoubleWithinToleranceValue(array, value, ppmTolerance)).Within(1e-6));

            value = 5.0 + 1e-6;
            Assert.That(value, Is.EqualTo(c.DoubleWithinToleranceValue(array, value, ppmTolerance)).Within(1e-6));
        }
    }
}
