using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Algorithms;
using MassSpectrometry.Deconvolution.Parameters;
using MzLibUtil;
using NUnit.Framework;
using Test.FileReadingTests;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class FlashDeconvTests
    {

        [Test]
        public void FlashDeconvWholeFileTest()
        {
            //txt file, not mgf, because it's an MS1. Most intense proteoform has mass of ~14037.9 Da
            string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "Temp");
            if(!Directory.Exists(tempDir))
            {
                Directory.CreateDirectory(tempDir);
            }
            string mzMlFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\SmallCalibratibleYeast.mzml");
            string outputFilePath = Path.Combine(tempDir, "output.mzML");

            DeconvolutionParameters deconParameters = new FlashDeconvDeconvolutionParamters();

            FlashDeconvOpenMsAlgorithm alg = new FlashDeconvOpenMsAlgorithm(deconParameters);

            alg.Bubba(mzMlFilePath, outputFilePath);

            Assert.That(File.Exists(outputFilePath), Is.True, "Output file was created.");

            // Clean up: delete temp directory and all contents
            if (Directory.Exists(tempDir))
            {
                Directory.Delete(tempDir, true);
            }
        }
    }
}
