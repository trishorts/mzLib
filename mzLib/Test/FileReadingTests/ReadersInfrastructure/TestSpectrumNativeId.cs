using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.ReadersInfrastructure
{
    /// <summary>
    /// <see cref="SpectrumNativeId"/> called directly. MzIdentMLResultFile turns every spectrumID through it, so
    /// TestMzIdentMLResultFile covers it in context; these cases pin the shared entry point itself.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestSpectrumNativeId
    {
        [TestCase(null, -1)]
        [TestCase("", -1)]
        [TestCase("scan=42", 42)]
        [TestCase("controllerType=0 controllerNumber=1 scan=7", 7)]
        [TestCase("scan=81409 file=454", 81409)]
        [TestCase("index=0", 1)]
        [TestCase("index=16332", 16333)]
        [TestCase("mzMLid=controllerType=0", -1)]
        [TestCase("index=12x", -1)]
        [TestCase("scan=1,scan=2", 1)]
        [TestCase("spectrum=2458", 2458)]                 // OpenMS: the mzML spectrum whose id is spectrum=2458
        [TestCase("spectrum=0", 0)]
        [TestCase("spectrum=12x", -1)]
        [TestCase("frame=1 scan=5", -1)]                  // Bruker: the scan counts within the frame
        [TestCase("merged=0 frame=1 scanStart=1 scanEnd=10", -1)]
        [TestCase("sample=1 period=1 cycle=12 experiment=3", -1)] // SCIEX: no single scan number
        [TestCase("123", -1)]                             // a bare number (Mascot through jmzTab): meaning unknown
        public void OneBasedScanNumber_FromTheNativeId(string nativeId, int expected)
        {
            Assert.That(SpectrumNativeId.OneBasedScanNumberOf(nativeId), Is.EqualTo(expected));
        }
    }
}
