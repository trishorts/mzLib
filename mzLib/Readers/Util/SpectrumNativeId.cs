using System.Globalization;
using System.Text.RegularExpressions;

namespace Readers
{
    /// <summary>
    /// Reads the scan number out of a spectrum's nativeID, the PSI-MS identifier that mzIdentML, mzTab and mzML
    /// all use to name a spectrum within its file.
    /// </summary>
    internal static class SpectrumNativeId
    {
        private static readonly Regex Index = new(@"(^|\s)index=(\d+)($|\s)", RegexOptions.Compiled);
        private static readonly Regex Spectrum = new(@"(^|\s)spectrum=(\d+)($|\s)", RegexOptions.Compiled);
        private static readonly Regex Frame = new(@"(^|\s)frame=", RegexOptions.Compiled);

        /// <summary>
        /// The nativeID's scan number:
        /// <list type="bullet">
        /// <item>"scan=N" gives N;</item>
        /// <item>"spectrum=N" gives N. It is the spectrum's id in its file, not a position: OpenMS's cut-down mzMLs
        /// number spectra 1011, 1012, ... from their first entry;</item>
        /// <item>"index=N", a zero-based position in a peak list, gives N + 1, as CasanovoMzTabFile does;</item>
        /// <item>a Bruker "frame=F scan=S" gives -1, because S counts within the frame and repeats in every frame;</item>
        /// <item>anything else gives -1.</item>
        /// </list>
        /// CAREFUL: these are the producer's numbers. mzLib's Mzml reader numbers a file's spectra by position and
        /// ignores the nativeID, so they equal its OneBasedScanNumber only for a complete file in its original order.
        /// </summary>
        internal static int OneBasedScanNumberOf(string? nativeId)
        {
            if (string.IsNullOrEmpty(nativeId))
            {
                return -1;
            }

            if (Frame.IsMatch(nativeId))
            {
                return -1;
            }

            var scan = Mzml.nativeIdScanNumberParser.Match(nativeId);
            if (scan.Success && int.TryParse(scan.Groups[2].Value, NumberStyles.None, CultureInfo.InvariantCulture, out int scanNumber))
            {
                return scanNumber;
            }

            var spectrum = Spectrum.Match(nativeId);
            if (spectrum.Success && int.TryParse(spectrum.Groups[2].Value, NumberStyles.None, CultureInfo.InvariantCulture, out int spectrumNumber))
            {
                return spectrumNumber;
            }

            var index = Index.Match(nativeId);
            if (index.Success && int.TryParse(index.Groups[2].Value, NumberStyles.None, CultureInfo.InvariantCulture, out int position))
            {
                return position + 1;
            }

            return -1;
        }
    }
}
