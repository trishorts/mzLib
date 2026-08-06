using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// Interprets an SDRF cell's key=value grammar, e.g.
    /// "NT=Oxidation;AC=UNIMOD:35;TA=M;MT=Variable".
    ///
    /// This is the opt-in projection that <see cref="SdrfDocument"/> deliberately does not perform
    /// at read time. Cells are stored raw so that documents round-trip losslessly; anything that
    /// wants to reason about a cell's meaning asks for it here.
    ///
    /// Detection is by LEADING KEY, never by shape. A cell is only treated as a term if it starts
    /// with one of the recognised SDRF keys. Testing for "contains '=' and ';'" instead would
    /// misfire on real data: comment[file uri] cells in the curated corpus hold pre-signed download
    /// URLs whose query strings carry "Signature=", "Expires=" and "Id=" -- 5,750 occurrences each
    /// -- and a shape-based parser decodes those as controlled-vocabulary descriptors.
    /// </summary>
    public static class SdrfCell
    {
        /// <summary>
        /// Keys observed across the whole 1,236-file curated corpus, most frequent first:
        /// NT name, AC accession, MT modification type, TA target amino acid, VV version value,
        /// PP position, MM monoisotopic mass, SN specificity name, CS cleavage site,
        /// CF chemical formula, CT compound type, QY quantity, SP species, CN common name,
        /// ML/MH mass low/high, A generic attribute.
        /// </summary>
        private static readonly HashSet<string> KnownKeys = new(StringComparer.Ordinal)
        {
            "NT", "AC", "MT", "TA", "VV", "PP", "MM", "SN", "CS", "CF", "CT", "QY", "SP", "CN", "ML", "MH", "A"
        };

        /// <summary>
        /// True when the cell is written in the key=value grammar rather than being free text.
        /// </summary>
        public static bool IsTerm(string cell)
        {
            if (string.IsNullOrEmpty(cell)) return false;
            int equals = cell.IndexOf('=');
            if (equals <= 0) return false;
            int semicolon = cell.IndexOf(';');
            if (semicolon >= 0 && semicolon < equals) return false;
            return KnownKeys.Contains(cell.Substring(0, equals).Trim());
        }

        /// <summary>
        /// Splits a cell into its key=value pairs, preserving order and duplicates-last-wins.
        /// Returns an empty dictionary for free text. Keys are upper-cased for lookup; values keep
        /// their original casing and inner whitespace.
        /// </summary>
        public static IReadOnlyDictionary<string, string> ParseKeyValues(string cell)
        {
            var pairs = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            if (!IsTerm(cell)) return pairs;

            foreach (var part in cell.Split(';'))
            {
                int equals = part.IndexOf('=');
                if (equals <= 0) continue;
                string key = part.Substring(0, equals).Trim();
                if (key.Length == 0) continue;
                pairs[key] = part.Substring(equals + 1).Trim();
            }
            return pairs;
        }

        /// <summary>
        /// Reads the cell as a controlled-vocabulary term. False for free text.
        ///
        /// The CV label is derived from the accession prefix ("MS:1001911" -> "MS"), because SDRF
        /// cells do not carry one separately the way an mzML cvParam does.
        /// </summary>
        public static bool TryParseTerm(string cell, out CvParam term)
        {
            term = null;
            var pairs = ParseKeyValues(cell);
            if (pairs.Count == 0) return false;

            pairs.TryGetValue("NT", out string name);
            pairs.TryGetValue("AC", out string accession);
            if (string.IsNullOrEmpty(name) && string.IsNullOrEmpty(accession)) return false;

            accession ??= "";
            int colon = accession.IndexOf(':');
            string cvLabel = colon > 0 ? accession.Substring(0, colon) : "";

            term = new CvParam(cvLabel, accession, name ?? "", "");
            return true;
        }
    }
}
