using MzLibUtil;
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Text;
using System.Threading.Tasks;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// Retrieves protein sequence data from UniProt's REST API (https://rest.uniprot.org) and writes it to
    /// disk for <see cref="ProteinDbLoader"/> to read: a whole proteome (<see cref="RetrieveProteome"/>), a
    /// single UniProtKB entry (<see cref="RetrieveEntry"/>), or the catalogue of proteomes UniProt offers
    /// (<see cref="DownloadAvailableUniProtProteomes"/>, or <see cref="SearchUniProtProteomes"/> to find one
    /// organism's without downloading all of them).
    /// </summary>
    /// <remarks>
    /// <para>
    /// The usual sequence for "give me the complete database for this organism" is to find the proteome ID
    /// with <see cref="SearchUniProtProteomes"/> and then download it with <see cref="RetrieveProteome"/>
    /// asking for <see cref="Reviewed.all"/>, in <see cref="ProteomeFormat.fasta"/> or
    /// <see cref="ProteomeFormat.xml"/>. The resulting file is what
    /// <see cref="ProteinDbLoader.LoadProteinFasta"/> and <see cref="ProteinDbLoader.LoadProteinXML"/> read.
    /// </para>
    /// <para>
    /// Failures are told apart by exception type, because they call for opposite responses. Following the
    /// contract <see cref="PrideArchiveClient"/> already uses for the same problem:
    /// </para>
    /// <list type="bullet">
    /// <item><description>
    /// <see cref="ArgumentException"/>, <see cref="DirectoryNotFoundException"/>,
    /// <see cref="FileNotFoundException"/> — the call itself is wrong: a blank identifier, an unsupported
    /// format, a storage directory that does not exist. <b>Fix the call.</b> No request is made.
    /// </description></item>
    /// <item><description>
    /// <see cref="MzLibException"/> — UniProt answered, and there is no such proteome or entry (or it
    /// answered with something the contract forbids). <b>Fix the identifier.</b> This is deliberately NOT
    /// an <see cref="HttpRequestException"/>: that type means "the service is unavailable" and
    /// <c>ExternalServiceTestHelper</c> turns it into a skipped test, which would let a withdrawn accession
    /// — or a URL that has quietly started returning 404 — pass unnoticed.
    /// </description></item>
    /// <item><description>
    /// <see cref="HttpRequestException"/> — UniProt was unreachable, timed out, rate-limited the caller, or
    /// returned 5xx. <b>Nobody's mistake; try again later.</b>
    /// </description></item>
    /// </list>
    /// <para>
    /// Every one of these used to be reported as <c>null</c>, so a caller could not tell a mistake from an
    /// outage, and the costliest guess — reading an outage as a bad proteome ID — sent users "correcting" an
    /// identifier that was right all along. Worse, a UniProt error page was written to disk and its path
    /// returned as a success; that is how the XML branch of <see cref="RetrieveProteome"/> came to save a
    /// 153-byte "404 Not Found" page under a <c>.xml</c> name with nothing reporting a problem.
    /// </para>
    /// <para>
    /// <see cref="RetrieveProteome"/> and <see cref="DownloadAvailableUniProtProteomes"/> stream to a
    /// sibling ".partial" file and move it into place only on success, so an interrupted transfer never
    /// leaves a truncated file at the destination path; re-running a retrieval overwrites rather than
    /// throwing on an existing file. <see cref="RetrieveEntry"/> instead buffers the whole (small) entry in
    /// memory, because it has to inspect the payload before deciding whether it is worth saving at all.
    /// Either way these methods do not use <see cref="Loaders.DownloadContent"/>, which ignores the response
    /// status and writes whatever came back.
    /// </para>
    /// <para>
    /// One place this class deliberately parts company with <see cref="PrideArchiveClient"/>: a storage
    /// directory that does not exist is a <see cref="DirectoryNotFoundException"/> here, where
    /// <see cref="PrideArchiveClient.DownloadFileAsync"/> creates it. A caller naming a directory that is
    /// not there has usually misconfigured a path, and silently creating it buries that — which is the
    /// distinction this class exists to preserve.
    /// </para>
    /// </remarks>
    public static class ProteinDbRetriever
    {
        /// <summary>The base address of the UniProt REST API.</summary>
        public const string UniProtRestBaseAddress = "https://rest.uniprot.org/";

        /// <summary>
        /// The first column heading of UniProt's tab-separated proteome catalogue. Recognised so the header
        /// line can be told from a data line, rather than being returned as a proteome whose ID is
        /// "Proteome Id" and whose organism is "Organism".
        /// </summary>
        private const string ProteomeCatalogueIdHeader = "Proteome Id";

        /// <summary>
        /// One reused client for every retrieval. The timeout is generous because a proteome download is
        /// large and slow; exceeding it is reported as an <see cref="HttpRequestException"/> like any other
        /// availability failure (see <see cref="Get"/>).
        /// </summary>
        private static readonly HttpClient SharedHttpClient = new() { Timeout = TimeSpan.FromMinutes(10) };

        /// <summary>
        /// Downloads a UniProt proteome — every protein belonging to one organism's proteome ID — and
        /// returns the full path of the file written.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <see cref="Reviewed"/> selects a subset rather than a sort order: <see cref="Reviewed.yes"/>
        /// asks for the Swiss-Prot (reviewed) entries only, <see cref="Reviewed.no"/> for the TrEMBL
        /// (unreviewed) entries only, and <see cref="Reviewed.all"/> for the complete proteome — which is
        /// why the file is named "_reviewed", "_unreviewed" or "_all". A proteome that has none of the kind
        /// asked for is an <see cref="MzLibException"/>, not a zero-byte file reported as a success.
        /// </para>
        /// <para>
        /// This costs two requests: one to ask how many entries match, and one to download them. The count
        /// is what makes "no such proteome" reportable, and it is asked first so that a mistaken ID costs
        /// nothing. The download uses UniProt's "stream" endpoint, which returns the whole result set —
        /// "search" would return only its first page of 25.
        /// </para>
        /// <para>
        /// A complete proteome is a large download and this method is synchronous: the human proteome as
        /// gzipped xml is around 400 MB and takes minutes. There is no progress reporting or cancellation,
        /// so a caller with a user interface should run it on a background thread and say that it will take
        /// a while.
        /// </para>
        /// <para>
        /// One file name now means something different: asking for fasta isoforms has always produced
        /// "..._reviewed_isoform.fasta", but until the parameter spelling was fixed that file did not
        /// actually hold the isoforms. Callers caching by path will find the same name now carries more
        /// sequences than it used to.
        /// </para>
        /// </remarks>
        /// <param name="proteomeID">A UniProt proteome ID, e.g. "UP000005640" (Homo sapiens).</param>
        /// <param name="absolutePathToStorageDirectory">An existing directory to write the file into.</param>
        /// <param name="format">The download format. Only <see cref="ProteomeFormat.fasta"/> and <see cref="ProteomeFormat.xml"/> are supported.</param>
        /// <param name="reviewed">
        /// Whether to take the reviewed (Swiss-Prot) entries, the unreviewed (TrEMBL) entries, or
        /// <see cref="Reviewed.all"/> of them — the last being the complete organism database.
        /// </param>
        /// <param name="compress">
        /// Whether to save the download gzipped, under a ".gz" name. Note that
        /// <see cref="ProteinDbLoader"/> decompresses a ".gz" database to a fixed "temp" name beside it, so
        /// two gzipped databases in one directory cannot be loaded at the same time; download them to
        /// separate directories, or uncompressed, if they are to be read concurrently.
        /// </param>
        /// <param name="include">
        /// Whether to request extra isoform sequences. They are only delivered where they exist: in fasta,
        /// and among reviewed entries — asking for them alongside <see cref="Reviewed.no"/> is ignored, and
        /// the file is not named "_isoform" when it could not contain any.
        /// </param>
        /// <returns>The full path of the file written. Never null.</returns>
        /// <exception cref="ArgumentException">The proteome ID or storage directory is blank, or the proteome ID cannot be used as a file name.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The format is not fasta or xml, or the review status is not a defined <see cref="Reviewed"/> value.</exception>
        /// <exception cref="DirectoryNotFoundException">The storage directory does not exist.</exception>
        /// <exception cref="MzLibException">UniProt answered, but no protein matched — an unknown proteome ID, or none of the requested review status — or it rejected the request, or returned an empty body for a non-empty proteome.</exception>
        /// <exception cref="HttpRequestException">UniProt was unreachable, timed out, rate-limited the caller, or returned 5xx.</exception>
        public static string RetrieveProteome(string proteomeID, string absolutePathToStorageDirectory, ProteomeFormat format,
            Reviewed reviewed, Compress compress, IncludeIsoforms include) =>
            RetrieveProteome(proteomeID, absolutePathToStorageDirectory, format, reviewed, compress, include, SharedHttpClient);

        /// <summary>
        /// <see cref="RetrieveProteome(string, string, ProteomeFormat, Reviewed, Compress, IncludeIsoforms)"/>
        /// over a caller-supplied <see cref="HttpClient"/>, so tests can drive the response classification
        /// without a live service.
        /// </summary>
        internal static string RetrieveProteome(string proteomeID, string absolutePathToStorageDirectory, ProteomeFormat format,
            Reviewed reviewed, Compress compress, IncludeIsoforms include, HttpClient httpClient)
        {
            if (string.IsNullOrWhiteSpace(proteomeID))
                throw new ArgumentException("A UniProt proteome ID is required, e.g. \"UP000005640\".", nameof(proteomeID));
            if (string.IsNullOrWhiteSpace(absolutePathToStorageDirectory))
                throw new ArgumentException("A storage directory is required.", nameof(absolutePathToStorageDirectory));
            if (format != ProteomeFormat.fasta && format != ProteomeFormat.xml)
                throw new ArgumentOutOfRangeException(nameof(format), format,
                    $"Proteome format '{format}' is not supported; use {nameof(ProteomeFormat.fasta)} or {nameof(ProteomeFormat.xml)}.");
            // Checked here beside the format, not further down: both are "this argument is not a value this
            // method accepts", and they should not be reported in a different order depending on which
            // other argument happens to be wrong as well.
            if (!Enum.IsDefined(typeof(Reviewed), reviewed))
                throw new ArgumentOutOfRangeException(nameof(reviewed), reviewed,
                    $"Review status '{reviewed}' is not a defined {nameof(Reviewed)} value.");

            // The proteome ID becomes part of a file name, so refuse anything that could escape the storage
            // directory or, on NTFS, be read as an alternate data stream (a ':' in the name).
            if (proteomeID.Any(c => Path.GetInvalidFileNameChars().Contains(c)))
                throw new ArgumentException(
                    $"'{proteomeID}' cannot be used as a file name; proteome IDs hold no path or reserved characters.",
                    nameof(proteomeID));

            if (!Directory.Exists(absolutePathToStorageDirectory))
                throw new DirectoryNotFoundException(
                    $"The storage directory '{absolutePathToStorageDirectory}' does not exist; create it before retrieving a proteome.");

            bool compressBool = compress == Compress.yes;

            // Isoforms are only ever delivered when they can actually arrive, so that a file named
            // "_isoform" is never one without any in it. Two conditions rule them out: only fasta carries
            // the extra sequences at all, and isoform ("alternative products") annotation is Swiss-Prot
            // curation, so an unreviewed-only query has none to give. Verified live across three organisms:
            // the unreviewed halves of UP000001450, UP000000803 and UP000000589 return 5043, 18080 and
            // 37597 sequences whether or not isoforms are asked for.
            bool isoformBool = format == ProteomeFormat.fasta
                && include == IncludeIsoforms.yes
                && reviewed != Reviewed.no;

            string reviewedSuffix = reviewed switch
            {
                Reviewed.yes => "_reviewed",
                Reviewed.no => "_unreviewed",
                _ => "_all",
            };

            string filename = proteomeID
                + reviewedSuffix
                + (isoformBool ? "_isoform" : "")
                + "." + format
                + (compressBool ? ".gz" : "");
            string destinationPath = Path.Combine(absolutePathToStorageDirectory, filename);

            // The query both requests share, field-qualified as "(proteome:UP...)". It used to be the bare
            // ID as free text, which happens to match today only because UniProt's default search fields
            // include the proteome ID; naming the field is what the API documents, and it stops an ID that
            // also occurs in some entry's text from dragging in proteins from other organisms.
            //
            // Reviewed.all omits the clause entirely rather than asking for both values, because
            // "reviewed:true OR reviewed:false" is not how UniProt spells "everything" — the absence of the
            // clause is. For Homo sapiens (UP000005640) that is 147,506 entries, exactly the 20,416
            // reviewed plus the 127,090 unreviewed, which is the arithmetic the live test asserts.
            string query = Uri.EscapeDataString("(proteome:" + proteomeID + ")");
            if (reviewed != Reviewed.all)
                query += "+AND+" + Uri.EscapeDataString("(reviewed:" + (reviewed == Reviewed.yes ? "true" : "false") + ")");

            // UniProt's parameter is "includeIsoform"; what this class used to send was "&includeIsoforms:"
            // — a plural name, joined with a colon instead of '=', so it was a parameter name with no value
            // and UniProt ignored it. Isoforms were therefore never included, in any release, even though
            // the file was still named "_isoform". Verified against the live API: streaming P02768 returns
            // one sequence with the old spelling and three with this one. It is only sent when isoforms are
            // wanted, since "false" is already the server's default.
            string isoformFragment = isoformBool ? "&includeIsoform=true" : "";

            // Ask how many entries match BEFORE downloading any of them. "size=0" answers with the
            // X-Total-Results header and an empty body, so this costs one small round trip and settles the
            // question the download itself cannot answer: a proteome ID that matches nothing is a mistake
            // the caller must hear about, not a zero-byte file whose path is returned as a success.
            string countUrl = UniProtRestBaseAddress + "uniprotkb/search?query=" + query + "&format=list&size=0";
            long totalResults = GetMatchCount(httpClient, countUrl, proteomeID);

            if (totalResults == 0)
                throw new MzLibException(reviewed == Reviewed.all
                    ? $"UniProt returned no proteins at all for proteome '{proteomeID}'; the proteome ID may not exist."
                    : $"UniProt returned no {(reviewed == Reviewed.yes ? "reviewed" : "unreviewed")} proteins for proteome '{proteomeID}'. " +
                      "The proteome ID may not exist, or it may have no entries of that review status.");

            // The download itself uses "stream", NOT "search". Both were wrong before: the fasta branch
            // pointed at "/uniprot/search" (a 301 hop to uniprotkb) and the xml branch at "/proteome/search"
            // (which does not exist, so its 404 page was saved as the proteome). But "search" is also PAGED —
            // it answers with the first 25 entries and a "Link: rel=next" cursor — so a human proteome came
            // back as 25 of 20416 proteins with nothing reporting a problem. "stream" returns the complete
            // result set in one response, which is what a caller asking for a proteome means.
            string url = UniProtRestBaseAddress + "uniprotkb/stream?query=" + query
                + "&compressed=" + (compressBool ? "true" : "false")
                + "&format=" + format
                + isoformFragment;

            using HttpResponseMessage response = Get(httpClient, url);
            ThrowIfServiceUnavailable(response, url);

            // The count probe above has already established that the query matches something, so a
            // non-success status here is never "your ID is wrong" — it is a rejected or misdirected request,
            // and must fail loudly rather than be mistaken for an outage and skipped.
            if (!response.IsSuccessStatusCode)
                throw new MzLibException(
                    $"UniProt rejected the proteome request for '{proteomeID}' with status " +
                    $"{(int)response.StatusCode} {response.ReasonPhrase} ('{url}').");

            WriteResponseToFile(response, destinationPath);

            // A 200 can still carry nothing — UniProt is known to answer streaming requests with an error
            // message in the body, and a proxy can truncate one. The count probe said there were entries, so
            // an empty file here is a broken download, not an absence.
            if (new FileInfo(destinationPath).Length == 0)
            {
                File.Delete(destinationPath);
                throw new MzLibException(
                    $"UniProt reported {totalResults} proteins for proteome '{proteomeID}' but returned an empty body ('{url}').");
            }

            return destinationPath;
        }

        /// <summary>
        /// Asks UniProt how many entries a query matches, without downloading any of them. Used to tell "no
        /// such proteome" apart from a successful download before the download is attempted.
        /// </summary>
        /// <exception cref="MzLibException">UniProt rejected the probe, or answered without a usable count.</exception>
        /// <exception cref="HttpRequestException">UniProt was unreachable, timed out, rate-limited the caller, or returned 5xx.</exception>
        private static long GetMatchCount(HttpClient httpClient, string countUrl, string proteomeID)
        {
            using HttpResponseMessage response = Get(httpClient, countUrl);
            ThrowIfServiceUnavailable(response, countUrl);

            if (!response.IsSuccessStatusCode)
                throw new MzLibException(
                    $"UniProt rejected the proteome count request for '{proteomeID}' with status " +
                    $"{(int)response.StatusCode} {response.ReasonPhrase} ('{countUrl}').");

            // Without the header there is no count, and guessing one either way is worse than saying so: a
            // silent "assume it's fine" is how the truncated and empty downloads got through in the first place.
            if (!TryGetTotalResults(response, out long totalResults))
                throw new MzLibException(
                    $"UniProt answered the proteome count request for '{proteomeID}' without an X-Total-Results header ('{countUrl}'); " +
                    "the API contract may have changed.");

            return totalResults;
        }

        /// <summary>
        /// Downloads a single UniProtKB entry — one protein with its annotations — and returns the full path
        /// of the file written. This is the single-accession counterpart to <see cref="RetrieveProteome"/>,
        /// for checking one protein, building a small test database, or answering a question about one
        /// target; <see cref="ProteinDbLoader.LoadProteinXML"/> reads the resulting one-entry XML happily.
        /// </summary>
        /// <param name="accession">
        /// A UniProtKB accession, e.g. "P02768", optionally with an isoform suffix, e.g. "P02768-2".
        /// Matching is case-insensitive at UniProt.
        /// </param>
        /// <param name="absolutePathToStorageDirectory">An existing directory to write the file into.</param>
        /// <param name="format">The download format. Only <see cref="ProteomeFormat.xml"/> (the default) and <see cref="ProteomeFormat.fasta"/> are supported.</param>
        /// <returns>The full path of the file written, named for the accession. Never null.</returns>
        /// <exception cref="ArgumentException">
        /// The accession is blank, holds a character outside letters, digits and '-' (which would also make
        /// it an unsafe file name), the storage directory is blank, or the format is not xml or fasta.
        /// </exception>
        /// <exception cref="DirectoryNotFoundException">The storage directory does not exist.</exception>
        /// <exception cref="MzLibException">
        /// UniProt has no such entry. It answers a malformed accession with <b>400</b> rather than 404, a
        /// well-formed but unknown one with 404, and a deleted one with 200 carrying no entry at all — all
        /// three mean the same thing to a caller and are reported the same way.
        /// </exception>
        /// <exception cref="HttpRequestException">UniProt was unreachable, timed out, rate-limited the caller, or returned 5xx.</exception>
        public static string RetrieveEntry(string accession, string absolutePathToStorageDirectory,
            ProteomeFormat format = ProteomeFormat.xml) =>
            RetrieveEntry(accession, absolutePathToStorageDirectory, format, SharedHttpClient);

        /// <summary>
        /// <see cref="RetrieveEntry(string, string, ProteomeFormat)"/> over a caller-supplied
        /// <see cref="HttpClient"/>, so tests can drive the response classification without a live service.
        /// </summary>
        internal static string RetrieveEntry(string accession, string absolutePathToStorageDirectory,
            ProteomeFormat format, HttpClient httpClient)
        {
            (bool found, string filePath) = TryRetrieveEntry(accession, absolutePathToStorageDirectory, format, httpClient);

            if (!found)
                throw new MzLibException(
                    $"UniProt has no entry with accession '{accession}'; it may be misspelled, or the entry may have been deleted or demerged.");

            return filePath;
        }

        /// <summary>
        /// Attempts to download a single UniProtKB entry, reporting a non-existent accession as a value
        /// rather than an exception — the cheap way to validate an accession a user typed, without catching
        /// an exception for control flow.
        /// </summary>
        /// <remarks>
        /// <c>Found</c> is false for one reason only: UniProt answered, and there is no such entry. That
        /// covers all three ways it says so — <b>400</b> for a malformed accession (it does NOT use 404
        /// here), 404 for a well-formed but unknown one, and 200 carrying a payload with no entry in it for
        /// one that has been deleted. Every other failure still throws, because collapsing an outage into
        /// "no such entry" would send a caller hunting for a typo in a perfectly good accession.
        /// </remarks>
        /// <param name="accession">A UniProtKB accession, e.g. "P02768", optionally with an isoform suffix.</param>
        /// <param name="absolutePathToStorageDirectory">An existing directory to write the file into.</param>
        /// <param name="format">The download format. Only <see cref="ProteomeFormat.xml"/> (the default) and <see cref="ProteomeFormat.fasta"/> are supported.</param>
        /// <returns>
        /// <c>Found</c> and the full path of the file written when the entry exists; <c>false</c> and null
        /// when UniProt reports no such entry, in which case nothing is written. The path is never null when
        /// <c>Found</c> is true.
        /// </returns>
        /// <exception cref="ArgumentException">The accession or storage directory is blank, the accession holds a character outside letters, digits and '-', or the format is not xml or fasta.</exception>
        /// <exception cref="DirectoryNotFoundException">The storage directory does not exist.</exception>
        /// <exception cref="HttpRequestException">UniProt was unreachable, timed out, rate-limited the caller, or returned 5xx.</exception>
        public static (bool Found, string FilePath) TryRetrieveEntry(string accession, string absolutePathToStorageDirectory,
            ProteomeFormat format = ProteomeFormat.xml) =>
            TryRetrieveEntry(accession, absolutePathToStorageDirectory, format, SharedHttpClient);

        /// <summary>
        /// <see cref="TryRetrieveEntry(string, string, ProteomeFormat)"/> over a caller-supplied
        /// <see cref="HttpClient"/>, so tests can drive the response classification without a live service.
        /// </summary>
        internal static (bool Found, string FilePath) TryRetrieveEntry(string accession, string absolutePathToStorageDirectory,
            ProteomeFormat format, HttpClient httpClient)
        {
            if (string.IsNullOrWhiteSpace(accession))
                throw new ArgumentException("A UniProtKB accession is required, e.g. \"P02768\".", nameof(accession));
            if (string.IsNullOrWhiteSpace(absolutePathToStorageDirectory))
                throw new ArgumentException("A storage directory is required.", nameof(absolutePathToStorageDirectory));
            if (format != ProteomeFormat.xml && format != ProteomeFormat.fasta)
                throw new ArgumentOutOfRangeException(nameof(format), format,
                    $"Entry format '{format}' is not supported; use {nameof(ProteomeFormat.xml)} or {nameof(ProteomeFormat.fasta)}.");

            // The accession becomes part of a file name, so refuse anything that could escape the storage
            // directory. This is deliberately a character guard and NOT the full UniProt accession pattern:
            // a well-formed-looking-but-wrong accession should reach UniProt and come back as a definite
            // "no such entry", rather than be second-guessed here against a pattern UniProt may widen.
            if (!accession.All(c => char.IsLetterOrDigit(c) || c == '-'))
                throw new ArgumentException(
                    $"'{accession}' is not a valid UniProtKB accession; accessions hold only letters, digits and '-'.",
                    nameof(accession));

            if (!Directory.Exists(absolutePathToStorageDirectory))
                throw new DirectoryNotFoundException(
                    $"The storage directory '{absolutePathToStorageDirectory}' does not exist; create it before retrieving an entry.");

            string destinationPath = Path.Combine(absolutePathToStorageDirectory, accession + "." + format);
            string url = UniProtRestBaseAddress + "uniprotkb/" + Uri.EscapeDataString(accession) + "." + format;

            using HttpResponseMessage response = Get(httpClient, url);
            ThrowIfServiceUnavailable(response, url);

            // 400 (malformed accession) and 404 (well-formed, unknown) are the two ways UniProt reports "no
            // such entry" by status. They are checked before the general guard below so that every OTHER
            // non-success status still throws rather than being reported as a plain absence.
            if (response.StatusCode is HttpStatusCode.BadRequest or HttpStatusCode.NotFound)
                return (false, null);

            // Anything else non-success is a rejected or misdirected request, not an outage —
            // ThrowIfServiceUnavailable above has already claimed every genuinely transient status. It must
            // be MzLibException for the same reason as everywhere else in this class: HttpRequestException
            // means "unavailable", and a live test skips on it, so a 401/403/410 here would silently stop
            // testing the contract instead of reporting that it broke.
            if (!response.IsSuccessStatusCode)
                throw new MzLibException(
                    $"UniProt rejected the entry request for '{accession}' with status " +
                    $"{(int)response.StatusCode} {response.ReasonPhrase} ('{url}').");

            // The third way: a deleted accession is answered with 200 and a payload holding no entry at all —
            // an XML document with only the <uniprot> wrapper and its copyright notice, or an empty fasta.
            // Written out as-is that file loads as zero proteins, so hold the body (one entry is small) and
            // refuse to save an empty answer. The BYTES are what gets written — the text is only read to
            // inspect it, so the file on disk is UniProt's own, byte for byte, whatever its encoding.
            byte[] body = ReadBody(response);
            string text = Encoding.UTF8.GetString(body);
            bool holdsEntry = format == ProteomeFormat.xml
                ? text.Contains("<entry", StringComparison.Ordinal)
                : text.TrimStart().StartsWith(">", StringComparison.Ordinal);
            if (!holdsEntry)
                return (false, null);

            File.WriteAllBytes(destinationPath, body);
            return (true, destinationPath);
        }

        /// <summary>
        /// Downloads and then returns the filepath to a compressed (.gz), tab-delimited text file of the
        /// available proteomes. Line one is the header; the columns are Proteome Id, Organism, Organism Id
        /// and Protein count.
        /// </summary>
        /// <remarks>
        /// <para>
        /// This is how a caller who knows an organism but not its proteome ID finds one: download the
        /// catalogue, read it with <see cref="UniprotProteomesList"/>, and pass the ID to
        /// <see cref="RetrieveProteome(string, string, ProteomeFormat, Reviewed, Compress, IncludeIsoforms)"/>.
        /// </para>
        /// <para>
        /// It is the complete catalogue, and it is large: over a million proteomes, some 50 MB of text, and
        /// minutes to transfer. Most of that is redundant strain-level entries. For anything interactive,
        /// <see cref="SearchUniProtProteomes(string)"/> answers the same question about one organism in a
        /// single small request and is the better choice; this method is for taking a copy of the whole
        /// list. UniProt's curated subset can be had by searching "(reference:true)" instead.
        /// </para>
        /// </remarks>
        /// <param name="destinationFolder">An existing directory to write the file into.</param>
        /// <returns>The full path of the file written. Never null.</returns>
        /// <exception cref="ArgumentException">The destination folder is blank.</exception>
        /// <exception cref="DirectoryNotFoundException">The destination folder does not exist.</exception>
        /// <exception cref="MzLibException">UniProt rejected the request or returned an empty catalogue.</exception>
        /// <exception cref="HttpRequestException">UniProt was unreachable, timed out, rate-limited the caller, or returned 5xx.</exception>
        public static string DownloadAvailableUniProtProteomes(string destinationFolder) =>
            DownloadAvailableUniProtProteomes(destinationFolder, SharedHttpClient);

        /// <summary>
        /// <see cref="DownloadAvailableUniProtProteomes(string)"/> over a caller-supplied
        /// <see cref="HttpClient"/>, so tests can drive the response classification without a live service.
        /// </summary>
        internal static string DownloadAvailableUniProtProteomes(string destinationFolder, HttpClient httpClient)
        {
            if (string.IsNullOrWhiteSpace(destinationFolder))
                throw new ArgumentException("A destination folder is required.", nameof(destinationFolder));
            if (!Directory.Exists(destinationFolder))
                throw new DirectoryNotFoundException(
                    $"The destination folder '{destinationFolder}' does not exist; create it before downloading the proteome list.");

            // "stream", not "search", for the same reason the proteome download uses it: "search" is paged,
            // so this returned the first 25 of roughly 100,000 proteomes and reported no problem. That is
            // why looking Homo sapiens up in the catalogue could not find it — UP000005640 was never in the
            // 25 rows that arrived.
            string url = UniProtRestBaseAddress + "proteomes/stream?query=*&format=tsv&compressed=true";
            string filepath = Path.Combine(destinationFolder, "availableUniProtProteomes.txt.gz");

            using HttpResponseMessage response = Get(httpClient, url);
            ThrowIfServiceUnavailable(response, url);

            if (!response.IsSuccessStatusCode)
                throw new MzLibException(
                    $"UniProt rejected the proteome-list request with status {(int)response.StatusCode} " +
                    $"{response.ReasonPhrase} ('{url}').");

            WriteResponseToFile(response, filepath);

            // "stream" sends no X-Total-Results header, so an empty catalogue cannot be spotted from the
            // headers the way RetrieveProteome spots an empty proteome; and because the body is gzipped, an
            // empty one is not a zero-byte file either. Read it back instead: a usable catalogue is the
            // header line plus at least one proteome.
            try
            {
                if (ReadAllGZippedLines(filepath).Take(2).Count() < 2)
                {
                    File.Delete(filepath);
                    throw new MzLibException($"UniProt returned an empty list of available proteomes ('{url}').");
                }
            }
            catch (InvalidDataException e)
            {
                File.Delete(filepath);
                throw new MzLibException(
                    $"UniProt answered the proteome-list request with a body that is not gzip ('{url}'), " +
                    "though gzip was requested; the API contract may have changed.", e);
            }

            return filepath;
        }

        /// <summary>
        /// Asks UniProt which proteomes match a search term — an organism name, a taxonomy ID, or a
        /// proteome ID — and returns what it knows about each.
        /// </summary>
        /// <remarks>
        /// <para>
        /// This is the short way round the problem that
        /// <see cref="RetrieveProteome(string, string, ProteomeFormat, Reviewed, Compress, IncludeIsoforms)"/>
        /// is keyed by proteome ID while people know organisms by name. Searching "Homo sapiens" returns
        /// UP000005640 among its hits, which is the argument the download actually wants; the alternative is
        /// downloading the whole catalogue with <see cref="DownloadAvailableUniProtProteomes(string)"/> —
        /// over a million rows — to look up one of them.
        /// </para>
        /// <para>
        /// The result carries the protein count and taxonomy ID, not just the name, because the name is not
        /// enough to choose by: UniProt lists <b>two</b> proteomes for Homo sapiens, UP000005640 and
        /// UP001055169, with the same organism string and the same taxonomy ID 9606, differing only in that
        /// one holds 147,506 proteins and the other 54,825. A caller picking by name alone picks at random
        /// between a complete human database and a third of one, with nothing to say which it got. That is
        /// why this returns a list of records rather than a name-keyed dictionary, and why it is ordered:
        /// duplicates are real and must survive to the caller.
        /// </para>
        /// <para>
        /// Matching is UniProt's, not this method's: a bare term is a free-text search over the proteomes
        /// collection, and a field-qualified one such as "(organism_id:9606)" or "(upid:UP000005640)" is
        /// passed through as written. Free text matches broadly — "human" returns close to 16,000 proteomes
        /// whose description merely mentions it, against 44 for "Homo sapiens" and 2 for "(organism_id:9606)"
        /// — so a field-qualified term is worth preferring.
        /// </para>
        /// <para>
        /// Like the download, this uses "stream" rather than the paged "search", so the answer is complete
        /// and is read into memory whole. That is a few hundred kilobytes for an organism and about 50 MB
        /// for "*", which is a reason to give it something specific to look for.
        /// </para>
        /// </remarks>
        /// <param name="searchTerm">An organism name, taxonomy ID, proteome ID, or UniProt query expression.</param>
        /// <returns>
        /// The matching proteomes, in the order UniProt returned them — which puts the reference proteome
        /// first for an organism query. Empty when nothing matches: an organism UniProt has no proteome for
        /// is an answer, not a fault. Never null.
        /// </returns>
        /// <exception cref="ArgumentException">The search term is blank.</exception>
        /// <exception cref="MzLibException">UniProt rejected the query, or answered in a shape this cannot read.</exception>
        /// <exception cref="HttpRequestException">UniProt was unreachable, timed out, rate-limited the caller, or returned 5xx.</exception>
        public static IReadOnlyList<UniProtProteome> SearchUniProtProteomes(string searchTerm) =>
            SearchUniProtProteomes(searchTerm, SharedHttpClient);

        /// <summary>
        /// <see cref="SearchUniProtProteomes(string)"/> over a caller-supplied <see cref="HttpClient"/>, so
        /// tests can drive the response classification without a live service.
        /// </summary>
        internal static IReadOnlyList<UniProtProteome> SearchUniProtProteomes(string searchTerm, HttpClient httpClient)
        {
            if (string.IsNullOrWhiteSpace(searchTerm))
                throw new ArgumentException(
                    "A search term is required, e.g. \"Homo sapiens\" or \"(organism_id:9606)\".", nameof(searchTerm));

            // The columns are named rather than left to the server default, so the meaning of each position
            // is this method's contract and not something UniProt can quietly reorder underneath it. A field
            // name UniProt does not know is a 400, so a rename would be reported rather than mis-parsed.
            string url = UniProtRestBaseAddress + "proteomes/stream?query=" + Uri.EscapeDataString(searchTerm)
                + "&format=tsv&compressed=false&fields=upid,organism,organism_id,protein_count";

            using HttpResponseMessage response = Get(httpClient, url);
            ThrowIfServiceUnavailable(response, url);

            // UniProt answers an unparseable query expression with 400. That is a mistake in the call, but
            // it is only discoverable from UniProt's reply, so it arrives here as MzLibException like every
            // other "UniProt answered and said no" in this class rather than as an argument exception.
            if (!response.IsSuccessStatusCode)
                throw new MzLibException(
                    $"UniProt rejected the proteome search for '{searchTerm}' with status " +
                    $"{(int)response.StatusCode} {response.ReasonPhrase} ('{url}').");

            string body = Encoding.UTF8.GetString(ReadBody(response));
            List<UniProtProteome> proteomes = new();

            foreach (string line in body.Split('\n'))
            {
                string[] columns = line.TrimEnd('\r').Split('\t');
                if (columns.Length < 2 || columns[0] == ProteomeCatalogueIdHeader)
                    continue;

                // The two numbers are parsed leniently and default to 0. They are supporting detail — the
                // ID and the organism are what a caller cannot proceed without — and losing the whole row
                // because UniProt sent "" for a protein count would be the worse failure.
                int organismId = columns.Length > 2 && int.TryParse(columns[2], out int parsedOrganismId) ? parsedOrganismId : 0;
                int proteinCount = columns.Length > 3 && int.TryParse(columns[3], out int parsedCount) ? parsedCount : 0;

                proteomes.Add(new UniProtProteome(columns[0], columns[1], organismId, proteinCount));
            }

            return proteomes;
        }

        /// <summary>
        /// One proteome UniProt offers, as <see cref="SearchUniProtProteomes(string)"/> reports it.
        /// </summary>
        /// <param name="ProteomeId">
        /// The UniProt proteome ID, e.g. "UP000005640" — the argument
        /// <see cref="RetrieveProteome(string, string, ProteomeFormat, Reviewed, Compress, IncludeIsoforms)"/> takes.
        /// </param>
        /// <param name="Organism">The organism name, e.g. "Homo sapiens (Human)". Not unique: an organism may have several proteomes.</param>
        /// <param name="OrganismId">The NCBI taxonomy ID, e.g. 9606. Also not unique across proteomes. 0 if UniProt did not supply it.</param>
        /// <param name="ProteinCount">
        /// How many proteins the proteome holds. Usually the only thing telling two proteomes of the same
        /// organism apart, so it is what a caller choosing between them should choose on. 0 if UniProt did
        /// not supply it.
        /// </param>
        public record UniProtProteome(string ProteomeId, string Organism, int OrganismId, int ProteinCount);

        /// <summary>
        /// Reads a downloaded UniProt proteome catalogue into a dictionary keyed by Proteome ID, valued by
        /// organism. Accepts the file gzipped (".gz"), zipped (".zip") or plain (".txt").
        /// </summary>
        /// <param name="completePathToAvailableUniProtProteomes">Path to a file written by <see cref="DownloadAvailableUniProtProteomes(string)"/>.</param>
        /// <returns>Proteome ID to organism. Never null.</returns>
        /// <exception cref="ArgumentException">The path is blank, or names a file whose extension is not .gz, .zip or .txt.</exception>
        /// <exception cref="FileNotFoundException">No file exists at that path.</exception>
        public static Dictionary<string, string> UniprotProteomesList(string completePathToAvailableUniProtProteomes)
        {
            if (string.IsNullOrWhiteSpace(completePathToAvailableUniProtProteomes))
                throw new ArgumentException("A path to a downloaded proteome list is required.",
                    nameof(completePathToAvailableUniProtProteomes));
            if (!File.Exists(completePathToAvailableUniProtProteomes))
                throw new FileNotFoundException(
                    "No downloaded UniProt proteome list exists at that path.", completePathToAvailableUniProtProteomes);

            string fileExtension = Path.GetExtension(completePathToAvailableUniProtProteomes);
            IEnumerable<string> lines = fileExtension switch
            {
                ".gz" => ReadAllGZippedLines(completePathToAvailableUniProtProteomes),
                ".zip" => ReadAllZippedLines(completePathToAvailableUniProtProteomes),
                ".txt" => File.ReadAllLines(completePathToAvailableUniProtProteomes),
                _ => throw new ArgumentException(
                    $"'{fileExtension}' is not a proteome-list extension; expected .gz, .zip or .txt.",
                    nameof(completePathToAvailableUniProtProteomes)),
            };

            Dictionary<string, string> dictionaryOfAvailableProteomes = new();
            foreach (string item in lines)
            {
                // A line with no tab (the trailing blank one, or a truncated download) used to raise
                // IndexOutOfRangeException, and a repeated proteome ID an ArgumentException that read
                // exactly like this method's own "you passed a bad path" — sending the caller to look in
                // the wrong place. Skip the first, let the last value win for the second.
                var lineValuesArray = item.Split("\t");
                if (lineValuesArray.Length < 2)
                    continue;

                // The catalogue's first line is its header, which otherwise arrived here as a proteome
                // whose ID was "Proteome Id" — harmless to iterate past, but it is offered to users as a
                // list of organisms to choose from, and it is not one.
                if (lineValuesArray[0] == ProteomeCatalogueIdHeader)
                    continue;

                dictionaryOfAvailableProteomes[lineValuesArray[0]] = lineValuesArray[1];
            }
            return dictionaryOfAvailableProteomes;
        }

        /// <summary>
        /// One can retrieve a table of information about specific proteins.
        /// This method returns the names of columns that could be included in the table
        /// </summary>
        /// <returns></returns>
        public static Dictionary<string, string> UniprotColumnsList()
        {
            string currentDirectory = Directory.GetCurrentDirectory();
            string filePath = Path.Combine(currentDirectory, "UniProtKB_columnNamesForProgrammaticAccess.txt");
            Dictionary<string, string> d = new();
            string[] idNameList = File.ReadAllLines(filePath);
            foreach (string item in idNameList)
            {
                if (!item.StartsWith('#'))
                {
                    var j = item.Split("\t");
                    d.Add(j[0], j[1]);
                }
            }
            return d;
        }

        /// <summary>
        /// Issues the request, reporting a client-side timeout as the same transient type as any other
        /// availability failure so callers have ONE exception to catch for "the service is having a bad day".
        /// </summary>
        /// <remarks>
        /// This class's public API is synchronous, so the async call has to be bridged. The bridge is
        /// <c>Task.Run(...).GetAwaiter().GetResult()</c> rather than a bare <c>.GetAwaiter().GetResult()</c>
        /// or <c>.Result</c>, matching the reasoning written out at
        /// <c>PredictionClients\Koina\AbstractClasses\FragmentIntensityModel.cs</c>: Task.Run schedules the
        /// work on a ThreadPool thread carrying no SynchronizationContext, so awaited continuations resume
        /// freely instead of deadlocking on a single-threaded context. MetaMorpheus's WPF "Download UniProt
        /// Database" window is exactly that kind of caller — it currently hand-rolls its own download, and
        /// adopting this API must not deadlock it.
        /// </remarks>
        private static HttpResponseMessage Get(HttpClient httpClient, string url)
        {
            try
            {
                return Task.Run(() => httpClient.GetAsync(url, HttpCompletionOption.ResponseHeadersRead))
                    .GetAwaiter().GetResult();
            }
            catch (TaskCanceledException e)
            {
                throw new HttpRequestException($"The UniProt request to '{url}' timed out after {httpClient.Timeout}.", e);
            }
        }

        /// <summary>Reads the whole response body, bridged off the caller's context — see <see cref="Get"/>.</summary>
        /// <remarks>
        /// The body arrives after the headers <see cref="Get"/> has already classified, so a stall or a
        /// dropped connection surfaces here rather than at the request. It has to reach the caller as
        /// <see cref="HttpRequestException"/> for the same reason <c>WriteResponseToFile</c> translates it:
        /// the documented contract is that one type covers "try again later", and a raw
        /// <see cref="IOException"/> escaping from here would be a hole in it.
        /// </remarks>
        private static byte[] ReadBody(HttpResponseMessage response)
        {
            try
            {
                return Task.Run(() => response.Content.ReadAsByteArrayAsync()).GetAwaiter().GetResult();
            }
            catch (Exception e) when (e is IOException or OperationCanceledException)
            {
                throw new HttpRequestException("The UniProt download failed while reading the response body.", e);
            }
        }

        /// <summary>
        /// Throws <see cref="HttpRequestException"/> if the response says the service — rather than the
        /// request — is at fault: a timeout (408), rate limiting (429), or any server error (5xx). This is
        /// the ONLY classification that means "try again later", and the only one a live test may skip on.
        /// </summary>
        /// <remarks>
        /// Named to stay distinct from <c>Test.ExternalServiceTestHelper.ThrowIfUnavailable</c>, which makes
        /// the same 408/429/5xx judgement but throws its own marker type for NUnit. Both exist because the
        /// test helper cannot reach into production code; keep the two in step if either changes.
        /// </remarks>
        internal static void ThrowIfServiceUnavailable(HttpResponseMessage response, string url)
        {
            int status = (int)response.StatusCode;
            if (status is 408 or 429 || status >= 500)
                throw new HttpRequestException(
                    $"UniProt is unavailable: HTTP {status} {response.ReasonPhrase} for '{url}'.");
        }

        /// <summary>
        /// Reads UniProt's "X-Total-Results" response header — how many records the query matched, sent
        /// alongside the first page. Returns false when the header is absent or unparseable, in which case
        /// the caller must not conclude anything about the result count.
        /// </summary>
        internal static bool TryGetTotalResults(HttpResponseMessage response, out long totalResults)
        {
            totalResults = 0;
            if (!response.Headers.TryGetValues("X-Total-Results", out IEnumerable<string> values))
                return false;

            foreach (string value in values)
            {
                if (long.TryParse(value, out totalResults))
                    return true;
            }
            return false;
        }

        /// <summary>
        /// Streams a response body to <paramref name="destinationPath"/> through a sibling ".partial" file
        /// that is moved into place only on success, so an interrupted transfer never leaves a truncated
        /// file where a complete one is expected.
        /// </summary>
        private static void WriteResponseToFile(HttpResponseMessage response, string destinationPath)
        {
            // The scratch name carries a unique token rather than being derived from the destination alone.
            // Two concurrent retrievals of the same proteome — this is a static class, reachable from a
            // Parallel.For — would otherwise collide on one scratch file and delete each other's download.
            string partialPath = destinationPath + "." + Guid.NewGuid().ToString("N") + ".partial";
            try
            {
                try
                {
                    using var fileStream = new FileStream(partialPath, FileMode.Create, FileAccess.Write, FileShare.None);
                    using Stream httpStream = response.Content.ReadAsStream();
                    httpStream.CopyTo(fileStream);
                }
                catch (Exception e) when (e is IOException or OperationCanceledException)
                {
                    // The body arrives AFTER the headers this class already classified, so a stall or a
                    // dropped connection surfaces here rather than at the request. It means the same thing as
                    // a 504 and must reach the caller as the same type, or the documented contract — catch
                    // HttpRequestException for "try again later" — silently fails to cover the phase where a
                    // large download actually breaks.
                    throw new HttpRequestException(
                        $"The UniProt download failed while reading the response body into '{destinationPath}'.", e);
                }

                File.Move(partialPath, destinationPath, overwrite: true);
            }
            finally
            {
                // Cleanup must never replace the exception that caused it: on Windows a locked scratch file
                // makes Delete throw, and an exception from a finally block discards the real one.
                try
                {
                    if (File.Exists(partialPath))
                        File.Delete(partialPath);
                }
                catch (IOException) { }
                catch (UnauthorizedAccessException) { }
            }
        }

        /// <summary>
        /// This is an enumerated list of a file types available for download at UniProt
        /// </summary>
        public enum ProteomeFormat
        {
            html,
            tab,
            xls,
            fasta,
            gff,
            txt,
            xml,
            rdf,
            list,
            rss
        }

        /// <summary>
        /// Which review status of protein to take from a proteome. This is a filter, not a flag: the two
        /// original members are mutually exclusive halves, so neither of them downloads a whole proteome.
        /// </summary>
        public enum Reviewed
        {
            /// <summary>Swiss-Prot (manually reviewed) entries only.</summary>
            yes,

            /// <summary>TrEMBL (unreviewed) entries only.</summary>
            no,

            /// <summary>
            /// Every entry in the proteome, Swiss-Prot and TrEMBL together — the complete organism
            /// database. Added because the two members above cannot express it between them: they are
            /// alternatives, and there was previously no way to ask for the whole thing at once.
            /// </summary>
            all
        }

        /// <summary>
        /// Compress will cause the table to be downloaded in .gz format
        /// </summary>
        public enum Compress
        {
            yes,
            no
        }

        /// <summary>
        /// Include isoform sequences when the format parameter is set to fasta.
        /// Include description of referenced data when the format parameter is set to rdf.
        /// </summary>
        public enum IncludeIsoforms
        {
            yes,
            no
        }

        /// <summary>
        /// Columns to select for retrieving results in tab or xls format.
        /// https://www.uniprot.org/help/return_fields
        /// </summary>
        public enum Columns
        {
        }

        public static IEnumerable<string> ReadAllGZippedLines(string filename)
        {
            using (var fileStream = File.OpenRead(filename))
            {
                using (var gzipStream = new GZipStream(fileStream, CompressionMode.Decompress))
                {
                    using (var reader = new StreamReader(gzipStream))
                    {
                        string currentLine;
                        while ((currentLine = reader.ReadLine()) != null)
                        {
                            yield return currentLine;
                        }
                    }
                }
            }
        }

        public static IEnumerable<string> ReadAllZippedLines(string filename)
        {
            using (ZipArchive archive = ZipFile.OpenRead(filename))
            {
                foreach (ZipArchiveEntry entry in archive.Entries)
                {
                    using (var reader = new StreamReader(entry.Open()))
                    {
                        string currentLine;
                        while ((currentLine = reader.ReadLine()) != null)
                        {
                            yield return currentLine;
                        }
                    }
                }
            }
        }
    }
}
