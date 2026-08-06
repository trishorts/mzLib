using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for the SDRF-Proteomics reader/writer.
    ///
    /// Each fixture here pins a property measured across the 1,236-file curated corpus at
    /// bigbio/sdrf-annotated-datasets rather than a property of the specification, because the two
    /// disagree in ways that would break a reader written from the spec alone:
    ///
    ///   PXD000070 -- typical file; repeats comment[modification parameters] 8 times, and uses a
    ///                PRIDE CV term (not PSI-MS) for the acquisition method.
    ///   PXD026824 -- contains literal double-quote characters inside cell values, which a
    ///                CSV parser with default quoting would strip.
    ///   PXD059974 -- ragged: a 46-column header with 17 of its 23 rows carrying only 42 cells.
    ///
    /// The load-bearing test is <see cref="RoundTrip_IsByteIdentical"/>. A reader that normalizes
    /// anything -- line endings, short rows, quoting, header case -- fails it, and normalizing a
    /// curated annotation silently is the failure mode this whole component exists to avoid.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrf
    {
        private static string FixtureDir => Path.Combine(
            TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes");

        private static string Fixture(string name) => Path.Combine(FixtureDir, name);

        private static IEnumerable<string> AllFixtures()
        {
            yield return "PXD000070.sdrf.tsv";
            yield return "PXD026824.sdrf.tsv";
            yield return "PXD059974.sdrf.tsv";
        }

        [Test]
        [TestCaseSource(nameof(AllFixtures))]
        public void RoundTrip_IsByteIdentical(string fixtureName)
        {
            string source = Fixture(fixtureName);
            var document = new SdrfDocument(source);

            string outputPath = Path.Combine(
                TestContext.CurrentContext.TestDirectory, $"roundtrip_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                document.WriteResults(outputPath);

                var expected = File.ReadAllBytes(source);
                var actual = File.ReadAllBytes(outputPath);

                Assert.That(actual, Is.EqualTo(expected),
                    $"round trip of {fixtureName} was not byte-identical");
            }
            finally
            {
                if (File.Exists(outputPath))
                    File.Delete(outputPath);
            }
        }

        [Test]
        public void Read_ParsesHeaderAndRows()
        {
            var document = new SdrfDocument(Fixture("PXD000070.sdrf.tsv"));

            Assert.That(document.Header[0], Is.EqualTo("source name"));
            Assert.That(document.Header.Contains("characteristics[organism]"), Is.True);
            Assert.That(document.Results, Is.Not.Empty);
            Assert.That(document.Results.Count, Is.EqualTo(6));
        }

        [Test]
        public void RepeatedColumns_AreAllPreserved()
        {
            // The house [Name(...)] pattern cannot express this: the column name repeats, and each
            // occurrence carries a different modification.
            var document = new SdrfDocument(Fixture("PXD000070.sdrf.tsv"));

            var indexes = document.Header.IndexesOf("comment[modification parameters]");
            Assert.That(indexes.Count, Is.EqualTo(8));

            var mods = document.Results[0].All("comment[modification parameters]");
            Assert.That(mods.Count, Is.EqualTo(8));
            Assert.That(mods[0], Does.StartWith("NT=Carbamidomethyl"));
            Assert.That(mods.Any(m => m.Contains("UNIMOD:35")), Is.True);

            // The single-value indexer returns the FIRST occurrence, not a join of them.
            Assert.That(document.Results[0]["comment[modification parameters]"],
                Is.EqualTo(mods[0]));
        }

        [Test]
        public void LiteralQuotes_AreNotTreatedAsCsvQuoting()
        {
            // SDRF defines no quoting or escaping. These double quotes are data.
            var document = new SdrfDocument(Fixture("PXD026824.sdrf.tsv"));

            bool anyQuoted = document.Results
                .SelectMany(r => r.Cells)
                .Any(c => c.StartsWith("\"") && c.EndsWith("\"") && c.Length > 1);

            Assert.That(anyQuoted, Is.True,
                "expected at least one cell whose literal double quotes survived parsing");
        }

        [Test]
        public void ShortRows_ArePreservedNotPadded()
        {
            // PXD059974 is genuinely malformed. The reader preserves it so the file round-trips;
            // reporting it is the validator's job (PR-3).
            var document = new SdrfDocument(Fixture("PXD059974.sdrf.tsv"));

            Assert.That(document.Header.Count, Is.EqualTo(46));

            var widths = document.Results.Select(r => r.Cells.Count).Distinct().OrderBy(x => x).ToList();
            Assert.That(widths, Is.EqualTo(new[] { 42, 46 }));

            // Reaching past the end of a short row yields null rather than throwing.
            var shortRow = document.Results.First(r => r.Cells.Count == 42);
            Assert.That(shortRow["comment[isolation window width]"], Is.Null);
        }

        [Test]
        public void HeaderCase_IsNotNormalized()
        {
            // The spec says lowercase and case-sensitive, but the corpus contains mixed-case
            // headers. Rewriting a user's column name is not the reader's call.
            var document = new SdrfDocument(Fixture("PXD059974.sdrf.tsv"));

            Assert.That(document.Header.Any(h => h.Any(char.IsUpper)), Is.True);
            Assert.That(document.Header.Contains("comment[MS max charge]"), Is.True);
            Assert.That(document.Header.Contains("comment[ms max charge]"), Is.False);
        }

        [Test]
        public void Write_IsDeterministic()
        {
            // Same document, written twice, must produce identical bytes. This is what makes the
            // accumulated corpus joinable: a document that serialized differently between runs
            // would show as a spurious change on every re-search.
            var document = new SdrfDocument(Fixture("PXD000070.sdrf.tsv"));

            string first = Path.Combine(TestContext.CurrentContext.TestDirectory, $"det1_{Guid.NewGuid():N}.sdrf.tsv");
            string second = Path.Combine(TestContext.CurrentContext.TestDirectory, $"det2_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                document.WriteResults(first);
                document.WriteResults(second);
                Assert.That(File.ReadAllBytes(second), Is.EqualTo(File.ReadAllBytes(first)));
            }
            finally
            {
                foreach (var p in new[] { first, second })
                    if (File.Exists(p)) File.Delete(p);
            }
        }

        [Test]
        public void InMemoryDocument_Writes()
        {
            var header = new SdrfHeader(new[] { "source name", "assay name", "comment[data file]" });
            var rows = new List<SdrfRow>
            {
                new SdrfRow(header, new[] { "Sample 1", "run 1", "a.raw" }),
                new SdrfRow(header, new[] { "Sample 2", "run 2", "b.raw" })
            };
            var document = new SdrfDocument(header, rows);

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $"mem_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                document.WriteResults(outputPath);

                var text = File.ReadAllText(outputPath);
                Assert.That(text, Is.EqualTo(
                    "source name\tassay name\tcomment[data file]\r\n" +
                    "Sample 1\trun 1\ta.raw\r\n" +
                    "Sample 2\trun 2\tb.raw\r\n"));

                // No byte-order mark -- not one corpus file has one.
                var bytes = File.ReadAllBytes(outputPath);
                Assert.That(bytes.Take(3).ToArray(), Is.Not.EqualTo(new byte[] { 0xEF, 0xBB, 0xBF }));
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }

        [Test]
        public void Write_AppendsExtensionWhenMissing()
        {
            var header = new SdrfHeader(new[] { "source name" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "Sample 1" }) });

            string stem = Path.Combine(TestContext.CurrentContext.TestDirectory, $"stem_{Guid.NewGuid():N}");
            string expected = stem + ".sdrf.tsv";
            try
            {
                document.WriteResults(stem);
                Assert.That(File.Exists(expected), Is.True);
            }
            finally
            {
                if (File.Exists(expected)) File.Delete(expected);
            }
        }

        [Test]
        public void Write_RejectsCellContainingTab()
        {
            // There is no escape mechanism, so a tab inside a cell would silently shift every
            // following column. Fail loudly rather than emit a file that reads back wrong.
            var header = new SdrfHeader(new[] { "source name", "assay name" });
            var document = new SdrfDocument(header,
                new[] { new SdrfRow(header, new[] { "Sample\t1", "run 1" }) });

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $"bad_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                Assert.That(() => document.WriteResults(outputPath), Throws.TypeOf<MzLibException>());
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }

        [Test]
        public void FileReader_ResolvesSdrfByExtension()
        {
            Assert.That(Fixture("PXD000070.sdrf.tsv").ParseFileType(), Is.EqualTo(SupportedFileType.Sdrf));
            Assert.That(SupportedFileType.Sdrf.GetResultFileType(), Is.EqualTo(typeof(SdrfDocument)));

            var loaded = FileReader.ReadFile<SdrfDocument>(Fixture("PXD000070.sdrf.tsv"));
            Assert.That(loaded.Results, Is.Not.Empty);
        }

        /// <summary>
        /// Round-trips EVERY file in the curated corpus. [Explicit] so it never runs in CI: it needs
        /// a local clone of bigbio/sdrf-annotated-datasets (~1,236 files), which is not committed.
        ///
        ///     git clone https://github.com/bigbio/sdrf-annotated-datasets
        ///     set MZLIB_SDRF_CORPUS=&lt;clone&gt;\datasets
        ///     dotnet test --filter "FullyQualifiedName~RoundTrip_EntireCorpus" -- NUnit.Explicit=true
        ///
        /// This is the regression that actually guards the format layer. The three committed
        /// fixtures pin the failure modes already found; this catches the ones that have not been.
        /// </summary>
        [Test]
        [Explicit("Requires a local clone of bigbio/sdrf-annotated-datasets; set MZLIB_SDRF_CORPUS.")]
        public void RoundTrip_EntireCorpus()
        {
            string? corpus = Environment.GetEnvironmentVariable("MZLIB_SDRF_CORPUS");
            if (string.IsNullOrWhiteSpace(corpus) || !Directory.Exists(corpus))
                Assert.Ignore($"MZLIB_SDRF_CORPUS not set or not found: '{corpus}'");

            var files = Directory.GetFiles(corpus, "*.sdrf.tsv", SearchOption.AllDirectories);
            Assert.That(files, Is.Not.Empty, "corpus directory contained no .sdrf.tsv files");

            var failures = new List<string>();
            string scratch = Path.Combine(Path.GetTempPath(), $"sdrf_corpus_{Guid.NewGuid():N}");
            Directory.CreateDirectory(scratch);
            try
            {
                foreach (var file in files)
                {
                    string outputPath = Path.Combine(scratch, "rt.sdrf.tsv");
                    try
                    {
                        new SdrfDocument(file).WriteResults(outputPath);
                        if (!File.ReadAllBytes(outputPath).SequenceEqual(File.ReadAllBytes(file)))
                            failures.Add($"{Path.GetFileName(file)}: bytes differ");
                    }
                    catch (Exception e)
                    {
                        failures.Add($"{Path.GetFileName(file)}: {e.GetType().Name}: {e.Message}");
                    }
                }
            }
            finally
            {
                if (Directory.Exists(scratch)) Directory.Delete(scratch, true);
            }

            TestContext.Progress.WriteLine($"round-tripped {files.Length} files, {failures.Count} failures");
            foreach (var f in failures.Take(25))
                TestContext.Progress.WriteLine("  " + f);

            Assert.That(failures, Is.Empty,
                $"{failures.Count} of {files.Length} corpus files did not round-trip byte-identically");
        }

        [Test]
        public void EmptyFile_ThrowsMzLibException()
        {
            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $"empty_{Guid.NewGuid():N}.sdrf.tsv");
            File.WriteAllText(outputPath, string.Empty);
            try
            {
                var document = new SdrfDocument(outputPath);
                Assert.That(() => document.LoadResults(), Throws.TypeOf<MzLibException>());
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }
    }
}
