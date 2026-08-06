using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for <see cref="SdrfValidator"/>.
    ///
    /// Two halves. The unit tests below build a minimal synthetic document and check that each rule
    /// fires exactly when it should. <see cref="CorpusCalibration"/> then runs the validator over
    /// the whole curated corpus and pins the observed rates, because a rule that condemns most of
    /// the community reference set is a wrong rule, not 1,236 wrong files.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfValidator
    {
        private static string FixtureDir => Path.Combine(
            TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes");

        /// <summary>A minimal document that should validate with no errors.</summary>
        private static SdrfDocument MinimalValid()
        {
            var header = new SdrfHeader(new[]
            {
                "source name", "characteristics[organism]", "characteristics[organism part]",
                "characteristics[disease]", "characteristics[biological replicate]",
                "assay name", "technology type",
                "comment[proteomics data acquisition method]", "comment[label]", "comment[instrument]",
                "comment[cleavage agent details]", "comment[fraction identifier]",
                "comment[technical replicate]", "comment[data file]",
                "factor value[disease]"
            });

            SdrfRow Row(string sample, string file) => new(header, new[]
            {
                sample, "homo sapiens", "liver", "normal", "1",
                "run " + sample, "proteomic profiling by mass spectrometry",
                "NT=Data-dependent acquisition;AC=PRIDE:0000627", "label free sample",
                "NT=Q Exactive;AC=MS:1001911", "NT=Trypsin;AC=MS:1001251", "1", "1", file,
                "normal"
            });

            return new SdrfDocument(header, new[] { Row("Sample 1", "a.raw"), Row("Sample 2", "b.raw") });
        }

        [Test]
        public void MinimalDocument_IsValidWithNoWarnings()
        {
            var result = SdrfValidator.Validate(MinimalValid());
            Assert.That(result.IsValid, Is.True, result.ToString());
            Assert.That(result.Messages, Is.Empty,
                "unexpected: " + string.Join(" | ", result.Messages.Select(m => m.ToString())));
        }

        [Test]
        public void MissingRequiredColumn_IsError()
        {
            var header = new SdrfHeader(new[] { "assay name", "technology type" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "run 1", "x" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.IsValid, Is.False);
            Assert.That(result.Errors.Any(e => e.Rule == "RequiredColumn" && e.ColumnName == "source name"));
        }

        [Test]
        public void MissingRecommendedColumn_IsWarningNotError()
        {
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "s", "a", "t" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.IsValid, Is.True, "a missing recommended column must not invalidate");
            Assert.That(result.Warnings.Any(w => w.Rule == "RecommendedColumn"));
        }

        [Test]
        public void RaggedRow_IsError()
        {
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "s", "a" }) });

            var result = SdrfValidator.Validate(document);
            var finding = result.Errors.FirstOrDefault(e => e.Rule == "RowWidth");
            Assert.That(finding, Is.Not.Null);
            Assert.That(finding!.RowIndex, Is.EqualTo(0));
            Assert.That(finding.LineNumber, Is.EqualTo(2), "line number is row index + 2 (header is line 1)");
        }

        [Test]
        public void DuplicateRowKey_IsError()
        {
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type", "comment[label]" });
            var document = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "Sample 1", "run 1", "t", "label free sample" }),
                new SdrfRow(header, new[] { "Sample 1", "run 1", "t", "label free sample" })
            });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Errors.Any(e => e.Rule == "RowKeyUniqueness"));
        }

        [Test]
        public void SameSourceDifferentLabel_IsNotDuplicate()
        {
            // The whole point of including comment[label] in the key: one sample, two TMT channels.
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type", "comment[label]" });
            var document = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "Sample 1", "run 1", "t", "TMT126" }),
                new SdrfRow(header, new[] { "Sample 1", "run 1", "t", "TMT127N" })
            });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Errors.Any(e => e.Rule == "RowKeyUniqueness"), Is.False);
        }

        [Test]
        public void UppercaseColumnName_IsWarning()
        {
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type", "comment[MS min charge]" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "s", "a", "t", "2" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.IsValid, Is.True);
            Assert.That(result.Warnings.Any(w => w.Rule == "ColumnNameCase"));
        }

        [Test]
        public void SpaceBeforeBracket_IsWarning()
        {
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type", "characteristics [organism]" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "s", "a", "t", "homo sapiens" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "MalformedColumnName"));
        }

        [Test]
        public void CasedReservedWord_IsWarning()
        {
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "s", "a", "Not Available" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "ReservedWordCase"));
        }

        [Test]
        public void NonIntegerReplicate_IsWarning()
        {
            var header = new SdrfHeader(new[]
                { "source name", "assay name", "technology type", "comment[technical replicate]" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "s", "a", "t", "one" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "NonIntegerValue"));
        }

        [Test]
        public void ReservedWordInIntegerColumn_IsAccepted()
        {
            var header = new SdrfHeader(new[]
                { "source name", "assay name", "technology type", "characteristics[biological replicate]" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "s", "a", "t", "pooled" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "NonIntegerValue"), Is.False,
                "'pooled' is the spec's own way to say the replicate is not a number");
        }

        [Test]
        public void FactorValueBeforeComment_IsOrderingWarning()
        {
            var header = new SdrfHeader(new[]
                { "source name", "factor value[disease]", "assay name", "technology type", "comment[data file]" });
            var document = new SdrfDocument(header,
                new[] { new SdrfRow(header, new[] { "s", "normal", "a", "t", "x.raw" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "ColumnOrdering"));
        }

        [Test]
        public void RaggedCorpusFile_IsReportedAsError()
        {
            // PXD059974 is the known-malformed corpus file. The reader loads it happily; the
            // validator is what says so.
            var document = new SdrfDocument(Path.Combine(FixtureDir, "PXD059974.sdrf.tsv"));
            var result = SdrfValidator.Validate(document);

            Assert.That(result.IsValid, Is.False);
            Assert.That(result.Errors.Count(e => e.Rule == "RowWidth"), Is.EqualTo(17),
                "17 of its 23 rows are short");
        }

        [Test]
        public void TypicalCorpusFile_HasNoErrors()
        {
            var document = new SdrfDocument(Path.Combine(FixtureDir, "PXD000070.sdrf.tsv"));
            var result = SdrfValidator.Validate(document);

            Assert.That(result.IsValid, Is.True,
                "PXD000070 is a well-formed curated file: " +
                string.Join(" | ", result.Errors.Select(e => e.ToString())));
        }

        /// <summary>
        /// Runs the validator over the entire curated corpus and reports how often each rule fires.
        /// [Explicit]; needs MZLIB_SDRF_CORPUS. See TestSdrf.RoundTrip_EntireCorpus.
        ///
        /// This is rule calibration, not a pass/fail check on the corpus. The assertion is that the
        /// great majority of curated files carry no ERROR -- if that stops being true, a rule has
        /// been miscategorised and the validator has become useless against real data.
        /// </summary>
        [Test]
        [Explicit("Requires a local clone of bigbio/sdrf-annotated-datasets; set MZLIB_SDRF_CORPUS.")]
        public void CorpusCalibration()
        {
            string? corpus = Environment.GetEnvironmentVariable("MZLIB_SDRF_CORPUS");
            if (string.IsNullOrWhiteSpace(corpus) || !Directory.Exists(corpus))
                Assert.Ignore($"MZLIB_SDRF_CORPUS not set or not found: '{corpus}'");

            var files = Directory.GetFiles(corpus, "*.sdrf.tsv", SearchOption.AllDirectories);

            var filesByRule = new Dictionary<string, int>(StringComparer.Ordinal);
            var severityByRule = new Dictionary<string, SdrfValidationSeverity>(StringComparer.Ordinal);
            int filesWithErrors = 0, filesClean = 0;

            foreach (var file in files)
            {
                var result = SdrfValidator.Validate(new SdrfDocument(file));
                if (!result.IsValid) filesWithErrors++;
                if (result.Messages.Count == 0) filesClean++;

                foreach (var rule in result.Messages.Select(m => m.Rule).Distinct(StringComparer.Ordinal))
                {
                    filesByRule.TryGetValue(rule, out int n);
                    filesByRule[rule] = n + 1;
                    severityByRule[rule] = result.Messages.First(m => m.Rule == rule).Severity;
                }
            }

            TestContext.Progress.WriteLine($"corpus files          : {files.Length}");
            TestContext.Progress.WriteLine($"entirely clean        : {filesClean}");
            TestContext.Progress.WriteLine($"with >=1 ERROR        : {filesWithErrors} " +
                                           $"({100.0 * filesWithErrors / files.Length:F1}%)");
            TestContext.Progress.WriteLine("rule                     severity  files   %");
            foreach (var kv in filesByRule.OrderByDescending(k => k.Value))
            {
                TestContext.Progress.WriteLine(
                    $"  {kv.Key,-22} {severityByRule[kv.Key],-9} {kv.Value,5}  {100.0 * kv.Value / files.Length,5:F1}");
            }

            // Observed 2026-08-06 against corpus @ 4f823dcd: 1,236 files, 1,143 entirely clean
            // (92.5%), exactly 1 with an error (PXD059974, genuinely ragged). Highest-firing
            // warning was ColumnNameCase at 2.8%.
            //
            // The bounds are loose enough to survive corpus growth but tight enough that promoting
            // any warning to an error, or adding an over-eager rule, fails here rather than quietly
            // condemning the reference set.
            double errorRate = 100.0 * filesWithErrors / files.Length;
            double cleanRate = 100.0 * filesClean / files.Length;

            Assert.That(errorRate, Is.LessThan(1.0),
                $"{errorRate:F1}% of the curated corpus is reported invalid; a rule is miscategorised " +
                "as Error. Only genuinely unreadable files should error.");
            Assert.That(cleanRate, Is.GreaterThan(85.0),
                $"only {cleanRate:F1}% of the curated corpus is warning-free; a rule is over-eager.");
        }
    }
}
