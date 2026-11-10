using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.ExternalFileReading
{
    /// <summary>
    /// Spectronaut19_ConditionSetup.tsv is the whole condition setup of MSV000096226: 10 runs, conditions Old (4) and RF (6),
    /// replicates 1-based, Fraction NA, Reference True for the Old condition.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestSpectronautConditionSetupFile
    {
        private static string FixturePath =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes", "Spectronaut19_ConditionSetup.tsv");

        [Test]
        public void ReadsEveryRun()
        {
            var file = new SpectronautConditionSetupFile(FixturePath);
            Assert.That(file.Count(), Is.EqualTo(10));
            SpectronautConditionSetupRow first = file.First();
            Assert.That(first.Number, Is.EqualTo(1));
            Assert.That(first.IsReferenceCondition, Is.True);
            Assert.That(first.RunLabel, Is.EqualTo("240418_AM_O1.raw"));
            Assert.That(first.RunStem, Is.EqualTo("240418_AM_O1"));
            Assert.That(first.Condition, Is.EqualTo("Old"));
            Assert.That(first.Fraction, Is.Null);
            Assert.That(first.Replicate, Is.EqualTo(1));
            Assert.That(first.QuantityCorrectionFactor, Is.EqualTo(1));
            Assert.That(first.FileName, Is.EqualTo("240418_AM_O1"));
            Assert.That(file.Count(r => r.Condition == "RF"), Is.EqualTo(6));
        }

        [Test]
        public void TheDesignIsZeroBasedAndNeverMarksAPooledReference()
        {
            SampleExperimentalDesign design = new SpectronautConditionSetupFile(FixturePath).ToExperimentalDesign();
            var files = design.FileNameSampleInfoDictionary.Values.SelectMany(s => s).Cast<SpectraFileInfo>().ToList();
            Assert.That(files, Has.Count.EqualTo(10));

            var o1 = files.Single(f => f.FilenameWithoutExtension == "240418_AM_O1");
            Assert.That(o1.Condition, Is.EqualTo("Old"));
            Assert.That(o1.BiologicalReplicate, Is.EqualTo(0));
            Assert.That(o1.Fraction, Is.EqualTo(0));
            Assert.That(o1.TechnicalReplicate, Is.EqualTo(0));

            var rf4 = files.Single(f => f.FilenameWithoutExtension == "240418_AM_ODR_RF_4");
            Assert.That(rf4.Condition, Is.EqualTo("RF"));
            Assert.That(rf4.BiologicalReplicate, Is.EqualTo(5));

            // Spectronaut's Reference marks the ratio denominator condition, not a pooled reference sample
            Assert.That(files.All(f => !f.IsReference));
        }

        [Test]
        public void RunsAreKeyedOnTheGivenSpectraFiles()
        {
            var paths = new SpectronautConditionSetupFile(FixturePath).Select(r => $@"D:\Data\{r.RunStem}.raw").ToList();
            var design = new SpectronautConditionSetupFile(FixturePath).ToExperimentalDesign(paths);
            var files = design.FileNameSampleInfoDictionary.Values.SelectMany(s => s).Cast<SpectraFileInfo>().ToList();
            Assert.That(files.Select(f => f.FullFilePathWithExtension), Is.EquivalentTo(paths));
        }

        [Test]
        public void ARunWithNoFileIsRefusedByName()
        {
            var e = Assert.Throws<MzLibException>(() =>
                new SpectronautConditionSetupFile(FixturePath).ToExperimentalDesign(new List<string> { @"D:\Data\240418_AM_O1.raw" }));
            Assert.That(e!.Message, Does.Contain("240418_AM_O2.raw"));
        }

        [Test]
        public void ARunWithNoConditionIsRefused()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"nocondition_{Guid.NewGuid():N}_ConditionSetup.tsv");
            File.WriteAllText(path, "Run Label\tCondition\tReplicate\r\nrun1.raw\tNot Defined\t1\r\n");
            try
            {
                var e = Assert.Throws<MzLibException>(() => new SpectronautConditionSetupFile(path).ToExperimentalDesign());
                Assert.That(e!.Message, Does.Contain("run1.raw"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void FactoryBuiltFileReportsSpectronautAndWritingIsNotSupported()
        {
            Assert.That(new SpectronautConditionSetupFile().Software, Is.EqualTo(Software.Spectronaut));
            Assert.Throws<NotSupportedException>(() => new SpectronautConditionSetupFile(FixturePath).WriteResults("x_ConditionSetup.tsv"));
            Assert.Throws<FileNotFoundException>(() => new SpectronautConditionSetupFile("missing_ConditionSetup.tsv").LoadResults());
        }
    }
}
