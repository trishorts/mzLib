using System;
using System.IO;
using System.Threading.Tasks;
using NUnit.Framework;
using PredictionClients.AI;  // Changed from Proteomics.AI

namespace Test
{
    [TestFixture]
    [Category("AI")]
    [Category("Koina")]
    public class TestPeptideOptimization
    {
        private string OutputDir => Path.Combine(@"C:\Users\trish\Downloads", "PeptideOptimization");

        [OneTimeSetUp]
        public void Setup()
        {
            if (!Directory.Exists(OutputDir))
                Directory.CreateDirectory(OutputDir);
        }

        [Test]
        [Explicit("Requires Koina API access")]
        public async Task OptimizePeptideFragmentation_Basic()
        {
            var optimizer = new PeptideFragmentationOptimizer(
                peptideLength: 12,
                populationSize: 30,
                collisionEnergy: 25,
                precursorCharge: 2,
                randomSeed: 42);

            var best = await optimizer.OptimizeAsync(
                generations: 10,
                onGenerationComplete: (gen, bestFitness) =>
                {
                    Console.WriteLine($"Gen {gen}: {bestFitness}");
                });

            // Save results
            var report = optimizer.GetProgressReport();
            var patterns = optimizer.GetLearnedPatternsSummary();

            File.WriteAllText(Path.Combine(OutputDir, "optimization_report.txt"), report);
            File.WriteAllText(Path.Combine(OutputDir, "learned_patterns.txt"), patterns);

            Console.WriteLine(report);
            Console.WriteLine(patterns);
            Console.WriteLine($"\nBEST PEPTIDE: {best.Sequence}");
            Console.WriteLine($"Fitness: {best.OverallFitness:F4}");

            Assert.That(best.Sequence.Length, Is.EqualTo(12));
            Assert.That(best.Sequence.Distinct().Count(), Is.EqualTo(12));
        }

        [Test]
        [Explicit("Requires Anthropic API key and Koina access")]
        public async Task OptimizePeptideFragmentation_WithClaude()
        {
            string apiKey = Environment.GetEnvironmentVariable("ANTHROPIC_API_KEY");
            if (string.IsNullOrEmpty(apiKey))
            {
                Assert.Ignore("ANTHROPIC_API_KEY not set");
            }

            var optimizer = new ClaudeGuidedPeptideOptimizer(
                anthropicApiKey: apiKey,
                peptideLength: 12,
                populationSize: 30);

            var best = await optimizer.OptimizeWithClaudeGuidanceAsync(
                generations: 15,
                claudeAnalysisInterval: 5,
                onProgress: (gen, fitness, insight) =>
                {
                    Console.WriteLine($"Gen {gen}: {fitness}");
                    if (insight != null)
                    {
                        Console.WriteLine("\n=== CLAUDE ANALYSIS ===");
                        Console.WriteLine(insight);
                        Console.WriteLine("=======================\n");
                    }
                });

            // Save everything
            File.WriteAllText(
                Path.Combine(OutputDir, "claude_guided_report.txt"),
                optimizer.GetProgressReport());
            
            File.WriteAllText(
                Path.Combine(OutputDir, "claude_insights.txt"),
                optimizer.GetAllClaudeInsights());
            
            File.WriteAllText(
                Path.Combine(OutputDir, "final_patterns.txt"),
                optimizer.GetLearnedPatternsSummary());

            Console.WriteLine($"\nBEST PEPTIDE: {best.Sequence}");
            Console.WriteLine($"Fitness: {best.OverallFitness:F4}");
        }
    }
}