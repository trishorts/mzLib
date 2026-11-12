namespace Readers
{
    /// <summary>A run named in a Spectronaut Run Pivot report's header: its 1-based index and its spectra file name, with extension.</summary>
    public sealed record SpectronautPivotRun(int Index, string FileName)
    {
        /// <summary>The file name without directory or extension.</summary>
        public string Stem => RunFileMatcher.StripPathAndExtension(FileName);
    }

    /// <summary>
    /// The columns of a Spectronaut Run Pivot report. Row-label columns name what each row is (PG.ProteinGroups, PG.Genes,
    /// GO terms...). Every other column is one cell value for one run, headed <c>[&lt;index&gt;] &lt;run file&gt;.&lt;column&gt;</c>.
    /// A schema can export several cell values (MSV000090777's custom pivot has NrOfPrecursorsIdentified,
    /// NrOfPrecursorsUsedForQuantification and Quantity for each run).
    /// </summary>
    public sealed class SpectronautPivotLayout
    {
        internal SpectronautPivotLayout(IReadOnlyList<string> rowLabelColumns, IReadOnlyList<SpectronautPivotRun> runs,
            IReadOnlyList<string> cellColumns, Dictionary<string, int> labelIndex, Dictionary<(int Run, string Column), int> cellIndex)
        {
            RowLabelColumns = rowLabelColumns;
            Runs = runs;
            CellColumns = cellColumns;
            LabelIndex = labelIndex;
            CellIndex = cellIndex;
        }

        /// <summary>The row-label columns, in header order.</summary>
        public IReadOnlyList<string> RowLabelColumns { get; }

        /// <summary>The runs, by index.</summary>
        public IReadOnlyList<SpectronautPivotRun> Runs { get; }

        /// <summary>The cell-value columns (<c>PG.Quantity</c>...), in header order.</summary>
        public IReadOnlyList<string> CellColumns { get; }

        internal Dictionary<string, int> LabelIndex { get; }
        internal Dictionary<(int Run, string Column), int> CellIndex { get; }
    }

    /// <summary>
    /// One row of a Spectronaut Run Pivot (wide) report: its row labels, and one cell per run and cell-value column.
    /// Cells hold numbers, <c>NaN</c> (nothing measured) or <c>Filtered</c> (removed by the q-value filter);
    /// <see cref="State"/> tells these apart.
    /// </summary>
    public sealed class SpectronautPivotRow
    {
        private readonly string[] _fields;

        internal SpectronautPivotRow(SpectronautPivotLayout layout, string[] fields)
        {
            Layout = layout;
            _fields = fields;
        }

        public SpectronautPivotLayout Layout { get; }

        /// <summary>A row label's text, or null when the report has no such row-label column.</summary>
        public string? Label(string column) => Layout.LabelIndex.TryGetValue(column, out int i) ? _fields[i] : null;

        /// <summary>A cell's text, or null when the report has no such run or cell-value column.</summary>
        public string? Text(int runIndex, string column) =>
            Layout.CellIndex.TryGetValue((runIndex, column), out int i) ? _fields[i] : null;

        /// <summary>A cell's number, or null for <c>NaN</c>, <c>Filtered</c>, a blank, or a cell the report does not have.</summary>
        /// <exception cref="FormatException">The cell holds text that is neither a number nor a Spectronaut token.</exception>
        public double? Value(int runIndex, string column) => SpectronautReports.ParseDouble(Text(runIndex, column));

        /// <summary>Whether a cell holds a number, nothing measured, or a filtered value.</summary>
        public SpectronautValueState State(int runIndex, string column) => SpectronautReports.StateOf(Text(runIndex, column));
    }

    /// <summary>One cell of a Run Pivot report in long form: a row, a run and a cell-value column.</summary>
    public readonly record struct SpectronautPivotCell(SpectronautPivotRow Row, SpectronautPivotRun Run, string Column)
    {
        public string? Text => Row.Text(Run.Index, Column);
        public double? Value => Row.Value(Run.Index, Column);
        public SpectronautValueState State => Row.State(Run.Index, Column);
    }
}
