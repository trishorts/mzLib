using System.Diagnostics.CodeAnalysis;
using System.Text;
using MzLibUtil;

namespace Readers.ExternalResults.ResultFiles;

/// <summary>
/// Reads an mzTab parameter, "[cvLabel, accession, name, value]", into a <see cref="CvParam"/>.
/// </summary>
/// <remarks>
/// Producers bend the grammar, and every bend below comes from a real file:
/// <list type="bullet">
/// <item>a user parameter leaves the label and accession empty: MaxQuant's "[, , CHEMMOD:57.0214637236,]" puts the
/// modification in the name;</item>
/// <item>the value may be missing, giving three fields (Proteome Discoverer);</item>
/// <item>a field may be wrapped in double quotes, and a comma inside the quotes is part of it (the spec's escape);</item>
/// <item>an unquoted value that itself holds commas gives more than four fields; the extra fields are kept as part
/// of the value rather than refusing the parameter.</item>
/// </list>
/// </remarks>
internal static class MzTabCvParam
{
    /// <summary>
    /// False when <paramref name="text"/> is not a bracketed parameter with at least three fields.
    /// </summary>
    internal static bool TryParse(string? text, [MaybeNullWhen(false)] out CvParam param)
    {
        param = null;
        if (text == null)
        {
            return false;
        }

        string trimmed = text.Trim();
        if (trimmed.Length < 2 || trimmed[0] != '[' || trimmed[^1] != ']')
        {
            return false;
        }

        var fields = SplitOutsideQuotes(trimmed[1..^1]);
        if (fields.Count < 3)
        {
            return false;
        }

        string value = fields.Count == 3 ? "" : Unquote(string.Join(",", fields.Skip(3)));
        param = new CvParam(Unquote(fields[0]), Unquote(fields[1]), Unquote(fields[2]), value);
        return true;
    }

    /// <summary>
    /// Every top-level parameter in a '|'-separated cell, such as search_engine's
    /// "[MS, MS:1001207, Mascot, ]|[MS, MS:1001583, MaxQuant, ]". Parts that are not parameters are left out.
    /// </summary>
    internal static List<CvParam> ParseList(string? cell)
    {
        var parameters = new List<CvParam>();
        if (string.IsNullOrWhiteSpace(cell))
        {
            return parameters;
        }

        foreach (var part in SplitTopLevel(cell, '|'))
        {
            if (TryParse(part, out var parameter))
            {
                parameters.Add(parameter);
            }
        }

        return parameters;
    }

    /// <summary>
    /// Splits <paramref name="text"/> on <paramref name="separator"/> where it is outside square brackets and double
    /// quotes. Used for the modifications cell and for '|'-separated parameter lists.
    /// </summary>
    internal static List<string> SplitTopLevel(string text, char separator)
    {
        var parts = new List<string>();
        var current = new StringBuilder();
        int depth = 0;
        bool quoted = false;
        foreach (char c in text)
        {
            if (c == '"')
            {
                quoted = !quoted;
            }
            else if (!quoted && c == '[')
            {
                depth++;
            }
            else if (!quoted && c == ']' && depth > 0)
            {
                depth--;
            }
            else if (!quoted && depth == 0 && c == separator)
            {
                parts.Add(current.ToString());
                current.Clear();
                continue;
            }

            current.Append(c);
        }

        parts.Add(current.ToString());
        return parts;
    }

    private static List<string> SplitOutsideQuotes(string inner)
    {
        var fields = new List<string>();
        var current = new StringBuilder();
        bool quoted = false;
        foreach (char c in inner)
        {
            if (c == '"')
            {
                quoted = !quoted;
            }
            else if (c == ',' && !quoted)
            {
                fields.Add(current.ToString());
                current.Clear();
                continue;
            }

            current.Append(c);
        }

        fields.Add(current.ToString());
        return fields;
    }

    private static string Unquote(string field)
    {
        string trimmed = field.Trim();
        return trimmed.Length >= 2 && trimmed[0] == '"' && trimmed[^1] == '"' ? trimmed[1..^1] : trimmed;
    }
}
