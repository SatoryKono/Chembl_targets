# Target Name Normalization

A Python-based utility for performing heavy but non-destructive normalization of protein target names. This tool is designed to clean, standardize, and enrich target names to aid in downstream matching, database lookups, and data integration.

## Purpose

The primary goal of this repository is to provide a robust and configurable pipeline for normalizing protein target names, which are often inconsistent and messy in biological datasets. The normalization process includes:

-   **Sanitization**: Removing control characters and standardizing whitespace.
-   **Unicode Normalization**: Converting text to a standard Unicode form (NFKC) and handling special characters like Greek letters and superscripts.
-   **Tokenization**: Splitting names into meaningful tokens.
-   **Variant Generation**: Creating alternative forms for hyphenated tokens and adjacent letter-digit pairs to improve matching.
-   **Mutation Detection**: Identifying and optionally stripping mutation-like substrings (e.g., `V600E`).
-   **Gene Candidate Inference**: Generating potential gene symbols from the normalized name using a comprehensive set of regex-based rules.

The output is a structured dataset that preserves the original name while providing cleaned versions, queryable tokens, and rich metadata.

## Installation

1.  **Create and activate a virtual environment:**
    ```bash
    python -m venv .venv
    source .venv/bin/activate
    ```

2.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

3.  **Optional quality tools:**
    For development, you can install the following tools:
    ```bash
    pip install black ruff mypy pandas-stubs
    ```

## Usage

The script is run from the command line, with options to specify input/output files and control the normalization process.

```bash
python main.py --input target_validation_new.csv --output normalized.csv
```

### Command-Line Arguments

-   `--input`: (Required) Path to the input CSV file.
-   `--output`: (Required) Path to the output CSV file.
-   `--column`: The name of the column containing the target names (default: `target_name`).
-   `--delimiter`: Overrides the auto-detected CSV delimiter.
-   `--encoding`: Overrides the auto-detected file encoding.
-   `--keep-mutations`: If set, mutation-like tokens are retained in the output.
-   `--no-mutations`: If set, mutation detection is disabled for faster processing.
-   `--mutation-whitelist`: Path to a text file with additional mutation tokens to preserve (one per line).
-   `--taxon`: The NCBI taxon ID to associate with the results (default: `9606` for *Homo sapiens*).
-   `--json-columns`: A comma-separated list of columns to serialize as JSON in the output (default: `hints,rules_applied`).

## Normalization Pipeline

The normalization process follows these steps:

1.  **Sanitization**: Removes control characters and normalizes all whitespace to single spaces.
2.  **Unicode Normalization**: Applies NFKC normalization, converts to lowercase, and standardizes special characters (e.g., Greek letters, smart quotes).
3.  **Mutation Detection**: Identifies mutation-like substrings based on a set of regex patterns and, if available, the `hgvs` library.
4.  **Roman Numeral Replacement**: Converts standalone Roman numerals (II-XX) to digits.
5.  **Parenthetical Extraction**: Extracts content from brackets `()` `[]` `{}`. Short, index-like tokens are kept, while others are moved to the `hints` field.
6.  **Tokenization**: The cleaned string is tokenized based on spaces and punctuation.
7.  **Variant Generation**:
    -   **Hyphenated Tokens**: For tokens like `beta-2`, `beta-2` and `beta2` are generated.
    -   **Letter-Digit Pairs**: For tokens like `h 3`, `h3` and `h-3` are generated.
8.  **Stop Word Removal**: Common, non-specific words (e.g., "protein", "receptor") are removed to create `clean_text`. The version with stop words is kept as `clean_text_alt`.
9.  **Gene Candidate Inference**: A final set of rules is applied to the normalized tokens to generate potential gene symbols.

## Output Columns

The output CSV file contains the original data along with the following new columns:

-   `clean_text`: The primary normalized name, generated after removing stop words. Multiple variants are joined by a pipe (`|`).
-   `clean_text_alt`: An alternative normalized name that retains stop words.
-   `query_tokens`: A pipe-separated list of all generated tokens, suitable for database queries.
-   `gene_like_candidates`: A space-separated list of inferred gene symbols.
-   `hint_taxon`: The NCBI taxon ID provided during execution.
-   `hints`: A JSON object containing metadata about the normalization process:
    -   `parenthetical`: A list of strings extracted from brackets.
    -   `dropped`: A list of stop words that were removed.
    -   `mutations`: A list of detected mutation substrings.
    -   `mutations_only`: A boolean flag that is `true` if the name consisted only of mutations.
-   `rules_applied`: A JSON list of the regex rule patterns that were successfully applied.

## Gene-like Candidate Rules

The script uses an extensive set of rules to infer gene symbols. Here are some examples:

| Input Text Pattern                  | Inferred Candidate(s)                                           |
| ----------------------------------- | --------------------------------------------------------------- |
| `histamine h3`                      | `hrh3`                                                          |
| `dopamine d2`                       | `drd2`                                                          |
| `adrenergic beta1`                  | `adrb1`                                                         |
| `p2x3`                              | `p2rx3`                                                         |
| `5-ht1a`                            | `htr1a`                                                         |
| `gaba a alpha2`                     | `gabra2`                                                        |
| `ampa glua2`                        | `gria2` (or `gria1-4` if subtype is absent)                     |
| `nmda nr2b`                         | `grin2b` (or `grin` family if subtype is absent)                |
| `kainate gluk3`                     | `grik3` (or `grik` family if subtype is absent)                 |
| `mglur5`                            | `grm5` (or `grm` family if subtype is absent)                   |
| `adenosine a2a receptor`            | `a2a`, `adora2a`                                                |
| `nociceptin receptor`, `orl1`, `nop` | `nop`, `orl1`, `oprl1`                                          |
| `melanocortin 4 receptor`           | `mc4r`                                                          |
| `prostaglandin ep3 receptor`        | `ep3`, `ptger3`                                                 |

The rules cover a wide range of protein families, including GPCRs, ion channels, and nuclear receptors.

## Development

To maintain code quality, run the following checks:

```bash
# Check for linting issues
ruff check .

# Check for formatting issues
black --check .

# Run unit tests
pytest

# Run static type checking
mypy .
```
