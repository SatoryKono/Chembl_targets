# Target Name Normalization

Utilities for heavy but non-destructive normalization of protein target names.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Optional quality tools:

```bash
pip install black ruff mypy
```

## Usage

```bash
python main.py --input target_validation_new.csv --output normalized.csv
```

The output CSV includes columns:

- `clean_text` – normalized string without stop words
- `clean_text_alt` – version retaining stop words
- `query_tokens`, `gene_like_candidates`, `hints`, `rules_applied`, `hint_taxon`

Hyphenated tokens (e.g. `beta2-adrenergic`) and letter–digit pairs with spaces
(`h 3`) produce additional variants in `clean_text` such as `beta2adrenergic`
and `h-3` to aid downstream matching.

## Development

Run formatting and tests:

```bash
ruff check .
black --check .
pytest
mypy .
```
