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
(`h 3`) emit both dashed and undashed variants in `clean_text` and
`query_tokens` (`beta2-adrenergic`/`beta2adrenergic`, `h-3`/`h3`) to aid
downstream matching.

Gene-like candidates are inferred via regex rules:

- `histamine h3` → `hrh3`
- `dopamine d2` → `drd2`
- `adrenergic beta1` → `adrb1`
- `p2x3` → `p2rx3`
- `5-ht1a` → `htr1a`
- `gaba a alpha2` → `gabra2`

## Development

Run formatting and tests:

```bash
ruff check .
black --check .
pytest
mypy .
```
