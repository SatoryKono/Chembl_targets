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

## Development

Run formatting and tests:

```bash
ruff check .
black --check .
pytest
mypy .
```
