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

- `clean_text` – normalized name without stop words
- `clean_text_alt` – version retaining stop words
- `query_tokens`, `gene_like_candidates`, `hints`, `rules_applied`, `hint_taxon`

`clean_text` and `clean_text_alt` collect all distinct textual variants
(hyphenated/concatenated forms, retained bracket indices, etc.) joined with a
pipe (`|`) only when more than one unique form exists. `query_tokens` uses the
same separator for multiple tokens. `hints` and `rules_applied` are stored as
JSON structures in the output CSV, preserving the original dictionaries and
lists used during processing.


Hyphenated tokens (e.g. `beta2-adrenergic`) and letter–digit pairs with spaces
(`h 3`) emit both dashed and undashed variants in `clean_text` and
`query_tokens` (`beta2-adrenergic`/`beta2adrenergic`, `h-3`/`h3`) to aid
downstream matching. Original split tokens (e.g. `h` and `3`) are also
preserved in `query_tokens`.

Short indices enclosed in brackets such as `h3`, `p2x7`, or `5-ht1a` are
kept; the bracketed content is stored under `hints.parenthetical`.

Mutation-like substrings (`V600E`, `p.Gly12Asp`, `ΔF508`, `F508del`, etc.) are
detected via regex rules, with optional assistance from the `hgvs` parser when
installed. By default these tokens are removed from the normalized output and
recorded under `hints.mutations`. Use `--keep-mutations` to retain them. Tokens
that resemble mutations but are valid receptor names are preserved through
`MUTATION_WHITELIST`; additional tokens can be supplied with
`--mutation-whitelist` (a text file with one token per line). If removing
mutations leaves no core tokens, the original tokens are restored and
`hints.mutations_only` is set to `true`.

Gene-like candidates are inferred via regex rules:

- `histamine h3` → `hrh3`
- `dopamine d2` → `drd2`
- `adrenergic beta1` → `adrb1`
- `p2x3` → `p2rx3`
- `5-ht1a` → `htr1a`
- `gaba a alpha2` → `gabra2`
- `trp v 1` → `trpv1`
- `ampa glua2` → `gria2` (or `gria1`–`gria4` if subtype absent)
- `nmda nr2b` → `grin2b` (family fallback to `grin1`, `grin2a`–`grin2d`, `grin3a`, `grin3b`)
- `kainate gluk3` → `grik3` (family fallback to `grik1`–`grik5`)
- `mglur5` / `metabotropic glutamate receptor` → `grm1`–`grm8`
- `chemokine cc receptor 5` / `cxcr4` → corresponding `ccr`/`cxcr` genes
- ligand aliases such as `sdf-1` → `cxcr4`, `il-8` → `cxcr1|cxcr2`,
  `rantes` → `ccr1|ccr3|ccr5`, `fractalkine` → `cx3cr1`
- `adenosine a2a receptor` → `a2a|adora2a`; `adenosine receptor` → `adora1`–`adora3`
- `nociceptin receptor` / `orl1` / `nop` → `nop|orl1|oprl1`
- `neuropeptide y1 receptor` → `y1|npy1r` (family `neuropeptide y receptor` → `npy1r`–`npy5r`)
- `melanocortin 4 receptor` → `mc4r` (family `melanocortin receptor` → `mc1r`–`mc5r`)
- `prostaglandin ep3 receptor` → `ep3|ptger3` (family `prostaglandin receptor` → `ptger1`–`ptger4`, `ptgdr`, `ptgdr2`, `ptgfr`, `ptgir`, `tbxa2r`)
- Additional GPCR families such as calcitonin/CGRP/amylin (CALCR/CALCRL + RAMP1–3), parathyroid hormone (PTH1R/2R), neuropeptide S/FF/B/W, neuromedin U, kisspeptin (KISS1R/GPR54), ghrelin (GHSR), motilin (MLNR/GPR38), prolactin-releasing peptide (PRLHR/GPR10), melanin-concentrating hormone (MCHR1/2), fractalkine (CX3CR1) and XCR1, platelet-activating factor (PTAFR), formyl peptide receptors (FPR1–3/ALX), free fatty acid and hydroxycarboxylic acid receptors (FFAR1–4/GPR40/41/43/120/84 and HCAR1–3/GPR81/109A/B), trace amine-associated receptors (TAAR1–9), bile acid (GPBAR1/TGR5), urotensin II (UTS2R) and apelin (APLNR) receptors.

## Development

Run formatting and tests:

```bash
ruff check .
black --check .
pytest
mypy .
```
