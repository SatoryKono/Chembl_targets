"""CLI entry point for target name normalization."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List

import pandas as pd

from mylib.io_utils import read_target_names, write_with_new_columns
from mylib.transforms import NormalizationResult, normalize_target_name


def configure_logging(level: str) -> None:
    logging.basicConfig(level=getattr(logging, level.upper(), logging.INFO))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Normalize target names")
    parser.add_argument("--input", required=True, help="Path to input CSV")
    parser.add_argument("--output", required=True, help="Path to output CSV")
    parser.add_argument("--column", default="target_name", help="Name column")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    return parser.parse_args()


def normalize_dataframe(df: pd.DataFrame, column: str) -> pd.DataFrame:
    results: List[NormalizationResult] = [
        normalize_target_name(str(val)) for val in df[column]
    ]
    df = df.copy()
    df["clean_text"] = [r.clean_text for r in results]

    df["clean_text_alt"] = [r.clean_text_alt for r in results]
    df["query_tokens"] = ["|".join(r.query_tokens) for r in results]
    df["gene_like_candidates"] = [" ".join(r.gene_like_candidates) for r in results]
    df["hints"] = [r.hints for r in results]
    df["rules_applied"] = [r.rules_applied for r in results]
    df["hint_taxon"] = [r.hint_taxon for r in results]
    return df


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level)
    inp = Path(args.input)
    out = Path(args.output)
    df = read_target_names(inp, column=args.column)
    df = normalize_dataframe(df, args.column)
    write_with_new_columns(df, out)


if __name__ == "__main__":
    main()
