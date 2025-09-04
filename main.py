"""CLI entry point for target name normalization."""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path
from typing import List, Sequence

import pandas as pd

from library.io_utils import read_target_names, write_with_new_columns
from library.transforms import NormalizationResult, normalize_target_name


def configure_logging(level: str) -> None:
    """Configure root logger."""
    logging.basicConfig(level=getattr(logging, level.upper(), logging.INFO))


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Normalize target names")
    parser.add_argument("--input", required=True, help="Path to input CSV")
    parser.add_argument("--output", required=True, help="Path to output CSV")
    parser.add_argument("--column", default="target_name", help="Name column")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    parser.add_argument("--delimiter", help="CSV delimiter override")
    parser.add_argument("--encoding", help="File encoding override")
    parser.add_argument(
        "--keep-mutations",
        action="store_true",
        help="Retain mutation-like tokens instead of stripping them",
    )
    parser.add_argument(
        "--no-mutations",
        action="store_true",
        help="Disable mutation detection for faster processing",
    )
    parser.add_argument(
        "--mutation-whitelist",
        help="File with additional mutation tokens to preserve (one per line)",
    )
    parser.add_argument("--taxon", type=int, default=9606, help="NCBI taxon ID")
    parser.add_argument(
        "--json-columns",
        help="Comma-separated columns to serialize as JSON when writing",
    )
    return parser.parse_args()


def normalize_dataframe(
    df: pd.DataFrame,
    column: str,
    *,
    strip_mutations: bool,
    mutation_whitelist: Sequence[str] | None,
    detect_mutations: bool,
    taxon: int,
) -> pd.DataFrame:
    """Apply target name normalization to a DataFrame column.

    This function iterates over the specified column, normalizes each target
    name using :func:`library.transforms.normalize_target_name`, and appends
    the structured results as new columns to the DataFrame.

    Args:
        df: The input DataFrame.
        column: The name of the column containing target names to normalize.
        strip_mutations: If True, mutation-like tokens are removed.
        mutation_whitelist: A sequence of tokens to exempt from mutation stripping.
        detect_mutations: If True, mutation detection is performed.
        taxon: The NCBI taxon ID to associate with the normalization.

    Returns:
        A new DataFrame with added columns containing the normalization results.
    """
    results: List[NormalizationResult] = [
        normalize_target_name(
            str(val),
            strip_mutations=strip_mutations,
            mutation_whitelist=mutation_whitelist,
            detect_mutations=detect_mutations,
            taxon=taxon,
        )
        for val in df[column]
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
    """Entry point for the CLI interface."""
    args = parse_args()
    configure_logging(args.log_level)
    inp = Path(args.input)
    out = Path(args.output)
    df = read_target_names(
        inp, column=args.column, encoding=args.encoding, delimiter=args.delimiter
    )
    extra_whitelist = None
    if args.mutation_whitelist:
        tokens = Path(args.mutation_whitelist).read_text().splitlines()
        extra_whitelist = [t.strip() for t in tokens if t.strip()]
    logging.info("Loaded %d records", len(df))
    start = time.perf_counter()
    df = normalize_dataframe(
        df,
        args.column,
        strip_mutations=not args.keep_mutations,
        mutation_whitelist=extra_whitelist,
        detect_mutations=not args.no_mutations,
        taxon=args.taxon,
    )
    elapsed = time.perf_counter() - start
    logging.info("Normalized %d records in %.2f s", len(df), elapsed)
    logging.info("Sample output: %s", df.head().to_dict(orient="records"))
    json_cols = (
        args.json_columns.split(",")
        if args.json_columns
        else ["hints", "rules_applied"]
    )
    write_with_new_columns(
        df,
        out,
        encoding=args.encoding or "utf-8",
        delimiter=args.delimiter or ",",
        json_columns=json_cols,
    )
    logging.info("Wrote results to %s", out)


if __name__ == "__main__":
    main()
