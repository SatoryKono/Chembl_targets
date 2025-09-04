"""I/O utilities for target normalization.

This module contains functions for reading the source CSV with
automatic encoding and delimiter detection, as well as helpers for
loading external mapping dictionaries.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Dict, Sequence, Tuple

import pandas as pd  # type: ignore[import-untyped]

# Default encodings and delimiters to try when reading CSV files.
ENCODINGS: Tuple[str, ...] = ("utf-8-sig", "cp1251", "latin1")
DELIMITERS: Tuple[str, ...] = (",", ";", "\t", "|")


def detect_csv_format(path: Path) -> Tuple[str, str]:
    """Detect encoding and delimiter of a CSV file.

    Args:
        path: Path to the CSV file.

    Returns:
        A tuple containing the detected encoding and delimiter.

    Raises:
        ValueError: If the file format could not be detected. The message
            includes a snippet of the file to aid debugging.
    """
    snippet = ""
    try:
        snippet = (
            path.read_bytes()[:200]
            .decode("utf-8", errors="replace")
            .replace("\n", "\\n")
        )
    except Exception:
        pass

    for enc in ENCODINGS:
        try:
            with path.open("r", encoding=enc) as fh:
                sample = fh.read(4096)
                sniffer = csv.Sniffer()
                try:
                    dialect = sniffer.sniff(sample, delimiters="".join(DELIMITERS))
                    return enc, dialect.delimiter
                except csv.Error:
                    first_line = sample.splitlines()[0] if sample else ""
                    for delim in DELIMITERS:
                        if delim in first_line:
                            return enc, delim
                    # fallback to comma when single column
                    return enc, ","
        except Exception:
            continue
    raise ValueError(
        "Could not detect CSV encoding or delimiter. "
        f"Snippet: {snippet!r}. "
        "Consider verifying the file or specifying the column with --column."
    )


def read_target_names(
    path: Path,
    column: str = "target_name",
    *,
    encoding: str | None = None,
    delimiter: str | None = None,
) -> pd.DataFrame:
    """Read the CSV and return DataFrame with target names.

    Args:
        path: Path to the input CSV file.
        column: Name of the column containing target names.
        encoding: Optional file encoding. If None, it's auto-detected.
        delimiter: Optional delimiter. If None, it's auto-detected.

    Returns:
        A DataFrame containing at least the raw target name column.
    """
    if encoding is None or delimiter is None:
        det_enc, det_delim = detect_csv_format(path)
        encoding = encoding or det_enc
        delimiter = delimiter or det_delim
    df = pd.read_csv(path, encoding=encoding, sep=delimiter)
    if column not in df.columns:
        sample = df.head().to_string(index=False)
        cols = ", ".join(df.columns)
        raise KeyError(
            f"Column '{column}' not found in {path}. "
            f"Available columns: {cols}. Sample:\n{sample}"
        )
    return df


def write_with_new_columns(
    df: pd.DataFrame,
    path: Path,
    *,
    encoding: str = "utf-8",
    delimiter: str = ",",
    json_columns: Sequence[str] | None = ("hints", "rules_applied"),
) -> None:
    """Write DataFrame to CSV with given encoding.

    Args:
        df: DataFrame to save.
        path: Output path for the CSV file.
        encoding: Target encoding (default is UTF-8).
        delimiter: Field delimiter for the output CSV.
        json_columns: A sequence of column names to serialize as JSON.
    """
    df = df.copy()
    if json_columns:
        for col in json_columns:
            if col in df.columns:
                df[col] = df[col].apply(lambda x: json.dumps(x, ensure_ascii=False))
    df.to_csv(path, index=False, encoding=encoding, sep=delimiter)


def load_mapping(path: Path) -> Dict[str, str]:
    """Load a mapping dictionary from a JSON or TSV file.

    Args:
        path: Path to the JSON or TSV file.

    Returns:
        A dictionary containing the loaded mapping.

    Raises:
        ValueError: If the file format is not supported.
    """
    if path.suffix.lower() == ".json":
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    if path.suffix.lower() in {".tsv", ".txt"}:
        mapping: Dict[str, str] = {}
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                key, val = line.rstrip().split("\t", 1)
                mapping[key] = val
        return mapping
    raise ValueError(f"Unsupported mapping format: {path.suffix}")
