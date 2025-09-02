"""I/O utilities for target normalization.

This module contains functions for reading the source CSV with
automatic encoding and delimiter detection, as well as helpers for
loading external mapping dictionaries.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple
import csv
import json

import pandas as pd

# Default encodings and delimiters to try when reading CSV files.
ENCODINGS: Tuple[str, ...] = ("utf-8-sig", "cp1251", "latin1")
DELIMITERS: Tuple[str, ...] = (",", ";", "\t", "|")


def detect_csv_format(path: Path) -> Tuple[str, str]:
    """Detect encoding and delimiter of a CSV file.

    Parameters
    ----------
    path:
        Path to the CSV file.

    Returns
    -------
    Tuple[str, str]
        Detected encoding and delimiter.

    Raises
    ------
    ValueError
        If the file format could not be detected.
    """
    for enc in ENCODINGS:
        try:
            with path.open("r", encoding=enc) as fh:
                sample = fh.read(4096)
                sniffer = csv.Sniffer()
                try:
                    dialect = sniffer.sniff(sample, delimiters="".join(DELIMITERS))
                    return enc, dialect.delimiter
                except csv.Error:
                    first_line = sample.splitlines()[0]
                    for delim in DELIMITERS:
                        if delim in first_line:
                            return enc, delim
                    # fallback to comma when single column
                    return enc, ","
        except Exception:
            continue
    raise ValueError("Could not detect CSV encoding or delimiter")


def read_target_names(path: Path, column: str = "target_name") -> pd.DataFrame:
    """Read the CSV and return DataFrame with target names.

    Parameters
    ----------
    path:
        Path to the input CSV file.
    column:
        Name of the column containing target names.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing at least the raw column.
    """
    enc, delim = detect_csv_format(path)
    df = pd.read_csv(path, encoding=enc, sep=delim)
    if column not in df.columns:
        raise KeyError(f"Column '{column}' not found in {path}")
    return df


def write_with_new_columns(
    df: pd.DataFrame, path: Path, encoding: str = "utf-8"
) -> None:
    """Write DataFrame to CSV with given encoding.

    Parameters
    ----------
    df:
        DataFrame to save.
    path:
        Output path.
    encoding:
        Target encoding, default UTF-8.
    """
    df.to_csv(path, index=False, encoding=encoding)


def load_mapping(path: Path) -> Dict[str, str]:
    """Load mapping dictionary from JSON or TSV.

    Parameters
    ----------
    path:
        Path to JSON or TSV file.

    Returns
    -------
    dict
        Mapping loaded from file.
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
