"""Validation helpers for target normalization."""

from __future__ import annotations

import pandas as pd


def ensure_column(df: pd.DataFrame, column: str) -> None:
    """Ensure DataFrame has required column.

    Raises
    ------
    KeyError
        If column is missing.
    """
    if column not in df.columns:
        raise KeyError(f"Column '{column}' is required")
