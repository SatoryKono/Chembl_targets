"""Validation helpers for target normalization."""

from __future__ import annotations

import pandas as pd  # type: ignore[import-untyped]


def ensure_column(df: pd.DataFrame, column: str) -> None:
    """Ensure DataFrame has a required column.

    Args:
        df: The DataFrame to check.
        column: The name of the column to ensure exists.

    Raises:
        KeyError: If the specified column is missing from the DataFrame.
    """
    if column not in df.columns:
        raise KeyError(f"Column '{column}' is required")
