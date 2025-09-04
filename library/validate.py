"""Validation helpers for target normalization."""

from __future__ import annotations

import logging

import pandas as pd

from .uniprot import extract_names, fetch_uniprot_record

logger = logging.getLogger(__name__)


def ensure_column(df: pd.DataFrame, column: str) -> None:
    """Ensure DataFrame has required column.

    Raises
    ------
    KeyError
        If column is missing.
    """
    if column not in df.columns:
        raise KeyError(f"Column '{column}' is required")


def validate_uniprot_name(uniprot_id: str, name: str) -> bool:
    """Check whether ``name`` matches the UniProt record for ``uniprot_id``.

    Parameters
    ----------
    uniprot_id:
        UniProt accession.
    name:
        Protein or gene name to validate.

    Returns
    -------
    bool
        ``True`` if the provided name matches the canonical protein name or any
        gene names/synonyms in UniProt, otherwise ``False``.
    """
    record = fetch_uniprot_record(uniprot_id)
    protein_name, gene_names = extract_names(record)
    candidates = {protein_name.lower(), *(g.lower() for g in gene_names)}
    return name.lower() in candidates


def validate_uniprot_dataframe(
    df: pd.DataFrame, id_column: str, name_column: str
) -> pd.DataFrame:
    """Validate UniProt mappings for an entire DataFrame.

    Parameters
    ----------
    df:
        Input DataFrame containing UniProt IDs and protein names.
    id_column:
        Column with UniProt accessions.
    name_column:
        Column with names to validate.

    Returns
    -------
    pandas.DataFrame
        Copy of ``df`` with an added boolean ``uniprot_match`` column.
    """
    ensure_column(df, id_column)
    ensure_column(df, name_column)
    matches: list[bool] = []
    for row in df.itertuples(index=False):
        uniprot_id = getattr(row, id_column)
        name = getattr(row, name_column)
        try:
            matches.append(validate_uniprot_name(str(uniprot_id), str(name)))
        except Exception as exc:  # pragma: no cover - network/parse issues
            logger.warning("Validation failed for %s: %s", uniprot_id, exc)
            matches.append(False)
    out = df.copy()
    out["uniprot_match"] = matches
    return out
