"""UniProt client utilities.

This module provides helper functions for retrieving protein information
from the UniProt REST API.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from typing import Dict, List, Tuple

import requests

logger = logging.getLogger(__name__)


@lru_cache(maxsize=128)
def fetch_uniprot_record(uniprot_id: str) -> Dict:
    """Fetch UniProt entry JSON by accession.

    Parameters
    ----------
    uniprot_id:
        UniProt accession identifier (e.g., ``"P68871"``).

    Returns
    -------
    dict
        Parsed JSON payload returned by the UniProt REST API.

    Raises
    ------
    requests.HTTPError
        If the request fails.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    logger.debug("Fetching UniProt record for %s", uniprot_id)
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    return response.json()


def extract_names(record: Dict) -> Tuple[str, List[str]]:
    """Extract canonical protein and gene names from a UniProt record.

    Parameters
    ----------
    record:
        JSON object representing the UniProt entry.

    Returns
    -------
    tuple[str, list[str]]
        Canonical protein name and a list of gene names including synonyms.
    """
    protein_name = record["proteinDescription"]["recommendedName"]["fullName"]["value"]
    gene_names: List[str] = []
    for gene in record.get("genes", []):
        gene_name = gene.get("geneName", {}).get("value")
        if gene_name:
            gene_names.append(gene_name)
        for syn in gene.get("synonyms", []):
            syn_val = syn.get("value")
            if syn_val:
                gene_names.append(syn_val)
    return protein_name, gene_names
