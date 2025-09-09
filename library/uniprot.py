"""UniProt client utilities.

This module provides helper functions for retrieving protein information
from the UniProt REST API.
"""

from __future__ import annotations

import logging
import asyncio
from typing import Dict, Iterable, List, Tuple

import httpx
import requests

logger = logging.getLogger(__name__)


_CACHE: Dict[str, Dict] = {}


def fetch_uniprot_record(uniprot_id: str) -> Dict:
    """Fetch UniProt entry JSON by accession.

    The result is cached to avoid repeated network calls when the same
    identifier is requested multiple times.

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
    if uniprot_id in _CACHE:
        return _CACHE[uniprot_id]
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    logger.debug("Fetching UniProt record for %s", uniprot_id)
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    _CACHE[uniprot_id] = response.json()
    return _CACHE[uniprot_id]


async def _fetch_single(
    client: httpx.AsyncClient, uniprot_id: str, retries: int = 3
) -> None:
    """Fetch a single UniProt record with retries and cache it."""

    if uniprot_id in _CACHE:
        return
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    for attempt in range(retries):
        try:
            resp = await client.get(url, timeout=10)
            resp.raise_for_status()
            _CACHE[uniprot_id] = resp.json()
            return
        except httpx.HTTPError as exc:  # pragma: no cover - network issues
            if attempt < retries - 1:
                await asyncio.sleep(2**attempt)
            else:
                logger.warning("Failed to fetch %s: %s", uniprot_id, exc)
                _CACHE[uniprot_id] = {}


async def fetch_uniprot_records(
    uniprot_ids: Iterable[str], *, concurrency: int = 10
) -> Dict[str, Dict]:
    """Fetch multiple UniProt records concurrently.

    Parameters
    ----------
    uniprot_ids:
        Iterable of UniProt accession identifiers.
    concurrency:
        Maximum number of simultaneous connections.

    Returns
    -------
    dict[str, dict]
        Mapping of UniProt ID to the retrieved record. Failed requests return
        an empty dict.
    """

    ids = [uid for uid in dict.fromkeys(uniprot_ids) if uid not in _CACHE]
    if not ids:
        return {uid: _CACHE.get(uid, {}) for uid in uniprot_ids}
    semaphore = asyncio.Semaphore(concurrency)

    async with httpx.AsyncClient() as client:

        async def worker(uid: str) -> None:
            async with semaphore:
                await _fetch_single(client, uid)

        await asyncio.gather(*(worker(uid) for uid in ids))
    return {uid: _CACHE.get(uid, {}) for uid in uniprot_ids}


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
