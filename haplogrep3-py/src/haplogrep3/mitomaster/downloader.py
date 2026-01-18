"""NCBI genome stream downloader for MitoMaster."""

import logging
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional, Callable

from Bio import Entrez

logger = logging.getLogger(__name__)


@dataclass
class DownloadConfig:
    """Configuration for NCBI downloads."""

    email: str
    api_key: Optional[str] = None
    batch_size: int = 100
    retry_attempts: int = 3
    retry_delay: float = 1.0
    rate_limit_delay: float = 0.34  # ~3 requests/sec without API key


class NCBIDownloader:
    """Stream downloader for NCBI mitochondrial genomes."""

    # Search query for complete mtDNA genomes
    HUMAN_MTDNA_QUERY = (
        "Homo sapiens[Organism] AND mitochondrion[Title] AND "
        "complete genome[Title] AND 16000:17000[Sequence Length] "
        "NOT ancient[Title] NOT Neanderthal[Title] NOT Denisova[Title]"
    )

    ALL_MTDNA_QUERY = (
        "mitochondrion[Title] AND complete genome[Title] AND "
        "16000:17000[Sequence Length]"
    )

    def __init__(self, config: DownloadConfig):
        """Initialize downloader with configuration.

        Args:
            config: Download configuration
        """
        self.config = config
        Entrez.email = config.email
        if config.api_key:
            Entrez.api_key = config.api_key

    def search_genomes(
        self,
        query: Optional[str] = None,
        organism: Optional[str] = None,
        max_results: Optional[int] = None,
    ) -> list[str]:
        """Search NCBI for genome accessions.

        Args:
            query: Custom search query (overrides organism)
            organism: Organism filter (e.g., "Homo sapiens")
            max_results: Maximum number of results

        Returns:
            List of accession IDs
        """
        if query is None:
            if organism and "homo sapiens" in organism.lower():
                query = self.HUMAN_MTDNA_QUERY
            else:
                query = self.ALL_MTDNA_QUERY
                if organism:
                    query = f"{organism}[Organism] AND {query}"

        logger.info(f"Searching NCBI with query: {query[:100]}...")

        # First, get count
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()

        total_count = int(record["Count"])
        logger.info(f"Found {total_count} genomes")

        if max_results:
            total_count = min(total_count, max_results)

        # Fetch all IDs in batches
        all_ids = []
        for start in range(0, total_count, self.config.batch_size):
            retmax = min(self.config.batch_size, total_count - start)
            handle = Entrez.esearch(
                db="nucleotide", term=query, retstart=start, retmax=retmax
            )
            record = Entrez.read(handle)
            handle.close()

            all_ids.extend(record["IdList"])
            time.sleep(self.config.rate_limit_delay)

            if len(all_ids) % 1000 == 0:
                logger.info(f"Retrieved {len(all_ids)} IDs...")

        logger.info(f"Total IDs retrieved: {len(all_ids)}")
        return all_ids

    def download_genbank(
        self, accession: str, output_path: Optional[Path] = None
    ) -> Path:
        """Download a single GenBank record.

        Args:
            accession: NCBI accession ID
            output_path: Optional output path (uses temp file if None)

        Returns:
            Path to downloaded GenBank file
        """
        for attempt in range(self.config.retry_attempts):
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=accession,
                    rettype="gb",
                    retmode="text",
                )
                content = handle.read()
                handle.close()

                if output_path is None:
                    fd, temp_path = tempfile.mkstemp(suffix=".gb")
                    output_path = Path(temp_path)
                    with open(fd, "w") as f:
                        f.write(content)
                else:
                    with open(output_path, "w") as f:
                        f.write(content)

                time.sleep(self.config.rate_limit_delay)
                return output_path

            except Exception as e:
                logger.warning(
                    f"Download attempt {attempt + 1} failed for {accession}: {e}"
                )
                if attempt < self.config.retry_attempts - 1:
                    time.sleep(self.config.retry_delay * (attempt + 1))
                else:
                    raise

        raise RuntimeError(f"Failed to download {accession} after all retries")

    def stream_genomes(
        self,
        accessions: list[str],
        progress_callback: Optional[Callable[[int, int], None]] = None,
        skip_existing: Optional[Callable[[str], bool]] = None,
    ) -> Iterator[tuple[str, Path]]:
        """Stream download genomes one at a time.

        Downloads to temp files that should be deleted after processing.

        Args:
            accessions: List of accession IDs to download
            progress_callback: Optional callback(current, total) for progress
            skip_existing: Optional callback to check if accession should be skipped

        Yields:
            Tuples of (accession, temp_file_path)
        """
        total = len(accessions)
        downloaded = 0
        skipped = 0

        for i, accession in enumerate(accessions):
            # Check if should skip
            if skip_existing and skip_existing(accession):
                skipped += 1
                if progress_callback:
                    progress_callback(i + 1, total)
                continue

            try:
                temp_path = self.download_genbank(accession)
                downloaded += 1
                yield (accession, temp_path)

            except Exception as e:
                logger.error(f"Failed to download {accession}: {e}")
                continue

            finally:
                if progress_callback:
                    progress_callback(i + 1, total)

        logger.info(
            f"Download complete: {downloaded} downloaded, {skipped} skipped, "
            f"{total - downloaded - skipped} failed"
        )

    def get_genome_count(
        self, query: Optional[str] = None, organism: Optional[str] = None
    ) -> int:
        """Get count of available genomes without downloading.

        Args:
            query: Custom search query
            organism: Organism filter

        Returns:
            Number of matching genomes
        """
        if query is None:
            if organism and "homo sapiens" in organism.lower():
                query = self.HUMAN_MTDNA_QUERY
            else:
                query = self.ALL_MTDNA_QUERY
                if organism:
                    query = f"{organism}[Organism] AND {query}"

        handle = Entrez.esearch(db="nucleotide", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()

        return int(record["Count"])
