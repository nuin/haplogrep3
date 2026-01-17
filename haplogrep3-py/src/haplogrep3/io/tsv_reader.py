"""TSV variant file reader for mtDNA variants."""

import csv
import re
from pathlib import Path
from typing import Optional

from haplogrep3.models import Sample, Polymorphism


class TsvReader:
    """Reader for TSV variant files with HGVS notation.

    Supports formats with columns containing mtDNA variant information,
    particularly the 'c. HGVS' column with notation like m.73A>G.
    """

    # Regex patterns for HGVS notation
    HGVS_SUB = re.compile(r'm\.(\d+)([ACGT])>([ACGT])')  # m.73A>G
    HGVS_INS = re.compile(r'm\.(\d+)ins([ACGT]+)')  # m.311insTC
    HGVS_DUP = re.compile(r'm\.(\d+)dup([ACGT]+)')  # m.316dupC
    HGVS_DEL = re.compile(r'm\.(\d+)del([ACGT]*)')  # m.3107delC

    def __init__(
        self,
        hgvs_column: str = "c. HGVS",
        sample_column: Optional[str] = None,
        het_threshold: float = 0.1,
        skip_hotspots: bool = True,
    ):
        """Initialize TSV reader.

        Args:
            hgvs_column: Column name containing HGVS notation
            sample_column: Column for sample ID (uses filename if None)
            het_threshold: Minimum heteroplasmy level to include variant
            skip_hotspots: Skip variants marked as HotSpot
        """
        self.hgvs_column = hgvs_column
        self.sample_column = sample_column
        self.het_threshold = het_threshold
        self.skip_hotspots = skip_hotspots

    def read(self, path: Path) -> list[Sample]:
        """Read samples from a TSV file.

        Args:
            path: Path to TSV file

        Returns:
            List of Sample objects
        """
        polymorphisms = []
        sample_id = path.stem

        with open(path, newline='', encoding='utf-8') as f:
            # Detect delimiter
            first_line = f.readline()
            f.seek(0)

            if '\t' in first_line:
                delimiter = '\t'
            else:
                delimiter = ','

            reader = csv.DictReader(f, delimiter=delimiter)

            # Find the HGVS column (handle variations in column names)
            hgvs_col = None
            filter_col = None
            coverage_col = None

            if reader.fieldnames:
                for col in reader.fieldnames:
                    col_lower = col.lower().strip()
                    if 'hgvs' in col_lower and hgvs_col is None:
                        hgvs_col = col
                    if col_lower in ('filter step', 'filter', 'hint'):
                        filter_col = col
                    if col_lower == 'coverage':
                        coverage_col = col

            if not hgvs_col:
                # Try to find by exact name
                hgvs_col = self.hgvs_column

            for row in reader:
                # Skip hotspots if configured
                if self.skip_hotspots and filter_col and row.get(filter_col):
                    if 'hotspot' in row[filter_col].lower():
                        continue

                # Get HGVS notation
                hgvs = row.get(hgvs_col, '').strip()
                if not hgvs:
                    continue

                # Parse coverage/heteroplasmy if available
                if coverage_col and row.get(coverage_col):
                    cov = row[coverage_col]
                    # Extract percentage like "100% (31979)" or "88% (4511)"
                    match = re.match(r'([\d.]+)%', cov)
                    if match:
                        het_level = float(match.group(1)) / 100
                        if het_level < self.het_threshold:
                            continue

                # Parse HGVS to polymorphism
                poly = self._parse_hgvs(hgvs)
                if poly:
                    polymorphisms.append(poly)

        # Calculate range from polymorphisms
        range_start = 1
        range_end = 16569
        if polymorphisms:
            positions = [p.position for p in polymorphisms]
            range_start = min(positions)
            range_end = max(positions)

        # Create sample
        sample = Sample(
            id=sample_id,
            polymorphisms=polymorphisms,
            range_start=range_start,
            range_end=range_end,
        )

        return [sample]

    def _parse_hgvs(self, hgvs: str) -> Optional[Polymorphism]:
        """Parse HGVS notation to Polymorphism.

        Args:
            hgvs: HGVS string like m.73A>G

        Returns:
            Polymorphism or None if parsing fails
        """
        hgvs = hgvs.strip()

        # Substitution: m.73A>G
        match = self.HGVS_SUB.match(hgvs)
        if match:
            pos, ref, alt = match.groups()
            return Polymorphism(position=int(pos), mutation=alt, reference=ref)

        # Insertion: m.311insTC
        match = self.HGVS_INS.match(hgvs)
        if match:
            pos, inserted = match.groups()
            # For insertions, we use .1, .2 notation for each inserted base
            # Return first inserted base with .1 position
            return Polymorphism(
                position=int(pos),
                mutation=inserted[0],
                insertion_position=1,
            )

        # Duplication: m.316dupC (treated as insertion)
        match = self.HGVS_DUP.match(hgvs)
        if match:
            pos, dup = match.groups()
            return Polymorphism(
                position=int(pos),
                mutation=dup[0],
                insertion_position=1,
            )

        # Deletion: m.3107delC
        match = self.HGVS_DEL.match(hgvs)
        if match:
            pos, deleted = match.groups()
            return Polymorphism(position=int(pos), mutation='d')

        return None


def read_tsv(path: Path, **kwargs) -> list[Sample]:
    """Convenience function to read TSV file.

    Args:
        path: Path to TSV file
        **kwargs: Arguments passed to TsvReader

    Returns:
        List of Sample objects
    """
    reader = TsvReader(**kwargs)
    return reader.read(path)
