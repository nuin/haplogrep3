"""DuckDB database management for MitoMaster."""

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterator, Optional

import duckdb


@dataclass
class GenomeRecord:
    """Genome record for storage."""
    accession: str
    organism: str
    taxonomy: str
    length: int
    haplogroup: Optional[str] = None
    haplogroup_quality: Optional[float] = None
    is_human: bool = False
    country: Optional[str] = None
    collection_date: Optional[str] = None


@dataclass
class VariantRecord:
    """Variant record for storage."""
    accession: str
    position: int
    ref: str
    alt: str
    mutation_type: Optional[str] = None
    locus: Optional[str] = None
    gene_type: Optional[str] = None
    amino_acid_change: Optional[str] = None


@dataclass
class VariantFrequency:
    """Variant frequency result."""
    position: int
    ref: str
    alt: str
    total_count: int
    total_frequency: float
    human_count: int
    human_frequency: float
    locus: Optional[str] = None
    gene_type: Optional[str] = None
    amino_acid_change: Optional[str] = None


@dataclass
class HaplogroupFrequency:
    """Haplogroup-specific variant frequency."""
    haplogroup: str
    count: int
    frequency: float


class MitoMasterDB:
    """DuckDB database for MitoMaster variant data."""

    # Standard mtDNA gene map
    GENE_MAP = [
        ("D-loop-1", 1, 576, "control_region", "H", "Control region (HVS2)"),
        ("TRNF", 577, 647, "tRNA", "H", "tRNA-Phe"),
        ("RNR1", 648, 1601, "rRNA", "H", "12S ribosomal RNA"),
        ("TRNV", 1602, 1670, "tRNA", "H", "tRNA-Val"),
        ("RNR2", 1671, 3229, "rRNA", "H", "16S ribosomal RNA"),
        ("TRNL1", 3230, 3304, "tRNA", "H", "tRNA-Leu (UUR)"),
        ("ND1", 3307, 4262, "protein_coding", "H", "NADH dehydrogenase subunit 1"),
        ("TRNI", 4263, 4331, "tRNA", "H", "tRNA-Ile"),
        ("TRNQ", 4329, 4400, "tRNA", "L", "tRNA-Gln"),
        ("TRNM", 4402, 4469, "tRNA", "H", "tRNA-Met"),
        ("ND2", 4470, 5511, "protein_coding", "H", "NADH dehydrogenase subunit 2"),
        ("TRNW", 5512, 5579, "tRNA", "H", "tRNA-Trp"),
        ("TRNA", 5587, 5655, "tRNA", "L", "tRNA-Ala"),
        ("TRNN", 5657, 5729, "tRNA", "L", "tRNA-Asn"),
        ("TRNC", 5761, 5826, "tRNA", "L", "tRNA-Cys"),
        ("TRNY", 5826, 5891, "tRNA", "L", "tRNA-Tyr"),
        ("COX1", 5904, 7445, "protein_coding", "H", "Cytochrome c oxidase subunit 1"),
        ("TRNS1", 7446, 7514, "tRNA", "L", "tRNA-Ser (UCN)"),
        ("TRND", 7518, 7585, "tRNA", "H", "tRNA-Asp"),
        ("COX2", 7586, 8269, "protein_coding", "H", "Cytochrome c oxidase subunit 2"),
        ("TRNK", 8295, 8364, "tRNA", "H", "tRNA-Lys"),
        ("ATP8", 8366, 8572, "protein_coding", "H", "ATP synthase F0 subunit 8"),
        ("ATP6", 8527, 9207, "protein_coding", "H", "ATP synthase F0 subunit 6"),
        ("COX3", 9207, 9990, "protein_coding", "H", "Cytochrome c oxidase subunit 3"),
        ("TRNG", 9991, 10058, "tRNA", "H", "tRNA-Gly"),
        ("ND3", 10059, 10404, "protein_coding", "H", "NADH dehydrogenase subunit 3"),
        ("TRNR", 10405, 10469, "tRNA", "H", "tRNA-Arg"),
        ("ND4L", 10470, 10766, "protein_coding", "H", "NADH dehydrogenase subunit 4L"),
        ("ND4", 10760, 12137, "protein_coding", "H", "NADH dehydrogenase subunit 4"),
        ("TRNH", 12138, 12206, "tRNA", "H", "tRNA-His"),
        ("TRNS2", 12207, 12265, "tRNA", "H", "tRNA-Ser (AGY)"),
        ("TRNL2", 12266, 12336, "tRNA", "H", "tRNA-Leu (CUN)"),
        ("ND5", 12337, 14148, "protein_coding", "H", "NADH dehydrogenase subunit 5"),
        ("ND6", 14149, 14673, "protein_coding", "L", "NADH dehydrogenase subunit 6"),
        ("TRNE", 14674, 14742, "tRNA", "L", "tRNA-Glu"),
        ("CYTB", 14747, 15887, "protein_coding", "H", "Cytochrome b"),
        ("TRNT", 15888, 15953, "tRNA", "H", "tRNA-Thr"),
        ("TRNP", 15956, 16023, "tRNA", "L", "tRNA-Pro"),
        ("D-loop-2", 16024, 16569, "control_region", "H", "Control region (HVS1)"),
    ]

    def __init__(self, db_path: Path):
        """Initialize database connection.

        Args:
            db_path: Path to DuckDB file
        """
        self.db_path = Path(db_path)
        self._conn: Optional[duckdb.DuckDBPyConnection] = None

    @property
    def conn(self) -> duckdb.DuckDBPyConnection:
        """Get or create database connection."""
        if self._conn is None:
            self._conn = duckdb.connect(str(self.db_path))
        return self._conn

    def close(self) -> None:
        """Close database connection."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def initialize(self) -> None:
        """Initialize database schema."""
        schema_path = Path(__file__).parent / "schema.sql"
        with open(schema_path) as f:
            schema_sql = f.read()

        # Execute schema
        for statement in schema_sql.split(";"):
            statement = statement.strip()
            if statement:
                self.conn.execute(statement)

        # Populate gene map if empty
        count = self.conn.execute("SELECT COUNT(*) FROM gene_map").fetchone()[0]
        if count == 0:
            self._populate_gene_map()

    def _populate_gene_map(self) -> None:
        """Populate the gene map table."""
        self.conn.executemany(
            """
            INSERT INTO gene_map (name, start_pos, end_pos, gene_type, strand, product)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            self.GENE_MAP,
        )

    # -------------------------------------------------------------------------
    # Insert operations
    # -------------------------------------------------------------------------

    def insert_genome(self, genome: GenomeRecord) -> None:
        """Insert a genome record."""
        self.conn.execute(
            """
            INSERT OR REPLACE INTO genomes
            (accession, organism, taxonomy, length, haplogroup, haplogroup_quality,
             is_human, country, collection_date, processed_date)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                genome.accession,
                genome.organism,
                genome.taxonomy,
                genome.length,
                genome.haplogroup,
                genome.haplogroup_quality,
                genome.is_human,
                genome.country,
                genome.collection_date,
                datetime.now(),
            ),
        )

    def insert_variants(self, variants: list[VariantRecord]) -> None:
        """Insert multiple variant records."""
        if not variants:
            return

        self.conn.executemany(
            """
            INSERT INTO variants
            (id, accession, position, ref, alt, mutation_type, locus, gene_type, amino_acid_change)
            VALUES (nextval('variant_id_seq'), ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    v.accession,
                    v.position,
                    v.ref,
                    v.alt,
                    v.mutation_type,
                    v.locus,
                    v.gene_type,
                    v.amino_acid_change,
                )
                for v in variants
            ],
        )

    def genome_exists(self, accession: str) -> bool:
        """Check if a genome already exists in the database."""
        result = self.conn.execute(
            "SELECT 1 FROM genomes WHERE accession = ?", (accession,)
        ).fetchone()
        return result is not None

    # -------------------------------------------------------------------------
    # Query operations
    # -------------------------------------------------------------------------

    def get_variant_at_position(
        self, position: int, ref: Optional[str] = None, alt: Optional[str] = None
    ) -> list[VariantFrequency]:
        """Get variant frequencies at a position.

        Args:
            position: mtDNA position (1-16569)
            ref: Optional reference base filter
            alt: Optional alternate base filter

        Returns:
            List of variant frequencies
        """
        query = """
            SELECT
                vf.position, vf.ref, vf.alt,
                vf.total_count, vf.total_frequency,
                vf.human_count, vf.human_frequency,
                v.locus, v.gene_type, v.amino_acid_change
            FROM variant_frequencies vf
            LEFT JOIN (
                SELECT DISTINCT position, ref, alt, locus, gene_type, amino_acid_change
                FROM variants
                WHERE position = ?
            ) v ON vf.position = v.position AND vf.ref = v.ref AND vf.alt = v.alt
            WHERE vf.position = ?
        """
        params = [position, position]

        if ref:
            query += " AND vf.ref = ?"
            params.append(ref)
        if alt:
            query += " AND vf.alt = ?"
            params.append(alt)

        query += " ORDER BY vf.total_count DESC"

        results = self.conn.execute(query, params).fetchall()
        return [
            VariantFrequency(
                position=r[0],
                ref=r[1],
                alt=r[2],
                total_count=r[3],
                total_frequency=r[4],
                human_count=r[5],
                human_frequency=r[6],
                locus=r[7],
                gene_type=r[8],
                amino_acid_change=r[9],
            )
            for r in results
        ]

    def get_haplogroup_frequencies(
        self, position: int, ref: str, alt: str
    ) -> list[HaplogroupFrequency]:
        """Get haplogroup-specific frequencies for a variant.

        Args:
            position: mtDNA position
            ref: Reference base
            alt: Alternate base

        Returns:
            List of haplogroup frequencies
        """
        results = self.conn.execute(
            """
            SELECT haplogroup, count, frequency
            FROM variant_haplogroup_freq
            WHERE position = ? AND ref = ? AND alt = ?
            ORDER BY count DESC
            """,
            (position, ref, alt),
        ).fetchall()

        return [
            HaplogroupFrequency(haplogroup=r[0], count=r[1], frequency=r[2])
            for r in results
        ]

    def search_variants(
        self,
        position_start: Optional[int] = None,
        position_end: Optional[int] = None,
        gene: Optional[str] = None,
        haplogroup: Optional[str] = None,
        min_frequency: Optional[float] = None,
        limit: int = 100,
    ) -> list[VariantFrequency]:
        """Search variants with filters.

        Args:
            position_start: Start position filter
            position_end: End position filter
            gene: Gene name filter
            haplogroup: Haplogroup filter (for haplogroup-specific frequency)
            min_frequency: Minimum frequency filter
            limit: Maximum results

        Returns:
            List of matching variant frequencies
        """
        conditions = []
        params = []

        if position_start is not None:
            conditions.append("vf.position >= ?")
            params.append(position_start)
        if position_end is not None:
            conditions.append("vf.position <= ?")
            params.append(position_end)
        if min_frequency is not None:
            conditions.append("vf.total_frequency >= ?")
            params.append(min_frequency)

        # Handle gene filter
        if gene:
            gene_row = self.conn.execute(
                "SELECT start_pos, end_pos FROM gene_map WHERE name = ?", (gene,)
            ).fetchone()
            if gene_row:
                conditions.append("vf.position BETWEEN ? AND ?")
                params.extend([gene_row[0], gene_row[1]])

        where_clause = " AND ".join(conditions) if conditions else "1=1"

        query = f"""
            SELECT DISTINCT
                vf.position, vf.ref, vf.alt,
                vf.total_count, vf.total_frequency,
                vf.human_count, vf.human_frequency,
                v.locus, v.gene_type, v.amino_acid_change
            FROM variant_frequencies vf
            LEFT JOIN (
                SELECT DISTINCT position, ref, alt, locus, gene_type, amino_acid_change
                FROM variants
            ) v ON vf.position = v.position AND vf.ref = v.ref AND vf.alt = v.alt
            WHERE {where_clause}
            ORDER BY vf.total_count DESC
            LIMIT ?
        """
        params.append(limit)

        results = self.conn.execute(query, params).fetchall()
        return [
            VariantFrequency(
                position=r[0],
                ref=r[1],
                alt=r[2],
                total_count=r[3],
                total_frequency=r[4],
                human_count=r[5],
                human_frequency=r[6],
                locus=r[7],
                gene_type=r[8],
                amino_acid_change=r[9],
            )
            for r in results
        ]

    def get_gene_for_position(self, position: int) -> Optional[tuple[str, str]]:
        """Get gene name and type for a position.

        Args:
            position: mtDNA position

        Returns:
            Tuple of (gene_name, gene_type) or None
        """
        result = self.conn.execute(
            """
            SELECT name, gene_type FROM gene_map
            WHERE start_pos <= ? AND end_pos >= ?
            """,
            (position, position),
        ).fetchone()
        return (result[0], result[1]) if result else None

    def get_stats(self) -> dict:
        """Get database statistics.

        Returns:
            Dictionary with counts and metadata
        """
        genome_count = self.conn.execute("SELECT COUNT(*) FROM genomes").fetchone()[0]
        human_count = self.conn.execute(
            "SELECT COUNT(*) FROM genomes WHERE is_human = true"
        ).fetchone()[0]
        variant_count = self.conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
        unique_variants = self.conn.execute(
            "SELECT COUNT(*) FROM variant_frequencies"
        ).fetchone()[0]

        haplogroups = self.conn.execute(
            """
            SELECT haplogroup, COUNT(*) as cnt
            FROM genomes
            WHERE haplogroup IS NOT NULL
            GROUP BY haplogroup
            ORDER BY cnt DESC
            LIMIT 20
            """
        ).fetchall()

        return {
            "genome_count": genome_count,
            "human_genome_count": human_count,
            "variant_count": variant_count,
            "unique_variant_count": unique_variants,
            "top_haplogroups": [{"haplogroup": h[0], "count": h[1]} for h in haplogroups],
        }

    def get_all_genes(self) -> list[dict]:
        """Get all genes from gene map.

        Returns:
            List of gene dictionaries
        """
        results = self.conn.execute(
            "SELECT name, start_pos, end_pos, gene_type, strand, product FROM gene_map ORDER BY start_pos"
        ).fetchall()
        return [
            {
                "name": r[0],
                "start": r[1],
                "end": r[2],
                "type": r[3],
                "strand": r[4],
                "product": r[5],
            }
            for r in results
        ]
