"""Frequency calculation for MitoMaster."""

import logging
from typing import Optional

from .database import MitoMasterDB

logger = logging.getLogger(__name__)


class FrequencyCalculator:
    """Calculates and updates variant frequency statistics."""

    def __init__(self, db: MitoMasterDB):
        """Initialize calculator.

        Args:
            db: MitoMaster database connection
        """
        self.db = db

    def compute_all_frequencies(self) -> None:
        """Compute all frequency statistics from scratch.

        This rebuilds the variant_frequencies and variant_haplogroup_freq tables.
        """
        logger.info("Computing variant frequencies...")

        # Get total counts
        total_genomes = self.db.conn.execute(
            "SELECT COUNT(*) FROM genomes"
        ).fetchone()[0]
        human_genomes = self.db.conn.execute(
            "SELECT COUNT(*) FROM genomes WHERE is_human = true"
        ).fetchone()[0]

        if total_genomes == 0:
            logger.warning("No genomes in database, skipping frequency calculation")
            return

        logger.info(f"Total genomes: {total_genomes}, Human genomes: {human_genomes}")

        # Clear existing frequencies
        self.db.conn.execute("DELETE FROM variant_frequencies")
        self.db.conn.execute("DELETE FROM variant_haplogroup_freq")

        # Compute overall frequencies
        logger.info("Computing overall variant frequencies...")
        self.db.conn.execute(
            f"""
            INSERT INTO variant_frequencies (position, ref, alt, total_count, total_frequency, human_count, human_frequency)
            SELECT
                v.position,
                v.ref,
                v.alt,
                COUNT(*) as total_count,
                CAST(COUNT(*) AS DOUBLE) / {total_genomes} as total_frequency,
                SUM(CASE WHEN g.is_human THEN 1 ELSE 0 END) as human_count,
                CASE WHEN {human_genomes} > 0
                     THEN CAST(SUM(CASE WHEN g.is_human THEN 1 ELSE 0 END) AS DOUBLE) / {human_genomes}
                     ELSE 0 END as human_frequency
            FROM variants v
            JOIN genomes g ON v.accession = g.accession
            GROUP BY v.position, v.ref, v.alt
            """
        )

        freq_count = self.db.conn.execute(
            "SELECT COUNT(*) FROM variant_frequencies"
        ).fetchone()[0]
        logger.info(f"Computed frequencies for {freq_count} unique variants")

        # Compute haplogroup-specific frequencies
        logger.info("Computing haplogroup-specific frequencies...")
        self._compute_haplogroup_frequencies()

        logger.info("Frequency calculation complete")

    def _compute_haplogroup_frequencies(self) -> None:
        """Compute frequencies by top-level haplogroup."""
        # Get haplogroup counts (using first letter as top-level)
        haplogroup_counts = self.db.conn.execute(
            """
            SELECT
                CASE
                    WHEN haplogroup LIKE 'L%' THEN
                        CASE
                            WHEN haplogroup LIKE 'L0%' THEN 'L0'
                            WHEN haplogroup LIKE 'L1%' THEN 'L1'
                            WHEN haplogroup LIKE 'L2%' THEN 'L2'
                            WHEN haplogroup LIKE 'L3%' THEN 'L3'
                            WHEN haplogroup LIKE 'L4%' THEN 'L4'
                            WHEN haplogroup LIKE 'L5%' THEN 'L5'
                            WHEN haplogroup LIKE 'L6%' THEN 'L6'
                            ELSE 'L'
                        END
                    ELSE SUBSTRING(haplogroup, 1, 1)
                END as top_haplogroup,
                COUNT(*) as cnt
            FROM genomes
            WHERE haplogroup IS NOT NULL
            GROUP BY top_haplogroup
            """
        ).fetchall()

        hg_count_map = {row[0]: row[1] for row in haplogroup_counts}
        logger.info(f"Found {len(hg_count_map)} top-level haplogroups")

        # Compute frequencies for each haplogroup
        self.db.conn.execute(
            """
            INSERT INTO variant_haplogroup_freq (position, ref, alt, haplogroup, count, frequency)
            SELECT
                v.position,
                v.ref,
                v.alt,
                CASE
                    WHEN g.haplogroup LIKE 'L%' THEN
                        CASE
                            WHEN g.haplogroup LIKE 'L0%' THEN 'L0'
                            WHEN g.haplogroup LIKE 'L1%' THEN 'L1'
                            WHEN g.haplogroup LIKE 'L2%' THEN 'L2'
                            WHEN g.haplogroup LIKE 'L3%' THEN 'L3'
                            WHEN g.haplogroup LIKE 'L4%' THEN 'L4'
                            WHEN g.haplogroup LIKE 'L5%' THEN 'L5'
                            WHEN g.haplogroup LIKE 'L6%' THEN 'L6'
                            ELSE 'L'
                        END
                    ELSE SUBSTRING(g.haplogroup, 1, 1)
                END as top_haplogroup,
                COUNT(*) as cnt,
                0.0 as frequency  -- Will update below
            FROM variants v
            JOIN genomes g ON v.accession = g.accession
            WHERE g.haplogroup IS NOT NULL
            GROUP BY v.position, v.ref, v.alt, top_haplogroup
            """
        )

        # Update frequencies based on haplogroup counts
        for hg, total in hg_count_map.items():
            if total > 0:
                self.db.conn.execute(
                    """
                    UPDATE variant_haplogroup_freq
                    SET frequency = CAST(count AS DOUBLE) / ?
                    WHERE haplogroup = ?
                    """,
                    (total, hg),
                )

        hg_freq_count = self.db.conn.execute(
            "SELECT COUNT(*) FROM variant_haplogroup_freq"
        ).fetchone()[0]
        logger.info(f"Computed {hg_freq_count} haplogroup-specific frequency entries")

    def update_frequencies_incremental(self, accession: str) -> None:
        """Update frequencies incrementally after adding a genome.

        This is more efficient than recomputing all frequencies.

        Args:
            accession: Accession ID of newly added genome
        """
        # Get genome info
        genome = self.db.conn.execute(
            "SELECT is_human, haplogroup FROM genomes WHERE accession = ?",
            (accession,),
        ).fetchone()

        if not genome:
            return

        is_human, haplogroup = genome

        # Get new variant counts
        new_variants = self.db.conn.execute(
            "SELECT position, ref, alt FROM variants WHERE accession = ?",
            (accession,),
        ).fetchall()

        # Get current totals
        total_genomes = self.db.conn.execute(
            "SELECT COUNT(*) FROM genomes"
        ).fetchone()[0]
        human_genomes = self.db.conn.execute(
            "SELECT COUNT(*) FROM genomes WHERE is_human = true"
        ).fetchone()[0]

        for pos, ref, alt in new_variants:
            # Update or insert overall frequency
            existing = self.db.conn.execute(
                """
                SELECT total_count, human_count
                FROM variant_frequencies
                WHERE position = ? AND ref = ? AND alt = ?
                """,
                (pos, ref, alt),
            ).fetchone()

            if existing:
                new_total = existing[0] + 1
                new_human = existing[1] + (1 if is_human else 0)
                self.db.conn.execute(
                    """
                    UPDATE variant_frequencies
                    SET total_count = ?,
                        total_frequency = CAST(? AS DOUBLE) / ?,
                        human_count = ?,
                        human_frequency = CASE WHEN ? > 0 THEN CAST(? AS DOUBLE) / ? ELSE 0 END
                    WHERE position = ? AND ref = ? AND alt = ?
                    """,
                    (
                        new_total,
                        new_total,
                        total_genomes,
                        new_human,
                        human_genomes,
                        new_human,
                        human_genomes,
                        pos,
                        ref,
                        alt,
                    ),
                )
            else:
                self.db.conn.execute(
                    """
                    INSERT INTO variant_frequencies
                    (position, ref, alt, total_count, total_frequency, human_count, human_frequency)
                    VALUES (?, ?, ?, 1, CAST(1 AS DOUBLE) / ?, ?, CASE WHEN ? > 0 THEN CAST(? AS DOUBLE) / ? ELSE 0 END)
                    """,
                    (
                        pos,
                        ref,
                        alt,
                        total_genomes,
                        1 if is_human else 0,
                        human_genomes,
                        1 if is_human else 0,
                        human_genomes,
                    ),
                )

            # Update haplogroup frequency if haplogroup is known
            if haplogroup:
                # Get top-level haplogroup
                if haplogroup.startswith("L"):
                    for i in range(7):
                        if haplogroup.startswith(f"L{i}"):
                            top_hg = f"L{i}"
                            break
                    else:
                        top_hg = "L"
                else:
                    top_hg = haplogroup[0] if haplogroup else None

                if top_hg:
                    # Get haplogroup total
                    hg_total = self.db.conn.execute(
                        """
                        SELECT COUNT(*) FROM genomes
                        WHERE haplogroup LIKE ? || '%'
                        """,
                        (top_hg,),
                    ).fetchone()[0]

                    existing_hg = self.db.conn.execute(
                        """
                        SELECT count FROM variant_haplogroup_freq
                        WHERE position = ? AND ref = ? AND alt = ? AND haplogroup = ?
                        """,
                        (pos, ref, alt, top_hg),
                    ).fetchone()

                    if existing_hg:
                        new_count = existing_hg[0] + 1
                        self.db.conn.execute(
                            """
                            UPDATE variant_haplogroup_freq
                            SET count = ?, frequency = CAST(? AS DOUBLE) / ?
                            WHERE position = ? AND ref = ? AND alt = ? AND haplogroup = ?
                            """,
                            (new_count, new_count, hg_total, pos, ref, alt, top_hg),
                        )
                    else:
                        self.db.conn.execute(
                            """
                            INSERT INTO variant_haplogroup_freq
                            (position, ref, alt, haplogroup, count, frequency)
                            VALUES (?, ?, ?, ?, 1, CAST(1 AS DOUBLE) / ?)
                            """,
                            (pos, ref, alt, top_hg, hg_total),
                        )
