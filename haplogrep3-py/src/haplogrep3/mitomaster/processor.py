"""GenBank to variants processor for MitoMaster."""

import logging
import re
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq

from haplogrep3.io import PhylotreeLoader
from haplogrep3.models import Sample, Polymorphism, Mutation
from haplogrep3.tasks import ClassificationTask
from haplogrep3.distance import Distance

from .database import MitoMasterDB, GenomeRecord, VariantRecord

logger = logging.getLogger(__name__)

# rCRS reference sequence (standard mtDNA reference)
RCRS_PATH = Path(__file__).parent.parent.parent.parent / "data" / "rCRS.fasta"


class GenomeProcessor:
    """Processes GenBank files into variant records."""

    # Transitions
    TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    def __init__(
        self,
        db: MitoMasterDB,
        tree_id: str = "phylotree-rcrs@17.2",
        rcrs_path: Optional[Path] = None,
    ):
        """Initialize processor.

        Args:
            db: MitoMaster database connection
            tree_id: Phylogenetic tree to use for classification
            rcrs_path: Path to rCRS reference FASTA
        """
        self.db = db
        self.tree_id = tree_id

        # Load phylotree
        loader = PhylotreeLoader()
        self.phylotree = loader.load(tree_id)

        # Load rCRS reference
        self.rcrs = self._load_rcrs(rcrs_path)

    def _load_rcrs(self, rcrs_path: Optional[Path] = None) -> str:
        """Load rCRS reference sequence.

        Args:
            rcrs_path: Path to rCRS FASTA file

        Returns:
            rCRS sequence string
        """
        if rcrs_path is None:
            # Try to find rCRS in tree data
            if self.phylotree.reference_fasta:
                rcrs_path = Path(self.phylotree.reference_fasta)
            else:
                rcrs_path = RCRS_PATH

        if not rcrs_path.exists():
            raise FileNotFoundError(f"rCRS reference not found: {rcrs_path}")

        with open(rcrs_path) as f:
            record = next(SeqIO.parse(f, "fasta"))
            return str(record.seq).upper()

    def process_genbank(self, gb_path: Path) -> tuple[GenomeRecord, list[VariantRecord]]:
        """Process a GenBank file into genome and variant records.

        Args:
            gb_path: Path to GenBank file

        Returns:
            Tuple of (GenomeRecord, list of VariantRecords)
        """
        # Parse GenBank
        record = SeqIO.read(str(gb_path), "genbank")
        sequence = str(record.seq).upper()

        # Extract metadata
        accession = record.id
        organism = record.annotations.get("organism", "Unknown")
        taxonomy = ";".join(record.annotations.get("taxonomy", []))
        is_human = "homo sapiens" in organism.lower()

        # Extract country and collection date from features
        country = None
        collection_date = None
        for feature in record.features:
            if feature.type == "source":
                country = feature.qualifiers.get("country", [None])[0]
                collection_date = feature.qualifiers.get("collection_date", [None])[0]
                break

        # Classify with haplogrep3
        haplogroup = None
        quality = None

        if is_human and len(sequence) >= 16000:
            try:
                # Convert sequence to polymorphisms vs rCRS
                polys = self._sequence_to_polymorphisms(sequence)
                sample = Sample(
                    id=accession,
                    polymorphisms=polys,
                    range_start=1,
                    range_end=len(sequence),
                )

                # Run classification
                task = ClassificationTask(self.phylotree, Distance.KULCZYNSKI, 1)
                classified = task.classify([sample])

                if classified and classified[0].classification:
                    haplogroup = classified[0].classification.haplogroup.name
                    quality = classified[0].classification.quality

            except Exception as e:
                logger.warning(f"Classification failed for {accession}: {e}")

        # Create genome record
        genome = GenomeRecord(
            accession=accession,
            organism=organism,
            taxonomy=taxonomy,
            length=len(sequence),
            haplogroup=haplogroup,
            haplogroup_quality=quality,
            is_human=is_human,
            country=country,
            collection_date=collection_date,
        )

        # Extract variants vs rCRS
        variants = self._extract_variants(accession, sequence)

        return genome, variants

    def _sequence_to_polymorphisms(self, sequence: str) -> list[Polymorphism]:
        """Convert sequence to list of polymorphisms vs rCRS.

        Args:
            sequence: Query sequence

        Returns:
            List of Polymorphism objects
        """
        polys = []
        min_len = min(len(sequence), len(self.rcrs))

        for i in range(min_len):
            ref = self.rcrs[i]
            alt = sequence[i]

            if ref != alt and alt in "ACGT":
                try:
                    poly = Polymorphism(
                        position=i + 1,  # 1-based
                        mutation=Mutation(alt),
                        reference=ref,
                    )
                    polys.append(poly)
                except ValueError:
                    pass  # Skip invalid mutations

        return polys

    def _extract_variants(
        self, accession: str, sequence: str
    ) -> list[VariantRecord]:
        """Extract all variants from sequence compared to rCRS.

        Args:
            accession: Genome accession ID
            sequence: Query sequence

        Returns:
            List of VariantRecord objects
        """
        variants = []
        min_len = min(len(sequence), len(self.rcrs))

        for i in range(min_len):
            ref = self.rcrs[i]
            alt = sequence[i]
            position = i + 1  # 1-based

            if ref != alt:
                # Determine mutation type
                if alt == "-" or alt == "N":
                    mutation_type = "deletion"
                elif ref == "-":
                    mutation_type = "insertion"
                elif (ref, alt) in self.TRANSITIONS:
                    mutation_type = "transition"
                else:
                    mutation_type = "transversion"

                # Get gene info
                gene_info = self.db.get_gene_for_position(position)
                locus = gene_info[0] if gene_info else None
                gene_type = gene_info[1] if gene_info else None

                # Calculate amino acid change for coding regions
                aa_change = None
                if gene_type == "protein_coding" and alt in "ACGT":
                    aa_change = self._get_amino_acid_change(
                        position, ref, alt, locus, sequence
                    )

                variant = VariantRecord(
                    accession=accession,
                    position=position,
                    ref=ref,
                    alt=alt,
                    mutation_type=mutation_type,
                    locus=locus,
                    gene_type=gene_type,
                    amino_acid_change=aa_change,
                )
                variants.append(variant)

        return variants

    def _get_amino_acid_change(
        self,
        position: int,
        ref: str,
        alt: str,
        gene: str,
        sequence: str,
    ) -> Optional[str]:
        """Calculate amino acid change for a variant in coding region.

        Args:
            position: 1-based position
            ref: Reference base
            alt: Alternate base
            gene: Gene name
            sequence: Full sequence

        Returns:
            Amino acid change string (e.g., "M1T") or None
        """
        # Get gene coordinates
        gene_info = None
        for g in self.db.GENE_MAP:
            if g[0] == gene:
                gene_info = g
                break

        if not gene_info or gene_info[3] != "protein_coding":
            return None

        start, end = gene_info[1], gene_info[2]
        strand = gene_info[4]

        # Calculate codon position
        gene_pos = position - start  # 0-based position within gene
        codon_num = gene_pos // 3 + 1
        codon_pos = gene_pos % 3

        # Extract codon
        codon_start = start + (codon_num - 1) * 3 - 1  # 0-based
        if codon_start < 0 or codon_start + 3 > len(sequence):
            return None

        ref_codon = list(self.rcrs[codon_start : codon_start + 3])
        alt_codon = list(sequence[codon_start : codon_start + 3])

        # Handle reverse complement for L-strand genes
        if strand == "L":
            ref_codon = [self._complement(b) for b in reversed(ref_codon)]
            alt_codon = [self._complement(b) for b in reversed(alt_codon)]

        try:
            ref_codon_str = "".join(ref_codon)
            alt_codon_str = "".join(alt_codon)

            ref_aa = str(Seq(ref_codon_str).translate())
            alt_aa = str(Seq(alt_codon_str).translate())

            if ref_aa != alt_aa:
                return f"{ref_aa}{codon_num}{alt_aa}"
            else:
                return f"syn:{ref_aa}{codon_num}"  # Synonymous
        except Exception:
            return None

    def _complement(self, base: str) -> str:
        """Get complement of a base."""
        complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        return complements.get(base, "N")

    def process_and_store(self, gb_path: Path, delete_after: bool = True) -> bool:
        """Process a GenBank file and store results in database.

        Args:
            gb_path: Path to GenBank file
            delete_after: Whether to delete the file after processing

        Returns:
            True if successful
        """
        try:
            genome, variants = self.process_genbank(gb_path)

            # Store in database
            self.db.insert_genome(genome)
            self.db.insert_variants(variants)

            logger.info(
                f"Processed {genome.accession}: {len(variants)} variants, "
                f"haplogroup={genome.haplogroup}"
            )

            return True

        except Exception as e:
            logger.error(f"Failed to process {gb_path}: {e}")
            return False

        finally:
            if delete_after and gb_path.exists():
                gb_path.unlink()
