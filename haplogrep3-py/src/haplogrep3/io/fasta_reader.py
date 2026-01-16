"""FASTA file reader for mtDNA sequences."""

from pathlib import Path
from typing import Optional

from haplogrep3.models import Sample, Polymorphism, Mutation


# rCRS (revised Cambridge Reference Sequence) - mtDNA reference
# This is a simplified version; full implementation would load from file
RCRS_LENGTH = 16569


class FastaReader:
    """Reads mtDNA sequences from FASTA files and extracts polymorphisms.

    Compares sequences against the rCRS reference to identify variants.
    """

    def __init__(
        self,
        reference_path: Optional[Path] = None,
        skip_alignment_rules: bool = False,
    ):
        """Initialize the FASTA reader.

        Args:
            reference_path: Path to reference FASTA (rCRS)
            skip_alignment_rules: Skip nomenclature fixing rules
        """
        self.reference_path = reference_path
        self.skip_alignment_rules = skip_alignment_rules
        self._reference: Optional[str] = None

    @property
    def reference(self) -> str:
        """Get the reference sequence."""
        if self._reference is None:
            if self.reference_path and self.reference_path.exists():
                self._reference = self._load_fasta_sequence(self.reference_path)
            else:
                raise ValueError(
                    "Reference sequence not loaded. Provide reference_path."
                )
        return self._reference

    def read(self, path: Path) -> list[Sample]:
        """Read samples from a FASTA file.

        Args:
            path: Path to FASTA file

        Returns:
            List of Sample objects with polymorphisms
        """
        try:
            from Bio import SeqIO
        except ImportError:
            raise ImportError(
                "biopython is required for FASTA parsing. Install with: pip install biopython"
            )

        samples = []

        for record in SeqIO.parse(str(path), "fasta"):
            sample_id = record.id
            sequence = str(record.seq).upper()

            # Extract polymorphisms by comparing to reference
            polys = self._extract_polymorphisms(sequence)

            sample = Sample(
                id=sample_id,
                polymorphisms=polys,
                range_start=1,
                range_end=len(sequence),
            )
            samples.append(sample)

        return samples

    def _extract_polymorphisms(self, sequence: str) -> list[Polymorphism]:
        """Extract polymorphisms by comparing sequence to reference.

        Args:
            sequence: Sample mtDNA sequence

        Returns:
            List of polymorphisms (differences from reference)
        """
        polys = []
        ref = self.reference

        # Compare position by position
        for i, (ref_base, sample_base) in enumerate(zip(ref, sequence)):
            pos = i + 1  # 1-based position

            # Skip if same as reference
            if ref_base.upper() == sample_base.upper():
                continue

            # Skip N's (unknown bases)
            if sample_base.upper() == 'N':
                continue

            # Skip gaps in sample
            if sample_base == '-':
                continue

            try:
                poly = Polymorphism(
                    position=pos,
                    mutation=Mutation(sample_base.upper()),
                    reference=ref_base.upper(),
                )
                polys.append(poly)
            except ValueError:
                pass  # Skip invalid mutations

        return polys

    def _load_fasta_sequence(self, path: Path) -> str:
        """Load a single sequence from a FASTA file.

        Args:
            path: Path to FASTA file

        Returns:
            Sequence string
        """
        try:
            from Bio import SeqIO
        except ImportError:
            raise ImportError(
                "biopython is required for FASTA parsing. Install with: pip install biopython"
            )

        record = next(SeqIO.parse(str(path), "fasta"))
        return str(record.seq).upper()


def read_fasta(
    path: Path,
    reference_path: Optional[Path] = None,
    skip_alignment_rules: bool = False,
) -> list[Sample]:
    """Convenience function to read a FASTA file.

    Args:
        path: Path to FASTA file
        reference_path: Path to reference FASTA
        skip_alignment_rules: Skip nomenclature rules

    Returns:
        List of Sample objects
    """
    reader = FastaReader(
        reference_path=reference_path,
        skip_alignment_rules=skip_alignment_rules,
    )
    return reader.read(path)
