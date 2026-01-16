"""VCF file reader for mtDNA variants."""

from pathlib import Path
from typing import Iterator, Optional

from haplogrep3.models import Sample, Polymorphism, Mutation


class VcfReader:
    """Reads mtDNA variants from VCF files.

    Supports both uncompressed (.vcf) and compressed (.vcf.gz) files.
    """

    # Valid mtDNA chromosome names
    MT_CHROMS = {"MT", "chrM", "chrMT", "M", "mitochondria", "mitochondrion"}

    def __init__(
        self,
        het_level: float = 0.9,
        chip: bool = False,
    ):
        """Initialize the VCF reader.

        Args:
            het_level: Minimum heteroplasmy level to include (default 0.9)
            chip: Whether data is from genotyping chip (limits to chip positions)
        """
        self.het_level = het_level
        self.chip = chip

    def read(self, path: Path) -> list[Sample]:
        """Read samples from a VCF file.

        Args:
            path: Path to VCF file

        Returns:
            List of Sample objects with polymorphisms
        """
        try:
            from cyvcf2 import VCF
        except ImportError:
            raise ImportError(
                "cyvcf2 is required for VCF parsing. Install with: pip install cyvcf2"
            )

        vcf = VCF(str(path))
        sample_names = vcf.samples

        # Initialize samples
        samples = {name: [] for name in sample_names}
        min_pos = 16569
        max_pos = 1

        for variant in vcf:
            # Check if this is mtDNA
            if variant.CHROM not in self.MT_CHROMS:
                continue

            pos = variant.POS
            ref = variant.REF
            alts = variant.ALT

            # Track range
            min_pos = min(min_pos, pos)
            max_pos = max(max_pos, pos)

            # Skip multiallelic variants
            if len(alts) != 1:
                continue

            alt = alts[0]

            # Skip indels for now (complex handling required)
            if len(ref) != 1 or len(alt) != 1:
                continue

            # Process each sample
            for i, sample_name in enumerate(sample_names):
                gt = variant.genotypes[i]

                # gt is [allele1, allele2, phased]
                allele1, allele2 = gt[0], gt[1]

                # Skip if no call or reference
                if allele1 == -1 or (allele1 == 0 and allele2 == 0):
                    continue

                # Check heteroplasmy level if available
                if allele1 != allele2:
                    # Heterozygous - check AF if available
                    af = self._get_af(variant, i)
                    if af is not None and af < self.het_level:
                        continue

                # Add polymorphism
                try:
                    poly = Polymorphism(
                        position=pos,
                        mutation=Mutation(alt.upper()),
                        reference=ref.upper(),
                    )
                    samples[sample_name].append(poly)
                except ValueError:
                    pass  # Skip invalid mutations

        vcf.close()

        # Create Sample objects
        result = []
        for name, polys in samples.items():
            sample = Sample(
                id=name,
                polymorphisms=polys,
                range_start=min_pos if polys else 1,
                range_end=max_pos if polys else 16569,
            )
            result.append(sample)

        return result

    def _get_af(self, variant, sample_idx: int) -> Optional[float]:
        """Get allele frequency for a sample.

        Args:
            variant: cyvcf2 Variant object
            sample_idx: Sample index

        Returns:
            Allele frequency or None if not available
        """
        # Try to get AF from FORMAT field
        try:
            af = variant.format("AF")
            if af is not None and len(af) > sample_idx:
                return float(af[sample_idx][0])
        except (KeyError, TypeError, IndexError):
            pass

        return None


def read_vcf(path: Path, het_level: float = 0.9, chip: bool = False) -> list[Sample]:
    """Convenience function to read a VCF file.

    Args:
        path: Path to VCF file
        het_level: Minimum heteroplasmy level
        chip: Whether data is from genotyping chip

    Returns:
        List of Sample objects
    """
    reader = VcfReader(het_level=het_level, chip=chip)
    return reader.read(path)
