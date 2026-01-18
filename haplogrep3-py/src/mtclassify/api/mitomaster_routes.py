"""MitoMaster API routes for variant database queries."""

import logging
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/mitomaster", tags=["mitomaster"])

# Default database path - can be overridden
_db_path: Optional[Path] = None


def set_database_path(path: Path) -> None:
    """Set the MitoMaster database path."""
    global _db_path
    _db_path = path


def get_database_path() -> Optional[Path]:
    """Get the MitoMaster database path."""
    return _db_path


def _get_db():
    """Get database connection, raising error if not configured."""
    from mtclassify.mitomaster import MitoMasterDB

    if _db_path is None or not _db_path.exists():
        raise HTTPException(
            status_code=503,
            detail="MitoMaster database not configured or not found. "
                   "Use CLI 'mtclassify mitomaster-build' to create the database.",
        )
    return MitoMasterDB(_db_path)


# Response models
class VariantFrequencyResponse(BaseModel):
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


class HaplogroupFrequencyResponse(BaseModel):
    haplogroup: str
    count: int
    frequency: float


class VariantDetailResponse(BaseModel):
    variant: VariantFrequencyResponse
    haplogroup_frequencies: list[HaplogroupFrequencyResponse]


class GeneInfo(BaseModel):
    name: str
    start: int
    end: int
    type: str
    strand: str
    product: str


class DatabaseStatsResponse(BaseModel):
    genome_count: int
    human_genome_count: int
    variant_count: int
    unique_variant_count: int
    top_haplogroups: list[dict]


class SequenceVariant(BaseModel):
    position: int
    ref: str
    alt: str
    mutation_type: str
    locus: Optional[str] = None
    gene_type: Optional[str] = None
    amino_acid_change: Optional[str] = None
    frequency: Optional[VariantFrequencyResponse] = None


class AnalyzeSequenceResponse(BaseModel):
    variant_count: int
    variants: list[SequenceVariant]


class SearchVariantsResponse(BaseModel):
    count: int
    variants: list[VariantFrequencyResponse]


# API endpoints
@router.get("/status")
async def mitomaster_status():
    """Check MitoMaster database status."""
    if _db_path is None:
        return {
            "status": "not_configured",
            "message": "Database path not set",
            "database_path": None,
        }
    if not _db_path.exists():
        return {
            "status": "missing",
            "message": "Database file does not exist",
            "database_path": str(_db_path),
        }
    try:
        with _get_db() as db:
            stats = db.get_stats()
        return {
            "status": "ok",
            "database_path": str(_db_path),
            "genome_count": stats["genome_count"],
            "variant_count": stats["variant_count"],
        }
    except Exception as e:
        return {
            "status": "error",
            "message": str(e),
            "database_path": str(_db_path),
        }


@router.get("/variant/{position}", response_model=VariantDetailResponse)
async def get_variant(
    position: int,
    ref: Optional[str] = Query(None, description="Reference base filter"),
    alt: Optional[str] = Query(None, description="Alternate base filter"),
):
    """Get variant frequency and annotations at a position.

    Args:
        position: mtDNA position (1-16569)
        ref: Optional reference base filter
        alt: Optional alternate base filter

    Returns:
        Variant frequency data with haplogroup breakdown
    """
    if position < 1 or position > 16569:
        raise HTTPException(
            status_code=400,
            detail="Position must be between 1 and 16569",
        )

    with _get_db() as db:
        variants = db.get_variant_at_position(position, ref, alt)

        if not variants:
            raise HTTPException(
                status_code=404,
                detail=f"No variants found at position {position}",
            )

        # Get the first (most common) variant
        variant = variants[0]

        # Get haplogroup frequencies
        hg_freqs = db.get_haplogroup_frequencies(position, variant.ref, variant.alt)

        return VariantDetailResponse(
            variant=VariantFrequencyResponse(
                position=variant.position,
                ref=variant.ref,
                alt=variant.alt,
                total_count=variant.total_count,
                total_frequency=variant.total_frequency,
                human_count=variant.human_count,
                human_frequency=variant.human_frequency,
                locus=variant.locus,
                gene_type=variant.gene_type,
                amino_acid_change=variant.amino_acid_change,
            ),
            haplogroup_frequencies=[
                HaplogroupFrequencyResponse(
                    haplogroup=hf.haplogroup,
                    count=hf.count,
                    frequency=hf.frequency,
                )
                for hf in hg_freqs
            ],
        )


@router.get("/variants", response_model=list[VariantFrequencyResponse])
async def list_variants_at_position(
    position: int,
    ref: Optional[str] = Query(None, description="Reference base filter"),
):
    """List all variant alleles at a position.

    Args:
        position: mtDNA position (1-16569)
        ref: Optional reference base filter

    Returns:
        List of all variants at the position
    """
    if position < 1 or position > 16569:
        raise HTTPException(
            status_code=400,
            detail="Position must be between 1 and 16569",
        )

    with _get_db() as db:
        variants = db.get_variant_at_position(position, ref, None)

        return [
            VariantFrequencyResponse(
                position=v.position,
                ref=v.ref,
                alt=v.alt,
                total_count=v.total_count,
                total_frequency=v.total_frequency,
                human_count=v.human_count,
                human_frequency=v.human_frequency,
                locus=v.locus,
                gene_type=v.gene_type,
                amino_acid_change=v.amino_acid_change,
            )
            for v in variants
        ]


@router.get("/search", response_model=SearchVariantsResponse)
async def search_variants(
    position_start: Optional[int] = Query(None, description="Start position"),
    position_end: Optional[int] = Query(None, description="End position"),
    gene: Optional[str] = Query(None, description="Gene name filter"),
    haplogroup: Optional[str] = Query(None, description="Haplogroup filter"),
    min_frequency: Optional[float] = Query(None, description="Minimum frequency"),
    limit: int = Query(100, description="Maximum results", le=1000),
):
    """Search variants with filters.

    Args:
        position_start: Start position filter
        position_end: End position filter
        gene: Gene name filter (e.g., "ND1", "CYTB")
        haplogroup: Haplogroup filter (not yet implemented)
        min_frequency: Minimum frequency filter (0-1)
        limit: Maximum results (default 100, max 1000)

    Returns:
        List of matching variants
    """
    with _get_db() as db:
        variants = db.search_variants(
            position_start=position_start,
            position_end=position_end,
            gene=gene,
            haplogroup=haplogroup,
            min_frequency=min_frequency,
            limit=limit,
        )

        return SearchVariantsResponse(
            count=len(variants),
            variants=[
                VariantFrequencyResponse(
                    position=v.position,
                    ref=v.ref,
                    alt=v.alt,
                    total_count=v.total_count,
                    total_frequency=v.total_frequency,
                    human_count=v.human_count,
                    human_frequency=v.human_frequency,
                    locus=v.locus,
                    gene_type=v.gene_type,
                    amino_acid_change=v.amino_acid_change,
                )
                for v in variants
            ],
        )


@router.get("/genes", response_model=list[GeneInfo])
async def list_genes():
    """List all mtDNA genes.

    Returns:
        List of gene information
    """
    with _get_db() as db:
        genes = db.get_all_genes()
        return [
            GeneInfo(
                name=g["name"],
                start=g["start"],
                end=g["end"],
                type=g["type"],
                strand=g["strand"],
                product=g["product"],
            )
            for g in genes
        ]


@router.get("/stats", response_model=DatabaseStatsResponse)
async def get_stats():
    """Get database statistics.

    Returns:
        Database statistics including genome and variant counts
    """
    with _get_db() as db:
        stats = db.get_stats()
        return DatabaseStatsResponse(
            genome_count=stats["genome_count"],
            human_genome_count=stats["human_genome_count"],
            variant_count=stats["variant_count"],
            unique_variant_count=stats["unique_variant_count"],
            top_haplogroups=stats["top_haplogroups"],
        )


@router.post("/analyze", response_model=AnalyzeSequenceResponse)
async def analyze_sequence(
    sequence: str,
    include_frequencies: bool = Query(True, description="Include frequency data"),
):
    """Analyze a sequence and return all variants vs rCRS.

    Args:
        sequence: mtDNA sequence to analyze
        include_frequencies: Whether to include frequency data for each variant

    Returns:
        List of variants found in the sequence
    """
    from Bio.Seq import Seq

    # Validate sequence
    sequence = sequence.upper().replace(" ", "").replace("\n", "")
    if not all(c in "ACGTN-" for c in sequence):
        raise HTTPException(
            status_code=400,
            detail="Invalid sequence characters. Only A, C, G, T, N, - are allowed.",
        )

    if len(sequence) < 100:
        raise HTTPException(
            status_code=400,
            detail="Sequence too short. Minimum 100 bases required.",
        )

    with _get_db() as db:
        # Get rCRS reference from processor
        from mtclassify.mitomaster.processor import GenomeProcessor

        # Use a temporary processor just to get rCRS
        try:
            processor = GenomeProcessor(db)
            rcrs = processor.rcrs
        except Exception as e:
            raise HTTPException(
                status_code=500,
                detail=f"Failed to load reference sequence: {e}",
            )

        # Compare sequence to rCRS
        variants = []
        min_len = min(len(sequence), len(rcrs))

        for i in range(min_len):
            ref = rcrs[i]
            alt = sequence[i]
            position = i + 1

            if ref != alt and alt in "ACGT":
                # Determine mutation type
                transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
                if (ref, alt) in transitions:
                    mutation_type = "transition"
                else:
                    mutation_type = "transversion"

                # Get gene info
                gene_info = db.get_gene_for_position(position)
                locus = gene_info[0] if gene_info else None
                gene_type = gene_info[1] if gene_info else None

                variant = SequenceVariant(
                    position=position,
                    ref=ref,
                    alt=alt,
                    mutation_type=mutation_type,
                    locus=locus,
                    gene_type=gene_type,
                )

                # Get frequency data if requested
                if include_frequencies:
                    freq_list = db.get_variant_at_position(position, ref, alt)
                    if freq_list:
                        freq = freq_list[0]
                        variant.frequency = VariantFrequencyResponse(
                            position=freq.position,
                            ref=freq.ref,
                            alt=freq.alt,
                            total_count=freq.total_count,
                            total_frequency=freq.total_frequency,
                            human_count=freq.human_count,
                            human_frequency=freq.human_frequency,
                            locus=freq.locus,
                            gene_type=freq.gene_type,
                            amino_acid_change=freq.amino_acid_change,
                        )

                variants.append(variant)

        return AnalyzeSequenceResponse(
            variant_count=len(variants),
            variants=variants,
        )


@router.get("/gene/{gene_name}")
async def get_gene_variants(
    gene_name: str,
    min_frequency: Optional[float] = Query(None, description="Minimum frequency"),
    limit: int = Query(100, description="Maximum results", le=1000),
):
    """Get all variants within a specific gene.

    Args:
        gene_name: Gene name (e.g., "ND1", "CYTB", "D-loop-1")
        min_frequency: Minimum frequency filter
        limit: Maximum results

    Returns:
        List of variants in the gene
    """
    with _get_db() as db:
        variants = db.search_variants(
            gene=gene_name,
            min_frequency=min_frequency,
            limit=limit,
        )

        if not variants:
            raise HTTPException(
                status_code=404,
                detail=f"No variants found for gene '{gene_name}'. "
                       "Use /api/mitomaster/genes to list valid gene names.",
            )

        return {
            "gene": gene_name,
            "count": len(variants),
            "variants": [
                VariantFrequencyResponse(
                    position=v.position,
                    ref=v.ref,
                    alt=v.alt,
                    total_count=v.total_count,
                    total_frequency=v.total_frequency,
                    human_count=v.human_count,
                    human_frequency=v.human_frequency,
                    locus=v.locus,
                    gene_type=v.gene_type,
                    amino_acid_change=v.amino_acid_change,
                )
                for v in variants
            ],
        }
