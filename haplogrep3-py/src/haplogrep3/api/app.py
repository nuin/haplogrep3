"""FastAPI application for haplogrep3 web service."""

import logging
import os
import tempfile
import uuid
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from haplogrep3 import __version__
from haplogrep3.distance import Distance
from haplogrep3.io import PhylotreeLoader, VcfReader, FastaReader, TsvReader
from haplogrep3.tasks import ClassificationTask
from haplogrep3.api.mitomaster_routes import router as mitomaster_router, set_database_path

logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Lifespan context manager for startup/shutdown events."""
    # Startup: preload all trees into cache
    logger.info("Preloading phylogenetic trees...")
    loader = PhylotreeLoader()
    available = loader.list_available()
    for tree_id in available:
        try:
            loader.load(tree_id)
            logger.info(f"  Loaded: {tree_id}")
        except Exception as e:
            logger.warning(f"  Failed to load {tree_id}: {e}")
    logger.info(f"Preloaded {len(available)} trees")

    # Configure MitoMaster database path
    mitomaster_db = os.environ.get("MITOMASTER_DB")
    if mitomaster_db:
        db_path = Path(mitomaster_db)
        if db_path.exists():
            set_database_path(db_path)
            logger.info(f"MitoMaster database configured: {db_path}")
        else:
            logger.warning(f"MitoMaster database not found: {db_path}")
    else:
        # Try default location
        default_db = Path.home() / ".haplogrep3" / "mitomaster.db"
        if default_db.exists():
            set_database_path(default_db)
            logger.info(f"MitoMaster database found at default location: {default_db}")
        else:
            logger.info("MitoMaster database not configured (set MITOMASTER_DB env var)")

    yield

    # Shutdown: nothing to clean up


def create_app() -> FastAPI:
    """Create and configure the FastAPI application."""
    app = FastAPI(
        title="Haplogrep3 API",
        description="mtDNA haplogroup classification service",
        version=__version__,
        lifespan=lifespan,
    )

    # CORS middleware for frontend
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],  # Configure appropriately for production
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # Include MitoMaster routes
    app.include_router(mitomaster_router)

    return app


app = create_app()


# Response models
class PolymorphismResponse(BaseModel):
    position: int
    mutation: str
    insertion_position: Optional[int] = None


class ClassificationResultResponse(BaseModel):
    haplogroup: str
    quality: float
    found_count: int
    missing_count: int
    remaining_count: int
    found_polymorphisms: list[str]
    missing_polymorphisms: list[str]
    remaining_polymorphisms: list[str]


class SampleResultResponse(BaseModel):
    id: str
    haplogroup: str
    quality: float
    range: str
    polymorphisms: list[str]
    top_result: ClassificationResultResponse
    other_results: list[ClassificationResultResponse]


class ClassifyResponse(BaseModel):
    job_id: str
    tree: str
    distance: str
    samples_count: int
    results: list[SampleResultResponse]


class TreeInfo(BaseModel):
    id: str
    name: str
    version: str


class TreesResponse(BaseModel):
    trees: list[TreeInfo]


class HealthResponse(BaseModel):
    status: str
    version: str


# API endpoints
@app.get("/api/health", response_model=HealthResponse)
async def health_check():
    """Health check endpoint."""
    return HealthResponse(status="ok", version=__version__)


@app.get("/api/trees", response_model=TreesResponse)
async def list_trees():
    """List available phylogenetic trees."""
    loader = PhylotreeLoader()
    available = loader.list_available()

    trees = []
    for tree_id in available:
        try:
            tree = loader.load(tree_id)
            trees.append(TreeInfo(
                id=tree.id,
                name=tree.name,
                version=tree.version,
            ))
        except Exception:
            trees.append(TreeInfo(id=tree_id, name=tree_id, version=""))

    return TreesResponse(trees=trees)


@app.post("/api/classify", response_model=ClassifyResponse)
async def classify(
    file: UploadFile = File(..., description="VCF or FASTA file to classify"),
    tree: str = Form(..., description="Tree ID to use for classification"),
    distance: str = Form("kulczynski", description="Distance metric"),
    hits: int = Form(1, description="Number of top hits to return"),
    het_level: float = Form(0.9, description="Heteroplasmy level threshold"),
    chip: bool = Form(False, description="VCF from genotyping chip"),
):
    """Classify samples from uploaded file."""
    job_id = str(uuid.uuid4())[:8]

    # Validate distance metric
    try:
        distance_enum = Distance(distance.lower())
    except ValueError:
        raise HTTPException(
            status_code=400,
            detail=f"Invalid distance metric: {distance}. Valid options: kulczynski, hamming, jaccard, kimura"
        )

    # Load phylotree
    loader = PhylotreeLoader()
    try:
        phylotree = loader.load(tree)
    except FileNotFoundError:
        raise HTTPException(
            status_code=404,
            detail=f"Tree '{tree}' not found. Use /api/trees to list available trees."
        )

    # Save uploaded file to temp location
    suffix = Path(file.filename).suffix if file.filename else ".vcf"
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        content = await file.read()
        tmp.write(content)
        tmp_path = Path(tmp.name)

    try:
        # Determine file type and read samples
        if suffix.lower() in (".vcf", ".vcf.gz"):
            reader = VcfReader(het_level=het_level, chip=chip)
            samples = reader.read(tmp_path)
        elif suffix.lower() in (".fasta", ".fa", ".fasta.gz", ".fa.gz"):
            if not phylotree.reference_fasta:
                raise HTTPException(
                    status_code=400,
                    detail="Tree does not have a reference FASTA for FASTA input"
                )
            reader = FastaReader(reference_path=Path(phylotree.reference_fasta))
            samples = reader.read(tmp_path)
        elif suffix.lower() in (".tsv", ".txt", ".csv"):
            reader = TsvReader(het_threshold=het_level)
            samples = reader.read(tmp_path)
        else:
            raise HTTPException(
                status_code=400,
                detail=f"Unsupported file format: {suffix}. Use .vcf, .fasta, .tsv, or .txt"
            )

        # Run classification
        task = ClassificationTask(phylotree, distance_enum, hits)
        samples = task.classify(samples)

        # Build response
        results = []
        for sample in samples:
            if sample.classification:
                top = sample.classification.top_result
                top_result = ClassificationResultResponse(
                    haplogroup=top.haplogroup.name,
                    quality=top.quality,
                    found_count=len(top.found_polymorphisms),
                    missing_count=len(top.missing_polymorphisms),
                    remaining_count=len(top.remaining_polymorphisms),
                    found_polymorphisms=[str(p) for p in top.found_polymorphisms],
                    missing_polymorphisms=[str(p) for p in top.missing_polymorphisms],
                    remaining_polymorphisms=[str(p) for p in top.remaining_polymorphisms],
                )

                other_results = []
                for other in sample.classification.other_results:
                    other_results.append(ClassificationResultResponse(
                        haplogroup=other.haplogroup.name,
                        quality=other.quality,
                        found_count=len(other.found_polymorphisms),
                        missing_count=len(other.missing_polymorphisms),
                        remaining_count=len(other.remaining_polymorphisms),
                        found_polymorphisms=[str(p) for p in other.found_polymorphisms],
                        missing_polymorphisms=[str(p) for p in other.missing_polymorphisms],
                        remaining_polymorphisms=[str(p) for p in other.remaining_polymorphisms],
                    ))

                results.append(SampleResultResponse(
                    id=sample.id,
                    haplogroup=top.haplogroup.name,
                    quality=top.quality,
                    range=sample.range_str,
                    polymorphisms=sample.get_polymorphism_strings(),
                    top_result=top_result,
                    other_results=other_results,
                ))
            else:
                results.append(SampleResultResponse(
                    id=sample.id,
                    haplogroup="?",
                    quality=0.0,
                    range=sample.range_str,
                    polymorphisms=sample.get_polymorphism_strings(),
                    top_result=ClassificationResultResponse(
                        haplogroup="?",
                        quality=0.0,
                        found_count=0,
                        missing_count=0,
                        remaining_count=len(sample.polymorphisms),
                        found_polymorphisms=[],
                        missing_polymorphisms=[],
                        remaining_polymorphisms=sample.get_polymorphism_strings(),
                    ),
                    other_results=[],
                ))

        return ClassifyResponse(
            job_id=job_id,
            tree=tree,
            distance=distance_enum.value,
            samples_count=len(results),
            results=results,
        )

    finally:
        # Clean up temp file
        tmp_path.unlink(missing_ok=True)


@app.get("/api/distances")
async def list_distances():
    """List available distance metrics."""
    return {
        "distances": [
            {"id": "kulczynski", "name": "Kulczynski", "default": True},
            {"id": "hamming", "name": "Hamming", "default": False},
            {"id": "jaccard", "name": "Jaccard", "default": False},
            {"id": "kimura", "name": "Kimura", "default": False},
        ]
    }
