# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

mtclassify is a Python mtDNA haplogroup classification tool. It classifies mitochondrial DNA samples against phylogenetic trees (PhyloTree) to determine haplogroups.

## Build and Development Commands

```bash
# Install dependencies
uv sync

# Run tests
uv run pytest tests/ -v

# Run a single test file
uv run pytest tests/test_polymorphism.py -v

# Run with coverage
uv run pytest tests/ --cov=mtclassify

# Run the CLI
uv run mtclassify --help
uv run mtclassify run --help
uv run mtclassify trees
```

## Project Structure

```
src/mtclassify/
├── __init__.py
├── cli.py              # Typer CLI entry point
├── models/             # Pydantic data models
│   ├── polymorphism.py # Polymorphism (mtDNA variant)
│   ├── haplogroup.py   # Haplogroup class
│   ├── sample.py       # Sample, RankedResult, ClassificationResult
│   └── phylotree.py    # PhyloTreeNode, Phylotree
├── distance/           # Distance metrics
│   ├── base.py         # DistanceMetric ABC
│   ├── kulczynski.py   # Default metric
│   ├── hamming.py
│   ├── jaccard.py
│   └── kimura.py
├── io/                 # Input/Output
│   ├── tree_loader.py  # Phylotree YAML loader
│   ├── vcf_reader.py   # VCF parser (cyvcf2)
│   └── fasta_reader.py # FASTA parser (biopython)
├── mitomaster/         # MitoMaster variant database
│   ├── database.py     # DuckDB database
│   ├── downloader.py   # NCBI downloader
│   ├── processor.py    # GenBank processor
│   └── frequency.py    # Frequency calculator
└── tasks/              # Classification pipeline
    ├── classify.py     # ClassificationTask
    └── export.py       # CSV/FASTA export
```

## Architecture

### Classification Pipeline
1. **Input** → VcfReader or FastaReader parses input file to Sample objects
2. **Tree Loading** → PhylotreeLoader loads phylogenetic tree from YAML
3. **Classification** → ClassificationTask scores all haplogroups using distance metric
4. **Export** → Results written to CSV, optionally FASTA/QC reports

### Key Classes
- `Polymorphism`: mtDNA variant (position + mutation)
- `Sample`: Sample with polymorphisms and classification results
- `Phylotree`: Tree structure with PhyloTreeNodes
- `ClassificationTask`: Orchestrates classification using distance metrics
- `DistanceMetric`: ABC for Kulczynski, Hamming, Jaccard, Kimura

### Dependencies
- `typer`: CLI framework
- `pydantic`: Data validation
- `cyvcf2`: VCF parsing
- `biopython`: FASTA parsing
- `rich`: CLI output formatting
- `fastapi`: Web API framework
- `uvicorn`: ASGI server
- `duckdb`: MitoMaster database

## Web Interface

### Running the Server

```bash
# Start the FastAPI server
uv run mtclassify server

# With custom host/port
uv run mtclassify server --host 0.0.0.0 --port 8000

# Development mode with auto-reload
uv run mtclassify server --reload
```

### API Endpoints
- `GET /api/health` - Health check
- `GET /api/trees` - List available phylogenetic trees
- `GET /api/distances` - List distance metrics
- `POST /api/classify` - Classify samples from uploaded file

### MitoMaster Endpoints
- `GET /api/mitomaster/status` - Database status
- `GET /api/mitomaster/stats` - Database statistics
- `GET /api/mitomaster/variant/{position}` - Get variant at position
- `GET /api/mitomaster/search` - Search variants

### Frontend Development

```bash
cd frontend

# Install dependencies
npm install

# Run development server (connects to backend at localhost:7001)
npm run dev

# Type check
npm run check

# Build for production
npm run build
```

### Project Structure (Web)
```
src/mtclassify/api/
├── __init__.py
├── app.py              # FastAPI application
└── mitomaster_routes.py # MitoMaster API routes

frontend/               # Svelte 5 + SvelteKit
├── src/
│   ├── lib/
│   │   └── api.ts      # API client
│   └── routes/
│       ├── +page.svelte      # Main classification UI
│       └── mitomaster/       # MitoMaster UI
│           └── +page.svelte
├── svelte.config.js
└── vite.config.ts      # Includes proxy to backend
```

## CLI Commands

### Classification
```bash
# Classify mtDNA samples
mtclassify run -i samples.vcf -o results.csv

# With options
mtclassify run -i samples.fasta -o results.csv -t phylotree-rcrs@17.2 -d kulczynski --hits 3
```

### MitoMaster Database
```bash
# Build database from NCBI
mtclassify mitomaster-build -e your@email.com

# Import local GenBank files
mtclassify mitomaster-import /path/to/genomes/ -o mitomaster.db

# Query a position
mtclassify mitomaster-query 3243

# Show database stats
mtclassify mitomaster-stats
```
