# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Haplogrep 3 is a free mtDNA haplogroup classification service. It can run as a web service or command-line tool, classifying mitochondrial DNA samples against phylogenetic trees (PhyloTree) to determine haplogroups.

## Build and Test Commands

```bash
# Build the project (creates haplogrep3.jar and platform-specific distributions)
mvn package

# Run tests
mvn test

# Run a single test class
mvn test -Dtest=ClassifyCommandTest

# Run a single test method
mvn test -Dtest=ClassifyCommandTest#testWithHsd
```

## Running Haplogrep

```bash
# Start local web server (default port 7000)
java -jar target/haplogrep3.jar server

# Classify samples via CLI
java -jar target/haplogrep3.jar classify --tree <tree-id> --in <input-file> --out <output-file>

# List available phylogenetic trees
java -jar target/haplogrep3.jar trees
```

## Architecture

### Entry Point and CLI
- `App.java` - Main entry point using picocli for CLI parsing
- `commands/` - CLI subcommands (classify, server, align, build-tree, install-tree, etc.)

### Core Classification Pipeline
1. **Input Readers** (`haplogrep/io/readers/`) - Parse VCF, FASTA, or HSD (text-based) input formats
2. **ClassificationTask** (`tasks/ClassificationTask.java`) - Orchestrates classification using phylotree and distance metrics
3. **Phylotree** (`model/Phylotree.java`) - Loads tree definitions from YAML, manages haplogroup lookups and classification via `haplogrep-core` library
4. **AnnotationTask** (`tasks/AnnotationTask.java`) - Enriches samples with annotation data from configured sources
5. **Export Tasks** (`tasks/Export*.java`) - Generate output in various formats (CSV, FASTA, QC reports, HTML)

### Web Application
- `web/WebApp.java` - Javalin-based web server with route definitions
- `web/handlers/` - Request handlers for jobs, phylogenies, mutations, clades
- Uses BasisTemplate for HTML rendering

### Key Model Classes
- `Phylotree` - Tree configuration and classification logic
- `PhylotreeRepository` - Manages installed phylogenetic trees
- `AnnotatedSample` - Sample with classification results and annotations
- `Distance` - Enum for distance metrics (Kulczynski, Hamming, Jaccard, Kimura)

### Configuration
- `haplogrep3.yaml` - Main configuration (port, examples, tree repositories, enabled phylotrees)
- Tree definitions loaded from remote repositories or local YAML files

### Dependencies
- `haplogrep-core` (genepi) - Core classification algorithms and phylotree parsing
- `genepi-annotate` - Annotation framework
- Javalin - Web framework
- picocli - CLI framework
- YamlBeans - YAML parsing

## Test Data

Test files in `test-data/` organized by format:
- `hsd/` - HSD format samples
- `fasta/` - FASTA format samples
- `vcf/` - VCF format samples
- `expected/` - Expected output files for test validation
