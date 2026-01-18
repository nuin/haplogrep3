"""Command-line interface for mtclassify."""

from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.table import Table

from mtclassify import __version__
from mtclassify.distance import Distance
from mtclassify.io import PhylotreeLoader, VcfReader, FastaReader, TsvReader
from mtclassify.tasks import ClassificationTask, export_csv, export_fasta
from mtclassify.tasks.export import export_qc_report

app = typer.Typer(
    name="mtclassify",
    help="mtDNA haplogroup classification tool",
    no_args_is_help=True,
)
console = Console()


def version_callback(value: bool):
    if value:
        console.print(f"mtclassify {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        Optional[bool],
        typer.Option("--version", "-v", callback=version_callback, is_eager=True),
    ] = None,
):
    """mtclassify - mtDNA haplogroup classification tool."""
    pass


@app.command()
def run(
    input_file: Annotated[
        Path,
        typer.Option("--input", "--in", "-i", help="Input file (VCF, FASTA, TSV, or TXT)"),
    ],
    output: Annotated[
        Path,
        typer.Option("--output", "--out", "-o", help="Output file path"),
    ],
    tree: Annotated[
        str,
        typer.Option("--tree", "-t", help="Tree ID or path to tree file"),
    ] = "phylotree-rcrs@17.2",
    distance: Annotated[
        Distance,
        typer.Option("--distance", "--metric", "-d", help="Distance metric"),
    ] = Distance.KULCZYNSKI,
    hits: Annotated[
        int,
        typer.Option("--hits", "-n", help="Number of top hits to return"),
    ] = 1,
    extend_report: Annotated[
        bool,
        typer.Option("--extend-report", help="Include extended polymorphism details"),
    ] = False,
    write_fasta: Annotated[
        bool,
        typer.Option("--write-fasta", help="Write results in FASTA format"),
    ] = False,
    write_qc: Annotated[
        bool,
        typer.Option("--write-qc", help="Write quality control report"),
    ] = False,
    het_level: Annotated[
        float,
        typer.Option("--het-level", help="Heteroplasmy level threshold"),
    ] = 0.9,
    chip: Annotated[
        bool,
        typer.Option("--chip", help="VCF data from genotyping chip"),
    ] = False,
    skip_alignment_rules: Annotated[
        bool,
        typer.Option("--skip-alignment-rules", help="Skip nomenclature fixes for FASTA"),
    ] = False,
    input_format: Annotated[
        Optional[str],
        typer.Option("--format", "-f", help="Input format: vcf, fasta, tsv, seqnext (auto-detect if not specified)"),
    ] = None,
    skip_hotspots: Annotated[
        bool,
        typer.Option("--skip-hotspots/--include-hotspots", help="Skip hotspot variants in TSV input"),
    ] = True,
):
    """Run mtDNA haplogroup classification on samples."""
    # Validate input file
    if not input_file.exists():
        console.print(f"[red]Error:[/red] Input file '{input_file}' not found.")
        raise typer.Exit(1)

    # Load phylotree
    console.print(f"Loading tree: {tree}")
    try:
        loader = PhylotreeLoader()
        phylotree = loader.load(tree)
    except FileNotFoundError as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)

    # Determine input format
    if input_format:
        fmt = input_format.lower()
    else:
        suffix = input_file.suffix.lower()
        if suffix == ".gz":
            stem_suffix = Path(input_file.stem).suffix.lower()
            suffix = stem_suffix + suffix

        if suffix in (".vcf", ".vcf.gz"):
            fmt = "vcf"
        elif suffix in (".fasta", ".fa", ".fasta.gz", ".fa.gz"):
            fmt = "fasta"
        elif suffix in (".tsv", ".txt", ".csv"):
            fmt = "tsv"
        else:
            # Try to detect from content
            with open(input_file) as f:
                first_line = f.readline()
            if first_line.startswith("##fileformat=VCF"):
                fmt = "vcf"
            elif first_line.startswith(">"):
                fmt = "fasta"
            elif "HGVS" in first_line or "Nuc Change" in first_line:
                fmt = "tsv"
            else:
                console.print(f"[red]Error:[/red] Could not detect file format. Use --format to specify.")
                raise typer.Exit(1)

    console.print(f"Reading input file: {input_file} (format: {fmt})")

    if fmt == "vcf":
        reader = VcfReader(het_level=het_level, chip=chip)
        samples = reader.read(input_file)
    elif fmt == "fasta":
        if phylotree.reference_fasta:
            ref_path = Path(phylotree.reference_fasta)
        else:
            console.print("[red]Error:[/red] Tree does not have a reference FASTA.")
            raise typer.Exit(1)

        reader = FastaReader(
            reference_path=ref_path,
            skip_alignment_rules=skip_alignment_rules,
        )
        samples = reader.read(input_file)
    elif fmt in ("tsv", "seqnext", "txt"):
        reader = TsvReader(
            het_threshold=het_level,
            skip_hotspots=skip_hotspots,
        )
        samples = reader.read(input_file)
    else:
        console.print(f"[red]Error:[/red] Unsupported format: {fmt}")
        raise typer.Exit(1)

    console.print(f"Loaded {len(samples)} samples")

    # Run classification
    console.print(f"Classifying with {distance.value} metric...")
    task = ClassificationTask(phylotree, distance, hits)
    samples = task.classify(samples)

    # Export results
    console.print(f"Writing results to: {output}")
    export_csv(samples, output, extended=extend_report, hits=hits)

    # Optional FASTA output
    if write_fasta and phylotree.reference_fasta:
        fasta_output = output.with_suffix(".fasta")
        with open(phylotree.reference_fasta) as f:
            from Bio import SeqIO
            ref_seq = str(next(SeqIO.parse(f, "fasta")).seq)
        export_fasta(samples, fasta_output, ref_seq)
        console.print(f"FASTA written to: {fasta_output}")

    # Optional QC report
    if write_qc:
        qc_output = output.with_suffix(".qc.txt")
        export_qc_report(samples, qc_output)
        console.print(f"QC report written to: {qc_output}")

    # Summary
    console.print(f"\n[green]Classification complete![/green]")
    console.print(f"  Samples processed: {len(samples)}")

    if samples and samples[0].classification:
        # Show top haplogroups
        haplogroups = {}
        for s in samples:
            if s.classification:
                hg = s.classification.haplogroup.name
                haplogroups[hg] = haplogroups.get(hg, 0) + 1

        console.print("\nTop haplogroups:")
        for hg, count in sorted(haplogroups.items(), key=lambda x: -x[1])[:5]:
            console.print(f"  {hg}: {count}")


@app.command()
def trees(
    trees_dir: Annotated[
        Optional[Path],
        typer.Option("--trees-dir", help="Directory containing tree files"),
    ] = None,
):
    """List available phylogenetic trees."""
    loader = PhylotreeLoader(trees_dir)
    available = loader.list_available()

    if not available:
        console.print("No trees installed.")
        console.print("\nInstall trees with:")
        console.print("  mtclassify install-tree <tree-id>")
        console.print("\nAvailable trees from repository:")
        console.print("  - phylotree-fu-rcrs@1.2")
        console.print("  - phylotree-fu-rcrs@1.0")
        console.print("  - phylotree-rcrs@17.2")
        return

    table = Table(title="Available Trees")
    table.add_column("Tree ID", style="cyan")

    for tree_id in available:
        table.add_row(tree_id)

    console.print(table)


@app.command()
def server(
    host: Annotated[
        str,
        typer.Option("--host", "-h", help="Host to bind to"),
    ] = "127.0.0.1",
    port: Annotated[
        int,
        typer.Option("--port", "-p", help="Port to bind to"),
    ] = 7001,
    reload: Annotated[
        bool,
        typer.Option("--reload", help="Enable auto-reload for development"),
    ] = False,
):
    """Start the web server."""
    import uvicorn

    console.print(f"Starting mtclassify server at http://{host}:{port}")
    console.print("API docs available at /docs")
    console.print("Press Ctrl+C to stop\n")

    uvicorn.run(
        "mtclassify.api:app",
        host=host,
        port=port,
        reload=reload,
    )


@app.command("install-tree")
def install_tree(
    tree_id: Annotated[
        str,
        typer.Argument(help="Tree ID to install (e.g., phylotree-fu-rcrs@1.2)"),
    ],
    trees_dir: Annotated[
        Optional[Path],
        typer.Option("--trees-dir", help="Directory to install trees"),
    ] = None,
):
    """Install a phylogenetic tree from the repository."""
    console.print(f"Installing tree: {tree_id}")

    # TODO: Implement tree download from repository
    console.print("[yellow]Note:[/yellow] Tree installation from remote repository not yet implemented.")
    console.print("\nFor now, download trees manually from:")
    console.print("  https://github.com/genepi/haplogrep-trees")
    console.print(f"\nAnd extract to: {PhylotreeLoader(trees_dir).trees_dir}")


# =============================================================================
# MitoMaster Commands
# =============================================================================

@app.command("mitomaster-build")
def mitomaster_build(
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Output database file path"),
    ] = Path.home() / ".mtclassify" / "mitomaster.db",
    email: Annotated[
        str,
        typer.Option("--email", "-e", help="NCBI Entrez email (required)"),
    ] = "",
    api_key: Annotated[
        Optional[str],
        typer.Option("--api-key", help="NCBI API key for faster downloads"),
    ] = None,
    organism: Annotated[
        Optional[str],
        typer.Option("--organism", help="Filter by organism (e.g., 'Homo sapiens')"),
    ] = None,
    max_genomes: Annotated[
        Optional[int],
        typer.Option("--max", "-m", help="Maximum number of genomes to process"),
    ] = None,
    tree: Annotated[
        str,
        typer.Option("--tree", "-t", help="Tree ID for haplogroup classification"),
    ] = "phylotree-rcrs@17.2",
    batch_size: Annotated[
        int,
        typer.Option("--batch-size", help="Genomes per batch for frequency updates"),
    ] = 1000,
):
    """Build MitoMaster database from NCBI mtDNA genomes.

    Downloads complete mitochondrial genomes from NCBI, classifies them with
    mtclassify, extracts variants vs rCRS, and stores everything in a portable
    DuckDB database.

    Example:
        mtclassify mitomaster-build -e your@email.com --max 1000
    """
    if not email:
        console.print("[red]Error:[/red] NCBI requires an email address. Use --email")
        raise typer.Exit(1)

    from rich.progress import Progress, TaskID

    from mtclassify.mitomaster import (
        MitoMasterDB,
        NCBIDownloader,
        DownloadConfig,
        GenomeProcessor,
        FrequencyCalculator,
    )

    # Create output directory
    output.parent.mkdir(parents=True, exist_ok=True)

    console.print(f"[bold]MitoMaster Database Builder[/bold]")
    console.print(f"Output: {output}")
    console.print(f"Tree: {tree}")
    if organism:
        console.print(f"Organism filter: {organism}")
    if max_genomes:
        console.print(f"Max genomes: {max_genomes}")
    console.print()

    # Initialize database
    console.print("Initializing database...")
    with MitoMasterDB(output) as db:
        db.initialize()

        # Setup downloader
        config = DownloadConfig(
            email=email,
            api_key=api_key,
            batch_size=100,
        )
        downloader = NCBIDownloader(config)

        # Search for genomes
        console.print("\nSearching NCBI for mitochondrial genomes...")
        accessions = downloader.search_genomes(organism=organism, max_results=max_genomes)
        console.print(f"Found {len(accessions)} genomes")

        # Filter existing
        new_accessions = [acc for acc in accessions if not db.genome_exists(acc)]
        console.print(f"New genomes to process: {len(new_accessions)}")

        if not new_accessions:
            console.print("[green]Database is up to date![/green]")
            return

        # Initialize processor
        processor = GenomeProcessor(db, tree_id=tree)

        # Process genomes with progress bar
        processed = 0
        failed = 0

        with Progress() as progress:
            task = progress.add_task("Processing genomes...", total=len(new_accessions))

            for accession, temp_path in downloader.stream_genomes(
                new_accessions,
                skip_existing=db.genome_exists,
            ):
                try:
                    success = processor.process_and_store(temp_path, delete_after=True)
                    if success:
                        processed += 1
                    else:
                        failed += 1
                except Exception as e:
                    console.print(f"[red]Error processing {accession}:[/red] {e}")
                    failed += 1
                    if temp_path.exists():
                        temp_path.unlink()

                progress.update(task, advance=1)

                # Periodic frequency updates
                if processed > 0 and processed % batch_size == 0:
                    console.print(f"\nUpdating frequencies (processed {processed})...")
                    calculator = FrequencyCalculator(db)
                    calculator.compute_all_frequencies()

        console.print(f"\n[green]Processing complete![/green]")
        console.print(f"  Processed: {processed}")
        console.print(f"  Failed: {failed}")

        # Final frequency calculation
        console.print("\nComputing final frequency statistics...")
        calculator = FrequencyCalculator(db)
        calculator.compute_all_frequencies()

        # Show stats
        stats = db.get_stats()
        console.print(f"\n[bold]Database Statistics:[/bold]")
        console.print(f"  Total genomes: {stats['genome_count']}")
        console.print(f"  Human genomes: {stats['human_genome_count']}")
        console.print(f"  Total variants: {stats['variant_count']}")
        console.print(f"  Unique variants: {stats['unique_variant_count']}")


@app.command("mitomaster-update")
def mitomaster_update(
    database: Annotated[
        Path,
        typer.Option("--database", "-d", help="Path to existing MitoMaster database"),
    ] = Path.home() / ".mtclassify" / "mitomaster.db",
    email: Annotated[
        str,
        typer.Option("--email", "-e", help="NCBI Entrez email (required)"),
    ] = "",
    api_key: Annotated[
        Optional[str],
        typer.Option("--api-key", help="NCBI API key for faster downloads"),
    ] = None,
    organism: Annotated[
        Optional[str],
        typer.Option("--organism", help="Filter by organism"),
    ] = None,
    tree: Annotated[
        str,
        typer.Option("--tree", "-t", help="Tree ID for haplogroup classification"),
    ] = "phylotree-rcrs@17.2",
):
    """Update existing MitoMaster database with new NCBI genomes.

    Checks NCBI for new genomes not already in the database and adds them.

    Example:
        mtclassify mitomaster-update -e your@email.com
    """
    if not database.exists():
        console.print(f"[red]Error:[/red] Database not found: {database}")
        console.print("Use 'mtclassify mitomaster-build' to create a new database.")
        raise typer.Exit(1)

    if not email:
        console.print("[red]Error:[/red] NCBI requires an email address. Use --email")
        raise typer.Exit(1)

    from rich.progress import Progress

    from mtclassify.mitomaster import (
        MitoMasterDB,
        NCBIDownloader,
        DownloadConfig,
        GenomeProcessor,
        FrequencyCalculator,
    )

    console.print(f"[bold]MitoMaster Database Update[/bold]")
    console.print(f"Database: {database}")

    with MitoMasterDB(database) as db:
        # Show current stats
        stats = db.get_stats()
        console.print(f"Current genomes: {stats['genome_count']}")

        # Setup downloader
        config = DownloadConfig(email=email, api_key=api_key)
        downloader = NCBIDownloader(config)

        # Search for genomes
        console.print("\nSearching NCBI for new genomes...")
        accessions = downloader.search_genomes(organism=organism)

        # Filter to new only
        new_accessions = [acc for acc in accessions if not db.genome_exists(acc)]
        console.print(f"New genomes found: {len(new_accessions)}")

        if not new_accessions:
            console.print("[green]Database is up to date![/green]")
            return

        # Process new genomes
        processor = GenomeProcessor(db, tree_id=tree)
        processed = 0

        with Progress() as progress:
            task = progress.add_task("Processing new genomes...", total=len(new_accessions))

            for accession, temp_path in downloader.stream_genomes(new_accessions):
                try:
                    if processor.process_and_store(temp_path, delete_after=True):
                        processed += 1
                except Exception as e:
                    console.print(f"[red]Error:[/red] {accession}: {e}")
                    if temp_path.exists():
                        temp_path.unlink()

                progress.update(task, advance=1)

        # Update frequencies
        console.print("\nUpdating frequency statistics...")
        calculator = FrequencyCalculator(db)
        calculator.compute_all_frequencies()

        # Show updated stats
        new_stats = db.get_stats()
        console.print(f"\n[green]Update complete![/green]")
        console.print(f"  Added: {new_stats['genome_count'] - stats['genome_count']} genomes")
        console.print(f"  Total: {new_stats['genome_count']} genomes")


@app.command("mitomaster-query")
def mitomaster_query(
    position: Annotated[
        int,
        typer.Argument(help="mtDNA position to query (1-16569)"),
    ],
    database: Annotated[
        Path,
        typer.Option("--database", "-d", help="Path to MitoMaster database"),
    ] = Path.home() / ".mtclassify" / "mitomaster.db",
    ref: Annotated[
        Optional[str],
        typer.Option("--ref", "-r", help="Reference base filter"),
    ] = None,
    alt: Annotated[
        Optional[str],
        typer.Option("--alt", "-a", help="Alternate base filter"),
    ] = None,
):
    """Query variant frequency at a specific position.

    Example:
        mtclassify mitomaster-query 3243
        mtclassify mitomaster-query 3243 --ref A --alt G
    """
    if not database.exists():
        console.print(f"[red]Error:[/red] Database not found: {database}")
        console.print("Use 'mtclassify mitomaster-build' to create a database.")
        raise typer.Exit(1)

    if position < 1 or position > 16569:
        console.print("[red]Error:[/red] Position must be between 1 and 16569")
        raise typer.Exit(1)

    from mtclassify.mitomaster import MitoMasterDB

    with MitoMasterDB(database) as db:
        # Get gene info
        gene_info = db.get_gene_for_position(position)
        if gene_info:
            console.print(f"[bold]Position {position}[/bold] ({gene_info[0]} - {gene_info[1]})")
        else:
            console.print(f"[bold]Position {position}[/bold]")

        # Get variants
        variants = db.get_variant_at_position(position, ref, alt)

        if not variants:
            console.print("No variants found at this position.")
            return

        # Display results
        table = Table(title="Variant Frequencies")
        table.add_column("Ref", style="cyan")
        table.add_column("Alt", style="magenta")
        table.add_column("Count", justify="right")
        table.add_column("Frequency", justify="right")
        table.add_column("Human Count", justify="right")
        table.add_column("Human Freq", justify="right")
        table.add_column("AA Change")

        for v in variants:
            table.add_row(
                v.ref,
                v.alt,
                str(v.total_count),
                f"{v.total_frequency:.4f}",
                str(v.human_count),
                f"{v.human_frequency:.4f}",
                v.amino_acid_change or "",
            )

        console.print(table)

        # Show haplogroup breakdown for top variant
        if variants:
            top_v = variants[0]
            hg_freqs = db.get_haplogroup_frequencies(position, top_v.ref, top_v.alt)

            if hg_freqs:
                console.print(f"\n[bold]Haplogroup breakdown for {top_v.ref}>{top_v.alt}:[/bold]")
                hg_table = Table()
                hg_table.add_column("Haplogroup", style="cyan")
                hg_table.add_column("Count", justify="right")
                hg_table.add_column("Frequency", justify="right")

                for hf in hg_freqs[:10]:  # Top 10
                    hg_table.add_row(
                        hf.haplogroup,
                        str(hf.count),
                        f"{hf.frequency:.4f}",
                    )

                console.print(hg_table)


@app.command("mitomaster-stats")
def mitomaster_stats(
    database: Annotated[
        Path,
        typer.Option("--database", "-d", help="Path to MitoMaster database"),
    ] = Path.home() / ".mtclassify" / "mitomaster.db",
):
    """Display MitoMaster database statistics.

    Example:
        mtclassify mitomaster-stats
    """
    if not database.exists():
        console.print(f"[red]Error:[/red] Database not found: {database}")
        console.print("Use 'mtclassify mitomaster-build' to create a database.")
        raise typer.Exit(1)

    from mtclassify.mitomaster import MitoMasterDB

    with MitoMasterDB(database) as db:
        stats = db.get_stats()

        console.print(f"[bold]MitoMaster Database Statistics[/bold]")
        console.print(f"Database: {database}")
        console.print()
        console.print(f"Total genomes: {stats['genome_count']:,}")
        console.print(f"Human genomes: {stats['human_genome_count']:,}")
        console.print(f"Total variants: {stats['variant_count']:,}")
        console.print(f"Unique variants: {stats['unique_variant_count']:,}")

        if stats['top_haplogroups']:
            console.print("\n[bold]Top Haplogroups:[/bold]")
            table = Table()
            table.add_column("Haplogroup", style="cyan")
            table.add_column("Count", justify="right")

            for hg in stats['top_haplogroups'][:15]:
                table.add_row(hg['haplogroup'], str(hg['count']))

            console.print(table)


@app.command("mitomaster-import")
def mitomaster_import(
    input_dir: Annotated[
        Path,
        typer.Argument(help="Directory containing GenBank (.gb) files"),
    ],
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Output database file path"),
    ] = Path.home() / ".mtclassify" / "mitomaster.db",
    tree: Annotated[
        str,
        typer.Option("--tree", "-t", help="Tree ID for haplogroup classification"),
    ] = "phylotree-rcrs@17.2",
    batch_size: Annotated[
        int,
        typer.Option("--batch-size", help="Genomes per batch for frequency updates"),
    ] = 500,
    pattern: Annotated[
        str,
        typer.Option("--pattern", "-p", help="File pattern to match"),
    ] = "*.gb",
    keep_files: Annotated[
        bool,
        typer.Option("--keep-files/--delete-files", help="Keep source files after import"),
    ] = True,
):
    """Import local GenBank files into MitoMaster database.

    Processes existing GenBank files from a directory, classifies them with
    mtclassify, extracts variants, and stores in the database.

    Example:
        mtclassify mitomaster-import /path/to/genomes/
        mtclassify mitomaster-import /path/to/genomes/ --pattern "*.gbk"
    """
    if not input_dir.exists():
        console.print(f"[red]Error:[/red] Directory not found: {input_dir}")
        raise typer.Exit(1)

    from rich.progress import Progress

    from mtclassify.mitomaster import (
        MitoMasterDB,
        GenomeProcessor,
        FrequencyCalculator,
    )

    # Find all GenBank files
    gb_files = sorted(input_dir.glob(pattern))
    if not gb_files:
        console.print(f"[red]Error:[/red] No files matching '{pattern}' found in {input_dir}")
        raise typer.Exit(1)

    console.print(f"[bold]MitoMaster Local Import[/bold]")
    console.print(f"Input: {input_dir}")
    console.print(f"Files found: {len(gb_files)}")
    console.print(f"Output: {output}")
    console.print(f"Tree: {tree}")
    console.print()

    # Create output directory
    output.parent.mkdir(parents=True, exist_ok=True)

    # Initialize database
    console.print("Initializing database...")
    with MitoMasterDB(output) as db:
        db.initialize()

        # Filter existing
        new_files = []
        for f in gb_files:
            # Use filename (without extension) as accession check
            accession = f.stem
            if not db.genome_exists(accession):
                new_files.append(f)

        console.print(f"New genomes to process: {len(new_files)}")

        if not new_files:
            console.print("[green]Database is up to date![/green]")
            return

        # Initialize processor
        processor = GenomeProcessor(db, tree_id=tree)

        # Process genomes with progress bar
        processed = 0
        failed = 0

        with Progress() as progress:
            task = progress.add_task("Processing genomes...", total=len(new_files))

            for gb_path in new_files:
                try:
                    success = processor.process_and_store(gb_path, delete_after=not keep_files)
                    if success:
                        processed += 1
                    else:
                        failed += 1
                except Exception as e:
                    console.print(f"[red]Error processing {gb_path.name}:[/red] {e}")
                    failed += 1

                progress.update(task, advance=1)

                # Periodic frequency updates
                if processed > 0 and processed % batch_size == 0:
                    console.print(f"\nUpdating frequencies (processed {processed})...")
                    calculator = FrequencyCalculator(db)
                    calculator.compute_all_frequencies()

        console.print(f"\n[green]Import complete![/green]")
        console.print(f"  Processed: {processed}")
        console.print(f"  Failed: {failed}")

        # Final frequency calculation
        console.print("\nComputing final frequency statistics...")
        calculator = FrequencyCalculator(db)
        calculator.compute_all_frequencies()

        # Show stats
        stats = db.get_stats()
        console.print(f"\n[bold]Database Statistics:[/bold]")
        console.print(f"  Total genomes: {stats['genome_count']}")
        console.print(f"  Human genomes: {stats['human_genome_count']}")
        console.print(f"  Total variants: {stats['variant_count']}")
        console.print(f"  Unique variants: {stats['unique_variant_count']}")


if __name__ == "__main__":
    app()
