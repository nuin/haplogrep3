"""Command-line interface for haplogrep3."""

from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.table import Table

from haplogrep3 import __version__
from haplogrep3.distance import Distance
from haplogrep3.io import PhylotreeLoader, VcfReader, FastaReader, TsvReader
from haplogrep3.tasks import ClassificationTask, export_csv, export_fasta
from haplogrep3.tasks.export import export_qc_report

app = typer.Typer(
    name="haplogrep3",
    help="mtDNA haplogroup classification tool",
    no_args_is_help=True,
)
console = Console()


def version_callback(value: bool):
    if value:
        console.print(f"Haplogrep 3 {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        Optional[bool],
        typer.Option("--version", "-v", callback=version_callback, is_eager=True),
    ] = None,
):
    """Haplogrep 3 - mtDNA haplogroup classification tool."""
    pass


@app.command()
def classify(
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
    """Classify mtDNA samples to determine haplogroups."""
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
        console.print("  haplogrep3 install-tree <tree-id>")
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

    console.print(f"Starting Haplogrep3 server at http://{host}:{port}")
    console.print("API docs available at /docs")
    console.print("Press Ctrl+C to stop\n")

    uvicorn.run(
        "haplogrep3.api:app",
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


if __name__ == "__main__":
    app()
