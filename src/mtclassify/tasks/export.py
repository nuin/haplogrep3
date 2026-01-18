"""Export functions for classification results."""

import csv
from pathlib import Path
from typing import TextIO

from mtclassify.models import Sample


def export_csv(
    samples: list[Sample],
    output: Path,
    extended: bool = False,
    hits: int = 1,
) -> None:
    """Export classification results to CSV.

    Args:
        samples: Classified samples
        output: Output file path
        extended: Include extended information (polymorphism details)
        hits: Number of hits to include per sample
    """
    with open(output, "w", newline="") as f:
        _write_csv(f, samples, extended, hits)


def _write_csv(
    f: TextIO,
    samples: list[Sample],
    extended: bool,
    hits: int,
) -> None:
    """Write CSV to file handle.

    Args:
        f: File handle
        samples: Classified samples
        extended: Include extended information
        hits: Number of hits to include
    """
    # Build header
    header = ["SampleID", "Haplogroup", "Quality", "Range"]

    if extended:
        header.extend([
            "Found_Polys",
            "Missing_Polys",
            "Remaining_Polys",
            "Input_Sample",
        ])

    if hits > 1:
        for i in range(2, hits + 1):
            header.append(f"Haplogroup_{i}")
            header.append(f"Quality_{i}")

    writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
    writer.writeheader()

    for sample in samples:
        row = {
            "SampleID": sample.id,
            "Range": sample.range_str,
        }

        if sample.classification:
            result = sample.classification.top_result
            row["Haplogroup"] = result.haplogroup.name
            row["Quality"] = f"{result.quality:.4f}"

            if extended:
                row["Found_Polys"] = " ".join(str(p) for p in result.found_polymorphisms)
                row["Missing_Polys"] = " ".join(str(p) for p in result.missing_polymorphisms)
                row["Remaining_Polys"] = " ".join(str(p) for p in result.remaining_polymorphisms)
                row["Input_Sample"] = " ".join(sample.get_polymorphism_strings())

            # Additional hits
            if hits > 1:
                for i, other in enumerate(sample.classification.other_results[:hits - 1], start=2):
                    row[f"Haplogroup_{i}"] = other.haplogroup.name
                    row[f"Quality_{i}"] = f"{other.quality:.4f}"
        else:
            row["Haplogroup"] = "?"
            row["Quality"] = "0.0000"

        writer.writerow(row)


def export_fasta(
    samples: list[Sample],
    output: Path,
    reference: str,
) -> None:
    """Export samples as FASTA sequences.

    Reconstructs sequences from reference + polymorphisms.

    Args:
        samples: Classified samples
        output: Output file path
        reference: Reference sequence (rCRS)
    """
    with open(output, "w") as f:
        for sample in samples:
            # Start with reference sequence
            seq = list(reference)

            # Apply polymorphisms
            for poly in sample.polymorphisms:
                if 1 <= poly.position <= len(seq):
                    seq[poly.position - 1] = poly.mutation.value

            # Write FASTA entry
            haplogroup = "?"
            if sample.classification:
                haplogroup = sample.classification.haplogroup.name

            f.write(f">{sample.id} {haplogroup}\n")

            # Write sequence in 70-character lines
            seq_str = "".join(seq)
            for i in range(0, len(seq_str), 70):
                f.write(seq_str[i:i + 70] + "\n")


def export_qc_report(
    samples: list[Sample],
    output: Path,
) -> None:
    """Export quality control report.

    Args:
        samples: Classified samples
        output: Output file path
    """
    with open(output, "w", newline="") as f:
        header = [
            "SampleID",
            "Haplogroup",
            "Quality",
            "Expected_Polys",
            "Found_Polys",
            "Missing_Polys",
            "Private_Polys",
        ]

        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()

        for sample in samples:
            row = {"SampleID": sample.id}

            if sample.classification:
                result = sample.classification.top_result
                row["Haplogroup"] = result.haplogroup.name
                row["Quality"] = f"{result.quality:.4f}"
                row["Expected_Polys"] = len(result.expected_polymorphisms)
                row["Found_Polys"] = len(result.found_polymorphisms)
                row["Missing_Polys"] = len(result.missing_polymorphisms)
                row["Private_Polys"] = len(result.remaining_polymorphisms)
            else:
                row["Haplogroup"] = "?"
                row["Quality"] = "0.0000"
                row["Expected_Polys"] = 0
                row["Found_Polys"] = 0
                row["Missing_Polys"] = 0
                row["Private_Polys"] = len(sample.polymorphisms)

            writer.writerow(row)
