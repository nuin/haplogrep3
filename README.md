# mtclassify

Python reimplementation of [Haplogrep3](https://github.com/genepi/haplogrep3) - a mtDNA haplogroup classification tool.

## Installation

```bash
# Using uv (recommended)
uv sync

# Or with pip
pip install -e .
```

## Usage

### CLI

```bash
# List available trees
mtclassify trees

# Classify samples from VCF
mtclassify run --input samples.vcf --tree phylotree-fu-rcrs@1.2 --output results.txt

# Classify with options
mtclassify run \
  --input samples.vcf \
  --tree phylotree-fu-rcrs@1.2 \
  --output results.txt \
  --distance kulczynski \
  --hits 10 \
  --extend-report \
  --write-fasta \
  --write-qc
```

### As a Library

```python
from pathlib import Path
from mtclassify.io import VcfReader, load_phylotree
from mtclassify.tasks import ClassificationTask
from mtclassify.distance import Distance

# Load samples
reader = VcfReader()
samples = reader.read(Path("samples.vcf"))

# Load phylotree
phylotree = load_phylotree("phylotree-fu-rcrs@1.2")

# Classify
task = ClassificationTask(phylotree, Distance.KULCZYNSKI, hits=1)
results = task.classify(samples)

# Print results
for sample in results:
    if sample.classification:
        print(f"{sample.id}: {sample.classification.haplogroup.name}")
```

## Supported Input Formats

- **VCF** (.vcf, .vcf.gz) - Variant Call Format
- **FASTA** (.fasta, .fa, .fasta.gz, .fa.gz) - Sequence format

## Distance Metrics

- **Kulczynski** (default) - Best for haplogroup classification
- **Hamming** - Simple position differences
- **Jaccard** - Set similarity
- **Kimura** - Transition/transversion weighted

## Development

```bash
# Run tests
uv run pytest tests/ -v

# Run with coverage
uv run pytest tests/ --cov=mtclassify
```

## License

MIT License - see original [Haplogrep3](https://github.com/genepi/haplogrep3) for details.
