# ScanEPIC

ScanEPIC (Scan **E**xitron **P**latform **I**ntegration **C**aller) is a bioinformatics tool for identifying and quantifying exitron splicing events across multiple RNA-seq platforms.

## Overview

Exitrons are novel intronic regions with both splice sites contained within a single annotated exon, effectively mimicking genomic deletions. ScanEPIC implements stringent filtering criteria to distinguish genuine exitron splicing events from technical artifacts such as those arising from reverse transcription slippage during library preparation.


## Installation

### Prerequisites

- Python 3.8+
- Conda (recommended)

### Install from source

```bash
# Create conda environment
conda create -n scanepic python=3.9
conda activate scanepic

# Install core dependencies
conda install -c conda-forge numpy pandas scipy statsmodels
conda install -c bioconda pysam gffutils
conda install -c conda-forge click biopython cython

# Install additional dependencies
pip install liqa

# Clone and install ScanEPIC
git clone https://github.com/yourusername/scanepic.git
cd scanepic
pip install -e .
```

## Quick Start

### Test installation
```bash
scanepic --help
```

### Short-read RNA-seq
```bash
scanepic extract short \
    -i input.bam \
    -g genome.fa \
    -r annotation.gtf.gz \
    -o output.exitron \
    --cores 4
```

### Long-read RNA-seq
```bash
scanepic extract long \
    -i input.bam \
    -g genome.fa \
    -r annotation.gtf.gz \
    -o output.exitron \
    --cores 4
```

### Single-cell RNA-seq
```bash
scanepic extract single \
    -i bam_list.txt \
    -t cell_types.tsv \
    -g genome.fa \
    -r annotation.gtf.gz \
    -o output_directory \
    --cores 4
```

## Input Requirements

### Required Files
- **BAM files**: Aligned RNA-seq reads with index (.bai)
- **Genome FASTA**: Reference genome sequence with index (.fai)
- **GTF annotation**: Gene annotation file (must be bgzip compressed and tabix indexed)

### File Preparation
```bash
# Compress and index GTF file
bgzip annotation.gtf
tabix -p gff annotation.gtf.gz

# Index genome FASTA
samtools faidx genome.fa

# Index BAM file
samtools index input.bam
```

## Output Format

ScanEPIC produces tab-separated files containing:

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| start | Exitron start position |
| end | Exitron end position |
| name | Exitron identifier |
| region | Genomic region (CDS, 3'UTR, 5'UTR) |
| ao | Alternative observations (supporting reads) |
| strand | Strand orientation |
| gene_symbol | Gene symbol |
| length | Exitron length |
| splice_site | Splice site motif |
| transcript_id | Transcript identifier |
| pso | Percent spliced-out |
| dp | Read depth |
| alignment50 | Alignment quality score |
| rt_repeat | RT-PCR repeat content |




## Citation

If you use ScanEPIC in your research, please cite:

> Fry J, Liu Q, Lu X, et al. 3'UTR exitron splicing reshapes the 3' regulatory landscape in cancer. *Manuscript in preparation.*

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

For questions, issues, or feature requests:
- Open an issue on GitHub
- Contact: joshua.fry@northwestern.edu
