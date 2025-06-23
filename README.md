[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15556565.svg)](https://doi.org/10.5281/zenodo.15556565) [![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/) [![License](https://img.shields.io/pypi/l/scanexitronlr)](https://opensource.org/licenses/MIT "License") 


# ScanEPIC

ScanEPIC (Scan **E**xitron **P**latform **I**ntegration **C**aller) is a bioinformatics tool for identifying and quantifying exitron splicing events across multiple RNA-seq platforms.

## Overview

Exitrons are novel intronic regions with both splice sites contained within a single annotated exon, effectively mimicking genomic deletions. ScanEPIC implements stringent filtering criteria to distinguish genuine exitron splicing events from technical artifacts such as those arising from reverse transcription slippage during library preparation.

Typical runtime is around ~10 minutes using 4 cores on an average sized RNA-seq experiment.


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
git clone git@github.com:ylab-hi/ScanEPIC.git
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

NOTE: The scRNA module takes as input a list of BAM files (path to BAM file on each line). If an exitron is detected in any BAM file, it will also be quantified in the BAM files where the exitron was not detected (for downstream analysis). Also required is a TSV file with three columns, e.g.:

| cells | sample_id | cell_group |
|--------|-------------|-------------|
| AACGGTACAAAAGC-1 | BT1249 | immune | 
| AACTCGGATCGCCT-1 | BT1301 | cancer |
| AAGAAGACACTGGT-1 | BT1300 | blood |


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
- **GTF annotation**: Gene annotation file (must be bgzip compressed and tabix indexed). the `gffutils` package will create a database file when ScanEPIC runs for the first time. This may take ~20 minutes but will only need to be done once. 

### File Preparation
```bash

# If needed, sort your GTF annotation
awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' in.gtf > out_sorted.gtf

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


## Test Data

ScanEPIC includes test datasets (generated using gencode v37) to verify installation and demonstrate usage:

### Short-read RNA-seq Test
```bash
scanepic extract short \
    -i test_data/test_data.bam \
    -g hg38.fa \
    -r annotation.gtf.gz \
    -o test_output_short.exitron \
    --cores 4
```

### Long-read RNA-seq Test
```bash
scanepic extract long \
    -i test_data/long_test_data.bam \
    -g hg38.fa \
    -r annotation.gtf.gz \
    -o test_output_long.exitron \
    --cores 4
```

### Single-cell RNA-seq Test
```bash
scanepic extract single \
    -i test_data/sc_test_bam_list.txt \
    -t test_data/sc_test_cell_types.tsv \
    -g hg38.fa \
    -r annotation.gtf.gz \
    -o test_output_sc \
    --cores 4
```

**Note**: The test BAM files are provided, but you must supply your own reference genome (hg38.fa) and annotation files (annotation.gtf.gz).

## Citation

If you use ScanEPIC in your research, please cite:

> Fry J, Liu Q, Lu X, et al. 3'UTR exitron splicing reshapes the 3' regulatory landscape in cancer. *Manuscript in preparation.*

## License

This project is licensed under the MIT License - see the LICENSE file for details.


## Support

For questions, issues, or feature requests:
- Open an issue on GitHub
- Contact: joshua.fry@northwestern.edu
