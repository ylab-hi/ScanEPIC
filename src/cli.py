#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Josh Fry, Northwestern University, YangLab

CLI Modules structure:
    extract ->
                short
                long
                singlecell
    tools ->
                exitron2vcf
                lrannotate
                ...

"""
import click
import sys
import os
import traceback
from shutil import rmtree






#=============================================================================
# command line interface
#=============================================================================
@click.group()
def cli():
    pass

#=============================================================================
# top level group structure
#=============================================================================
@cli.group()
def extract():
    """Quantify exitron splicing."""
    pass



@cli.group()
def tools():
    """Tools for analyzing exitron splicing data."""

    pass

#=============================================================================
# subcommand structure
#=============================================================================
## extract, short reads
@extract.command()
@click.option('-i', '--input', 'input_',
              help = 'Path to input BAM file.',
              required=True)
@click.option('-g', '--genome',
              help = 'Path to genome FASTA file.',
              required=True)
@click.option('-r', '--reference-transcriptome',
              help = 'Path to reference transcriptome annotation (GTF/GFF format). Database file will be created if not found.',
              required=True)
@click.option('-o', '--out',
              help = 'Output filename.')
@click.option('-c', '--cores',
              help = 'Number of cores for parallel processing. If 0 then no parallel processing is used',
              default = 0)
@click.option('-m', '--mapq',
              help = 'Consider only reads with MAPQ >= cutoff',
              default = 50)
@click.option('-a', '--ao',
              help = 'AO cutoff',
              default = 2)
@click.option('-p', '--pso',
              help = 'PSO cutoff',
              default = 0.01)
@click.option('-id', '--id', 'id_',
              help = 'Add ID column of string in output file.',
              default = None)
@click.option('-vcf', '--vcf',
              help = 'Format output as a VCF file.',
              default = None)
def short(input_,
          genome,
          reference_transcriptome,
          out,
          cores,
          mapq,
          ao,
          pso,
          id_,
          vcf):
    """Extract exitrons from short-read RNA-seq data."""
    from .extract.short import short as extract_short
    extract_short.main(input_,
                       genome,
                       reference_transcriptome,
                       out,
                       cores,
                       mapq,
                       ao,
                       pso,
                       id_,
                       vcf)

## extract, single-cell reads
@extract.command()
@click.option('-i', '--input', 'input_',
              help = 'Path to input BAM list.',
              required=True)
@click.option('-t', '--cell-types', 'cell_types',
              help = 'Input TSV file mapping cell barcodes to cell types.',
              required=True)
@click.option('-g', '--genome',
              help = 'Path to genome FASTA file.',
              required=True)
@click.option('-r', '--reference-transcriptome',
              help = 'Path to reference transcriptome annotation (GTF/GFF format). Database file will be created if not found.',
              required=True)
@click.option('-o', '--out',
              help = 'Output filename.')
@click.option('-c', '--cores',
              help = 'Number of cores for parallel processing. If 0 then no parallel processing is used',
              default = 0)
@click.option('-mq', '--mapq',
              help = 'Consider only reads with MAPQ >= cutoff',
              default = 50)
@click.option('-m', '--unique-exitron-mol',
              help = 'AO cutoff',
              default = 2)
@click.option('-p', '--pso',
              help = 'PSO cutoff',
              default = 0.01)
@click.option('-al', '--alignment50',
              help = 'PSO cutoff',
              default = 0.7)
@click.option('-id', '--id', 'id_',
              help = 'Add ID column of string in output file.',
              default = None)
def single(input_,
           cell_types,
           genome,
           reference_transcriptome,
           out,
           cores,
           mapq,
           unique_exitron_mol,
           pso,
           alignment50,
           id_):
    """Extract exitrons from single-cell RNA-seq data."""
    from .extract.single import single as extract_single
    extract_single.main(input_,
               cell_types,
               genome,
               reference_transcriptome,
               out,
               cores,
               mapq,
               unique_exitron_mol,
               pso,
               alignment50,
               id_)

## extract, long reads
@extract.command()
@click.option('-i', '--input', 'input_',
              help = 'Path to input BAM file.',
              required=True)
@click.option('-g', '--genome',
              help = 'Path to genome FASTA file.',
              required=True)
@click.option('-r', '--reference-transcriptome',
              help = 'Path to reference transcriptome annotation (GTF/GFF format). Database file will be created if not found.',
              required=True)
@click.option('-o', '--out',
              help = 'Output filename.')
@click.option('-c', '--cores',
              help = 'Number of cores for parallel processing. If 0 then no parallel processing is used',
              default = 0)
@click.option('-m', '--mapq',
              help = 'Consider only reads with MAPQ >= cutoff',
              default = 50)
@click.option('-a', '--ao_min',
              help = 'AO cutoff',
              default = 2)
@click.option('-p', '--pso_min',
              help = 'PSO cutoff',
              default = 0.01)
@click.option('-j', '--jitter',
              help = 'Treat splice-sites with fuzzy boundry of +/- INT',
              default = 10)
@click.option('-cp', '--cluster-purity',
              help = 'Cluster purity cutoff.',
              default = 0)
@click.option('-sr', '--skip-realign',
              help = 'Skip realignment step. We suggest skipping if alignment data is very clean, e.g. HiFi pacbio reads',
              is_flag = True)
@click.option('-sa', '--save-abundance',
              help = 'Save transcript abundance information for downstream processing',
              is_flag = True)
@click.option('-id', '--id', 'id_',
              help = 'Add ID column of string in output file. Files are of the form: input.isoform.exitrons, input.isoform.normals',
              default = None)
def long(input_,
          genome,
          reference_transcriptome,
          out,
          cores,
          mapq,
          ao_min,
          pso_min,
          jitter,
          cluster_purity,
          skip_realign,
          save_abundance,
          id_):
    """Extract exitrons from long-read RNA-seq data."""
    from .extract.long import long as extract_long

    # create tmp path
    this_dir = os.getcwd()
    tmp_path = os.path.join(this_dir, f'scanexitron_tmp{os.getpid()}')
    try:
        os.mkdir(tmp_path)
    except FileExistsError:
        pass
    try:
        extract_long.main(tmp_path, input_,
                           genome,
                           reference_transcriptome,
                           out,
                           cores,
                           mapq,
                           ao_min,
                           pso_min,
                           jitter,
                           cluster_purity,
                           skip_realign,
                           save_abundance,
                           id_)
    except Exception as e:
        if e.__class__.__name__ == 'InterruptedError':
            sys.stderr.write("User interrupt!")
        else:
            traceback.print_exc()
        rmtree(tmp_path)
        sys.exit(1)





## tools
@tools.command()
@click.option('-i', '--input', 'input_file',
              help='Path to input exitron TSV file.',
              required=True)
@click.option('-o', '--output', 'output_vcf',
              help='Path to output VCF file. If not specified, uses input filename with .vcf extension.')
@click.option('-g', '--genome',
              help='Path to reference genome FASTA file (optional, for extracting reference sequences).')
@click.option('-s', '--sample-name',
              help='Sample name for VCF header. If not specified, uses input filename.')
@click.option('--compress',
              help='Compress output VCF with bgzip.',
              is_flag=True)
@click.option('--index',
              help='Index compressed VCF with tabix (requires --compress).',
              is_flag=True)
def exitron2vcf(input_file, output_vcf, genome, sample_name, compress, index):
    """Convert exitron TSV files to VCF format."""
    from .tools.exitron2vcf import main as exitron2vcf_main
    exitron2vcf_main(input_file, output_vcf, genome, sample_name, compress, index)

@tools.command()
@click.option('--test3', default = 2, help = 'cube me')
def lrannotate(test3):
    click.echo(test3**3)

if __name__ == '__main__':
    cli()
