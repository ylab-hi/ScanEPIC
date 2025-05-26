#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exitron to VCF conversion tool for ScanEPIC

This module converts exitron TSV files to VCF format for compatibility with 
standard genomic analysis tools and pipelines.

@author: Josh Fry, Northwestern University, YangLab
"""

import os
import sys
import click
import pysam
import pandas as pd
from time import strftime, localtime


def pretty_print(text):
    """Pretty print with timestamp."""
    click.echo(click.style(f'[scanepic tools exitron2vcf -- {strftime("%Y-%m-%d %I:%M:%S %p", localtime())}]',
                           fg='blue') + '\t' + text)


def validate_exitron_file(exitron_file):
    """
    Validate that the exitron file has the required columns.
    
    Parameters
    ----------
    exitron_file : str
        Path to exitron TSV file
        
    Returns
    -------
    pd.DataFrame
        Validated exitron data
    """
    required_columns = ['chrom', 'start', 'end', 'name', 'ao', 'strand', 
                       'gene_symbol', 'length', 'splice_site', 'pso', 'dp']
    
    try:
        df = pd.read_csv(exitron_file, sep='\t')
    except Exception as e:
        pretty_print(f'ERROR: Could not read exitron file: {e}')
        sys.exit(1)
    
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        pretty_print(f'ERROR: Missing required columns: {missing_columns}')
        pretty_print(f'Available columns: {list(df.columns)}')
        sys.exit(1)
    
    pretty_print(f'Successfully loaded {len(df)} exitrons from {exitron_file}')
    return df


def generate_vcf_header(sample_name, reference_genome=None):
    """
    Generate VCF header with appropriate metadata.
    
    Parameters
    ----------
    sample_name : str
        Sample name for the VCF
    reference_genome : str, optional
        Reference genome file path
        
    Returns
    -------
    str
        VCF header string
    """
    header = f'''##fileformat=VCFv4.2
##source=ScanEPIC_exitron2vcf
##fileDate={strftime("%Y%m%d", localtime())}
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=DEL,Description="Deletion representing exitron splicing event">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the exitron">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the exitron (negative for deletion)">
##INFO=<ID=AO,Number=1,Type=Integer,Description="Number of reads supporting the exitron">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at exitron position">
##INFO=<ID=PSO,Number=1,Type=Float,Description="Percent spliced-out value">
##INFO=<ID=PSI,Number=1,Type=Float,Description="Percent spliced-in value (1-PSO)">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand of the exitron">
##INFO=<ID=SPLICE_SITE,Number=1,Type=String,Description="Splice site motif (e.g. GT-AG)">
##INFO=<ID=GENE_NAME,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=GENE_ID,Number=1,Type=String,Description="Gene identifier">
##INFO=<ID=REGION,Number=1,Type=String,Description="Genomic region (CDS, 3UTR, 5UTR)">
##INFO=<ID=TRANSCRIPT_ID,Number=1,Type=String,Description="Transcript identifier">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AO,Number=1,Type=Integer,Description="Alternate observation count">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=PSO,Number=1,Type=Float,Description="Percent spliced-out">'''
    
    if reference_genome:
        header += f'\n##reference={reference_genome}'
    
    header += f'\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n'
    
    return header


def exitron_to_vcf_record(exitron, genome_fasta=None, sample_name='SAMPLE'):
    """
    Convert a single exitron record to VCF format.
    
    Parameters
    ----------
    exitron : pd.Series
        Single exitron record
    genome_fasta : pysam.FastaFile, optional
        Reference genome for extracting sequences
    sample_name : str
        Sample name for genotype field
        
    Returns
    -------
    str
        VCF record line
    """
    chrom = exitron['chrom']
    pos = int(exitron['start'])
    end_pos = int(exitron['end'])
    exitron_id = exitron['name']
    ao = int(exitron['ao'])
    dp = int(exitron['dp']) if pd.notna(exitron['dp']) else ao
    pso = float(exitron['pso']) if pd.notna(exitron['pso']) else 0.0
    psi = 1.0 - pso
    strand = exitron['strand']
    length = int(exitron['length'])
    splice_site = exitron['splice_site']
    gene_symbol = exitron['gene_symbol']
    
    # Get reference and alternate sequences
    if genome_fasta:
        try:
            ref_seq = genome_fasta[chrom][pos-1:pos].upper()
            alt_seq = genome_fasta[chrom][pos-1:end_pos-1].upper()
        except Exception:
            # Fallback if genome access fails
            ref_seq = 'N'
            alt_seq = 'N' * (length + 1)
    else:
        ref_seq = 'N'
        alt_seq = 'N' * (length + 1)
    
    # Build INFO field
    info_fields = [
        f'SVTYPE=DEL',
        f'END={end_pos-1}',
        f'SVLEN=-{length}',
        f'AO={ao}',
        f'DP={dp}',
        f'PSO={pso:.4f}',
        f'PSI={psi:.4f}',
        f'STRAND={strand}',
        f'SPLICE_SITE={splice_site}',
        f'GENE_NAME={gene_symbol}'
    ]
    
    # Add optional fields if present
    if 'gene_id' in exitron and pd.notna(exitron['gene_id']):
        info_fields.append(f'GENE_ID={exitron["gene_id"]}')
    if 'region' in exitron and pd.notna(exitron['region']):
        info_fields.append(f'REGION={exitron["region"]}')
    if 'transcript_id' in exitron and pd.notna(exitron['transcript_id']):
        info_fields.append(f'TRANSCRIPT_ID={exitron["transcript_id"]}')
    
    info_string = ';'.join(info_fields)
    
    # Quality score based on supporting reads and PSO
    qual = min(60, max(10, int(ao * 5 + pso * 100)))
    
    # Filter - PASS if meets basic criteria
    filter_field = 'PASS' if ao >= 2 and pso >= 0.01 else 'LowQual'
    
    # Format and sample fields
    format_field = 'GT:AO:DP:PSO'
    sample_field = f'0/1:{ao}:{dp}:{pso:.4f}'
    
    # Construct VCF record
    vcf_record = f'{chrom}\t{pos}\t{exitron_id}\t{ref_seq}\t{alt_seq}\t{qual}\t{filter_field}\t{info_string}\t{format_field}\t{sample_field}'
    
    return vcf_record


def convert_exitrons_to_vcf(exitron_file, output_vcf, genome_file=None, sample_name=None):
    """
    Convert exitron TSV file to VCF format.
    
    Parameters
    ----------
    exitron_file : str
        Path to input exitron TSV file
    output_vcf : str
        Path to output VCF file
    genome_file : str, optional
        Path to reference genome FASTA file
    sample_name : str, optional
        Sample name for VCF (inferred from filename if not provided)
    """
    # Load and validate exitron data
    exitrons_df = validate_exitron_file(exitron_file)
    
    # Determine sample name
    if not sample_name:
        sample_name = os.path.splitext(os.path.basename(exitron_file))[0]
        sample_name = sample_name.replace('.exitron', '').replace('.tsv', '')
    
    # Load genome if provided
    genome_fasta = None
    if genome_file:
        try:
            genome_fasta = pysam.FastaFile(genome_file)
            pretty_print(f'Using reference genome: {genome_file}')
        except Exception as e:
            pretty_print(f'WARNING: Could not load genome file {genome_file}: {e}')
            pretty_print('Proceeding without reference sequences')
    
    # Sort exitrons by chromosome and position for proper VCF format
    exitrons_df = exitrons_df.sort_values(['chrom', 'start'])
    
    pretty_print(f'Converting {len(exitrons_df)} exitrons to VCF format')
    
    # Write VCF file
    with open(output_vcf, 'w') as vcf_out:
        # Write header
        header = generate_vcf_header(sample_name, genome_file)
        vcf_out.write(header)
        
        # Convert each exitron to VCF record
        for idx, exitron in exitrons_df.iterrows():
            vcf_record = exitron_to_vcf_record(exitron, genome_fasta, sample_name)
            vcf_out.write(vcf_record + '\n')
    
    # Close genome file if opened
    if genome_fasta:
        genome_fasta.close()
    
    pretty_print(f'Successfully converted exitrons to VCF: {output_vcf}')


def main(exitron_file, output_vcf, genome_file, sample_name, compress, index):
    """
    Main function for exitron2vcf conversion.
    
    Parameters
    ----------
    exitron_file : str
        Input exitron TSV file
    output_vcf : str
        Output VCF file
    genome_file : str
        Reference genome FASTA file
    sample_name : str
        Sample name for VCF
    compress : bool
        Whether to compress output VCF
    index : bool
        Whether to index output VCF
    """
    # Validate input file
    if not os.path.exists(exitron_file):
        pretty_print(f'ERROR: Input file does not exist: {exitron_file}')
        sys.exit(1)
    
    # Set default output filename if not provided
    if not output_vcf:
        base_name = os.path.splitext(exitron_file)[0]
        output_vcf = f'{base_name}.vcf'
    
    # Validate genome file if provided
    if genome_file and not os.path.exists(genome_file):
        pretty_print(f'ERROR: Genome file does not exist: {genome_file}')
        sys.exit(1)
    
    # Convert exitrons to VCF
    convert_exitrons_to_vcf(exitron_file, output_vcf, genome_file, sample_name)
    
    # Compress VCF if requested
    if compress:
        pretty_print('Compressing VCF file with bgzip')
        try:
            pysam.tabix_compress(output_vcf, output_vcf + '.gz', force=True)
            os.remove(output_vcf)  # Remove uncompressed version
            output_vcf = output_vcf + '.gz'
            pretty_print(f'Compressed VCF: {output_vcf}')
        except Exception as e:
            pretty_print(f'WARNING: Could not compress VCF: {e}')
    
    # Index VCF if requested and compressed
    if index and output_vcf.endswith('.gz'):
        pretty_print('Indexing VCF file with tabix')
        try:
            pysam.tabix_index(output_vcf, preset='vcf', force=True)
            pretty_print(f'Indexed VCF: {output_vcf}.tbi')
        except Exception as e:
            pretty_print(f'WARNING: Could not index VCF: {e}')
    elif index and not output_vcf.endswith('.gz'):
        pretty_print('WARNING: VCF indexing requires compressed file. Use --compress flag.')


if __name__ == '__main__':
    # Test function
    pass
