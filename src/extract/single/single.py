# -*- coding: utf-8 -*-
"""
Single-cell exitron extraction module for ScanEPIC

This module identifies and quantifies exitron splicing events from single-cell RNA-seq data.
It processes cell barcodes and UMIs to provide cell-type-specific exitron quantification.

@author: Josh Fry, Northwestern University, YangLab
@original_author: Tingyou Wang, Northwestern University, YangLab (v1)
"""
# ===============================================================================

__version__ = 'v2'
import sys
import os
import argparse
import pysam
import traceback
import gffutils
import multiprocessing as mp
import pandas as pd
import warnings
import click
import re
from collections import defaultdict
from shutil import rmtree
# TODO replace pairwise2 with Bio.Align.PairwiseAligner
from Bio import Align
from statsmodels.stats.proportion import proportions_ztest
from scipy.stats import combine_pvalues
from time import strftime, localtime
from .cython_helpers import find_introns, get_unique_mols





#=============================================================================
# Helpers
#=============================================================================

def pretty_print(text):
    click.echo(click.style(f'[exitrontools extract single -- {strftime("%Y-%m-%d %I:%M:%S %p", localtime())}]',
                           fg = 'green') + '\t' + text)

def repeat_test(seq_3, seq_5, kmer_min, kmer_max):
    for k in range(kmer_max, kmer_min - 1, -1):
        kmers = [seq_5[i:i+k] for i in range(len(seq_5) - (k - 1))]
        matches = [kmer for kmer in kmers if re.search(kmer, seq_3)]
        if matches:
            return matches[0]
        seq_5 = seq_5[1:-1]
        seq_3 = seq_3[1:-1]

    return ''

# ===============================================================================
# Modules
# ===============================================================================





def exitron_caller(bamfile, referencename, chrm, db, mapq=50):
    """


    Parameters
    ----------
    bamfile : pysam.AlignmentFile
    referencename : str
    chrms : list
    mapq : int, optional
        Only considers reads from bamfile with quality >= mapq. The default is 50.
    utr : bool
        if True, UTR exitrons will be called.

    Returns
    -------
    List of unfiltered exitrons.  Each exitron is a dictionary of features.

    """
    # call introns
    introns, reads = find_introns(
        (read for read in bamfile.fetch(chrm) if read.mapping_quality >= mapq)
    )

    gtf = pysam.TabixFile(referencename, parser=pysam.asGTF())
    exitrons = []
    exitrons_added = []
    known_splices = set()


    for intron in introns:
        intron_start = intron[0]
        intron_end = intron[1]
        intron_witnesses = introns[intron]
        intersection = gtf.fetch(chrm, intron_start - 1, intron_end + 2)
        for feature in intersection:
            region_type = feature.feature
            region_start = feature.start
            region_end = feature.end
            gene_name = feature.gene_name
            gene_id = feature.gene_id
            if (region_type == 'exon' and
                'nonsense_mediated_decay' not in feature.transcript_type):
                if intron_start == region_end:
                    # intron matches a known donor
                    known_splices.add(
                        (chrm, intron_start, intron_end + 1, 'D', feature.strand))
                if intron_end == region_start:
                    # intron matches a known acceptor
                    known_splices.add(
                        (chrm, intron_start, intron_end + 1, 'A', feature.strand))

            elif region_type == 'CDS' and region_start < intron_start \
                    and region_end > intron_end:
                if (intron_start, intron_end, region_type, feature.strand) not in exitrons_added:
                    exitrons.append({'chrom': chrm,
                                    'start': intron_start,
                                     'end': intron_end + 1,  # plus 1 because of bedtools conventions,
                                     'name': f'{gene_name}{intron_start}{intron_end + 1}',
                                     'region': region_type,
                                     'read_total': intron_witnesses,
                                     'strand': feature.strand,
                                     'gene_symbol': gene_name,
                                     'gene_id': gene_id,
                                     'length': intron_end + 1 - (intron_start) - 1,
                                     'splice_site': 'splice_site',
                                     'transcript_id': feature.transcript_id})
                    exitrons_added.append((intron_start, intron_end, region_type, feature.strand))

            elif region_type == 'UTR':
                if feature.transcript_type == 'protein_coding' and \
                    (intron_start, intron_end, region_type, feature.strand) not in exitrons_added:
                        # require that exitron is within exon
                        exon_intervals = [(x.start, x.end) for x in db.children(db[feature.transcript_id],
                                                                                limit = (chrm,
                                                                                         intron_start - 10,
                                                                                         intron_end + 10,
                                                                                         ),
                                                                                featuretype = 'exon')]
                        if any(map(lambda x: x[0] < intron_start and intron_end < x[1], exon_intervals)):

                            # determine 5' or 3'
                            cds_start = [(x.start, x.end) for x in db.children(db[feature.transcript_id], featuretype = 'CDS')]
                            utr_front = any(map(lambda x: region_end <= x[0], cds_start))
                            if utr_front:
                                region = '5\'UTR' if feature.strand == '+' else '3\'UTR'
                            else:
                                region = '3\'UTR' if feature.strand == '+' else '5\'UTR'

                            # determine overlap with CDS
                            if any(map(lambda x: x[0] <= intron_start <= x[1] or x[0] <= intron_end <= x[1], cds_start)):
                                region += ' + CDS'

                            exitrons.append({'chrom': chrm,
                                            'start': intron_start,
                                             'end': intron_end + 1,  # plus 1 because of bedtools conventions,
                                             'name': f'{gene_name}{intron_start}{intron_end + 1}',
                                             'region': region,
                                             'read_total': intron_witnesses,
                                             'strand': feature.strand,
                                             'gene_symbol': gene_name,
                                             'gene_id': gene_id,
                                             'length': intron_end + 1 - (intron_start) - 1,
                                             'splice_site': 'splice_site',
                                             'transcript_id': feature.transcript_id})
                            exitrons_added.append((intron_start, intron_end, region_type, feature.strand))







    return ([exitron for exitron in exitrons if ((exitron['chrom'], exitron['start'], exitron['end'], 'D', exitron['strand']) not in known_splices
                                                  and (exitron['chrom'], exitron['start'], exitron['end'], 'A', exitron['strand']) not in known_splices)],
            reads,
            )



def filter_exitrons(exitrons, reads, bamfile, genome, db, cell_types, mapq = 50, alignment50_cutoff = 0.7, pso_min = 0.01, m_min = 2):
    """
    Parameters
    ----------
    exitrons : list
        list of unfiltered exitrons from exitron_caller.
    reads : dict
        Each intron is a key, and the value is list of (read_seq, ref_seq, left_anchor, right_anchor)
    bamfile : pysam.AlignmentFile
    genome : pysam.libcfaidx.FastaFile
        Random access to reference genome.
    meta_data : dict
    verbose : bool
        If true, we report optional statistics
    mapq : TYPE, optional
        Only considers reads from bamfile with quality >= mapq. The default is 50.
        This is needed to calculate pso.
    pso_min : float, optional
        Number reads that witness the exitron over the number of reads within the
        spliced out region. The default is 0.05.
    ao_min : int, optional
        Minimum number of reads witnessing the exitron. The default is 3.

    pso_ufmin : float, optional
        Number reads that witness the exitron over the number of reads within the
        spliced out region, before filtering (used for backwards compatibility). The default is 0.05.
    ao_ufmin : int, optional
        Minimum number of reads witnessing the exitron, before filtering (used for
        backwards compatibility). The default is 3.
    anchor_min : int, optional
        Minimum anchor length.  The default is 5.

    Returns
    -------
    filtered exitrons and meta data
    """

    #=============================================================================
    # helpers
    #=============================================================================
    # Need to compute reverse complement for splice junctions.
    tab = str.maketrans("ACTG", "TGAC")
    def rc(seq):
        return seq.upper().translate(tab)[::-1]

    # define the aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1

    genome_fasta = pysam.FastaFile(genome)


    res = []

    cell_type_names = [x for x in cell_types.cell_group.unique() if not pd.isna(x)]

    # adjust exitrons for cell barcode and unique molecules
    # filter
    for exitron in exitrons:
        # Each unprocessed exitron needs to be cleaved into exitrons from
        #   different cell types.
        # filter by total AO
        if exitron['read_total'] < m_min:
            continue
        strand = exitron['strand']
        chrm = exitron['chrom']
        start = exitron['start']
        end = exitron['end']

        if strand == '+':
            pos = start - db[exitron['transcript_id']].start + 1
            genome_seq = db[exitron['transcript_id']].sequence(
                genome)[pos-10:pos + exitron['length'] + 10].upper()
            intron_seq = genome_seq[10:-10]
            splice_site = intron_seq[:2] + '-' + intron_seq[-2:]
        elif strand == '-':
            pos = db[exitron['transcript_id']].end - start
            genome_seq = db[exitron['transcript_id']].sequence(
                genome)[pos - exitron['length'] -10 :pos + 10].upper() # NB: db.sequence returns - strand seq if gene is -
            intron_seq = genome_seq[10:-10]
            left = intron_seq[:2]
            right = intron_seq[-2:]
            splice_site = left + '-' + right
            intron_seq = rc(intron_seq) # tests below assume all sequences are + stranded
        if splice_site not in ['GT-AG', 'GC-AG', 'AT-AC']:
            continue
        exitron['splice_site'] = splice_site

        # Now, examine reads supporting the exitron splicing event
        read_data = reads[(start, end - 1)]

        #=============================================================================
        # Require unique anchors
        #=============================================================================
        passed_reads = []
        for read in read_data:
            l_anchor_len = read[2]
            r_anchor_len = read[3]

            # check 5' intron and 3' exon similarity
            three_prime_exon_r = read[0][read[4]:read[4] + r_anchor_len]
            five_prime_intron_r = intron_seq[:len(three_prime_exon_r)]
            # check 5' exon and 3' intron similarity
            five_prime_exon_r = read[0][read[4] - l_anchor_len:read[4]]
            three_prime_intron_r = intron_seq[-len(five_prime_exon_r):]

            if three_prime_exon_r != five_prime_intron_r and five_prime_exon_r != three_prime_intron_r:
                passed_reads.append(read)

        if not passed_reads:
            continue

        #=============================================================================
        # Calculate alignment50 scores
        #=============================================================================
        try:
            alignment50_lr = aligner.score(intron_seq[:50],
                                           genome_seq[-50:]) / 50
        except:
            alignment50_lr = 'NA'
        try:
            alignment50_rl = aligner.score(genome_seq[:50],
                                           intron_seq[-50:]) / 50
        except:
            alignment50_rl = 'NA'

        # TODO: make this an argument
        if max(alignment50_lr, alignment50_rl) >= alignment50_cutoff:
            continue
        exitron['alignment50'] = max(alignment50_lr, alignment50_rl)

        #=============================================================================
        # Calculate RT repeats
        #=============================================================================
        # calculate repeats
        seq_5 = genome_seq[50 - 19: 50 + 21]
        seq_3 = genome_seq[-50 - 21: 19 - 50]
        exitron['rt_repeat'] = repeat_test(seq_3, seq_5, 4, 20)

        #=============================================================================
        # Calculate AO and PSO for each cell type
        #=============================================================================
        cells_filtered = len(set(x[6] for x in read_data if x[5] not in cell_types.cells.unique()))/len(read_data)
        exitron['exitron_cells_not_found'] = round(cells_filtered, 5)
        exitrons_to_test = []
        for cell_type in cell_type_names:
            exitron_cell_type = exitron.copy()
            exitron_cell_type['cell_type'] = cell_type
            cells_in_type = cell_types[cell_types['cell_group'] == cell_type].cells.unique()
            unique_exitron_mol = len(set(x[6] for x in read_data if x[5] in cells_in_type))


            unique_mols = get_unique_mols(chrm,
                                          exitron['start'],
                                          exitron['end'],
                                          bamfile,
                                          set(cells_in_type),
                                          mapq)

            exitron_cell_type['exitron_mols'] = unique_exitron_mol
            exitron_cell_type['unique_mols'] = unique_mols
            exitron_cell_type['pso'] = exitron_cell_type['exitron_mols']/exitron_cell_type['unique_mols'] \
                if exitron_cell_type['unique_mols'] else 0.0

            exitrons_to_test.append(exitron_cell_type)


        #=============================================================================
        # Filter by PSO
        #=============================================================================
        if not any(exitron['pso'] > pso_min for exitron in exitrons_to_test) \
           or sum(exitron['exitron_mols'] for exitron in exitrons_to_test) < m_min:
            continue


        #=============================================================================
        # Calculate pvals and add to final list
        #=============================================================================
        # divide by zero warnings occur when cells have dropout.
        #   in these cases, pval is forced to be 1.0
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        for e1 in exitrons_to_test:
            pvals = []
            for e2 in exitrons_to_test:
                count = (e1['exitron_mols'] + 1e-5, e2['exitron_mols'] + 1e-5)
                nobs = (e1['unique_mols'], e2['unique_mols'])
                pval = proportions_ztest(count = count, nobs = nobs)[1]
                if pd.isna(pval):
                    pval = 1.0
                pvals.append(pval)
            e1['pairwise_z_test_pvals'] = ','.join(map(lambda x: str(round(x, 5)), pvals))
            # remove diagonal comparison
            pvals.remove(max(pvals))
            e1['meta_pval'] = round(combine_pvalues(pvals, method = 'fisher')[1], 5)

            # append to final list
            res.append(e1)











    return res


# ===============================================================================
# Main
# ===============================================================================


def exitrons_in_chrm(bamfilename, referencename, genomename, chrm, cell_types, mapq, alignment50_cutoff, pso_min, m_min):
    """
    Wrapper that calls main functions *per chromosome*.
    """
    pretty_print(f'Finding exitrons in {chrm}')
    sys.stdout.flush()
    bamfile = pysam.AlignmentFile(bamfilename, 'rb', require_index=True)
    db = gffutils.FeatureDB(referencename + '.db')
    exitrons, reads = exitron_caller(bamfile,
                                      referencename,
                                      chrm,
                                      db,
                                      mapq,
                                      )
    exitrons = filter_exitrons(exitrons,
                                reads,
                                bamfile,
                                genomename,
                                db,
                                cell_types,
                                mapq,
                                alignment50_cutoff,
                                pso_min,
                                m_min,
                                )
    bamfile.close()
    del db

    return exitrons, chrm



def call_exitrons(bamfilename, chrms, sample_id,
                  input_,
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
    #=============================================================================
    # Given a bamfilename, call exitrons and write to file
    #=============================================================================

    try:
        bamfile = pysam.AlignmentFile(bamfilename, 'rb', require_index=True)
    except FileNotFoundError:
        try:
            pretty_print('Building bam index file')
            pysam.index(bamfilename)
            bamfile = pysam.AlignmentFile(bamfilename, 'rb', require_index=True)
        except FileNotFoundError:
            pretty_print(
                f'ERROR: There is a problem opening bam file at: {bamfilename}')
            sys.exit(1)
    bamfile.close()

    # Begin exitron calling
    global results
    results = {}

    def collect_result(output):
        exitrons = output[0]
        chrm = output[1]
        pretty_print(f'Collecting data from {chrm}')
        sys.stdout.flush()
        results[chrm] = exitrons

    if cores:
        pool = mp.Pool(int(cores))
        threads = []
        for chrm in chrms:
            threads.append(pool.apply_async(exitrons_in_chrm, args=(bamfilename,
                                                                    reference_transcriptome,
                                                                    genome,
                                                                    chrm,
                                                                    cell_types,
                                                                    mapq,
                                                                    alignment50,
                                                                    pso,
                                                                    unique_exitron_mol,
                                                                    ), callback=collect_result))
        pool.close()
        for t in threads:
            t.get()  # this gets any exceptions raised
        pool.join()
    else:
        for chrm in chrms:
            sys.stdout.flush()
            output = exitrons_in_chrm(bamfilename,
                                      reference_transcriptome,
                                      genome,
                                      chrm,
                                      cell_types,
                                      mapq,
                                      alignment50,
                                      pso,
                                      unique_exitron_mol,
                                      )
            collect_result(output)


    prefix_out = os.path.splitext(os.path.basename(
        bamfile.filename.decode('UTF-8')))[0]
    out_file_name = f'{out}/{prefix_out}.exitron'
    pretty_print(
        f'Finished exitron calling and filtering. Printing to {out_file_name}')
    sys.stdout.flush()
    with open(out_file_name, 'w') as out:
        header = ['chrom',
                  'start',
                  'end',
                  'strand',
                  'name',
                  'gene_symbol',
                  'region',
                  'splice_site',
                  'length',
                  'alignment50',
                  'rt_repeat',
                  'read_total',
                  'cell_type',
                  'exitron_mols',
                  'unique_mols',
                  'exitron_cells_not_found',
                  'pso',
                  'pairwise_z_test_pvals',
                  'meta_pval'
                  ]
        if id_:
            header += ['id']
        # write header
        for column in header:
            out.write(column + '\t')
        out.write('\n')
        out_res = []
        for chrm in chrms:
            # check if chromosome is empty or not
            try:
                if results[chrm]:
                    for exitron in results[chrm]:
                        if id_:
                            exitron['id'] = id_
                        out_res.append(exitron)
                        out.write('\t'.join([str(exitron[column])
                                  for column in header]) + '\n')

            except KeyError:
                pretty_print(
                    f'Thread most likely crashed on chromosome \'{chrm}\' without reporting exception.  Try fewer cores or allocate more memory.')
                sys.stdout.flush()
                sys.exit(1)

        return out_res


def main(input_,
           cell_types_fn,
           genome,
           reference_transcriptome,
           out,
           cores,
           mapq,
           unique_exitron_mol,
           pso,
           alignment50,
           id_):
    #=============================================================================
    # check for requirements and preprocess
    #=============================================================================
    # Check if annotation has a tabix index
    try:
        try:
            gtf = pysam.TabixFile(reference_transcriptome, parser=pysam.asGTF())
            gtf.close()
        except OSError:
            pretty_print('Building tabix index.')
            pysam.tabix_index(reference_transcriptome, preset='gff')
    except:
        pretty_print(
            f'ERROR: There is a problem reading the annotation file at: {reference_transcriptome}')
        pretty_print(f'Please make sure to use bgzip to compress your annotation file.')
        sys.exit(1)

    # Check for gziped annotation
    if reference_transcriptome[-2:] != 'gz':
        pretty_print('ERROR: Annotation file is required to be zipped in .gz format. Please compress your GTF file with a command such as: gzip -c in.gtf > out.gtf.gz')
        sys.exit(1)

    # Prepage gffutils database
    try:
        db = gffutils.FeatureDB(reference_transcriptome + '.db')
        pretty_print(f'Using annotation databse {reference_transcriptome + ".db"}')
    except ValueError:
        pretty_print('Preparing annotation database... This may take awhile ... ')
        db = gffutils.create_db(reference_transcriptome,
                                dbfn=reference_transcriptome + '.db',
                                disable_infer_transcripts=True,
                                disable_infer_genes=True)
        db = gffutils.FeatureDB(reference_transcriptome + '.db')
        pretty_print(f'Using annotation databse {reference_transcriptome + ".db"}')

    del db

    # Define chrms
    chrms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5',
              'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
              'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
              'chr21', 'chr22', 'chrX', 'chrY']



    #=============================================================================
    # For each bamfile in input_, find exitrons
    #=============================================================================
    bamlist = pd.read_csv(input_, delimiter = '\t') #TODO use input_
    n_bams = bamlist.shape[0]
    # let exitrons_all_samples be a dictionory where each key is an exitron
    exitrons_all_samples = defaultdict(list)
    try:
        os.mkdir(out)
    except FileExistsError:
        pass
    pretty_print(f'Processing {n_bams} samples')
    for i in range(n_bams):
        bam_name = bamlist.iloc[i, 0]
        bam = os.path.join(os.path.dirname(input_) if os.path.dirname(input_) else '.', bam_name)
        sample_id = bamlist.iloc[i, 1]
        group = bamlist.iloc[i, 2]

        # process cell type data
        cell_types = pd.read_csv(cell_types_fn, sep = '\t')
        cell_types = cell_types[cell_types['sample_id'] == sample_id]

        called_exitrons = call_exitrons(bam, chrms, sample_id, input_,
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
        for exitron in called_exitrons:
            exitron['bam_name'] = bam_name
            exitrons_all_samples[exitron['name']].append(exitron)


    #=============================================================================
    # For each exitron found in at least one sample, calculate exitron_mols and unique_mol for
    #   each sample (those who have the exitron have already been calculated)
    #=============================================================================
    pretty_print(f'Preping splicing summary between cell types.')

    # get cell type names
    cell_type_names = [x for x in cell_types.cell_group.unique() if not pd.isna(x)]

    # exitrons_mols_by_group is a dict with keys "group".
    # values are dicts with keys  "exitrons name"
    exitron_mols_by_group = defaultdict(dict)
    for group in bamlist.group.unique():
        # initialize group names
        exitron_mols_by_group[group] = defaultdict(list)
    res = []
    for exitron_name in exitrons_all_samples.keys():
        # for each exitron, find molecules mapping to locus in each sample (bamlist rows)
        for _, row in bamlist.iterrows():
            bam_name = row['bam']
            group = row['group']
            called_exitrons = [e for e in exitrons_all_samples[exitron_name] \
                               if e['bam_name'] == bam_name]
            if called_exitrons:
                # get exitron_mol and unique_mol from called exitron
                for exitron in called_exitrons:
                    entry = (exitron['cell_type'], exitron['exitron_mols'], exitron['unique_mols'])
                    exitron_mols_by_group[group][exitron_name].append(entry)
                    bam = os.path.dirname(input_) + '/' + bam_name
                    res.append([exitron_name, bam, group, entry[0], entry[1], entry[2]])

            else:
                chrm = exitrons_all_samples[exitron_name][0]['chrom']
                start = exitrons_all_samples[exitron_name][0]['start']
                end = exitrons_all_samples[exitron_name][0]['end']
                bam = os.path.join(os.path.dirname(input_) if os.path.dirname(input_) else '.', bam_name)
                bamfile = pysam.AlignmentFile(bam, 'rb')
                # for each cell type, get unique mols.
                for cell_type in cell_type_names:
                    cells_in_type = cell_types[cell_types['cell_group'] == cell_type].cells.unique()
                    unique_mols = get_unique_mols(chrm,
                                                  start,
                                                  end,
                                                  bamfile,
                                                  set(cells_in_type),
                                                  mapq)
                    entry = (cell_type, 0, unique_mols)
                    exitron_mols_by_group[group][exitron_name].append(entry)
                    res.append([exitron_name, bam, group, entry[0], entry[1], entry[2]])



    #=============================================================================
    # write to file for differential expression analysis using R script
    #=============================================================================

    # pretty_print(res)
    out_file_name = f'{out}/exitron_splicing_summary.tsv'
    pretty_print(f'Done. Writing to file {out_file_name}')

    with open(out_file_name, 'w') as out:
        header = ['exitron',
                  'bam',
                  'group',
                  'cell_type',
                  'spliced_cov',
                  'unspliced_cov']
        out.write('\t'.join(header) + '\n')
        out.write('\n'.join('\t'.join(map(str, x)) for x in res))



if __name__ == '__main__':
    main('/Users/jpfry/Documents/Projects/sesc/test_data/bam_list.txt',
               '/Users/jpfry/Documents/Projects/sesc/test_data/test_cell_types.tsv',
               '/Users/jpfry/Documents/Projects/scanexitron/hg38.fa',
               '/Users/jpfry/Documents/Projects/scanexitron/gencode.v37.annotation.sorted.gtf.gz',
               '/Users/jpfry/Documents/Projects/sesc/test_data/testing',
               1,
               50,
               2,
               0.01,
               0.7,
               'test')
