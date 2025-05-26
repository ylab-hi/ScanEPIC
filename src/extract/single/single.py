# -*- coding: utf-8 -*-
"""
Single-cell exitron extraction module for ScanEPIC

This module identifies and quantifies exitron splicing events from single-cell RNA-seq data.
It processes cell barcodes and UMIs to provide cell-type-specific exitron quantification.

@author: Josh Fry, Northwestern University, YangLab
@original_author: Tingyou Wang, Northwestern University, YangLab (v1)
"""
# ===============================================================================

__version__ = 'v1'
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



# from _cython_fi import find_introns


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


# def find_introns(read_iterator):
#     """Return a dictionary {(start, stop): count}
#     Listing the intronic sites in the reads (identified by 'N' in the cigar strings),
#     and their support ( = number of reads ).
#     read_iterator can be the result of a .fetch(...) call.
#     Or it can be a generator filtering such reads. Example
#     samfile.find_introns((read for read in samfile.fetch(...) if read.is_reverse)
#     """
#     # cdef:
#     #     int base_position, junc_start, nt, read_position, i
#     #     int op
#     #     AlignedSegment r
#     #     int BAM_CREF_SKIP = 3 #BAM_CREF_SKIP

#     BAM_CREF_SKIP = 3
#     res = collections.Counter()
#     reads = collections.defaultdict(list)
#     match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
#     for r in read_iterator:
#         base_position = r.pos
#         read_position = 0

#         for i, (op, nt) in enumerate(r.cigartuples):
#             if op in match_or_deletion:
#                 base_position += nt
#                 read_position += nt
#             elif op == BAM_CREF_SKIP:
#                 junc_start = base_position
#                 base_position += nt
#                 try:
#                     reads[(junc_start, base_position)].append((r.seq, '.', r.cigartuples[i-1][1], r.cigartuples[i+1][1], read_position,
#                                                                 r.get_tag('CB'), r.get_tag('UB')))
#                     res[(junc_start, base_position)] += 1
#                 except:
#                     continue

#     return res, reads

# def find_introns(read_iterator, stranded):
#     """

#     Parameters
#     ----------
#     read_iterator : iterator of reads from a pysam.AlignmentFile.
#         Expected that the iterator will be an entire chromosome. See exitron_caller


#     Returns
#     -------
#     introns -- counter of (intron_start, intron_end, strand)
#     reads -- dictionary of reads that support the junctions in introns.  This
#         will be used in later filtering steps.
#     meta_data -- dictionary consisting of metadata collected along the way.
#         This is used in later steps.

#     """
#     BAM_CREF_SKIP = 3

#     introns = Counter()
#     reads = defaultdict(list)

#     # only M/=/X (0/7/8) and D (2) are related to genome position
#     match_or_deletion = {0, 2, 7, 8}
#     for r in read_iterator:
#         base_position = r.pos
#         read_position = 0
#         # if cigarstring is * (r.cigartuples == None), unmatched, continue
#         if r.cigartuples == None:
#             continue
#         # iterate through cigar string looking for N
#         for i, (tag, nt) in enumerate(r.cigartuples):
#             # if (0, X), keep track of base_position.
#             # if (3, X), which corresponds to N,
#             # look at match before and after
#             if tag in match_or_deletion:
#                 base_position += nt
#                 read_position += nt
#             elif r.cigartuples[i][0] == BAM_CREF_SKIP:
#                 junc_start = base_position
#                 base_position += nt
#                 if stranded == 'no':
#                     try:
#                         introns[(junc_start, base_position,
#                                   r.get_tag('XS'))] += 1
#                         reads[(junc_start, base_position, r.get_tag('XS'))].append(
#                             (r.seq, '.', r.cigartuples[i-1][1], r.cigartuples[i+1][1], read_position))
#                     except KeyError:  # this ignores junctions without XS tags, usually because they are non-canonical
#                         pass

#                 else:
#                     if stranded == 'fr-firststrand':
#                         strand = '+' if (r.is_read2 and not r.is_reverse) or \
#                                         (r.is_read1 and r.is_reverse) else '-'
#                     elif stranded == 'fr-secondstrand':
#                         strand = '+' if (r.is_read1 and not r.is_reverse) or \
#                                         (r.is_read2 and r.is_reverse) else '-'
#                     introns[(junc_start, base_position, strand)] += 1
#                     reads[(junc_start, base_position, strand)].append((r.seq, r.get_reference_sequence(
#                     ), r.cigartuples[i-1][1], r.cigartuples[i+1][1], read_position))

#     return introns, reads


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
            # try:
            #     # some entries like 'gene' do not have transcript types
            #     # in this case, just continue
            #     gene_type = feature.gene_type
            # except:
            #     continue
            # if gene_type != 'protein_coding': continue
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






                    # if

                    #     # get transcript from which this UTR originates
                    #     transcript = data.transcript_by_id(
                    #         feature.attrs['transcript_id'])
                    #     exon_ranges = transcript.exon_intervals
                    #     for exon in exon_ranges:  # find the exon for which this UTR originates
                    #         if feature.start + 1 >= exon[0] and feature.end <= exon[1]:
                    #             UTR_exon = exon
                    #             break
                    #     # require that the exitron is an intron within the exon that contains the UTR
                    #     if UTR_exon[0] < intron_start + 1 and UTR_exon[1] > intron_end - 1:
                    #         # determine 5' UTR or 3' UTR
                    #         middle_exon_start = median(
                    #             exon[0] for exon in exon_ranges)
                    #         if feature.strand == '+':
                    #             region = '5\'UTR' if UTR_exon[0] < middle_exon_start else '3\'UTR'
                    #         elif feature.strand == '-':
                    #             region = '3\'UTR' if UTR_exon[0] < middle_exon_start else '5\'UTR'

                    #         # determine whether or not the intron overlaps a coding region.
                    #         coding_ranges = transcript.coding_sequence_position_ranges
                    #         for c_region in coding_ranges:
                    #             if (intron_start + 1 < c_region[0] and intron_end - 2 + 1 > c_region[0]) or \
                    #                     (intron_start + 1 < c_region[1] and intron_end - 2 + 1 > c_region[1]):
                    #                 region += ' + CDS'
                    #                 break
                    #         exitrons.append({'chrom': feature.chrom,
                    #                         'start': intron_start + 1,
                    #                          'end': intron_end - 2 + 1,  # plus 1 because of bedtools conventions,
                    #                          'name': f'{gene_name}{intron_start+1}{intron_end - 2 + 1}',
                    #                          'region': region,
                    #                          'ao': intron_witnesses,
                    #                          'strand': feature.strand,
                    #                          'gene_symbol': gene_name,
                    #                          'length': intron_end - 2 + 1 - (intron_start + 1) - 1,
                    #                          # this is calculated post filtering to save genome lookups
                    #                          'splice_site': 'splice_site',
                    #                          'transcript_id': feature.attrs['transcript_id']})
                    #         exitrons_added.append(
                    #             (intron_start, intron_end, region_type))

    # for feature in intersection:
    #     # Check for intron within coding exon.
    #     region_type = feature.fields[2]
    #     region_start = feature.start
    #     region_end = feature.end
    #     gene_name = feature.attrs['gene_name']

    #     intron_start = int(feature.fields[10])
    #     intron_end = int(feature.fields[11])
    #     intron_witnesses = int(feature.fields[12])

    #     # Use the ends to check for known donors or acceptors
    #     if (region_type == 'exon' or region_type == 'utr') and feature.attrs['transcript_type'] != 'nonsense_mediated_decay':
    #         if intron_start + 1 == region_end:
    #             # intron matches a known donor
    #             known_splices.add(
    #                 (feature.chrom, intron_start+1, intron_end - 2 + 1, 'D'))
    #         if intron_end - 2 == region_start:
    #             # intron matches a known acceptor
    #             known_splices.add(
    #                 (feature.chrom, intron_start+1, intron_end - 2 + 1, 'A'))

    #     # previous logic for reference:
    #     # if (intron_start + 1 == region_end and intron_end - 2 == region_start) and \
    #     #     (region_type == 'exon' or region_type == 'utr') and feature.attrs['transcript_type'] == 'protein_coding':
    #     #     # Intron matches a known acceptor or known donor (or both)
    #     #     # Add to list of known splices.  Some introns may be an exitron for
    #     #     # some annotated gene but a known splice for others.
    #     #     known_splices.add((feature.chrom, intron_start+1, intron_end - 2 + 1))

    #     elif region_type == 'CDS' and region_start < intron_start + 1 \
    #             and region_end > intron_end - 1:
    #         if (intron_start, intron_end, region_type) not in exitrons_added:
    #             exitrons.append({'chrom': feature.chrom,
    #                             'start': intron_start + 1,
    #                              'end': intron_end - 2 + 1,  # plus 1 because of bedtools conventions,
    #                              'name': f'{gene_name}{intron_start+1}{intron_end - 2 + 1}',
    #                              'region': region_type,
    #                              'ao': intron_witnesses,
    #                              'strand': feature.strand,
    #                              'gene_symbol': gene_name,
    #                              'length': intron_end - 2 + 1 - (intron_start + 1) - 1,
    #                              'splice_site': 'splice_site',
    #                              'transcript_id': feature.attrs['transcript_id']})
    #             exitrons_added.append((intron_start, intron_end, region_type))

    #     elif region_type == 'UTR':
    #         if feature.attrs['transcript_type'] == 'protein_coding':
    #             if (intron_start, intron_end, region_type) not in exitrons_added:
    #                 # get transcript from which this UTR originates
    #                 transcript = data.transcript_by_id(
    #                     feature.attrs['transcript_id'])
    #                 exon_ranges = transcript.exon_intervals
    #                 for exon in exon_ranges:  # find the exon for which this UTR originates
    #                     if feature.start + 1 >= exon[0] and feature.end <= exon[1]:
    #                         UTR_exon = exon
    #                         break
    #                 # require that the exitron is an intron within the exon that contains the UTR
    #                 if UTR_exon[0] < intron_start + 1 and UTR_exon[1] > intron_end - 1:
    #                     # determine 5' UTR or 3' UTR
    #                     middle_exon_start = median(
    #                         exon[0] for exon in exon_ranges)
    #                     if feature.strand == '+':
    #                         region = '5\'UTR' if UTR_exon[0] < middle_exon_start else '3\'UTR'
    #                     elif feature.strand == '-':
    #                         region = '3\'UTR' if UTR_exon[0] < middle_exon_start else '5\'UTR'

    #                     # determine whether or not the intron overlaps a coding region.
    #                     coding_ranges = transcript.coding_sequence_position_ranges
    #                     for c_region in coding_ranges:
    #                         if (intron_start + 1 < c_region[0] and intron_end - 2 + 1 > c_region[0]) or \
    #                                 (intron_start + 1 < c_region[1] and intron_end - 2 + 1 > c_region[1]):
    #                             region += ' + CDS'
    #                             break
    #                     exitrons.append({'chrom': feature.chrom,
    #                                     'start': intron_start + 1,
    #                                      'end': intron_end - 2 + 1,  # plus 1 because of bedtools conventions,
    #                                      'name': f'{gene_name}{intron_start+1}{intron_end - 2 + 1}',
    #                                      'region': region,
    #                                      'ao': intron_witnesses,
    #                                      'strand': feature.strand,
    #                                      'gene_symbol': gene_name,
    #                                      'length': intron_end - 2 + 1 - (intron_start + 1) - 1,
    #                                      # this is calculated post filtering to save genome lookups
    #                                      'splice_site': 'splice_site',
    #                                      'transcript_id': feature.attrs['transcript_id']})
    #                     exitrons_added.append(
    #                         (intron_start, intron_end, region_type))
    # del intersection
    # data.clear_cache()
    # del data
    return ([exitron for exitron in exitrons if ((exitron['chrom'], exitron['start'], exitron['end'], 'D', exitron['strand']) not in known_splices
                                                  and (exitron['chrom'], exitron['start'], exitron['end'], 'A', exitron['strand']) not in known_splices)],
            reads,
            )

    # return [], []


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
    # def get_unique_mol(chrm, start, stop, ct):
    #     umi = set()
    #     region_reads = [read for read in bamfile.fetch(chrm, start - 1, stop + 1)]
    #     for read in region_reads:
    #         try:
    #             if read.get_tag('CB') == cb:
    #                 umi.add(read.get_tag('UB') )
    #         except:
    #             continue
    #     return len(umi)

    # def get_unique_mols(chrm, pos, exitron_cell_type, cell_types, mapq):
    #     umi = set()
    #     region_reads = [read for read in bamfile.fetch(chrm, pos - 1, pos) if (read.mapq > mapq and pos - 1 in read.get_reference_positions())]
    #     for read in region_reads:
    #         try:
    #             # get cell type
    #             cb = read.get_tag('CB')
    #             cell_type = exitron_cell_types[exitron_cell_types['cells'].str.contains(cb)].cell_group.values[0]
    #             if cell_type == exitron_cell_type:
    #                 umi.add(read.get_tag('UB'))
    #         except:
    #             pass
    #     return len(umi)

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

            # mid = (exitron['start'] + exitron['end'])//2
            # a = len(set(read.get_tag('UB')  \
            #                    for read in bamfile.fetch(chrm, exitron['start'] - 1, exitron['start']) if \
            #                    read.has_tag('CB') and \
            #                    read.has_tag('UB') and \
            #                    (read.get_tag('CB') in cells_in_type) and read.mapping_quality >= mapq))
            # b = len(set(read.get_tag('UB')  \
            #                    for read in bamfile.fetch(chrm, exitron['end'] - 1, exitron['end']) if \
            #                    read.has_tag('CB') and \
            #                    read.has_tag('UB') and \
            #                    (read.get_tag('CB') in cells_in_type) and read.mapping_quality >= mapq))
            # c = len(set(read.get_tag('UB')  \
            #                    for read in bamfile.fetch(chrm, mid- 1, mid) if \
            #                    read.has_tag('CB') and \
            #                    read.has_tag('UB') and \
            #                    (read.get_tag('CB') in cells_in_type) and read.mapping_quality >= mapq))

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
            # print(f'{exitron_cell_type["unique_mols"]} vs {((uniqe_mols2 - unique_exitron_mol*exitron_cell_type["length"])//exitron_cell_type["length"] + unique_exitron_mol)})')
            #append exitron to penultimate list

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


            #     print(pval)
            # print('-----------------')

            # # pvals = ','.join(str(round(proportions_ztest(count = (e1['exitron_mols'], e2['exitron_mols']),
            # #                            nobs = (e1['unique_mols'], e2['unique_mols']))[1], 6)) \
            # #          for e2 in exitrons_to_test)
            # e1['z_test_pvals'] = pvals
            # res.append(e1)


        # unique_cells = set([r[5] for r in read_data])

        # # TODO: add warning if there is no annotation for a cell
        # # This determines the cell types present for this exitron
        # exitron_cell_types = cell_types[cell_types['cells'].isin(unique_cells)]

        # # For each cell type, calculate PSO, AO, etc, relative to all reads
        # #   within this cell type.
        # for _, cell_type in exitron_cell_types.iterrows():
        #     # copy because we are constructing a cell_type exitron
        #     #   perhaps this would be cleaner to do with classes and inheritance...
        #     exitron_cell = exitron.copy()
        #     exitron_cell['cell_type'] = cell_type.cell_group

        #     # calculate AO by finding those cell type specific reads supporting the exitron
        #     ao = 0
        #     unique_mol = set()
        #     for r in read_data:
        #         try:
        #             r_cell_type = exitron_cell_types[exitron_cell_types['cells'].str.contains(r[5])].cell_group.values[0]
        #             if r_cell_type == exitron_cell['cell_type']:
        #                 ao += 1
        #                 # r[6] is UMI of read
        #                 unique_mol.add(r[6])
        #         except IndexError:
        #             # there are no annotated cells with this barcode.
        #             # they may have been filtered when assigning barcodes to cell types
        #             continue
        #     exitron_cell['ao'] = len(unique_mol)
        #     a = get_unique_mols(chrm, exitron['start'], exitron_cell['cell_type'], cell_types, mapq)
        #     b = get_unique_mols(chrm, exitron['end'], exitron_cell['cell_type'], cell_types, mapq)
        #     c = get_unique_mols(chrm, (exitron['start']+exitron['end'])//2, exitron_cell['cell_type'], cell_types, mapq)
        #     exitron_cell['unique_mol'] = max(a, b, c)
        # # a = bamfile.count(chrm, start=start - 1, stop=start,
        # #                   read_callback=lambda x: x.mapq > mapq)
        # # b = bamfile.count(chrm, start=end - 1, stop=end,
        # #                   read_callback=lambda x: x.mapq > mapq)
        # # c = bamfile.count(chrm, start=mid - 1, stop=mid,
        # #                   read_callback=lambda x: x.mapq > mapq)





        # old way
        # read_data = reads[(start, end - 1)]
        # unique_cells = set([r[5] for r in read_data])
        # for cell in unique_cells:
        #     exitron_cell = exitron.copy()
        #     ao = 0
        #     for r in read_data:
        #         ao += r[5] == cell
        #     unique_mol = get_unique_mol(chrm, start, end, cell)
        #     exitron_cell['unique_mol'] = unique_mol
        #     exitron_cell['exitron_mol'] = len(set([r[6] for r in read_data if r[5] == cell]))
        #     exitron_cell['ao'] = ao
        #     exitron_cell['CB'] = cell
        #     exitron_cell['splice_site'] = splice_site
        #     exitron_cell['pso'] = exitron_cell['exitron_mol']/exitron_cell['unique_mol']
        #     if exitron_cell['pso'] >= pso_min and ao >= ao_min:
        #         res.append(exitron_cell)

    # # filter one exitron at a time
    # for exitron in processed_exitrons:
    #     ao = exitron['ao']
    #     if ao < ao_min: continue
    #     chrm = exitron['chrom']
    #     start = exitron['start']
    #     end = exitron['end']
    #     strand = exitron['strand']


    #     exitron['splice_site'] = splice_site
    #     exitron['pso'] = exitron['exitron_mol']/exitron['unique_mol']

    #     if exitron['pso'] > pso_min:
    #         res.append(exitron)

        # right_anchor = (0, '')
        # left_anchor = (0, '')
        # ao_true = 0
        # for read in read_data:

        #     l_anchor_len = read[2]
        #     r_anchor_len = read[3]

        #     # check 5' intron and 3' exon similarity
        #     three_prime_exon_r = read[0][read[4]:read[4] + r_anchor_len]
        #     five_prime_intron_r = intron_seq[:len(three_prime_exon_r)]
        #     # check 5' exon and 3' intron similarity
        #     five_prime_exon_r = read[0][read[4] - l_anchor_len:read[4]]
        #     three_prime_intron_r = intron_seq[-len(five_prime_exon_r):]

        #     if three_prime_exon_r != five_prime_intron_r and five_prime_exon_r != three_prime_intron_r:
        #         ao_true += 1

        # if ao_true == 0 and ao_min > 0:
        #     exitron['ao_unfiltered'] = ao

        # # We subtract 1 because these coords are in BED format.
        # mid = (start+end)/2
        # a = bamfile.count(chrm, start=start - 1, stop=start,
        #                   read_callback=lambda x: x.mapq > mapq)
        # b = bamfile.count(chrm, start=end - 1, stop=end,
        #                   read_callback=lambda x: x.mapq > mapq)
        # c = bamfile.count(chrm, start=mid - 1, stop=mid,
        #                   read_callback=lambda x: x.mapq > mapq)

        # pso = ao/((a + b + c - ao*3)/3.0 + ao)
        # psi = 1 - pso
        # dp = int(ao/pso) if pso > 0 else 0

        # pso_true = ao_true/((a + b + c - ao_true*3)/3.0 + ao_true)
        # psi_true = 1 - pso_true
        # dp_true = int(ao_true/pso_true) if pso_true > 0 else 0

        # # Check whether attributes exceed minimum values
        # if (pso >= pso_ufmin and
        #     pso_true >= pso_min and
        #         ao_true >= ao_min):

        #     exitron['ao'] = ao_true
        #     exitron['pso'] = pso_true
        #     exitron['psi'] = psi_true
        #     exitron['dp'] = dp_true

        #     exitron['ao_unfiltered'] = ao
        #     exitron['pso_unfiltered'] = pso
        #     exitron['psi_unfiltered'] = psi
        #     exitron['dp_unfiltered'] = dp

        #     exitron['delta_ao'] = ao - ao_true

        #     left_anchor = max(read_data, key=lambda x: x[2])
        #     right_anchor = max(read_data, key=lambda x: x[3])
        #     if not (left_anchor[2] < anchor_min or right_anchor[3] < anchor_min):
        #         res.append(exitron)
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

    # chrms = set()
    # query = db.execute("select seqid from features")
    # for x in query:
    #     chrms.add(x['seqid'])
    # chrms = list(chrms)
    # chrms.sort()
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
    # pretty_print('here')
    # # exitron_mols_by_group is a dict with keys group, whose values are dictionaries with keys exitron
    # exitron_mols_by_group = defaultdict(dict)
    # for group in bamlist.group.unique():
    #     exitron_mols_by_group[group] = defaultdict(dict)
    # for exitron_name in exitrons_all_samples:
    #     # get exitron_mol and unique_mol for each exitron
    #     for row in bamlist.iterrows():
    #         bam_name = row['bam']
    #         group = row['group']
    #         called_exitrons = [e for e in exitrons_all_samples[exitron_name] if e['bam_name'] == bam_name]
    #         if called_exitrons:
    #             # get exitron_mol and unique_mol from called exitron
    #             for exitron in called_exitrons:
    #                 exitron_mols_by_group[group][]

    #                 exitron_mols[exitron_name][exitron['cell_type']][group] = (exitron['exitron_mols'],
    #                                                                            exitron['unique_mols'])

    #         else:
    #             # get exitron_mol and unique_mol from bamfile
    #             bam = os.path.dirname(input_) + '/' + bam_name
    #             bamfile = pysam.AlignmentFile(bamfilename, 'rb', require_index=True)

    # get cell type names
    cell_type_names = [x for x in cell_types.cell_group.unique() if not pd.isna(x)]

    # exitrons_mols_by_group is a dict with keys "group".
    # values are dicts with keys  "exitrons name"
    exitron_mols_by_group = defaultdict(dict)
    for group in bamlist.group.unique():
        # initialize group names
        exitron_mols_by_group[group] = defaultdict(list)
    # for group in bamlist.group.unique():
    #     # initialize group names
    #     exitron_mols_by_group[group] = defaultdict(list)
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

    # #=============================================================================
    # # aggregate counts between groups
    # #=============================================================================


    # # stats_res is a dictionary with keys groups and values a dictionary with keys exitron_name
    # stats_res = defaultdict(lambda: defaultdict(dict))
    # for group in exitron_mols_by_group.keys():
    #     for exitron_name in exitron_mols_by_group[group]:
    #         cell_types = set(obs[0] for obs in exitron_mols_by_group[group][exitron_name])
    #         for cell_type in cell_types:
    #             aggregate_counts_e = sum(obs[1] for obs in exitron_mols_by_group[group][exitron_name] if obs[0] == cell_type)
    #             aggregate_counts_u = sum(obs[2] for obs in exitron_mols_by_group[group][exitron_name] if obs[0] == cell_type)
    #             stats_res[group][exitron_name][cell_type] = (aggregate_counts_e, aggregate_counts_u)

    # #=============================================================================
    # # proportionz test
    # #=============================================================================
    # pvals_res = defaultdict(lambda: defaultdict(dict))
    # for group in stats_res.keys():
    #     for exitron_name in stats_res[group]:
    #         for cell_type in stats_res[group][exitron_name].keys():
    #             pvals = []
    #             x1 = stats_res[group][exitron_name][cell_type]
    #             for comparison_type in stats_res[group][exitron_name].keys():
    #                 x2 = stats_res[group][exitron_name][comparison_type]
    #                 cell_type_count = (x1[0], x2[0])
    #                 nobs = (x1[1] + 1e-5, x2[1] + 1e-5)
    #                 pval = proportions_ztest(count = cell_type_count, nobs = nobs)[1]
    #                 if pd.isna(pval):
    #                     pval = 1.0
    #                 pvals.append(pval)
    #             # append pvals and aggregate PSO for later use
    #             pvals_res[group][exitron_name][cell_type] = (pvals, x1)

    # exitron_statistics = []
    # # fisher method for summarizing p vals
    # meta_pvals = defaultdict(lambda: defaultdict(dict))
    # for group in pvals_res.keys():
    #     for exitron_name in pvals_res[group].keys():
    #         for cell_type in pvals_res[group][exitron_name].keys():
    #             pvals, x1 = pvals_res[group][exitron_name][cell_type]
    #             # finally, construct results and apped to exitron_statistics
    #             res = {}
    #             res['name'] = exitron_name
    #             res['group'] = group
    #             res['chrom'] = exitrons_all_samples[exitron_name][0]['chrom']
    #             res['start'] = exitrons_all_samples[exitron_name][0]['start']
    #             res['end'] = exitrons_all_samples[exitron_name][0]['end']
    #             res['region'] = exitrons_all_samples[exitron_name][0]['region']
    #             res['gene_symbol'] = exitrons_all_samples[exitron_name][0]['gene_symbol']
    #             res['strand'] = exitrons_all_samples[exitron_name][0]['strand']
    #             res['length'] = exitrons_all_samples[exitron_name][0]['length']
    #             res['alignment50'] = exitrons_all_samples[exitron_name][0]['alignment50']
    #             res['cell_type'] = cell_type
    #             res['exitron_mols'] = x1[0]
    #             res['unique_mols'] = x1[1]
    #             res['pso'] = x1[0]/(x1[1] + 1e-5)
    #             res['pvals'] = ','.join(map(lambda x: str(round(x, 5)), pvals))
    #             # because one pval is identity comparison, which is kept in the result in order
    #             #   to keep things readable. But we need to remove it from meta pval
    #             pvals.remove(max(pvals))
    #             res['meta_pval'] = round(combine_pvalues(pvals, method = 'fisher')[1], 5)
    #             exitron_statistics.append(res)

    # #=============================================================================
    # # write statistics to file
    # #=============================================================================
    # out_file_name = f'{out}/exitron_differential_splicing.tsv'
    # pretty_print(f'Done. Writing to file {out_file_name}')

    # with open(out_file_name, 'w') as out:
    #     header = ['chrom',
    #               'start',
    #               'end',
    #               'name',
    #               'gene_symbol',
    #               'region',
    #               'length',
    #               'alignment50',
    #               'cell_type',
    #               'exitron_mols',
    #               'unique_mols',
    #               'pso',
    #               'pvals',
    #               'meta_pval']
    #     out.write('\t'.join(header) + '\n')
    #     for e in exitron_statistics:
    #         out.write('\t'.join([str(e[column])
    #                   for column in header]) + '\n')


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
