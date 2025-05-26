# -*- coding: utf-8 -*-
"""
Long-read exitron extraction module for ScanEPIC

This module identifies and quantifies exitron splicing events from long-read RNA-seq data.
It includes clustering, realignment, and transcript quantification capabilities.

@author: Josh Fry, Northwestern University, YangLab
"""
__version__ = 'v2.0.0'

import argparse
import multiprocessing as mp
import os
import re
import subprocess
import sys
import traceback
import gffutils
import pandas as pd
import pysam
import click
from collections import Counter
from collections import defaultdict
from itertools import islice
from time import strftime, localtime
from Bio import pairwise2, Align

# ===============================================================================
# Helper Methods
# ===============================================================================

def pretty_print(text):
    click.echo(click.style(f'[exitrontools extract short -- {strftime("%Y-%m-%d %I:%M:%S %p", localtime())}]',
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

def find_introns(read_iterator):
    """
    Identify introns from long-read RNA-seq data by parsing CIGAR strings.
    
    This function processes reads to identify splice junctions (N operations in CIGAR)
    and returns both the intron coordinates and supporting read information.

    Parameters
    ----------
    read_iterator : iterator
        Iterator of reads from a pysam.AlignmentFile for a specific chromosome

    Returns
    -------
    introns : Counter
        Counter of (intron_start, intron_end, strand) tuples
    reads : defaultdict
        Dictionary mapping intron coordinates to supporting read names
    """
    BAM_CREF_SKIP = 3

    introns = Counter()
    reads = defaultdict(list)

    # only M/=/X (0/7/8) and D (2) are related to genome position
    match = {0, 7, 8}
    for r in read_iterator:
        base_position = r.pos
        read_position = 0
        # if cigarstring is * (r.cigartuples == None), unmatched, continue
        if r.cigartuples == None:
            continue
        # iterate through cigar string looking for N
        for i, (tag, nt) in enumerate(r.cigartuples):
            # if (0, X), keep track of base_position.
            # if (3, X), which corresponds to N,
            # look at match before and after
            if tag in match or (tag == 2 and nt < 30):
                base_position += nt
                read_position += nt
            elif tag == BAM_CREF_SKIP or (tag == 2 and nt >= 30):
                junc_start = base_position
                base_position += nt
                junc_end = base_position
                try:
                    if r.get_tag('ts') == '+':
                        strand = '-' if r.is_reverse else '+'
                    elif r.get_tag('ts') == '-':
                        strand = '+' if r.is_reverse else '-'
                except KeyError:
                    # ts tag is not present, thus we don't know strand, continue
                    continue
                if r.cigartuples[i - 1][0] == 2:
                    junc_start -= r.cigartuples[i - 1][1]
                if r.cigartuples[i + 1][0] == 2:
                    junc_end += r.cigartuples[i + 1][1]
                introns[(junc_start, junc_end, strand)] += 1
                reads[(junc_start, junc_end, strand)].append(r.query_name)

    return introns, reads


def exitron_caller(bamfile, referencename, chrm, db, mapq=50, jitter=10):
    """


    Parameters
    ----------
    bamfile : pysam.AlignmentFile
    referencename : str
    chrms : list
    db : gffulits.database
    mapq : int, optional
        Only considers reads from bamfile with quality >= mapq. The default is 50.
    jitter : int

    Returns
    -------
    List of unfiltered exitrons.  Each exitron is a dictionary of features.

    """

    introns, reads = find_introns(
        read for read in bamfile.fetch(chrm) if read.mapping_quality >= mapq
    )

    gtf = pysam.TabixFile(referencename, parser=pysam.asGTF())
    exitrons = []
    exitrons_utr = []
    exitrons_added = []
    known_splices = set()

    # load blacklist of known false-positive hot spots
    this_dir = os.path.dirname(os.path.realpath(__file__))
    blacklist_path = os.path.join(this_dir, 'blacklist.tsv')
    try:
        with open(blacklist_path) as b:
            b.readline()
            blacklist = [l.split('\t')[1].rstrip() for l in b]
    except:
        blacklist = []

    for intron in introns:
        intron_start = intron[0] - 1 - jitter
        intron_end = intron[1] + 1 + jitter + 1
        intron_witnesses = introns[intron]
        intersection = gtf.fetch(chrm, intron_start, intron_end)
        for feature in intersection:
            # Check for intron within coding exon.
            if feature.strand != intron[2]:
                continue  # gtf.fetch is strand agnostic
            region_type = feature.feature
            region_start = feature.start
            region_end = feature.end
            try:
                gene_name = feature.gene_name
                gene_id = feature.gene_id
            except:
                try:
                    # Some Arabidopsis GTF filies have only a gene_id
                    gene_name = feature.gene_id
                    gene_id = feature.gene_id
                except:
                    # This is to catch cases where exon is not associated with any gene
                    continue
            # Use the ends to check for known donors or acceptors
            if region_type == 'exon':
                transcript_id = feature.transcript_id
                if intron_start in range(region_end - jitter*2, region_end + 1) and \
                    (chrm, intron_start + 1 + jitter, intron_end - 2 - jitter + 1, 'D') not in known_splices and \
                        feature.end != db[transcript_id].end:
                    # intron matches a known donor
                    known_splices.add(
                        (chrm, intron_start + 1 + jitter, intron_end - 2 - jitter + 1, 'D'))
                if intron_end in range(region_start, region_start + 1 + jitter*2) and \
                    (chrm, intron_start + 1 + jitter, intron_end - 2 - jitter + 1, 'A') not in known_splices and \
                        feature.start + 1 != db[transcript_id].start:
                    # intron matches a known acceptor
                    known_splices.add(
                        (chrm, intron_start + 1 + jitter, intron_end - 2 - jitter + 1, 'A'))

            # CDS branch
            elif region_type == 'CDS' and region_start < intron_start + 1 + jitter \
                    and region_end > intron_end - 1 - jitter:
                transcript_id = feature.transcript_id
                if gene_id in blacklist:
                    continue
                if (intron_start, intron_end, region_type) not in exitrons_added:
                    exitrons.append({'chrom': chrm,
                                    'start': intron_start + 1 + jitter,
                                     'end': intron_end - 2 - jitter + 1,  # plus 1 because of bedtools conventions,
                                     'name': f'{gene_name}d{intron_start + 1 + jitter}-{intron_end - 2 - jitter + 1}',
                                     'region': region_type,
                                     'ao': intron_witnesses,
                                     'strand': feature.strand,
                                     'gene_symbol': gene_name,
                                     'gene_id': gene_id,
                                     'length': intron_end - 2 - jitter + 1 - (intron_start + 1 + jitter) - 1,
                                     'splice_site': 'splice_site',
                                     'annotation': 'unannotated',
                                     'transcript_id': transcript_id})
                    exitrons_added.append(
                        (intron_start, intron_end, region_type))
                elif 'tag' in feature.keys():
                    if ('basic' in feature.tag or
                        'CCDS' in feature.tag or
                            'ccdsid' in feature.keys()):
                        # TODO pybedtools only grants access to ONE tag, we need to look at all the tags
                        # submit an issue with pybedtools
                        for e in exitrons:
                            if e['name'] == f'{gene_name}d{intron_start + 1 + jitter}-{intron_end - 2 - jitter + 1}':
                                e['transcript_id'] = transcript_id
                                break

            # UTR branch
            elif region_type == 'UTR':
                if feature.transcript_type == 'protein_coding' and \
                    (intron_start, intron_end, region_type, feature.strand) not in exitrons_added:
                        # require that exitron is within exon
                        exon_intervals = [(x.start, x.end) for x in db.children(db[feature.transcript_id], featuretype = 'exon')]
                        if any(map(lambda x: x[0] < intron_start + 1 + jitter and intron_end - 1 - jitter < x[1], exon_intervals)):

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
                                            'start': intron_start + 1 + jitter,
                                             'end': intron_end - 2 - jitter + 1,  # plus 1 because of bedtools conventions,
                                             'name': f'{gene_name}d{intron_start + 1 + jitter}-{intron_end - 2 - jitter + 1}',
                                             'region': region,
                                             'ao': intron_witnesses,
                                             'strand': feature.strand,
                                             'gene_symbol': gene_name,
                                             'gene_id': gene_id,
                                             'length': intron_end - 2 - jitter + 1 - (intron_start + 1 + jitter) - 1,
                                             'splice_site': 'splice_site',
                                             'transcript_id': feature.transcript_id})
                            exitrons_added.append((intron_start, intron_end, region_type, feature.strand))

    gtf.close()

    # Filter exitrons to remove those that match known splice sites
    res = [exitron for exitron in exitrons if ((exitron['chrom'], exitron['start'], exitron['end'], 'D') not in known_splices and
                                                 ((exitron['chrom'], exitron['start'], exitron['end'], 'A') not in known_splices))]
    # determine if UTR exitron is annotated
    for exitron in exitrons_utr:
        if ((exitron['chrom'], exitron['start'], exitron['end'], 'D') not in known_splices and
                                                     ((exitron['chrom'], exitron['start'], exitron['end'], 'A') not in known_splices)):
            res.extend(exitrons_utr)
        else:
            pass
    # res.extend(exitrons_utr)

    return res, reads


def filter_exitrons(exitrons, reads, bamfile, genome, db, skip_realign, mapq=50, pso_min=0.01, ao_min=2, cluster_purity=0, jitter=10):
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
    db : gffutils database
    skip_realign : bool
        If true, we skip realignment step
    mapq : TYPE, optional
        Only considers reads from bamfile with quality >= mapq. The default is 50.
        This is needed to calculate pso.
    pso_min : float, optional
        Number reads that witness the exitron over the number of reads within the
        spliced out region. The default is 0.05.
    ao_min : int, optional
        Minimum number of reads witnessing the exitron. The default is 3.
    jitter : int, optional
        treats intron borders as fuzzy with +/- jitter

    Returns
    -------
    List of filtered exitrons.  Each exitron is a dictionary of features.
    """

    if not exitrons:
        return []

    #=============================================================================
    # helpers
    #=============================================================================
    tab = str.maketrans("ACTG", "TGAC")
    def rc(seq):
        return seq.upper().translate(tab)[::-1]

    genome_fasta = pysam.FastaFile(genome)

    # define the aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1

    # define realigner
    realigner_a = Align.PairwiseAligner()
    realigner_a.mode = 'local'
    realigner_a.match = 2
    realigner_a.mismatch_score = -1
    realigner_a.open_gap_score = -2
    realigner_a.extend_gap_score = -1

    realigner_b = Align.PairwiseAligner()
    realigner_b.mode = 'local'
    realigner_b.match = 4
    realigner_b.mismatch_score = -2
    realigner_b.open_gap_score = -2
    realigner_b.extend_gap_score = 0

    #=============================================================================
    # main code
    #=============================================================================

    jitter = jitter*2
    # Need to compute reverse complement for splice junctions.
    res = []
    groups = {'+': [],
              '-': []}
    collection = {'+': [],
                  '-': []}

    # Need to sort exitrons for the clustering algorithm
    exitrons.sort(key=lambda x: (x['start'], x['end'], x['length']))

    # DEBUG: [(e['start'], e['end']) for e in exitrons]
    # filter one exitron at a time
    for exitron in exitrons:
        ao = exitron['ao']
        strand = exitron['strand']

        if collection[strand]:
            if (collection[strand][-1]['start'] - jitter <= exitron['start'] <= collection[strand][-1]['start'] + jitter and
                    collection[strand][-1]['end'] - jitter <= exitron['end'] <= collection[strand][-1]['end'] + jitter):
                collection[strand].append(exitron)
            else:
                groups[strand].append(collection[strand])
                collection[strand] = [exitron]
        else:
            collection[strand] = [exitron]

    groups['+'].append(collection['+'])
    groups['-'].append(collection['-'])

    for strand in groups:
        for group in groups[strand]:
            if not group:
                continue  # no exitrons found in this strand
            # calculate canonical spice sites and append read names.
            for e in group:
                start = e['start']
                end = e['end']
                # pos = db[e['transcript_id']].end - start if strand == '-' else start - db[e['transcript_id']].start
                # genome_seq = genome[e['chrom']][start:end - 1].upper()
                if strand == '+':
                    pos = start - db[e['transcript_id']].start + 1
                    genome_seq = db[e['transcript_id']].sequence(
                        genome)[pos:pos + e['length']]
                    e['splice_site'] = genome_seq[:2] + '-' + genome_seq[-2:]
                elif strand == '-':
                    pos = db[e['transcript_id']].end - start
                    genome_seq = db[e['transcript_id']].sequence(
                        genome)[pos - e['length']:pos]
                    right = genome_seq[:2]
                    left = genome_seq[-2:]
                    e['splice_site'] = genome_seq[:2] + '-' + genome_seq[-2:]
            try:
                consensus_e = max((e for e in group if e['splice_site'] in ['GT-AG', 'GC-AG', 'AT-AC']),
                                  key=lambda e: (e['ao'], e['length']))
                tot_ao = sum(e['ao'] for e in group)
                consensus_e['cluster_purity'] = round(
                    consensus_e['ao']/tot_ao, ndigits=2)
                if consensus_e['cluster_purity'] < cluster_purity:
                    continue
                consensus_e['ao'] = tot_ao
                consensus_reads = ''
                for e in group:
                    consensus_reads += ',' + ','.join(reads[(e['start'],
                                                             e['end'] - 1,
                                                             e['strand'])])
                consensus_e['reads'] = consensus_reads
            except ValueError:
                continue  # this occurs when there are no cannonical splice sites within the group

            ao = consensus_e['ao']
            start = consensus_e['start']
            end = consensus_e['end']
            chrm = consensus_e['chrom']

            # We subtract 1 because these coords are in BED format.
            mid = (start+end)/2
            a = bamfile.count(chrm, start=start - 1, stop=start,
                              read_callback=lambda x: x.mapq >= mapq)
            b = bamfile.count(chrm, start=end - 1, stop=end,
                              read_callback=lambda x: x.mapq >= mapq)
            c = bamfile.count(chrm, start=mid - 1, stop=mid,
                              read_callback=lambda x: x.mapq >= mapq)

            pso = ao/((a + b + c - ao*3)/3.0 + ao)
            dp = int(ao/pso) if pso > 0 else 0

            # Check whether attributes exceed minimum values
            if pso >= pso_min and ao >= ao_min:
                consensus_e['pso'] = round(pso, ndigits=4)
                consensus_e['dp'] = dp
                consensus_e['a'] = a
                consensus_e['b'] = b
                consensus_e['c'] = c

                # calculate alignment scores
                genome_seq = genome_fasta[chrm][start - 50:end + 50 - 1].upper()
                intron_seq = genome_seq[50:-50]

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
                consensus_e['alignment50'] = max(alignment50_lr, alignment50_rl)

                # calculate repeats
                seq_5 = genome_seq[50 - 19: 50 + 21]
                seq_3 = genome_seq[-50 - 21: 19 - 50]
                consensus_e['rt_repeat'] = repeat_test(seq_3, seq_5, 4, 20)
                res.append(consensus_e)

    # realignmnet may double count
    if not skip_realign and res:
        # realign
        pretty_print(f'Realigning exitrons in {res[0]["chrom"]}')
        called_reads = ''.join(e['reads'] for e in res)
        for exitron in res:
            e_start = exitron['start']
            e_end = exitron['end']
            chrm = exitron['chrom']
            e_length = exitron['length']

            for exon in db.children(db[exitron['transcript_id']], featuretype='exon', order_by='start'):
                if exon.start <= e_start <= exon.end:
                    exon_start = exon.start
                    exon_end = exon.end
                    break

            for read in bamfile.fetch(chrm, e_start, e_end):
                if read.query_name in called_reads or read.mapping_quality < mapq:
                    continue
                if (any([max(0, min(exon_end, i[1]) - max(exon_start, i[0])) > 0
                        for i in bamfile.find_introns([read]).keys()])):
                    pos = [p for p in read.get_aligned_pairs() if (
                        p[1] != None and exon_start <= p[1] <= exon_end)]
                    try:
                        start = min(p[0] for p in pos if p[0] != None)
                        end = max(p[0] for p in pos if p[0] != None)
                    except:
                        continue  # no nt overlap with exon

                    try:
                        r_seq = read.seq[start:end].upper()
                    except:
                        continue  # strangely, sometimes pysam returns None from read.seq
                        # it is a rare bug and I don't know why it happens.
                    if not r_seq:
                        continue
                    # exon.sequence(-) is much faster than using pysam FastaFile
                    g_seq = exon.sequence(genome).upper(
                    ) if exon.strand == '+' else rc(exon.sequence(genome).upper())

                    exitron_seq = g_seq[e_start - exon.start +
                                        1 - 2: e_end - exon.start + 2]

                    left = exitron_seq[:2]
                    right = exitron_seq[-2:]
                    exitron_seq = exitron_seq[2:-2]
                    similarity = realigner_a.score(
                        read.seq[max(0, start-e_length):end+e_length], exitron_seq)
                    # if no alignment, just continue
                    similarity = similarity/(e_length*2) if similarity else 1
                    if similarity <= 0.7:
                        alignments = realigner_b.align(
                            g_seq, r_seq)
                        #TODO: recheck that rewriting function from pairwise2 doesn't mess things up
                        alignments = islice(alignments, 10)
                        if (any(re.findall(f'{left}--*', aln[1]) and (e_length - 10 - len(r_seq)*0.05 <= aln[1].count('-') <= e_length + 10 + len(r_seq)*0.05) for aln in alignments) or
                                any(re.findall(f'--*{right}', aln[0]) and (e_length - 10 - len(r_seq)*0.05 <= aln[0].count('-') <= e_length + 10 + len(r_seq)*0.05) for aln in alignments)):
                            exitron['reads'] += f',REALIGNED_{read.query_name}'
                            exitron['ao'] += 1
                            exitron['pso'] = exitron['ao']/(
                                (exitron['a'] + exitron['b'] + exitron['c'] - exitron['ao']*3)/3.0 + exitron['ao'])
                            called_reads += f',{read.query_name}'

    return res


def identify_transcripts(exitrons, db, bamfilename, tmp_path, save_abundance, out_fn, cores):
    """


    Parameters
    ----------
    exitrons : list
        list of filtered exitrons.
    db : gffutils db file
    bamfilename : str
        bam filename.
    tmp_path : str
        temporary path
    save_abundance : bool
        if True, save abundance files
    out_fn : str
    cores : int
        number of cores for bamfile indexing/sorting

    Returns
    -------
    exitrons : list
        mutates exitrons with transcript features.

    """
    bamfile = pysam.AlignmentFile(bamfilename, 'rb', require_index=True)

    # construct new bamfile
    tmp_bamfile_exitrons = pysam.AlignmentFile(
        tmp_path + '/e_tmp.bam', 'wb', template=bamfile)
    if save_abundance:
        tmp_bamfile_normals = pysam.AlignmentFile(
            tmp_path + '/n_tmp.bam', 'wb', template=bamfile)

    for e in exitrons:
        e_reads = e['reads'].split(',')
        # fetch reads at exitron junction
        for read in bamfile.fetch(e['chrom'], start=int(e['start']), stop=int(e['end'])):
            if read.query_name in e_reads:
                tmp_bamfile_exitrons.write(read)
            elif save_abundance:
                tmp_bamfile_normals.write(read)
    tmp_bamfile_exitrons.close()
    if save_abundance:
        tmp_bamfile_normals.close()
    bamfile.close()

    # sort bamfiles and index
    pysam.sort('-@', str(cores), '-o', tmp_path +
               '/e_tmp_sorted.bam', tmp_path + '/e_tmp.bam')
    if save_abundance:
        pysam.sort('-o', tmp_path + '/n_tmp_sorted.bam',
                   tmp_path + '/n_tmp.bam')
    pysam.index(tmp_path + '/e_tmp_sorted.bam')
    if save_abundance:
        pysam.index(tmp_path + '/n_tmp_sorted.bam')

    # build a small gtf file of only thoes exitron spliced genes
    with open(tmp_path + '/tmp.gtf', 'w') as f:
        for e in exitrons:
            for region in db.children(db[e['gene_id']], order_by='start'):
                f.write(str(region) + '\n')

    # run liqa to create refgene
    subprocess.run(['liqa',
                    '-task',
                    'refgene',
                    '-ref',
                    f'{tmp_path + "/tmp.gtf"}',
                    '-format',
                    'gtf',
                    '-out',
                    f'{tmp_path + "/tmp.refgene"}'])

    # run liqa to quantify isoform expression
    jitter = 10  # TODO make this an argument
    subprocess.run(['liqa',
                    '-task',
                    'quantify',
                    '-refgene',
                    f'{tmp_path + "/tmp.refgene"}',
                    '-bam',
                    f'{tmp_path + "/e_tmp_sorted.bam"}',
                    '-out',
                    f'{tmp_path + "/isoform_estimates.out"}',
                    '-max_distance',
                    f'{jitter}',
                    '-f_weight',
                    '1'])
    if save_abundance:
        # run liqa to quantify isoform expression
        jitter = 10  # TODO make this an argument
        subprocess.run(['liqa',
                        '-task',
                        'quantify',
                        '-refgene',
                        f'{tmp_path + "/tmp.refgene"}',
                        '-bam',
                        f'{tmp_path + "/n_tmp_sorted.bam"}',
                        '-out',
                        f'{os.path.splitext(out_fn)[0] + ".isoform.normals"}',
                        '-max_distance',
                        f'{jitter}',
                        '-f_weight',
                        '1'])

    ie = pd.read_csv(f'{tmp_path}/isoform_estimates.out', sep='\t')
    if save_abundance:
        ie.to_csv(os.path.splitext(out_fn)[
                  0] + ".isoform.exitrons", sep='\t', index=False)
    for e in exitrons:
        gene = e['gene_symbol']
        ie_slice = ie[ie['GeneName'] == gene].sort_values(
            ascending=False, by='RelativeAbundance')
        transcripts = list(
            ie_slice[ie_slice['RelativeAbundance'] > 0.1]['IsoformName'])
        try:
            if not transcripts:
                transcripts = list(ie_slice[ie_slice['RelativeAbundance'] == max(
                    ie_slice['RelativeAbundance'])]['IsoformName'])
        except ValueError:
            # transcript could not be measured by liqa
            e['transcript_id'] += ',NA'  # revert back to default transcript_id
            continue
        t_str = ''
        for t in transcripts:
            t_str += f'{t},{round(float(ie_slice[ie_slice["IsoformName"] == t]["RelativeAbundance"]), 4)};'
        e['transcript_id'] = t_str[:-1]  # leave off trailing ;
    return exitrons


# ===============================================================================
# Main
# ===============================================================================


def exitrons_in_chrm(bamfilename, referencename, genomename, chrm, mapq, pso_min, ao_min, cluster_purity, jitter, skip_realign):
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
                                     jitter)
    exitrons = filter_exitrons(exitrons,
                               reads,
                               bamfile,
                               genomename,
                               db,
                               skip_realign,
                               mapq,
                               pso_min,
                               ao_min,
                               cluster_purity,
                               jitter)
    bamfile.close()
    del db

    return exitrons, chrm


def main(tmp_path, input_,
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
    # try:
    #     # Define chrms
    #     fa = pysam.FastaFile(genome)
    #     chrms = fa.references
    #     fa.close()
    # except:
    #     pretty_print('Building FASTA index file')
    #     pysam.faidx(genome)
    #     try:
    #         fa = pysam.FastaFile(genome)
    #         chrms = fa.references
    #         fa.close()
    #     except:
    #         pretty_print(f'ERROR: Unable to read FASTA file {genome}')
    #         rmtree(tmp_path)
    #         sys.exit(1)

    chrms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5',
              'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
              'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
              'chr21', 'chr22', 'chrX', 'chrY']

    # Check if bamfile can be opened and there is an index.
    try:
        bamfile = pysam.AlignmentFile(input_, 'rb', require_index=True)
    except FileNotFoundError:
        try:
            pretty_print('Building bam index file')
            pysam.index(input_)
            bamfile = pysam.AlignmentFile(input_, 'rb', require_index=True)
        except FileNotFoundError:
            pretty_print(
                f'ERROR: There is a problem opening bam file at: {input_}')
            sys.exit(1)
    bamfile.close()

    # Check if annotation has a tabix index
    try:
        try:
            gtf = pysam.TabixFile(reference_transcriptome, parser=pysam.asGTF())
            gtf.close()
        except OSError:
            pretty_print('Building tabix index.')
            pysam.tabix_index(reference_transcriptome, preset='gtf')
    except:
        pretty_print(
            f'ERROR: There is a problem reading the annotation file at: {reference_transcriptome}')
        pretty_print(f'Please make sure to use bgzip to compress your annotation file.')
        sys.exit(1)

    # Check if LIQA is available
    try:
        subprocess.run(['liqa'], capture_output=True)
    except FileNotFoundError:
        pretty_print('ERROR: Unable to locate LIQA. Please install with "pip install liqa"')
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

    # Define chrms
    # chrms = set()
    # query = db.execute("select seqid from features")
    # for x in query:
    #     chrms.add(x['seqid'])
    # chrms = list(chrms)
    # chrms.sort()

    # check for blacklist
    this_dir = os.path.dirname(os.path.realpath(__file__))
    blacklist_path = os.path.join(this_dir, 'blacklist.tsv')

    try:
        b = open(blacklist_path)
        b.close()
    except FileNotFoundError:
        pretty_print('No blacklist found -- continuing without it, beware of false positives')

    # Begin exitron calling
    global results
    results = {}

    def collect_result(output):
        exitrons = output[0]
        chrm = output[1]
        pretty_print(f'Collecting data from {chrm}')
        sys.stdout.flush()
        results[chrm] = exitrons

    if int(cores) > 1:
        pool = mp.Pool(int(cores))
        threads = []
        for chrm in chrms:
            threads.append(pool.apply_async(exitrons_in_chrm, args=(input_,
                                                                    reference_transcriptome,
                                                                    genome,
                                                                    chrm,
                                                                    mapq,
                                                                    pso_min,
                                                                    ao_min,
                                                                    cluster_purity,
                                                                    jitter,
                                                                    skip_realign,
                                                                    ), callback=collect_result))
        pool.close()
        for t in threads:
            t.get()  # this gets any exceptions raised
        pool.join()
    else:
        for chrm in chrms:
            sys.stdout.flush()
            output = exitrons_in_chrm(input_,
                                      reference_transcriptome,
                                      genome,
                                      chrm,
                                      mapq,
                                      pso_min,
                                      ao_min,
                                      cluster_purity,
                                      jitter,
                                      skip_realign,
                                      )
            collect_result(output)
    exitrons = []
    for chrm in results:
        exitrons.extend(results[chrm])

    out_file_name = out
    if not out_file_name:
        out_file_name = f'{os.path.splitext(input_)[0]}.exitron'

    pretty_print('Quantifying transcripts.')
    sys.stdout.flush()
    # update transcripts
    identify_transcripts(exitrons,
                         db,
                         input_,
                         tmp_path,
                         save_abundance,
                         out_file_name,
                         cores)
    pretty_print(
        f'Finished exitron calling and filtering. Printing to {out_file_name}')
    sys.stdout.flush()
    with open(out_file_name, 'w') as out:
        header = ['chrom',
                  'start',
                  'end',
                  'name',
                  'region',
                  'ao',
                  'strand',
                  'gene_symbol',
                  'length',
                  'splice_site',
                  'transcript_id',
                  'pso',
                  'dp',
                  'cluster_purity',
                  'alignment50',
                  'rt_repeat',
                  'reads']
        # write header
        out.write('\t'.join(header) + '\n')
        for exitron in exitrons:
            out.write('\t'.join([str(exitron[column])
                      for column in header]) + '\n')

    # Clear tmp directory
    # rmtree(tmp_path)


if __name__ == '__main__':
     input_ ='test_data/test_data_long.bam'
     genome = '../scanexitron/hg38.fa'
     reference_transcriptome = '../scanexitron/gencode.v37.annotation.sorted.gtf.gz'
     out = 'test_data/test_data_long.exitron'
     cores = 0
     mapq = 50
     ao_min = 2
     pso_min = 0.05
     id_ = 'test'
     tmp_path = 'tmpscanexitron/'
     main(tmp_path, input_,
                        genome,
                        reference_transcriptome,
                        out,
                        cores,
                        mapq,
                        ao_min,
                        pso_min,
                        10,
                        0,
                        False,
                        False,
                        id_)

