# -*- coding: utf-8 -*-
"""
Short-read exitron extraction module for ScanEPIC

This module identifies and quantifies exitron splicing events from short-read RNA-seq data.
Exitrons are internal splicing events within annotated exons that can affect gene regulation.

@author: Josh Fry, Northwestern University, YangLab
"""
__version__ = 'v2.0.0'
import sys
import os
import pysam
import gffutils
import pickle
import click
import re
import multiprocessing as mp
from Bio import Align
from time import strftime, localtime
from ._cython_fi import find_introns

#=============================================================================
# Helpers
#=============================================================================

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

#=============================================================================
# Main functions
#=============================================================================

def exitron_caller(bamfile, referencename, chrm, db, known_introns, mapq=50):
    """


    Parameters
    ----------
    bamfile : pysam alignment file
    referencename : str
        transcriptome reference location.
    chrm : str
    db : gffutils databse
    known_introns : list
        optional list of known intron locations
    mapq : int

    Returns
    -------
    list of exitrons
    list of read information for each exitron

    """
    # call introns
    introns_raw, reads = find_introns(
        (read for read in bamfile.fetch(chrm) if read.mapping_quality >= mapq)
    )


    if known_introns:
        introns = [intron for intron in introns_raw if (
            (chrm, intron[1]) not in known_introns['D'] and
            (chrm, intron[0]+1) not in known_introns['A'])]
    else:
        # if no preproccessing, just continue
        introns = introns_raw

    gtf = pysam.TabixFile(referencename, parser=pysam.asGTF())
    exitrons = []
    exitrons_added = []
    known_splices = set()


    for intron in introns:
        intron_start = intron[0]
        intron_end = intron[1]
        intron_witnesses = introns_raw[intron]
        intersection = gtf.fetch(chrm, intron_start - 1, intron_end + 2)
        for feature in intersection:
            try:
                # some entries like 'gene' do not have transcript types
                # in this case, just continue
                gene_type = feature.gene_type
            except:
                continue
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
                                     'name': f'{gene_name}d{intron_start}-{intron_end + 1}',
                                     'region': region_type,
                                     'ao': intron_witnesses,
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
                        exon_intervals = [(x.start, x.end) for x in db.children(db[feature.transcript_id], featuretype = 'exon')]
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
                                             'name': f'{gene_name}d{intron_start}-{intron_end + 1}',
                                             'region': region,
                                             'ao': intron_witnesses,
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


def filter_exitrons(exitrons, reads, bamfile, genome, db, mapq=50, pso_min=0.005, ao_min=1, pso_ufmin=0, ao_ufmin=0, anchor_min=5):
    """


    Parameters
    ----------
    exitrons : list
        list of exitrons, each of which is a dictionary
    reads : list
        list of read information
    bamfile : pysam alignment file
    genome : str
        location of reference genome
    db : gffutils database file
    mapq : int
    pso_min : TYPE, optional
    ao_min : TYPE, optional
    pso_ufmin : float, optional
    ao_ufmin : int, optional
    anchor_min : TYPE, optional

    Returns
    -------
    list of exitrons that pass filtering, each of which is a dictionary

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


    #=============================================================================
    # filtering
    #=============================================================================
    res = []
    # filter one exitron at a time
    for exitron in exitrons:
        ao = exitron['ao']
        if ao < ao_min: continue
        chrm = exitron['chrom']
        start = exitron['start']
        end = exitron['end']
        strand = exitron['strand']

        if ao < ao_ufmin:
            continue  # not enough unfiltered reads support this exitron
        read_data = reads[(start, end - 1)]

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

        max_r_anchor_len = 0
        max_l_anchor_len = 0
        ao_true = 0
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
                ao_true += 1
                if r_anchor_len > max_r_anchor_len: max_r_anchor_len = r_anchor_len
                if l_anchor_len > max_l_anchor_len: max_l_anchor_len = l_anchor_len

        if ao_true == 0 and ao_min > 0:
            exitron['ao_unfiltered'] = ao

        # We subtract 1 because these coords are in BED format.
        mid = (start+end)/2
        a = bamfile.count(chrm, start=start - 1, stop=start,
                          read_callback=lambda x: x.mapq > mapq)
        b = bamfile.count(chrm, start=end - 1, stop=end,
                          read_callback=lambda x: x.mapq > mapq)
        c = bamfile.count(chrm, start=mid - 1, stop=mid,
                          read_callback=lambda x: x.mapq > mapq)

        pso_true = ao_true/((a + b + c - ao_true*3)/3.0 + ao_true)
        psi_true = 1 - pso_true
        dp_true = int(ao_true/pso_true) if pso_true > 0 else 0

        # Check whether attributes exceed minimum values
        if (pso_true >= pso_min and
            max_l_anchor_len >= anchor_min and
            max_r_anchor_len >= anchor_min and
            ao_true >= ao_min):

            exitron['ao'] = ao_true
            exitron['pso'] = pso_true
            exitron['psi'] = psi_true
            exitron['dp'] = dp_true

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
            exitron['alignment50'] = max(alignment50_lr, alignment50_rl)

            # calculate repeats
            seq_5 = genome_seq[50 - 19: 50 + 21]
            seq_3 = genome_seq[-50 - 21: 19 - 50]
            exitron['rt_repeat'] = repeat_test(seq_3, seq_5, 4, 20)
            res.append(exitron)
    return res


# ===============================================================================
# Main
# ===============================================================================


def exitrons_in_chrm(bamfilename, referencename, genomename, chrm, known_introns, mapq, pso_min, ao_min, pso_ufmin, ao_ufmin, anchor_min):
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
                                      known_introns,
                                      mapq,
                                      )
    exitrons = filter_exitrons(exitrons,
                                reads,
                                bamfile,
                                genomename,
                                db,
                                mapq,
                                pso_min,
                                ao_min,
                                pso_ufmin,
                                ao_ufmin,
                                anchor_min)
    bamfile.close()
    del db

    return exitrons, chrm


def main(input_,
        genome,
        reference_transcriptome,
        out,
        cores,
        mapq,
        ao,
        pso,
        id_,
        vcf):

    #=============================================================================
    # set internal defaults:
    #   change these only if you want backwards compatibility scanexitron v1
    #=============================================================================
    pso_ufmin = 0
    ao_ufmin = 0
    anchor_min = 5


    #=============================================================================
    # check that all files are properly located
    #=============================================================================
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
            pretty_print('Building tabix index for GTF file.')
            pysam.tabix_index(reference_transcriptome, preset='gtf')
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

    try:
        # prune introns with known donors and acceptors
        with open('introns_v37.obj', 'rb') as p:
            known_introns = pickle.load(p)
    except:
        known_introns = []

    #=============================================================================
    # begin exitron calling
    #=============================================================================

    global results
    results = {}

    def collect_result(output):
        exitrons = output[0]
        chrm = output[1]
        pretty_print(f'Collecting data from {chrm}')
        sys.stdout.flush()
        results[chrm] = exitrons

    if cores > 1:
        pool = mp.Pool(int(cores))
        threads = []
        for chrm in chrms:
            threads.append(pool.apply_async(exitrons_in_chrm, args=(input_,
                                                                    reference_transcriptome,
                                                                    genome,
                                                                    chrm,
                                                                    known_introns,
                                                                    mapq,
                                                                    pso,
                                                                    ao,
                                                                    pso_ufmin,
                                                                    ao_ufmin,
                                                                    anchor_min,
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
                                      known_introns,
                                      mapq,
                                      pso,
                                      ao,
                                      pso_ufmin,
                                      ao_ufmin,
                                      anchor_min,
                                      )
            collect_result(output)

    out_file_name = out
    if not out_file_name:
        prefix = os.path.splitext(os.path.basename(
            bamfile.filename.decode('UTF-8')))[0]
        out_file_name = f'{prefix}.exitron'

    #=============================================================================
    # write exitron output
    #=============================================================================
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
                  'alignment50',
                  'rt_repeat']
        if id_:
            header += ['id']
        # write header
        for column in header:
            out.write(column + '\t')
        out.write('\n')
        for chrm in chrms:
            # check if chromosome is empty or not
            try:
                if results[chrm]:
                    for exitron in results[chrm]:
                        if id_:
                            exitron['id'] = id_
                        out.write('\t'.join([str(exitron[column])
                                  for column in header]) + '\n')

            except KeyError:
                pretty_print(
                    f'Thread most likely crashed on chromosome \'{chrm}\' without reporting exception.  Try fewer cores or allocate more memory.')
                sys.stdout.flush()
                sys.exit(1)

        #=============================================================================
        # write VCF file if requested
        #=============================================================================
        if vcf:
            out_file_name_vcf = out_file_name + '.vcf'
            pretty_print(
                f'Printing VCF to {out_file_name}')
            sys.stdout.flush()
            with open(out_file_name_vcf, 'w') as out:
                vcf_header ='''##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##FILTER=<ID=LowQual,Description="PE/SR support below 3 or mapping quality below 20.">
##INFO=<ID=DP,Number=2,Type=Integer,Description="Total depth of junction">
##INFO=<ID=SplicedSite,Number=2,Type=String,Description="Splieced site">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Junction strand">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=AO,Number=1,Type=Integer,Description="Reads support of the exitron">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=PSO,Number=1,Type=Integer,Description="Percent spliced-out">
##INFO=<ID=GeneName,Number=1,Type=String,Description="Gene name">
##INFO=<ID=GeneID,Number=1,Type=String,Description="Gene ID">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
##INFO=<ID=SRQ,Number=1,Type=Float,Description="Split-read consensus alignment quality">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Split-read consensus sequence">
##INFO=<ID=CE,Number=1,Type=Float,Description="Consensus sequence entropy">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Predicted length of the insertion">
##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Predicted microhomology length using a max. edit distance of 2">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}
                '''.format(os.path.splitext(os.path.basename(out_file_name))[0])
                # write header
                out.write(vcf_header)
                genome_seq = pysam.FastaFile(genome)
                for chrm in chrms:
                    # check if chromosome is empty or not
                    try:
                        if results[chrm]:
                            for exitron in results[chrm]:
                                out.write('{}\t{}\t{}\t{}\t{}\t.\t.\tSVTYPE=DEL;END={};AO={};DP={};STRAND={};SplicedSite={};GeneName={};GeneID={};SVLEN=-{};PSO={};PSI={}\tGT\t0/1\n'.
                                              format(exitron['chrom'],
                                                     exitron['start'],
                                                     exitron['name'],
                                                     genome_seq[exitron['chrom']][exitron['start'] - 1: exitron['start']],
                                                     genome_seq[exitron['chrom']][exitron['start'] - 1: exitron['end'] - 1],
                                                     exitron['end'] - 1,
                                                     exitron['ao'],
                                                     exitron['dp'],
                                                     exitron['strand'],
                                                     exitron['splice_site'],
                                                     exitron['gene_symbol'],
                                                     exitron['gene_id'],
                                                     exitron['length'],
                                                     exitron['pso'],
                                                     1 - exitron['pso']))
                    except KeyError:
                        pretty_print(
                            f'Thread most likely crashed on chromosome \'{chrm}\' without reporting exception.  Try fewer cores or allocate more memory.')
                        sys.stdout.flush()
                        sys.exit(1)


if __name__ == '__main__':
     input_ ='test_data/test_data.bam'
     genome = '../scanexitron/hg38.fa'
     reference_transcriptome = '../scanexitron/gencode.v37.annotation.sorted.gtf.gz'
     out = 'test_data/test_data.exitron'
     cores = 0
     mapq = 50
     ao = 3
     pso = 0.05
     id_ = 'test'
     vcf = True

     main(input_,
          genome,
          reference_transcriptome,
          out,
          cores,
          mapq,
          ao,
          pso,
          id_,
          vcf)








