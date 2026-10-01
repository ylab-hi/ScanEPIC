"""
Fast unit tests for ScanEPIC. These do not need a reference genome or GTF;
end-to-end runs on test_data/ are documented in the README.
"""
import pysam
import pytest
from click.testing import CliRunner

from src.cli import cli
from src.extract.short._cython_fi import find_introns
from src.extract.short.short import repeat_test
from src.tools.exitron2vcf import exitron_to_vcf_record


@pytest.mark.parametrize('args', [
    [],
    ['extract'],
    ['extract', 'short'],
    ['extract', 'long'],
    ['extract', 'single'],
    ['tools', 'exitron2vcf'],
])
def test_cli_help(args):
    result = CliRunner().invoke(cli, args + ['--help'])
    assert result.exit_code == 0, result.output


def test_repeat_test():
    # repeat_test is called on 40 bp windows around the two splice sites and
    # returns the longest shared k-mer (k = 4..20) near the junctions
    seq_5 = 'C' * 18 + 'GATTACA' + 'C' * 15
    seq_3 = 'T' * 18 + 'GATTACA' + 'T' * 15
    assert repeat_test(seq_3, seq_5, 4, 20) == 'GATTACA'
    assert repeat_test('A' * 40, 'C' * 40, 4, 20) == ''


def make_bam(path, reads):
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': 'chr1', 'LN': 10000}]}
    with pysam.AlignmentFile(path, 'wb', header=header) as out:
        for name, pos, cigar, seq in reads:
            r = pysam.AlignedSegment()
            r.query_name = name
            r.reference_id = 0
            r.reference_start = pos
            r.cigarstring = cigar
            r.query_sequence = seq
            r.mapping_quality = 60
            out.write(r)
    pysam.index(str(path))


def test_find_introns(tmp_path):
    bam_path = tmp_path / 'test.bam'
    make_bam(bam_path, [('r1', 100, '10M50N10M', 'A' * 10 + 'C' * 10),
                        ('r2', 105, '5M50N8M', 'A' * 5 + 'C' * 8),
                        ('r3', 200, '20M', 'G' * 20)])
    with pysam.AlignmentFile(str(bam_path)) as bam:
        introns, reads = find_introns(bam.fetch('chr1'))
    assert dict(introns) == {(110, 160): 2}
    # (seq, '.', left anchor, right anchor, read position of junction)
    assert reads[(110, 160)][0][2:] == (10, 10, 10)
    assert reads[(110, 160)][1][2:] == (5, 8, 5)


def test_exitron_to_vcf_record_is_deletion(tmp_path):
    fasta = tmp_path / 'genome.fa'
    fasta.write_text('>chr1\nACGTACGTACGTACGTACGT\n')
    pysam.faidx(str(fasta))
    exitron = {'chrom': 'chr1', 'start': 3, 'end': 9, 'name': 'Xd3-9',
               'ao': 5, 'dp': 10, 'pso': 0.5, 'strand': '+', 'length': 5,
               'splice_site': 'GT-AG', 'gene_symbol': 'X'}
    with pysam.FastaFile(str(fasta)) as genome:
        record = exitron_to_vcf_record(exitron, genome).split('\t')
    pos, ref, alt = record[1], record[3], record[4]
    # REF = anchor base + 5 deleted bases, ALT = anchor base
    assert pos == '3'
    assert ref == 'GTACGT'
    assert alt == 'G'
    assert 'END=8' in record[7] and 'SVLEN=-5' in record[7]
