#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 20:50:43 2022

@author: jpfry
"""
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment

import collections

def find_introns(read_iterator):
    """Return a dictionary {(start, stop): count}
    Listing the intronic sites in the reads (identified by 'N' in the cigar strings),
    and their support ( = number of reads ).
    read_iterator can be the result of a .fetch(...) call.
    Or it can be a generator filtering such reads. Example
    samfile.find_introns((read for read in samfile.fetch(...) if read.is_reverse)
    """
    cdef:
        int base_position, junc_start, nt, read_position, i
        int op
        AlignedSegment r
        int BAM_CREF_SKIP = 3 #BAM_CREF_SKIP


    res = collections.Counter()
    reads = collections.defaultdict(list)
    match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
    for r in read_iterator:
        base_position = r.pos
        read_position = 0

        for i, (op, nt) in enumerate(r.cigartuples):
            if op in match_or_deletion:
                base_position += nt
                read_position += nt
            elif op == BAM_CREF_SKIP:
                junc_start = base_position
                base_position += nt
                res[(junc_start, base_position)] += 1
                reads[(junc_start, base_position)].append((r.seq, '.', r.cigartuples[i-1][1], r.cigartuples[i+1][1], read_position))
    return res, reads