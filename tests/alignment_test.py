#!/usr/bin/env python

import pytest

from src import motif
from src import alignment

def test_find_spaced_motifs():
    A_dist = [0, 0.2, 0.4, 0.3, 0.3, 0.2, 0.5]
    C_dist = [0.3, 0.3, 0.2, 0.5, 0.3, 0.4, 0.1]
    G_dist = [0.6, 0.4, 0.2, 0.2, 0.2, 0.1, 0.2]
    T_dist = [0.1, 0.1, 0.2, 0, 0.2, 0.3, 0.2]
    profile_matrix = motif.Profile([A_dist, C_dist, G_dist, T_dist])

    dna = 'GCCATCA'
    k1 = 'GC'
    k2 = 'TCA'
    gap_lengths = [1,2,3]
    alignment_matrix = alignment.align_strings(k1, k2, profile_matrix)
    assert(alignment_matrix == 'test')
