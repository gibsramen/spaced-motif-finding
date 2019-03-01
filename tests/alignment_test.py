#!/usr/bin/env python

import pytest

from src import motif
from src import alignment

def test_find_spaced_motifs():
    A_dist = [0.2, 0.3, 0.1, 0.4]
    C_dist = [0.1, 0.5, 0.2, 0.1]
    G_dist = [0.4, 0.1, 0.3, 0.2]
    T_dist = [0.3, 0.1, 0.4, 0.3]
    profile_matrix = motif.Profile([A_dist, C_dist, G_dist, T_dist])

    A_dist = [0.2, 0.1, 0.1, 0.1]
    C_dist = [0.3, 0.1, 0.3, 0.7]
    G_dist = [0.1, 0.6, 0.5, 0.1]
    T_dist = [0.4, 0.2, 0.1, 0.1]
    profile_matrix2 = motif.Profile([A_dist, C_dist, G_dist, T_dist])

    dna = 'GCTATCATGGC'
    gap_lengths = [1,2,3,4]
    out = alignment.find_best_spaced_motif(profile_matrix, profile_matrix2,
            dna, gap_lengths)
    assert(out == [[0, 'GCTA'], 3, [7, 'TGGC']])
