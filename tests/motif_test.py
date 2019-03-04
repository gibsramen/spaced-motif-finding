#!/usr/bin/env python

import pytest

from src import motif

def test_motif_gap_length():
    dna = 'GCCATCAGTTCAGAGTCC'
    positions = [3, 7, 12, 14]
    test_motif = motif.SpacedMotif(dna, positions)
    assert(test_motif.gap_length == 5)

def test_motif_k1_mer():
    dna = 'GCCATCAGTTCAGAGTCC'
    positions = [3, 7, 12, 14]
    test_motif = motif.SpacedMotif(dna, positions)
    assert(test_motif.k1_mer == 'ATCA')

def test_motif_k2_mer():
    dna = 'GCCATCAGTTCAGAGTCC'
    positions = [3, 7, 12, 14]
    test_motif = motif.SpacedMotif(dna, positions)
    assert(test_motif.k2_mer == 'GAG')

def test_motif_k1_mer_pos():
    dna = 'GCCATCAGTTCAGAGTCC'
    positions = [3, 7, 12, 14]
    test_motif = motif.SpacedMotif(dna, positions)
    assert(test_motif.k1_mer_pos == 3)

def test_motif_k2_mer_pos():
    dna = 'GCCATCAGTTCAGAGTCC'
    positions = [3, 7, 12, 14]
    test_motif = motif.SpacedMotif(dna, positions)
    assert(test_motif.k2_mer_pos == 12)

def test_motif_distance():
    adna = 'GCCATCAGTTCAGAGTCC'
    apositions = [3, 7, 12, 14]
    atest_motif = motif.SpacedMotif(adna, apositions)

    bdna = 'GCCATCAGTTCAGAGTCC'
    bpositions = [4, 8, 12, 14]
    btest_motif = motif.SpacedMotif(bdna, bpositions)
    assert(atest_motif.distance(btest_motif) == 4)

def test_profile_A():
    A_dist = [0, 0.2, 0.4, 0.3]
    C_dist = [0.3, 0.3, 0.2, 0.5]
    G_dist = [0.6, 0.4, 0.2, 0.2]
    T_dist = [0.1, 0.1, 0.2, 0]
    profile_matrix = motif.Profile([A_dist, C_dist, G_dist, T_dist])
    assert(profile_matrix.get_prob('A', 2) == 0.4)

def test_profile_OOB():
    A_dist = [0, 0.2, 0.4, 0.3]
    C_dist = [0.3, 0.3, 0.2, 0.5]
    G_dist = [0.6, 0.4, 0.2, 0.2]
    T_dist = [0.1, 0.1, 0.2, 0]
    profile_matrix = motif.Profile([A_dist, C_dist, G_dist, T_dist])
    with pytest.raises(IndexError):
        profile_matrix.get_prob('A', 5) == 0
