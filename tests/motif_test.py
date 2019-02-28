#!/usr/bin/env python

from spacedmotif.src import motif

def test_motif_gap_length():
    k1 = 'ATCG'
    k2 = 'GAA'
    dna = 'GCCATCAGTTCAGAGTCC'
    positions = [3, 7, 12, 14]
    test_motif = motif.SpacedMotif(dna, positions)
    assert(test_motif.gap_length == 5)

def test_motif_k1_mer():
    k1 = 'ATCG'
    k2 = 'GAA'
    dna = 'GCCATCAGTTCAGAGTCC'
    positions = [3, 7, 12, 14]
    test_motif = motif.SpacedMotif(dna, positions)
    assert(test_motif.k1_mer == 'ATCA')

def test_motif_k2_mer():
    k1 = 'ATCG'
    k2 = 'GAA'
    dna = 'GCCATCAGTTCAGAGTCC'
    positions = [3, 7, 12, 14]
    test_motif = motif.SpacedMotif(dna, positions)
    assert(test_motif.k2_mer == 'GAG')

def test_motif_k1_mer_pos():
    k1 = 'ATCG'
    k2 = 'GAA'
    dna = 'GCCATCAGTTCAGAGTCC'
    positions = [3, 7, 12, 14]
    test_motif = motif.SpacedMotif(dna, positions)
    assert(test_motif.k1_mer_pos == 3)

def test_motif_k2_mer_pos():
    k1 = 'ATCG'
    k2 = 'GAA'
    dna = 'GCCATCAGTTCAGAGTCC'
    positions = [3, 7, 12, 14]
    test_motif = motif.SpacedMotif(dna, positions)
    assert(test_motif.k2_mer_pos == 7)
