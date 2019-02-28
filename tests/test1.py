#!/usr/bin/env python

from spacedmotif.src import alignment

def test1():
   k1 = 'ATCG'
   k2 = 'GAA'
   dna = 'GCCATCAGTTCAGAGTCC'
   assert (alignment.find_best_spaced_motif(k1, k2, dna, [4, 5, 6]) ==
           [[3, 'ATCA'], 5, [12, 'GAG']])

