#!/usr/bin/env python

"""This module contains the code for performing alignment."""

import math
import sys

def create_empty_matrix(y, x, val):
    """Initialize and return an empty array of y rows by x columns."""
    row = [val for i in range(x)]
    empty_mat = []
    for i in range(y): empty_mat.append(row[:])
    return(empty_mat)

def print_matrix(mat):
    """Print out matrix cleanly."""
    for row in mat: print(' '.join(list(map(str, row))))

def align_strings(profile, q):
    """Create the alignment graph of profile vs q.
    
    p(rofile) is the vertical string (rows)
    q is the horizontal string (columns)
    """
    p_len = profile.length
    alignment_matrix = create_empty_matrix(p_len+1, len(q)+1, 1)
    for x in range(1, len(q)+1): # columns
        for y in range(1, p_len+1): # rows
            probs = [x[y-1] for x in profile.matrix]
            nucls = {'A':0, 'C':1, 'G':2, 'T':3}
            q_char = q[x-1]

            align_score = profile.get_prob(q_char, y-1)
            alignment_matrix[y][x] = round(alignment_matrix[y-1][x-1]
                    * align_score, 6)
    return(alignment_matrix)

def find_best_spaced_motif(profile_1, profile_2, dna, gap_lengths):
    """Find the best spaced motif from dna aligning text_k1 and
    text_k2 with gap, g, in G.
    """
    k1_align = align_strings(profile_1, dna)
    k2_align = align_strings(profile_2.reverse(), dna[::-1])

    last_row_k1 = k1_align[-1]
    last_row_k2 = k2_align[-1][::-1]
    for i in range(profile_1.length): 
        last_row_k1[i] = 0
    for i in range(profile_2.length):
        last_row_k2[-i-1] = 0
    
    best_pos = -1
    best_gap = -1
    best_sum = -1
    # print("Allowed gap lengths are {}".format(gap_lengths))
    for i in range(profile_1.length, len(dna) + 1 - profile_2.length - 2):
        # print("Starting at new i of length {}".format(i))
        dists = []
        for gap in gap_lengths:
            # print("  DNA string is of length {}".format(len(dna)))
            # print("  gap plus index is of length {}".format(i + gap))
            if i + gap >= len(dna): 
                continue
            else:
                # [pos, gap, prod]
                # print("    Adding to the distance list")
                dists.append([i, gap, math.log(last_row_k1[i] * last_row_k2[i+gap] + 1)])
        if len(dists) > 0:
            max_dist = max(dists, key = lambda x: x[2])
            if max_dist[2] > best_sum:
                best_pos = max_dist[0]
                best_gap = max_dist[1]
                best_sum = max_dist[2]

    best_k1_pos = best_pos - profile_1.length 
    best_new_k1 = dna[best_k1_pos : best_k1_pos + profile_1.length]
    best_k2_pos = best_pos + best_gap
    best_new_k2 = dna[best_k2_pos : best_k2_pos + profile_2.length]
    if best_new_k1+best_new_k2 == "CATCCAGTA":
        print("Weird return from alignment")
        print(last_row_k1)
    return(best_k1_pos, best_new_k1+best_new_k2, best_gap)
    # return([[best_k1_pos, best_new_k1], best_gap, [best_k2_pos, best_new_k2]])
