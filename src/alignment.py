#!/usr/bin/env python

"""This module contains the code for performing alignment."""

import sys

def create_empty_matrix(y, x):
    """Initialize and return an empty array of y rows by x columns."""
    row = [0 for i in range(x)]
    empty_mat = []
    for i in range(y): empty_mat.append(row[:])
    return(empty_mat)

def print_matrix(mat):
    """Print out matrix cleanly."""
    for row in mat: print(' '.join(list(map(str, row))))

def create_scoring_dict():
    """Create dictionary to hold alignment scores."""
    alphabet = 'ACTG'
    score_dict = {}
    for nucl in alphabet:
        for nucl2 in alphabet:
            score_dict['%s%s' % (nucl, nucl2)] = int(nucl == nucl2)
    return(score_dict)

def align_strings(p, q, profile):
    """Create the alignment graph and backtrack matrix of p and q.
    
    p is the vertical string (rows)
    q is the horizontal string (columns)
    """
    score_dict = create_scoring_dict()
    alignment_matrix = create_empty_matrix(len(p)+1, len(q)+1)
    for x in range(1, len(q)+1): # columns
        for y in range(1, len(p)+1): # rows
            p_char = p[y-1]
            q_char = q[x-1]

            align_score = profile.get_prob(p_char, y-1)
            alignment_matrix[y][x] = alignment_matrix[y-1][x-1] * align_score
    return(alignment_matrix)

def find_best_spaced_motif(text_k1, text_k2, dna, gap_lengths):
    """Find the best spaced motif from dna aligning text_k1 and
    text_k2 with gap, g, in G.
    """
    k1_align = align_strings(text_k1, dna)
    k2_align = align_strings(text_k2[::-1], dna[::-1])

    last_row_k1 = k1_align[-1]
    last_row_k2 = k2_align[-1][::-1]
    for i in range(len(text_k1)): 
        last_row_k1[i] = -sys.maxsize
    
    best_pos = -1
    best_gap = -1
    best_sum = -1
    for i in range(len(last_row_k1) - len(text_k2) - 2):
        dists = []
        for gap in gap_lengths:
            # [pos, gap, sum]
            if i + gap >= len(dna): 
                continue
            else:
                dists.append([i, gap, last_row_k1[i] + last_row_k2[i+gap]])
        max_dist = max(dists, key = lambda x: x[2])
        if max_dist[2] > best_sum:
            best_pos = max_dist[0]
            best_gap = max_dist[1]
            best_sum = max_dist[2]

    best_k1_pos = best_pos - len(text_k1) 
    best_new_k1 = dna[best_k1_pos : best_k1_pos+len(text_k1)]
    best_k2_pos = best_pos + best_gap
    best_new_k2 = dna[best_k2_pos : best_k2_pos+len(text_k2)]
    return([[best_k1_pos, best_new_k1], best_gap, [best_k2_pos, best_new_k2]])
