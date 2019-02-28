#!/usr/bin/env python

"""This module contains the code for performing alignment."""

def create_empty_matrix(y, x):
    """Initialize and return an empty array of y rows by x columns."""
    row = [0 for i in range(x)]
    empty_mat = []
    for i in range(y): empty_mat.append(row[:])
    return(empty_mat)

def print_matrix(mat):
    for row in mat: print(' '.join(list(map(str, row))))

def create_scoring_dict():
    """Create dictionary to hold alignment scores."""
    alphabet = 'ACTG'
    score_dict = {}
    for nucl in alphabet:
        for nucl2 in alphabet:
            score_dict['%s%s' % (nucl, nucl2)] = int(nucl == nucl2)
    return(score_dict)

def align_strings(p, q):
    """Create the alignment graph and backtrack matrix of p and q.
    
    p is the vertical string (rows)
    q is the horizontal string (columns)
    """
    score_dict = create_scoring_dict()
    backtrack_matrix = create_empty_matrix(len(p)+1, len(q)+1)
    alignment_matrix = create_empty_matrix(len(p)+1, len(q)+1)

    for x in range(1, len(q)+1):
        for y in range(1, len(p)+1):
            p_char = p[y-1]
            q_char = q[x-1]

            align_score = score_dict.get('%s%s' % (p_char, q_char))
            alignment_matrix[y][x] = alignment_matrix[y-1][x-1] + align_score
    return(alignment_matrix)

def find_best_spaced_motif(text_k1, text_k2, dna, G):
    """Find the best spaced motif from dna aligning text_k1 and
    text_k2 with gap, g in G.
    """
    k1_align = align_strings(text_k1, dna)
    k2_align = align_strings(text_k2[::-1], dna[::-1])

    last_row_k1 = k1_align[-1]
    last_row_k2 = k2_align[-1][::-1]
    # need to set some of last row values to -inf
    # first k1 values can't be reached
    for i in range(len(text_k1)):
        last_row_k1[i] = -100
    print(last_row_k1)
    print(last_row_k2)
    print()
    
    best_pos = -1
    best_gap = -1
    best_sum = -1
    for i in range(len(last_row_k1) - len(text_k2) - 2):
        dists = []
        for g in G:
            # [pos, gap, sum]
            if i + g >= len(dna):
                continue
            dists.append([i, g, last_row_k1[i] + last_row_k2[i+g]])
            continue
        max_dist = max(dists, key = lambda x: x[2])
        if max_dist[2] > best_sum:
            best_pos = max_dist[0]
            best_gap = max_dist[1]
            best_sum = max_dist[2]

    best_k1_pos = best_pos - len(text_k1) # account for extra node, text
    best_new_k1 = dna[best_k1_pos : best_k1_pos+len(text_k1)]
    best_k2_pos = best_pos + best_gap
    best_new_k2 = dna[best_k2_pos : best_k2_pos+len(text_k2)]
    return([[best_k1_pos, best_new_k1], best_gap, [best_k2_pos, best_new_k2]])

if __name__ == '__main__':
    #k1 = 'ATCG'
    #k2 = 'GAA'
    #dna = 'GCCATCAGTTCAGAGTCC'
    k1 = 'AGCCT'
    k2 = 'TAAT'
    dna = 'GCCACGCTTTGAAAATCAGTAAT'
    out = find_best_spaced_motif(k1, k2, dna, [4, 5, 6])
    print(out)
