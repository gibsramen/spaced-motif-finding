#!/usr/bin/env python

"""Script to simulate strings, DNA, each with implanted motif separated
by variable gap length.
"""

import argparse
import collections
import os
import random
import sys

NUCLS = 'ACGT'

def mutate(text, mut_rate):
    new_text = ''
    for i in range(len(text)):
        rand_val = random.random()
        if rand_val > mut_rate:
            new_text += text[i]
            continue
        else:
            remaining_nucls = NUCLS.replace(text[i],'')
            new_text += random.choice(remaining_nucls)
    return new_text

def random_string(length):
    return ''.join([random.choice(NUCLS) for i in range(length)])

def create_dna(k1_mer, k2_mer, gaps, mut_rate, length):
    rand_gap_length = random.choice(gaps)
    rand_gap = random_string(rand_gap_length)

    mut_k1_mer = mutate(k1_mer, mut_rate)
    mut_k2_mer = mutate(k2_mer, mut_rate)
    motif = mut_k1_mer + rand_gap + mut_k2_mer

    remaining_nucls = length - len(motif)
    breakpt = random.choice(range(remaining_nucls))
    length_1 = breakpt
    length_2 = remaining_nucls - length_1
    rand_string_1 = random_string(length_1)
    rand_string_2 = random_string(length_2)

    dna = rand_string_1 + motif + rand_string_2
    k1_pos = length_1
    gap_pos = length_1 + len(k1_mer)
    k2_pos = gap_pos + rand_gap_length
    end_pos = k2_pos + len(k2_mer)
                
    return dna, (k1_pos, gap_pos, k2_pos, end_pos)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create simulated dataset')
    parser.add_argument('--k1', metavar = 'k1', type = int,
                        help = 'Length of k1 motif', required = True)
    parser.add_argument('--k2', metavar = 'k2', type = int,
                        help = 'Length of k2 motif', required = True)
    parser.add_argument('--gaps', metavar = 'gaps', type = int,
                        nargs = '+', help = 'Possible gap lengths',
                        required = True)
    parser.add_argument('--mut-rate', metavar = 'mut_rate', required = True, 
                        type = float, help = 'Percent mismatches in k1-k2-mer')
    parser.add_argument('--num-strings', metavar = 'num_strings',
                        required = True, type = int, 
                        help = 'Number of dna strings')
    parser.add_argument('--length', metavar = 'length', required = True,
                        type = int, help = 'Length of each dna string')
    parser.add_argument('--output', metavar = 'output', required = True,
                        type = str, help = 'Output directory')
    args = parser.parse_args()

    if args.k1 <= 0:
        print('Invalid value for k1')
        sys.exit()
    if args.k2 <= 0:
        print('Invalid value for k2')
        sys.exit()
    if args.mut_rate > 1 or args.mut_rate < 0:
        print('Invalid mutation rate')
        sys.exit()
    if args.num_strings <= 0:
        print('Invalid number of strings')
        sys.exit()
    if args.length <= 0:
        print('Length must be greater than 0')

    out_dir = args.output
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    all_dna = []
    k1_mer = random_string(args.k1)
    k2_mer = random_string(args.k2)
    for i in range(args.num_strings):
        dna = create_dna(k1_mer, k2_mer, args.gaps, args.mut_rate,
                args.length)
        all_dna.append(dna)

    gap_dist = {}
    for dna_string, positions in all_dna:
        gap = positions[2] - positions[1]
        if gap_dist.get(gap):
            gap_dist[gap] += 1
        else:
            gap_dist[gap] = 1
    for poss_gap in args.gaps:
        if not gap_dist.get(poss_gap):
            gap_dist[poss_gap] = 0

    gap_dist = collections.OrderedDict(sorted(gap_dist.items()))
    gap_dist_labels = gap_dist.keys()
    gap_dist_counts = gap_dist.values()

    # Files to be written:
    # - motif_lengths.txt
    # - motif_indices.txt
    # - motif_matrix.txt
    # - gap_distribution.txt

    motif_lengths_name = '%s/motif_lengths.txt' % out_dir
    motif_indices_name = '%s/motif_indices.txt' % out_dir
    motif_matrix_name = '%s/motif_matrix.txt' % out_dir
    gap_dist_name = '%s/gap_distribution.txt' % out_dir

    with open(motif_lengths_name, 'w+') as f:
        f.write('k1 length = %s\n' % args.k1)
        f.write('k2 length = %s' % args.k2)

    with open(motif_indices_name, 'w+') as f:
        f.write('k1_start\tk2_start\n')
        for entry in all_dna:
            f.write('%s\t%s\n' % (entry[1][0], entry[1][2]))

    with open(motif_matrix_name, 'w+') as f:
        f.write('\t'.join(list(map(str, range(args.k1+args.k2)))))
        f.write('\n')
        for entry in all_dna:
            k1_pos = entry[1][0]
            gap_pos = entry[1][1]
            k2_pos = entry[1][2]
            end_pos = entry[1][3]
            motif1 = entry[0][k1_pos:gap_pos]
            motif2 = entry[0][k2_pos:end_pos]
            f.write('\t'.join(list(motif1 + motif2)))
            f.write('\n')
    
    with open(gap_dist_name, 'w+') as f:
        f.write('%s\n' % '\t'.join(list(map(str, gap_dist_labels))))
        f.write('%s\n' % '\t'.join(list(map(str, gap_dist_counts))))

    print('Files created!')
