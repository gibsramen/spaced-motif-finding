#!/usr/bin/env python

"""Script to simulate strings, DNA, each with implanted motif separated
by variable gap length.
"""

import argparse
import itertools
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
    #rand_gap = random_string(random.choice(gaps))
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
                        type = str, help = 'Output location')
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

    out_file = args.output
    out_file_stem = args.output[:args.output.index('.')]
    out_file_truth = out_file_stem + '_truth.txt'

    all_dna = []
    for i in range(args.num_strings):
        k1_mer = random_string(args.k1)
        k2_mer = random_string(args.k2)
        dna = create_dna(k1_mer, k2_mer, args.gaps, args.mut_rate,
                args.length)
        all_dna.append(dna)

    with open(out_file, 'w+') as f:
        for dna in all_dna:
            f.write('%s\n' % dna[0])

    with open(out_file_truth, 'w+') as f:
        for dna in all_dna:
            pos_string = '\t'.join(list(map(str, dna[1])))
            f.write('%s\n' % pos_string)
    print('Files created!')
