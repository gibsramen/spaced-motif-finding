#!/usr/bin/env python

import argparse
import os
import sys

def create_parser():
    parser = argparse.ArgumentParser(
            description = 'Evaluate performance of algorithm')
    parser.add_argument('--truth', type = str, help='truth directory',
            required = True)
    parser.add_argument('--result', type = str, help='result directory',
            required = True)
    parser.add_argument('--output', type = str, help='output location',
            required = True)
    return parser

if __name__ == '__main__':
    args = create_parser().parse_args()
    data_dir = args.truth
    result_dir = args.result
    out_dir = args.output

    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    motif_lengths_file = '%s/motif_truth.txt' % data_dir
    with open(motif_lengths_file, 'r') as f:
        next(f)
        kmers = f.read().split('\t')

    k1 = len(kmers[0])
    k2 = len(kmers[1])
    total_nucl = k1 + k2

    truth_indices_file = '%s/motif_indices.txt' % data_dir
    result_indices_file = '%s/motif_indices.txt' % result_dir

    with open(truth_indices_file, 'r') as f:
        next(f)
        truth_indices = f.read().splitlines()
    truth_indices = [list(map(int, x.split('\t')[0:2])) for x in truth_indices]

    with open(result_indices_file, 'r') as f:
        next(f)
        result_indices = f.read().splitlines()
    result_indices = [list(map(int, x.split('\t')[0:2])) for x in result_indices]

    coverage = []
    for truth_index, result_index in zip(truth_indices, result_indices):
        truth_k1mer_range = set(range(truth_index[0], truth_index[0] + k1))
        truth_k2mer_range = set(range(truth_index[1], truth_index[1] + k2))
        result_k1mer_range = set(range(result_index[0], result_index[0] + k1))
        result_k2mer_range = set(range(result_index[1], result_index[1] + k2))
        k1_intersection = truth_k1mer_range.intersection(result_k1mer_range)
        k2_intersection = truth_k2mer_range.intersection(result_k2mer_range)
        intersecting_nucl = len(k1_intersection) + len(k2_intersection)
        coverage.append(intersecting_nucl / total_nucl)

    coverage_file = '%s/coverage.txt' % out_dir
    with open(coverage_file, 'w+') as f:
        for entry in coverage:
            f.write('%s\n' % round(entry, 3))
    print('Coverage calculated!')
