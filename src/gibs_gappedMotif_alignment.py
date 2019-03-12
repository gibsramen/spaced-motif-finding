import sys
import argparse
import numpy as np
import os
from numpy.random import choice
import alignment
from motif import Profile
from collections import Counter
import time


def profileWithPseudoCounts(Dna):
    '''returns positional weight matrix from inputted list of dna'''
    As, Cs, Gs, Ts = [], [], [], []
    # print([len(x) for x in Dna])

    for j in range(len(Dna[0])):  # for each nt position
        frequencies = {'A': 1, 'C': 1, 'G': 1, 'T': 1}
        for i in range(len(Dna)):
            try:
                nt = Dna[i][j]
            except BaseException:
                print(i, j)
            frequencies[nt] += 1
        As += [frequencies['A'] / sum(frequencies.values())]
        Cs += [frequencies['C'] / sum(frequencies.values())]
        Gs += [frequencies['G'] / sum(frequencies.values())]
        Ts += [frequencies['T'] / sum(frequencies.values())]

    return [As, Cs, Gs, Ts]


def scoreMotifs(profile, motifs):
    # get consensus sequence
    consensus = []
    for pA, pC, pG, pT in zip(profile[0], profile[1], profile[2], profile[3]):
        pMax = max([pA, pT, pG, pC])
        if pA == pMax:
            consensus += ['A']
        elif pC == pMax:
            consensus += ['C']
        elif pG == pMax:
            consensus += ['G']
        elif pT == pMax:
            consensus += ['T']

    # get score
    score = 0
    for m in motifs:
        for ntIdx in range(len(motifs[0])):
            if m[ntIdx] != consensus[ntIdx]:
                score += 1

    return score

# # Test scoreMotifs()
# print(scoreMotifs(profileWithPseudoCounts(['AAA','AAA','AAA','AAA']),['AAT','AAA','TAA']))


def kmers(text, k):
    '''returns list of kmers, a sliding window of length k across a string'''
    kmers = []
    for i in range(len(text) - k + 1):
        kmers += [text[i:i + k]]
    return kmers


def gappedMotifGibbsAlignment(Dna, G, k1, k2, n):
    '''
    Gibbs-sampler for gapped composite motif finding from t strings dna.

    Inputs:
            Dna - a collection of dna strings of length t.
            G - a collection of integers corresponding to gap length.
            k1 - integer length of k1 motif.
            k2 - integer length of k2 motif.
            N - an integer specifying number of iterations of sampling to perform.

    Returns (in order):
            bestMotifs - t x k1+k2 matrix. Contains unspaced motif of length k1+k2. Rows correspond to each dna in Dna.
            bestMotifs_indices - t x 2 matrix. col[0] k1 start in dna, col[1] k2 start in dna. Rows correspond to each dna in Dna.
            bestGaps - t x 1 matrix. Contains a gap in each row corresponding to gap observed in gap-free motif of bestMotifs matrix.
            gapDistribution - 2 x G matrix. row[0] contains gap lengths tested. row[1] contains frequency of each gap in G.
    '''

    t = len(Dna)

    # Choose arbitrary motifs to instantiate motifs.
    # Choose first k1-g-k2mer where k1 starts at 0 and gap = min(G).
    # Append unspaced motif to motifs.
    motifs = []
    motifs_indices = []
    gaps = []
    g_min = min(G)

    for d in Dna:
        k1_start_ind = choice(range(0, len(d) - k1 - g_min - k2 - 1))
        if len(G) > 1:
            k2_start_ind = choice(
                range(
                    k1_start_ind +
                    g_min,
                    k1_start_ind +
                    max(G)))
        else:
            k2_start_ind = k1_start_ind + G[0]
        motifs.append(d[k1_start_ind:k1_start_ind + k1] +
                      d[k2_start_ind: k2_start_ind + k2])
        motifs_indices.append(k1_start_ind)
        gaps.append(k2_start_ind - k1_start_ind)

    # Initialize best motifs.
    bestMotifs = motifs[:]
    bestMotif_indices = motifs_indices[:]
    bestGaps = gaps[:]
    score_bestMotifs = scoreMotifs(
        profileWithPseudoCounts(bestMotifs), bestMotifs)

    # Sample n times.
    for j in range(n):

        # Choose random d in Dna to sample.
        i = choice(range(t))

        # print('\n-------------\ni',i)
        # print('\nmotifs')
        # for z in motifs: print('',z)
        # print('\ndna[i]',Dna[i])

        # Create profile of all motifs except ith motif.
        assert all([len(x) == len(motifs[0])
                    for x in motifs]), print(i, motifs)
        motifs_exceptI = [
            m for m_index,
            m in enumerate(motifs) if m_index != i]
        profile_exceptI = profileWithPseudoCounts(motifs_exceptI)

        # print('\nmotifs except i')
        # for z in motifs_exceptI: print('',z)

        # Select the best gapped kmer for the next iteration.
        k1_profile_exceptI = Profile([x[:k1] for x in profile_exceptI])
        k2_profile_exceptI = Profile([x[-k2:] for x in profile_exceptI])
        motifs_indices[i], motifs[i], gaps[i] = alignment.find_best_spaced_motif(
            k1_profile_exceptI,
            k2_profile_exceptI,
            Dna[i],
            G
        )

        assert all([len(x) == len(motifs[0])
                    for x in motifs]), print(i, motifs)
        # print('\n>> New Motifs')
        # for z in motifs: print('',z)

        # Check if new motifs are better than old motifs.
        if scoreMotifs(profileWithPseudoCounts(motifs),
                       motifs) < score_bestMotifs:
            bestMotifs = motifs
            bestMotif_indices = motifs_indices
            bestGaps = gaps
            score_bestMotifs = scoreMotifs(
                profileWithPseudoCounts(bestMotifs), motifs)

    # Calculate frequency of each gap length observed
    gapDist = []
    for g in G:
        gapDist.append(bestGaps.count(g))

    # print(bestMotifs, bestMotif_indices, bestGaps)
    return bestMotifs, bestMotif_indices, bestGaps, [G, gapDist]


def pipeline_gibbs_sampler_with_alignment(Dna, G, k1, k2, num_iters, num_runs):

    run_results = list()
    for _ in range(num_runs):
        print(_)
        run_results.append(
            gappedMotifGibbsAlignment(
                Dna, G, k1, k2, num_iters))

    best_results = [[], [], []]

    for dna_ind in range(len(Dna)):

        curr_poss_matches = [x[0][dna_ind] for x in run_results]

        dna_res_counts = Counter(curr_poss_matches)

        best_string_in_dna = dna_res_counts.most_common(1)[0][0]

        for result in run_results:

            if result[0][dna_ind] == best_string_in_dna:
                best_results[0].append(result[0][dna_ind])
                best_results[1].append(result[1][dna_ind])
                best_results[2].append(result[2][dna_ind])
                break

    return best_results


def get_gap_distribution(gaps, allowed_gaps):

    counts = [0 for _ in allowed_gaps]

    for gap in gaps:
        counts[allowed_gaps.index(gap)] += 1

    return [counts, allowed_gaps]


if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser(
        description='Find gapped composite motifs in a collection of dna strings.')
    parser.add_argument(
        '--k1',
        metavar='k1',
        type=int,
        help='Integer length of k1.',
        required=True)
    parser.add_argument(
        '--k2',
        metavar='k2',
        type=int,
        help='Integer length of k2.',
        required=True)
    parser.add_argument(
        '--gaps',
        metavar='gaps',
        type=int,
        nargs='+',
        help='Valid gaps to search for. Collection of integers (min length 1) delimited by spaces.',
        required=True)
    parser.add_argument(
        '--num_iters',
        metavar='num_iters',
        type=int,
        help='Number of Gibs sampling iterations.',
        required=True)
    parser.add_argument(
        '--num_runs',
        metavar='num_runs',
        type=int,
        help='Number of random restarts for Gibbs sampling.',
        required=True)
    parser.add_argument(
        '--input',
        metavar='input',
        type=str,
        help='A .txt file. Inputted collection Dna of string dna. Each string dna should be on its own line.',
        required=True)
    parser.add_argument(
        '--outputdir',
        metavar='outputdir',
        type=str,
        help='Output /dir/ (do not append with .txt or other suffix).',
        required=True)
    args = parser.parse_args()

    # Read input File.
    Dna = [
        line.strip() for line in open(
            args.input,
            'r').readlines() if line != '\n']

    try:
        os.mkdir(args.outputdir)
    except:
        pass

    logFN = '%s/log.txt' % args.outputdir
    motifLengthsFN = '%s/motif_lengths.txt' % args.outputdir
    bestMotifsFN = '%s/motif_strings.txt' % args.outputdir
    bestMotif_indicesFN = '%s/motif_indices.txt' % args.outputdir
    bestGapsFN = '%s/gap_counts.txt' % args.outputdir
    gapDistributionFN = '%s/gap_distribution.txt' % args.outputdir

    with open(logFN, 'w+') as f:
        f.write('# Inputs and Parameters' + '\n')
        f.write('Dna input: ' + args.input + '\n')
        f.write('gaps: ' + str(args.gaps) + '\n')
        f.write('k1 length: ' + str(args.k1) + '\n')
        f.write('k2 length: ' + str(args.k2) + '\n')
        f.write('Sampling iterations per run: ' + str(args.num_iters) + '\n')
        f.write('Sampling runs used to define consensus: ' +
                str(args.num_runs) + '\n')

    bestMotifs, bestMotif_indices, bestGaps = pipeline_gibbs_sampler_with_alignment(
        Dna,
        args.gaps,
        args.k1,
        args.k2,
        args.num_iters,
        args.num_runs
    )

    gapDistribution = get_gap_distribution(bestGaps, args.gaps)

    # Write Motif Lengths
    with open(motifLengthsFN, 'w+') as f:
        f.write('k1_length\tk2_length\n')
        f.write(str(args.k1) + '\t' + str(args.k2) + '\n')

    # Write bestMotifs
    with open(bestMotifsFN, 'w+') as f:
        f.write('k1_motif\tk2_motif\n')
        out_motifs = [[unspaced[:args.k1], unspaced[args.k1:]] for 
            unspaced in bestMotifs]
        for entry in out_motifs:
            f.write('%s\t%s\n' % (entry[0], entry[1]))

    # Write bestMotif Indices
    with open(bestMotif_indicesFN, 'w+') as f:
        f.write('k1_start\tk2_start\tg\n')
        for k1gk2_index, g in zip(bestMotif_indices, bestGaps):
            f.write('%s\t%s\n' % ('\t'.join(list(map(str, k1gk2_index))), g))

    # Write bestGaps
    with open(bestGapsFN, 'w+') as f:
        f.write('gap' + '\n')
        f.write('\n'.join([str(i) for i in bestGaps][::-1]))

    with open(gapDistributionFN, 'w+') as f:
        f.write('gap\tcounts\n')
        for g, g_freq in zip(gapDistribution[0], gapDistribution[1]):
            f.write(str(g) + '\t' + str(g_freq) + '\n')
    # print('\ngapDistribution',gapDistribution)
    print(time.time() - start)
