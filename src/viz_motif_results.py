#!/usr/bin/env python

"""Script to simulate strings, DNA, each with implanted motif separated
by variable gap length.
"""

import argparse
import collections
import os
import random
import sys


import matplotlib as mpl
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
from motif import Profile
from gibs_gappedMotif_alignment import profileWithPseudoCounts
from math import log2
from matplotlib import gridspec


fp = FontProperties(family="Arial", weight="bold") 
globscale = 1.34
LETTERS = { "T" : TextPath((-0.305, 0), "T", size = 1, prop = fp),
            "G" : TextPath((-0.384, 0), "G", size = 1, prop = fp),
            "A" : TextPath((-0.35, 0), "A", size = 1, prop = fp),
            "C" : TextPath((-0.366, 0), "C", size = 1, prop = fp) }
COLOR_SCHEME = {'G': 'orange', 
                'A': 'red', 
                'C': 'blue', 
                'T': 'darkgreen'}

def letterAt(letter, x, y, yscale = 1, ax = None):
    text = LETTERS[letter]

    t = mpl.transforms.Affine2D().scale(1 * globscale, yscale * globscale) + \
        mpl.transforms.Affine2D().translate(x, y) + ax.transData
    p = PathPatch(text, lw = 0, fc = COLOR_SCHEME[letter],  transform = t)
    if ax != None:
        ax.add_artist(p)
    return p

def entropy_normalize(profile):
    pos_entropy = list()
    for col_ind in range(len(profile[0])):
        curr_entropy = 0
        for row_ind in range(len(profile)):
            val = profile[row_ind][col_ind]
            if val != 0:
                curr_entropy += - log2(val) * val
        pos_entropy.append(curr_entropy)
    
    output_profile = [[], [], [], []]
    for i, norm_val in enumerate(pos_entropy):
        for nuc_id in range(len(profile)):
            output_profile[nuc_id].append(profile[nuc_id][i] / norm_val)
    return output_profile


def plot_motif(kmer_motif, ax):
    
    x = 1
    for ind in range(len(kmer_motif[0])):
        y = 0
        for i, nuc in enumerate("ACGT"):
            score = kmer_motif[i][ind]
            if score < 0.05:
                continue
            letterAt(nuc, x, y, score, ax)
            y += score
        x += 1

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_xticklabels([])
    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    ax.set_xlim((.5, x)) 
    ax.set_ylim((0, 2.5))
    ax.set_ylabel("bits")
    
    return ax


def plot_gap_distribution(gap_dist, allowed_gaps, ax):
    
    ax.hist(gap_dist, align = "left", width = 0.14)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    x_range = list(range(min(allowed_gaps) - 1, max(allowed_gaps) + 2))
    ax.set_xticks([])
    ax.set_xlim((min(allowed_gaps) - 1.5, max(allowed_gaps) + 1.5))
    ax.set_ylabel("Gap Count")
    
    return ax


def plot_motif_dist_motif(str_fname, gap_counts_fname, log_fname):
    
    with open(str_fname) as f:
        lines = [x.strip().split("\t") for x in f.readlines()][1:]
        
    with open(gap_counts_fname) as f:
        gap_distribution = [int(x.strip()) for x in f.readlines()[1:]]

    with open(log_fname) as f:
        allowed_gaps = list(f.readlines()[2].split(": ")[1])
    
    k1mers = [line[0] for line in lines]
    k2mers = [line[1] for line in lines]

    k1mer_profile = profileWithPseudoCounts(k1mers)
    k2mer_profile = profileWithPseudoCounts(k2mers)

    k1mer_entropy_norm = entropy_normalize(k1mer_profile)
    k2mer_entropy_norm = entropy_normalize(k2mer_profile)

    gs = gridspec.GridSpec(1, 15)

    fig = plt.figure(figsize = (18, 4))
    
    ax1 = fig.add_subplot(gs[0, 0:5])

    ax1 = plot_motif(k1mer_entropy_norm, ax1)

    ax2 = fig.add_subplot(gs[0, 6:9])

    ax2 = plot_gap_distribution(gap_distribution, [12], ax2)

    ax3 = fig.add_subplot(gs[0, 10:])

    ax3 = plot_motif(k2mer_entropy_norm, ax3)

    plt.show()

    return



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Visualize output of gapped motif finding')
    parser.add_argument('--string', metavar = 'string', type = str,
                        help = 'Location of motif_strings.txt file', required = True)
    parser.add_argument('--gap_dist', metavar = 'gap_dist', type = str,
                        help = 'Location of gap_distribution.txt file', required = True)
    parser.add_argument('--log', metavar = 'log', type = str,
                        help = 'location of motif finding log file',
                        required = True)
    args = parser.parse_args()

    plot_motif_dist_motif(args.string, args.gap_dist, args.log)

    print("Done!")