#!/usr/bin/env python

"""Module containing class definitions for motif-related objects"""

class SpacedMotif:
    def __init__(self, text, positions=None):
        """Provide 4 positions:
            1) start of k1_mer in text
            2) start of gap in text (end of k1_mer)
            3) start of k2_mer in text
            4) end of k2_mer in text
        """
        if positions is None:
            raise Exception('No positions input')
        if len(positions) < 4:
            raise Exception('Not enough inputs given')
        self._positions = positions
        self._k1_mer = text[positions[0]:positions[1]]
        self._k2_mer = text[positions[2]:positions[3]+1]
        self._gap_length = positions[2] - positions[1]

    def distance(self, other):
        this_motif = self._k1_mer + self._k2_mer
        other_motif = other.k1_mer + other.k2_mer
        return sum([int(x != y) for x, y in zip(this_motif, other_motif)])

    def __repr__(self):
        return '%s-[%i]-%s' % (self._k1_mer, self._gap_length, self._k2_mer)

    @property
    def gap_length(self):
        return self._gap_length

    @property
    def k1_mer(self):
        return self._k1_mer

    @property
    def k2_mer(self):
        return self._k2_mer

    @property
    def k1_mer_pos(self):
        return self._positions[0]

    @property
    def k2_mer_pos(self):
        return self._positions[1]

class Profile:
    def __init__(self, profile_matrix=None):
        """Input a 4 x k matrix where each entry in profile_matrix
        is a list. Lists should correspond to A, C, G, T probability
        distributions.
        """
        if profile_matrix is None:
            raise ValueError('No profile matrix input')
        self._matrix = profile_matrix
        self._a_dist = profile_matrix[0]
        self._c_dist = profile_matrix[1]
        self._g_dist = profile_matrix[2]
        self._t_dist = profile_matrix[3]

    def get_prob(self, nucleotide, position):
        """Given a nucleotide and position in the profile, return
        associated probability.
        """
        nucl_dist_dict = {'A':self._a_dist, 'C':self._c_dist,
                          'G':self._g_dist, 'T':self._t_dist}
        return nucl_dist_dict[nucleotide][position]

    def print(self):
        print('\t'.join(list(map(str, self._a_dist))))
        print('\t'.join(list(map(str, self._c_dist))))
        print('\t'.join(list(map(str, self._g_dist))))
        print('\t'.join(list(map(str, self._t_dist))))

    def reverse(self):
        new_prof = [x[::-1] for x in self._matrix]
        return Profile(new_prof)

    @property
    def matrix(self):
        return self._matrix

    @property
    def length(self):
        return len(self._matrix[0])

    @property
    def a_dist(self):
        return self._a_dist

    @property
    def c_dist(self):
        return self._c_dist

    @property
    def g_dist(self):
        return self._g_dist

    @property
    def t_dist(self):
        return self._t_dist
