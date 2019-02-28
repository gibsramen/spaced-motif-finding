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
