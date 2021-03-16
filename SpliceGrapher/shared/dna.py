# Copyright (C) 2010 by Colorado State University
# Contact: Mark Rogers <rogersma@cs.colostate.edu>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307,
# USA.
"""
Data structures and methods used for handling DNA sequences.
"""
# Sequence codes:
#   R = G/A (purine)
#   Y = T/C (pyrimidine)
#   K = G/T (ketone)
#   M = A/C (amino group)
#   S = G/C (strong interaction)
#   W = A/T (weak interaction)
#   B = G/T/C (not A)
#   D = G/A/T (not C)
#   H = A/C/T (not G)
#   V = G/C/A (not T/U)
# Extended list of sequence characters copied from Sircah:
COMPLEMENT = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N', 'Y':'R',
              'R':'Y', 'U':'A', 'S':'S', 'W':'W', 'K':'M', 'M':'K',
              'B':'V', 'D':'H', 'H':'D', 'V':'B',
              'a':'t', 'c':'g', 't':'a', 'g':'c', 'n':'n', 'y':'r',
              'r':'y', 'u':'a', 's':'s', 'w':'w', 'k':'m', 'm':'k',
              'b':'v', 'd':'h', 'h':'d', 'v':'b' } 

def complement(s):
    """
    Returns the complement for the given sequence, retaining upper/lower case integrity.
    """
    try:
        return ''.join([COMPLEMENT[c] for c in s])
    except KeyError, e:
        raise ValueError("DNA sequence error")

def reverseComplement(s):
    """
    Returns the reverse complement for the given sequence.
    """
    return complement(s)[::-1]

def randomDNA(n):
    """
    Returns a string of random DNA characters.
    """
    import random
    DNA    = 'ACGT'
    result = ''
    for i in range(n):
        result += DNA[random.randint(0,3)]
    return result
