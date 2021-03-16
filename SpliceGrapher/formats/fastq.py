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
A parser for FASTQ files.
"""
from SpliceGrapher.shared.utils import ezopen, getAttribute
import os, sys, time

class MalformedInput :
    "Exception raised when the input file does not look like a fastq file."
    pass

class FastqRecord :
    "a fastq record."
    def __eq__(self, other) :
        return self.header == other.header and self.sequence == other.sequence

    def __init__(self, header, sequence, header2=None, quality=None):
        "Create a record with the given header, sequence and quality string"
        self.header   = header
        self.sequence = sequence
        self.header2  = header2
        self.quality  = quality

        if self.header2 is None :
            self.header2 = self.header

        # 'I' should be the maximum quality score for Illumina reads
        if self.quality is None :
            self.quality = 'I' * len(self.sequence)

    def __str__(self) :
        # FASTQ format:
        #   @<header>
        #   sequence
        #   +<header>
        #   quality scores
        return '@' + self.header + '\n' + self.sequence + '\n+' + self.header2 + '\n' + self.quality + '\n'
        
def _fastq_itr_from_file(file) :
    "Provide an iteration through the fastq records in file."

    # Process first record
    line = file.readline()
    if line[0] != '@':
        sys.stderr.write("INITIAL LINE STARTS WITH INVALID CHARACTER: '%s'\n" % line.rstrip())
        raise MalformedInput()

    # NB: assumes sequence is on a single line
    h       = line.rstrip()[1:]
    seq     = file.readline().rstrip()
    h2      = file.readline().rstrip()[1:]
    qual    = file.readline().rstrip()
    linectr = 4
    line    = file.readline()

    while line :
        linectr += 1
        if len(line) < 1 :
            raise Exception("ERROR: empty string at line %d" % linectr)

        if line[0] == '@':
            # Omit 2nd header to save space:
            yield FastqRecord(h, ''.join(seq), '', ''.join(qual))
            h        = line.rstrip()[1:]
            seq      = file.readline().rstrip()
            h2       = file.readline().rstrip()[1:]
            qual     = file.readline().rstrip()
            linectr += 3

        line = file.readline()

    # Omit 2nd header to save space:
    yield FastqRecord(h, ''.join(seq), '', ''.join(qual))

def _fastq_itr_from_name(fname):
    "Provide an iteration through the fastq records in the file named fname. "

    f = ezopen(fname)
    for rec in _fastq_itr_from_file(f) :
        yield rec

def _fastq_itr(src):
    """Provide an iteration through the fastq records in file `src'.
    
    Here `src' can be either a file object or the name of a file.
    """
    if type(src) == str :
        return _fastq_itr_from_name(src)
    elif type(src) == file :
        return _fastq_itr_from_file(src)
    else:
        raise TypeError("Cannot determine input stream type.")

def fastq_get_by_name(itr,name, byLength=False):
    "Return the record in itr with the given name."
    x = name.strip()
    n = len(x)
    for rec in itr:
        s = rec.header.strip()
        if byLength and n < len(s) :
            s = s[:n]
        if s == x:
            return rec
    return None

class fastq_itr (object) :
    "An iterator through a sequence of fastq records."

    def __init__(self,src) :
        "Create an iterator through the records in src."
        self.__itr = _fastq_itr(src)

    def __iter__(self) :
        return self

    def __next__(self) :
        return next(self.__itr)

    def __getitem__(self,name) :
        return fastq_get_by_name(iter(self),name)

def get_sequence(src, name):
    "Return the record in src with the given name."
    return fastq_itr(src)[name]

def fastq_count(src) :
    """
    count the number of records in a fastq file
    """
    return sum([1 for r in fastq_itr(src)])

def phred_prob(k) :
    """
    For a given phred score, returns the probability that a
    character is incorrect.
       k = -10 log(p) ---> p = 10**(-k/10), where 0 <= k <= 40
    """
    expon = -0.10 * k
    return 10**expon

def solexa_prob(k) :
    """
    For a given solexa score, returns the probability that a
    character is incorrect.
       k = -10 log(p/(1-p)) ---> p = 10**t/(10**t + 1), where t = -k/10 and -5 <= k <= 40 
    """
    t   = -0.10 * k
    val = 10.0**t
    return val/(val+1)

def quality_scores(s, **args) :
    """
    Converts a quality code stored in a string to a list of probabilities.
    Assumes phred-type scores unless 'scoretype' is specified.  Values are
    'phred' or 'solexa'.  For more information on quality score types see
    http://maq.sourceforge.net/qual.shtml
    """
    scoreType = getAttribute('scoretype', 'phred', **args)
    if scoreType == 'solexa' :
        return [solexa_prob(ord(c)-64) for c in s]
    elif scoreType == 'phred' :
        return [phred_prob(ord(c)-33) for c in s]
    else :
        raise ValueError('Unrecognized FASTQ quality score type: %s' % scoreType)
