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
Module that provides facilities for loading FASTA sequences into a structure
for easy retrieval.  Includes methods for formatting aligned sequences.
"""
from SpliceGrapher.formats.fasta import fasta_itr
from SpliceGrapher.shared        import dna
from SpliceGrapher.shared.utils  import getAttribute

import sys, os, warnings

class FastaLoader(object) :
    """
    Class that loads a complete set of FASTA sequences into a dictionary indexed
    by sequence header.  Provides access to subsequences and their reverse complements.
    """

    def __init__(self, fastaPath, **args) :
        """
        Loads a FASTA sequence file from the given path.
        
        :Parameters:
           'sequenceID'  - specifies a single sequence ID to load in lieu of the whole file
           'sequenceIDs' - specifies a subset of sequence IDs to load
        """
        self.verbose = getAttribute('verbose', False, **args)
        idParser     = getAttribute('parser', self.headerToSequenceID, **args)

        sequenceIDs = None
        if 'sequenceIDs' in args :
            tmpList     = list(args['sequenceIDs'])
            sequenceIDs = [s.lower() for s in tmpList]
        elif 'sequenceID' in args :
            sequenceIDs = [args['sequenceID'].lower()]

        self.sequences      = {}
        self.revcmp         = {}
        self.totalBasePairs = 0

        iter = fasta_itr(fastaPath)
        if sequenceIDs :
            found = {}.fromkeys(sequenceIDs,False)

        if self.verbose : sys.stderr.write('Loading FASTA records from %s\n' % fastaPath)

        for rec in iter :
            seqID = idParser(rec.header)
            seqID = seqID.lower()
            if sequenceIDs and seqID not in sequenceIDs : continue

            self.sequences[seqID.lower()]  = rec.sequence
            self.totalBasePairs   += len(rec.sequence)

            if sequenceIDs :
                found[seqID] = True
                if all(found) : break

        if sequenceIDs :
            if not any(found) :
                raise Exception('None of the sequence IDs (%s) found in %s.' % (','.join(sequenceIDs),fastaPath))
            elif not all(found) :
                missing = [k for k in found if not found[k]]
                warnings.warn('Not all of the sequence IDs found in %s: missing %s' % (fastaPath, ','.join(missing)))
        if self.verbose : sys.stderr.write('Found %d FASTA records in %s\n' % (len(self.sequences), fastaPath))

    def keys(self) :
        return self.sequences.keys()

    def __len__(self) :
        return len(self.sequences)

    def __getitem__(self, key) :
        return self.sequences[key.lower()]

    def formatAlignment(self, est, ref, reverse=False, offset=1, showMatches=True, width=100) :
        """
        Given a pair of aligned sequences, 'est' and 'ref', writes the alignment to
        stdout with position labels at the start of each line.
        """
        ## assert(len(ref) == len(est))
        i = 0
        while i < len(ref) :
            last = min(i+width, len(ref))
            print "%10d: %s" % (i+offset, ref[i:last])
            if showMatches :
                midstr = "%12s" % ' '
                for j in range(i,last) :
                    if ref[j] == est[j] :
                        midstr += "|"
                    else  :
                        midstr += " "
                print midstr

            print "%8d: %s" % (i+offset, est[i:last])
            print ""
            i += width

    def headerToSequenceID(self, header) :
        """Returns the sequence ID from within a FASTA header."""
        parts = header.split()
        return parts[0]

    def sequence(self, sequenceID, reverse=False) :
        """
        Returns the DNA sequence associated with the given sequence ID.
        If 'reverse' is true, it returns the reverse-complement.
        """
        if not sequenceID.lower() in self.sequences :
            raise KeyError("Header %s not found in FASTA sequences." % sequenceID)
        result = self.sequences[sequenceID.lower()]
        if reverse :
            try :
                result = self.revcmp[sequenceID.lower()]
            except KeyError :
                result = self.revcmp[sequenceID.lower()] = dna.reverseComplement(result)
        return result

    def sequenceLength(self, sequenceID) :
        """
        Returns the length of the DNA sequence associated with the given scaffold.
        """
        return len(self.sequence(sequenceID))

    def sequenceString(self, seq, offset=1, width=100) :
        """
        Returns a multi-line representation the given sequence.
        """
        ## assert(type(seq) == str)
        result = ''
        i = 0
        while i < len(seq) :
            last = min((i+width), len(seq))
            if offset > 0 :
                result += "%8d: %s\n" % (i+offset, seq[i:last])
            else :
                result += "%8d: %s\n" % (-offset-i, seq[i:last])
            i += width
        return result

    def showSequenceByID(self, sequenceID, startpos=0, endpos=None, reverse=False, width=100) :
        """
        Writes the sequence given by the ID to stdout, from 'startpos' to 'endpos'.
        """
        if not endpos :
            endpos = self.sequenceLength(sequenceID)
        seq = self.subsequence(sequenceID, startpos, endpos)
        ## assert(type(seq) == type(''))
        if reverse :
            self.showSequence(seq, offset=-endpos, width=width)
        else :
            self.showSequence(seq, offset=startpos, width=width)

    def showSequence(self, seq, offset=1, width=100) :
        """
        Writes a multi-line representation of the given sequence to stdout.
        """
        print self.sequenceString(seq, offset=offset, width=width)

    def subsequence(self, sequenceID, start, end, reverse=False) :
        """
        Returns the subsequence from 'start' to 'end' from within
        a given sequence.  Assumes that 'start' < 'end' unless it is a
        reverse-complement.
        """
        seq = self.sequence(sequenceID)
        if reverse :
            result = seq[start:end+1] if start < end else seq[end:start+1]
            return dna.reverseComplement(result)
        else :
            return seq[start:end+1]

    def totalLength(self) :
        """
        Returns the total length of all sequences currently stored.
        """
        return self.totalBasePairs
