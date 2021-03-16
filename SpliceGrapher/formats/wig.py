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
Module for encapsulating TopHat .wig files that contain
contiguous ranges of read depths.

Columns
    chromA  chromStartA  chromEndA  dataValueA
Example:
    track type=bedGraph name="TopHat - read coverage"
    Chr1	0	111925	0
    Chr1	111925	111961	5
    Chr1	111961	176472	0
"""
from SpliceGrapher.shared.utils import ezopen

HEADER_PFX = 'track'

# WIG Format Columns
CHROM     = 0
CHR_START = 1
CHR_END   = 2
VALUE     = 3
ALL_COLUMNS = [CHROM, CHR_START, CHR_END, VALUE]

def loadWigDepths(wig_file, minpos, maxpos) :
    """
    Loads read depth values from a Tophat .wig output file.
    """
    depths = {}
    wigStream = ezopen(wig_file)
    ignore    = wigStream.readline() # skip header
    for s in wigStream :
        rec = WIGRecord(s)
        if rec.endpos() < minpos : continue
        if rec.startpos() > maxpos : break
        for i in range(rec.startpos(), rec.endpos()+1) :
            if minpos <= i <= maxpos :
                depths[i] = (rec.value(),0)
    wigStream.close()
    return depths

def loadWIGRecords(f, header=True) :
    """
    Loads WIG records from a file and returns them in a list.
    Assumes first line is a header string unless 'header' is set to False.
    """
    hstring   = None
    result    = []
    wigStream = ezopen(f)
    for s in wigStream :
        if header and not hstring :
            hstring = s.strip()
        else :
            result.append(WIGRecord(s))
    return result

class WIGRecord(object) :
    def __init__(self,s) :
        parts      = s.strip().split('\t')
        self.attrs = {}
        for col in ALL_COLUMNS :
            try :
                self.attrs[col] = parts[col]
            except IndexError as ie :
                raise ie

    def chromosome(self) :
        return self.attrs[CHROM].lower()

    def __cmp__(self, other) :
        result = int(self.attrs[CHR_START]) - int(other.attrs[CHR_START])
        if result == 0 :
            result = int(self.attrs[CHR_END]) - int(other.attrs[CHR_END])
        return result

    def endpos(self) :
        return int(self.attrs[CHR_END])

    def __getitem__(self,key) :
        return self.attrs[key]

    # Required for creating dicts/sets
    def __hash__(self) :
        return str(self).__hash__()

    def startpos(self) :
        return int(self.attrs[CHR_START])

    def __str__(self) :
        return '\t'.join([str(self.attrs[x]) for x in ALL_COLUMNS])

    def value(self) :
        return int(self.attrs[VALUE])
