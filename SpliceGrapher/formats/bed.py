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
Module that encapsulates TopHat .bed files that contain
splice junction information.

Example:

  track name=junctions description="TopHat junctions"
  mitochondria	361439	363002	JUNC00000001	6	-	361439	363002	255,0,0	2	10,29	0,1534
  Chr2	4074	7764	JUNC00000002	1	-	4074	7764	255,0,0	2	28,8	0,3682
  Chr2	4080	8224	JUNC00000003	5	-	4080	8224	255,0,0	2	33,8	0,4136
  Chr2	4257	9081	JUNC00000004	2	+	4257	9081	255,0,0	2	28,10	0,4814
"""
from SpliceGrapher.shared.utils     import ezopen, getAttribute
from SpliceGrapher.shared.ShortRead import SpliceJunction

HEADER_PFX  = 'track'

# BED Format Columns
CHROM       = 0
CHR_START   = 1
CHR_END     = 2
NAME        = 3
SCORE       = 4
STRAND      = 5
# Display values for browsers
THICK_START = 6
THICK_END   = 7
RGB_VALS    = 8
BLOCKS      = 9
SIZES       = 10
STARTS      = 11

ALL_COLUMNS      = [CHROM, CHR_START, CHR_END, NAME, SCORE, STRAND, THICK_START, THICK_END, RGB_VALS, BLOCKS, SIZES, STARTS]
REQUIRED_COLUMNS = [CHROM, CHR_START, CHR_END]

def loadBedJunctions(bedFile, **args) :
    """Loads splice junction data from a BED file.  Relies on Tophat
    format that extends the BED format."""
    minpos = getAttribute('minpos', 0, **args)
    maxpos = getAttribute('maxpos', 0, **args)

    bedStream = ezopen(bedFile)
    ignore    = bedStream.readline()
    result    = {}
    for line in bedStream :
        s   = line.strip()
        if not s : continue
        rec = BEDRecord(s)
        if rec.startpos() < minpos : continue
        if rec.endpos() > maxpos : continue
        jct       = SpliceJunction(rec.chromosome(), rec.donor(), rec.acceptor(), rec.strand())
        jct.count = rec.count()
        keyStr    = '%s;%s;%d;%d' % (rec.chromosome(), rec.strand(), rec.donor(), rec.acceptor())
        result[keyStr] = jct
    return result

def loadBEDRecords(f, header=True) :
    """Loads BED records from a file and returns them in a list.
    Assumes first line is a header string unless 'header' is set to False."""
    from SpliceGrapher.shared.utils import ezopen
    hstring   = None
    result    = []
    bedStream = ezopen(f)
    for s in bedStream :
        if header and not hstring :
            hstring = s.strip()
        else :
            result.append(BEDRecord(s))
    return result

class BEDRecord(object) :
    """Encapsulates the information in a single line of a TopHat-style BED file."""
    def __init__(self,s) :
        parts       = s.strip().split('\t')
        self.attrs  = {}
        self.header = s.startswith(HEADER_PFX)
        if self.header : return

        for col in ALL_COLUMNS :
            try :
                self.attrs[col] = parts[col]
            except IndexError, ie :
                if col in REQUIRED_COLUMNS :
                    raise ie
                self.attrs[col] = None

        # extract upstream/downstream overlap sizes
        sizes = [int(x) for x in self.attrs[SIZES].split(',')]
        if self.strand() == '+' :
            self.us,self.ds = sizes
        else :
            self.ds,self.us = sizes

    def acceptor(self) :
        """Returns the acceptor site inferred by the record."""
        return self.endpos()-self.ds-1 if self.strand() == '+' else self.startpos()+self.ds

    def chromosome(self) :
        """Returns the chromosome given by the record."""
        return self.attrs[CHROM].lower()

    def __cmp__(self, other) :
        """Permits BED records to be sorted based on their start position in the chromosome."""
        result = int(self.attrs[CHR_START]) - int(other.attrs[CHR_START])
        if result == 0 :
            result = int(self.attrs[CHR_END]) - int(other.attrs[CHR_END])
        return result

    def count(self) :
        """Returns the read depth (coverage) associated with a record."""
        return int(self.attrs[SCORE])

    def donor(self) :
        """Returns the donor site inferred by the record."""
        return self.startpos()+self.us if self.strand() == '+' else self.endpos()- self.us-1

    def endpos(self) :
        """Convenience method that returns the last position for the record."""
        return int(self.attrs[CHR_END])

    def __getitem__(self,key) :
        return self.attrs[key]

    def __hash__(self) :
        return str(self).__hash__()

    def isHeader(self) :
        """Returns true if this is a header record; false otherwise"""
        return self.header

    def startpos(self) :
        """Convenience method that returns the first position for the record."""
        return int(self.attrs[CHR_START])

    def __str__(self) :
        vals = [self.attrs[x] for x in ALL_COLUMNS]
        return '\t'.join([x for x in vals if x is not None])

    def strand(self) :
        """Returns the strand given by the record."""
        return self.attrs[STRAND]
