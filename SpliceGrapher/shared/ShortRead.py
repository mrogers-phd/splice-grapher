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
Module containing classes and methods for handling short-read data.
"""
from SpliceGrapher.formats.fasta import *
from SpliceGrapher.shared.utils  import ezopen, idFactory, getAttribute, ProgressIndicator, commaFormat, MAXINT

import sys

#################################################
#  Constants
#################################################
CHROM_CODE  = 'C'
DEPTH_CODE  = 'D'
KNOWN_CODES = [READ_CODE, JCT_CODE, SPLICE_CODE] = ['R', 'J', 'S']
DEPTH_CODES = [CHROM_CODE, DEPTH_CODE, JCT_CODE]

# Splice junction codes:
SJ_CODES    = [KNOWN_JCT, UNKNOWN_JCT, PREDICTED_JCT, UNLABELED_JCT] = ['K', 'U', 'P', '']

#################################################
#  Methods
#################################################
def depthsHeader(path):
    """Returns the header (chromosome) information from a
    SpliceGrapher depths file.  Note that the chromosome
    information must all come at the start of a file."""
    inStream = ezopen(path)
    result   = {}
    ctr      = 0
    for line in inStream:
        ctr  += 1
        s     = line.strip()
        parts = s.split('\t')
        if parts[0] == CHROM_CODE:
            if len(parts) != 3:
                raise ValueError('** %s has invalid chromosome record at line %d:\n%s\n' % (path, ctr, s))
            result[parts[1]] = int(parts[2])
        else:
            break
    return result

def depthsToClusters(chromosome, depths, **args):
    """
    Given a list of read depths, returns a list of clusters where
    coverage exceeds a given threshold.  To get correct positions
    for a subset of depths within a chromosome, include a reference
    start position.  You may also specify a minimum average depth,
    in which case clusters having a lower average depth will be omitted.
    """
    minDepth  = getAttribute('minDepth', 1, **args)
    minpos    = getAttribute('minpos', 0, **args)
    maxpos    = getAttribute('maxpos', MAXINT, **args)
    reference = getAttribute('reference', 0, **args)
    threshold = getAttribute('threshold', 1, **args)
    verbose   = getAttribute('verbose', False, **args)

    result  = []
    current = None
    maxpos  = min(len(depths), maxpos)

    # Ensure that endpoint clusters are completely captured
    while maxpos > minpos > 0 and depths[minpos] >= threshold:
        minpos -= 1
    while minpos < maxpos < len(depths) and depths[maxpos] >= threshold:
        maxpos += 1

    for i in range(minpos,maxpos):
        if depths[i] >= threshold:
            if not current:
                # start with a 1-nt cluster at a single position
                current = Cluster(chromosome, i+reference, depths[i])
            else:
                current.addPosition(i+reference, depths[i])
        elif current:
            if current.avgDepth() >= minDepth:
                result.append(current)
            current = None

    if current and current.avgDepth() >= minDepth:
        result.append(current)
    return result

def isDepthsFile(f):
    """Returns True if a file is a SpliceGrapher depths file; False otherwise."""
    if type(f) == str:
        if not os.path.isfile(f):
            return False
        inStream  = ezopen(f)
        firstLine = inStream.readline()
        inStream.close()
    elif type(f) == file:
        inStream = f
        firstLine = inStream.readline()
        f.seek(0)
    else:
        return False

    parts = firstLine.strip().split('\t')
    return (parts and parts[0] in DEPTH_CODES)

def readDepths(f, **args):
    """Loads read depth data from a SpliceGrapher depth file.
    Returns a dictionary of read depth arrays, one per chromosome,
    along with a dictionary of splice junctions."""
    limit        = getAttribute('maxpos', MAXINT, **args)
    minanchor    = getAttribute('minanchor', 0, **args)
    minjct       = getAttribute('minjct', 1, **args)
    useDepths    = getAttribute('depths', True, **args)
    useJunctions = getAttribute('junctions', True, **args)
    verbose      = getAttribute('verbose', False, **args)

    d      = {}
    jcts   = {}

    # NB: limit is global position limit;
    #     maxpos is chromosome-specific
    maxpos = limit

    # Loading into buffer all at once is 18-30%
    # faster than doing millions of reads, and not
    # too memory-intensive even for large genomes
    lines     = ezopen(f).readlines()
    indicator = ProgressIndicator(1000000, verbose=verbose)
    for line in lines:
        indicator.update()
        s       = line.strip()
        parts   = s.split('\t')
        if len(parts) < 3:
            raise ValueError('Bad depths record at line %d:\n%s' % (indicator.ctr, line))
        recType = parts[0]
        c       = parts[1]
        if recType == CHROM_CODE:
            # C	chr1	123456789
            if not useDepths:
                continue
            # maxpos may change for each chromosome
            maxpos = min(limit, int(parts[2]))
            d[c]   = [0]*maxpos
        elif recType == JCT_CODE:
            if not useJunctions:
                continue
            jct = stringToJunction(s)
            if jct.count < minjct:
                continue
            if jct.minAnchor() < minanchor:
                continue
            if jct.minpos > maxpos:
                continue
            jcts.setdefault(c,[])
            jcts[c].append(jct)
        else:
            if not useDepths:
                continue
            if c not in d:
                raise ValueError('No chromosome information specified for %s depths' % c)

            # D	chr1	1234:0,66:1,...
            rlString = parts[2]

            # Initial implementation uses split();
            # later see if it's faster/less memory
            # to do direct string indexing:
            rparts = rlString.split(',')
            pos    = 0
            for rpair in rparts:
                [rlen,height] = [int(x) for x in rpair.split(':')]
                ubound = min(limit, pos+rlen)
                d[c][pos:ubound] = [height]*rlen
                pos += rlen
                if pos >= limit:
                    break

    indicator.finish()
    return (d,jcts)

def stringToJunction(s):
    """Converts a tab-delimited representation created by
    SpliceJunction.toString() into a SpliceJunction record."""
    parts = s.split('\t')
    if len(parts) != 9:
        raise ValueError('Invalid SpliceJunction record has %d columns (expected 9)' % len(parts))

    if parts[0] != JCT_CODE:
        raise ValueError('Invalid SpliceJunction code: %s' % parts[0])

    if parts[7] not in SJ_CODES:
        raise ValueError('Invalid SpliceJunction type: %s' % parts[7])

    chromosome = parts[1]
    strand     = parts[2]
    intParts   = [int(x) for x in parts[3:7]]
    p1         = min(intParts[0],intParts[1])
    p2         = max(intParts[0],intParts[1])
    anchors    = (intParts[2], intParts[3])
    sjCode     = parts[7]

    jct        = SpliceJunction(chromosome, p1, p2, anchors, sjCode, strand)
    jct.count  = int(parts[8])
    return jct

def writeDepths(ostr, depthDict, jctDict={}, verbose=False):
    """Writes read depths to a SpliceGrapher depth file.  Calling method
    must provide an output destination (file path or writeable stream),
    and read depths stored as a dictionary of chromosome ids mapped
    to lists of integer values."""
    if type(ostr) == file:
        outStream = ostr
    elif type(ostr) == str:
        outStream = open(ostr, 'w')
    else:
        raise ValueError('Unrecognized file type: %s' % type(ostr))

    indicator = ProgressIndicator(1000000, verbose=verbose)

    # Write chromosome details as header:
    for c in sorted(depthDict.keys()):
        indicator.update()
        depths = depthDict[c]
        outStream.write('%s\t%s\t%d\n' % (CHROM_CODE, c, len(depths)))

    for c in sorted(depthDict.keys()):
        outStream.write('%s\t%s\t' % (DEPTH_CODE, c))
        depths = depthDict[c]
        height = -1
        rlen   = 0
        # Run-length encoded depths
        for i in range(len(depths)):
            if depths[i] != height:
                if rlen > 0:
                    indicator.update()
                    outStream.write('%d:%d,' % (rlen,height))
                height = depths[i]
                rlen = 0
            rlen += 1
            i    += 1
        if rlen > 0:
            outStream.write('%d:%d\n' % (rlen,height))

        if c not in jctDict:
            continue
        junctions = jctDict[c]
        junctions.sort()
        for j in junctions:
            indicator.update()
            outStream.write('%s\n' % j.toString())
    indicator.finish()

#################################################
#  Classes
#################################################
class Read(object):
    """
    Class that encapsulates non-spliced reads.
    """
    # Static unique cluster id generator
    id_gen = idFactory('r_')

    def __init__(self, chromosome, p1, p2, strand):
        self.chromosome = chromosome.lower()

        # NB: enforces pi < pj whenever i < j
        self.p1     = min(p1,p2)
        self.p2     = max(p1,p2)
        self.minpos = self.p1
        self.maxpos = self.p2

        self.strand = strand
        self.count  = 1
        self.code   = READ_CODE
        self.id     = next(Read.id_gen)

    def __eq__(self, other):
        return self.chromosome == other.chromosome and self.strand == other.strand \
                and self.p1 == other.p1 and self.p2 == other.p2

    def __cmp__(self, other):
        if self.chromosome < other.chromosome:
            return -1
        elif self.minpos == other.minpos:
            return self.maxpos - other.maxpos
        else:
            return self.minpos - other.minpos

    def __hash__(self):
        return self.__str__().__hash__()

    def __len__(self):
        return self.p2 - self.p1 + 1

    def __repr__(self):
        return "%d reads on %s from %d to %d (%s)" % (self.count, self.chromosome, self.p1, self.p2, self.strand)

    def __str__(self):
        return "%s %d-%d (%s)" % (self.chromosome, self.p1, self.p2, self.strand)

    def contains(self, pos):
        return self.minpos <= pos <= self.maxpos

    def overlaps(self, other):
        """
        Returns true if the two objects overlap; false otherwise.
        """
        return self.chromosome == other.chromosome and self.minpos <= other.maxpos and other.minpos <= self.maxpos

    def sitesMatch(self, other):
        return (other.chromosome == self.chromosome and other.strand == self.strand \
                and other.accval == self.accval and other.donval == self.donval)

    def uniqueString(self):
        return self.id

class ReadPair(object):
    """
    Class that encapsulates paired-end reads.
    """
    def __init__(self, r1, r2):
        self.r1     = r1
        self.r2     = r2
        self.minpos = min(self.r1.minpos, self.r2.minpos)
        self.maxpos = max(self.r1.maxpos, self.r2.maxpos)
        self.name   = self.r1.id

    def __str__(self):
        return 'PE 1=%s, 2=%s' % (self.r1, self.r2)

# Spliced read has two "exon" segments
class SplicedRead(Read):
    """
    Class that encapsulates a spliced read.  This tracks the donor/acceptor sites as
    well as the read's position across the junction.  Used for tracking read depths.
    """
    def __init__(self, chromosome, p1, p2, p3, p4, sjCode, strand):
        Read.__init__(self, chromosome, p1, p2, strand)

        # NB: enforce pi < pj whenever i < j
        self.p3 = min(p3, p4)
        self.p4 = max(p3, p4)
        if self.p1 > self.p3:
            (self.p1, self.p2, self.p3, self.p4) = (self.p3, self.p4, self.p1, self.p2)
        self.minpos = self.p1
        self.maxpos = self.p4

        self.sjCode = sjCode
        self.code   = SPLICE_CODE
        self.accval = self.p3 if self.strand == '+' else self.p2
        self.donval = self.p2 if self.strand == '+' else self.p3
        self.hi     = self.p4-self.p3+1
        self.lo     = self.p2-self.p1+1

    def __eq__(self, other):
        return self.chromosome == other.chromosome and self.strand == other.strand \
                and self.p1 == other.p1 and self.p2 == other.p2 \
                and self.p3 == other.p3 and self.p4 == other.p4

    def __len__(self):
        return self.lo + self.hi

    def __repr__(self):
        return "%d spliced reads on %s from %d-%d to %d-%d (%s)" % (self.count, self.chromosome, self.p1, self.p2, self.p3, self.p4, self.strand)

    def __str__(self):
        return "%s [%d-%d,%d-%d] (%s)" % (self.chromosome, self.p1, self.p2, self.p3, self.p4, self.strand)
            
    def acceptor(self):
        return self.accval

    def donor(self):
        return self.donval

    def highOverlap(self):
        """Returns the overlap size on the high end."""
        return self.hi

    def lowOverlap(self):
        """Returns the overlap size on the low end."""
        return self.lo

# Splice junction has an acceptor and a donor site
class SpliceJunction(Read):
    """
    Class that encapsulates a splice junction.  This tracks only the donor/acceptor
    sites and not the read's position across the junction.  Used for tracking the number
    of reads that map to a junction.  (Note that the SpliceJunction count will be at
    least as high as the counts for the SplicedReads that cross it.)
    """
    def __init__(self, chromosome, p1, p2, anchors, sjCode, strand):
        Read.__init__(self, chromosome, p1, p2, strand)
        self.code    = JCT_CODE
        self.sjCode  = sjCode
        self.anchors = anchors
        self.accval  = self.maxpos if self.strand == '+' else self.minpos
        self.donval  = self.minpos if self.strand == '+' else self.maxpos

    def __repr__(self):
        return "Junction: %d reads on %s between donor %d and acceptor %d (%s)" % \
                (self.count, self.chromosome, self.donval, self.accval, self.strand)

    def __str__(self):
        return "%s don=%d, acc=%d (%s)" % (self.chromosome, self.donval, self.accval, self.strand)

    def acceptor(self, strand=None):
        """Returns the acceptor site for the junction based on the strand."""
        if not strand:
            return self.accval
        return self.maxpos if strand == '+' else self.minpos

    def donor(self, strand=None):
        """Returns the donor site for the junction based on the strand."""
        if not strand:
            return self.donval
        return self.minpos if strand == '+' else self.maxpos

    def isPredicted(self, requiredReads, requiredOverlap):
        """
        Returns true if the splice junction has the required read coverage
        and the required overlap on either side; false otherwise.
        """
        return self.count >= requiredReads and min(self.anchors) >= requiredOverlap

    def minAnchor(self):
        """Returns the smaller of the two anchor values."""
        return min(self.anchors)

    def toString(self):
        """Tab-delimited representation of a junction.  Columns are:
          code, chromosome, strand, minpos, maxpos, anchor1, anchor2, K/P/U
        """
        return '%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%d' % (JCT_CODE, self.chromosome, self.strand, self.minpos, self.maxpos, self.anchors[0], self.anchors[1], self.sjCode, self.count)

    def update(self, o):
        """Updates read count and anchor information using the given SplicedRead or SpliceJunction record."""
        if o.chromosome != self.chromosome or o.strand != self.strand \
                or o.accval != self.accval or o.donval != self.donval:
            raise Exception('Splice sites do not match: %s <--> %s' % (o,self))

        if o.code == SPLICE_CODE:
            self.count     += o.count
            self.anchors[0] = max(self.anchors[0], o.lo)
            self.anchors[1] = max(self.anchors[1], o.hi)
        elif o.code == JCT_CODE:
            self.count     += o.count
            self.anchors[0] = max(self.anchors[0], o.anchors[0])
            self.anchors[1] = max(self.anchors[1], o.anchors[1])
        else:
            raise ValueError('Attempted to update junction with object that was not a SplicedRead or SpliceJunction')


# Cluster is a collection of reads
class Cluster(object):
    """
    A cluster summarizes the information for overlapping reads that map to the same
    gene region.  Unlike reads, clusters make no distinction about strand: reads on
    both strands map to the same cluster.
    """
    # Static unique cluster id generator
    id_gen = idFactory('c_')

    def __init__(self, chromosome, pos, totDepth, id=None):
        self.id          = id if id else next(Cluster.id_gen)
        self.chromosome  = chromosome.lower()
        self.minpos      = pos
        self.maxpos      = pos
        self.totDepth    = int(totDepth)
        self.depths      = {}
        self.depths[pos] = totDepth
        #assert(self.minpos <= self.maxpos)

    def __cmp__(self, other):
        return self.minpos-other.minpos if self.minpos != other.minpos else self.maxpos-other.maxpos

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.__str__().__hash__()

    def __len__(self):
        return self.maxpos-self.minpos+1

    def __repr__(self):
        return "%s: %s %d-%d (%d nt) avg. depth %.2f" % (self.id, self.chromosome, self.minpos, self.maxpos, len(self), self.avgDepth())

    def __str__(self):
        result = "%s: %s %d-%d" % (self.id, self.chromosome, self.minpos, self.maxpos)
        return result

    def addPosition(self, pos, depth):
        """Allows read depth to be updated one nt at a time."""
        self.minpos      = min(self.minpos, pos)
        self.maxpos      = max(self.maxpos, pos)
        self.depths[pos] = self.depths.setdefault(pos,0) + depth
        #assert(self.minpos <= self.maxpos)

    def avgDepth(self, minpos=0, maxpos=MAXINT):
        """Returns the average read depth between the start and end of the cluster,
        or between two points within the cluster."""
        start = max(minpos, self.minpos)
        end   = min(maxpos, self.maxpos)
        dist  = end-start+1
        if dist <= 0:
            return 0
        total = sum([self.depths[k] for k in list(self.depths.keys()) if start <= k <= end])
        return float(total)/dist

    def contains(self, pos):
        return self.minpos <= pos <= self.maxpos

    def overlaps(self, other):
        """Returns true if the two clusters overlap; false otherwise."""
        return self.chromosome == other.chromosome and self.minpos <= other.maxpos and other.minpos <= self.maxpos

    def ssString(self, site):
        """Convenience method for formatting a splice site string. """
        return '%d' % site if site else 'None'
