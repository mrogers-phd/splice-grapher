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
from SpliceGrapher.shared.utils import *
from SpliceGrapher.shared.ShortRead import SpliceJunction
import math,copy

MIN_INTRON_SIZE = 30

class AdjustableRange(object):
    """Simple class that encapsulates a position range that can be
    used to shift values within the range."""
    def __init__(self, minpos, maxpos):
        self.minpos = minpos
        self.maxpos = maxpos
        self.delta  = 0

    def __cmp__(self,o):
        return self.minpos-o.minpos

    def __len__(self):
        return self.maxpos-self.minpos+1

    def contains(self,pos):
        """Returns true if the given position is within the range; false otherwise."""
        return self.minpos <= pos <= self.maxpos

    def containsFeature(self,o):
        """Returns true if the given feature is within the range; false otherwise."""
        return self.minpos <= o.minpos and o.maxpos <= self.maxpos

    def overlaps(self, o, slack=MIN_INTRON_SIZE):
        """Returns true if the given object overlaps the range; false otherwise."""
        return self.minpos <= o.maxpos+slack and o.minpos <= self.maxpos+slack

    def setShift(self, pos):
        """Use the given position as a reference for shifting values within the range."""
        self.delta = pos - self.minpos

    def shift(self, pos):
        """Return the shifted value for the given position."""
        return pos + self.delta

    def shiftedMax(self):
        """Convenience method that returns the shifted value for the range maximum position."""
        return self.shift(self.maxpos)

    def shiftedMin(self):
        """Convenience method that returns the shifted value for the range minimum position."""
        return self.shift(self.minpos)

    def __str__(self):
        """String representation for debugging."""
        return 'AdjustableRange %d-%d --> %d-%d (delta=%d)' % (self.minpos, self.maxpos, self.shiftedMin(), self.shiftedMax(), self.delta)

def adjustDepths(depths, ranges, **args):
    """Adjusts the values in an array of depths based on their positions
    relative to the given set of ranges.  Only returns depth values within
    those ranges."""
    if not ranges:
        return depths
    result    = [0]*len(depths)
    gaps      = getGaps(ranges)
    minpos    = ranges[0].minpos
    maxpos    = min(len(depths), ranges[-1].maxpos)
    fullRange = range(minpos,maxpos+1)

    # Easy part -- simply shift values within an exon:
    # Create map from positions to associated ranges
    rangeMap = {}.fromkeys(fullRange, None)
    for r in ranges:
        for i in range(r.minpos,r.maxpos+1):
            rangeMap[i] = r

    # Simply shift values that fall into exon ranges 
    for i in range(minpos,maxpos):
        if not rangeMap[i]:
            continue
        result[rangeMap[i].shift(i)] = depths[i]

    # Create map from positions to associated gaps
    gapMap = {}.fromkeys(fullRange, None)
    for g in gaps:
        # Associate gap with range that immediately precedes it:
        g.setShift(rangeMap[g.minpos-1].shift(g.minpos))
        for i in range(g.minpos,g.maxpos+1):
            gapMap[i] = g

    # Shift and scale values within an intron:
    for g in gaps:
        intronSize = adjustIntron(g, **args)
        # Ratio of original to shrunken intron:
        gamma      = float(len(g))/intronSize
        # Each value of the shrunken intron is the average of
        # the values in its corresponding range:
        for i in range(g.shiftedMin(), g.shiftedMin()+intronSize+1):
            a         = g.minpos + int(round(gamma*(i-g.shiftedMin())))
            b         = int(round(a+gamma))
            result[i] = float(sum(depths[a:b+1]))/(b-a+1)
    return result

def adjustGene(g, **args):
    """Given a gene and a list of ranges, returns a new gene
    with exons adjusted according to their overlapping ranges."""
    ranges   = getAttribute('ranges', [], **args)
    verbose  = getAttribute('verbose', False, **args)

    if not ranges:
        ranges = getGeneRanges(g, **args)
    result         = copy.deepcopy(g)
    result.origMin = g.minpos
    result.origMax = g.maxpos
    # Reset the gene boundaries based on min/max exon boundaries
    result.minpos  = g.maxpos
    result.maxpos  = g.minpos

    if verbose:
        sys.stderr.write('Adjusting gene %s (%d-%d) using %d ranges\n' % (g.id, g.minpos, g.maxpos, len(ranges)))
    for n in result.exons + result.cds:
        modified = False
        for r in ranges:
            if r.overlaps(n):
                n.minpos = r.shift(n.minpos)
                n.maxpos = r.shift(n.maxpos)
                modified = True
                break

        if not modified:
            loPos = adjustPosition(n.minpos, ranges)
            hiPos = adjustPosition(n.maxpos, ranges)
            if loPos and hiPos:
                n.minpos = loPos
                n.maxpos = hiPos
                modified = True

        if modified:
            result.minpos = min(result.minpos, n.minpos)
            result.maxpos = max(result.maxpos, n.maxpos)
        else:
            raise ValueError('Error shrinking introns: unable to adjust %s\nMake sure splice graph and gene model are both up to date.' % n)

    if verbose:
        sys.stderr.write('Adjusted gene %s (%d-%d) is now (%d-%d)\n' % (g.id, g.minpos, g.maxpos, result.minpos, result.maxpos))
    return result

def adjustGraph(G, **args):
    """Given a graph and a list of ranges, returns a new graph
    with nodes adjusted according to their overlapping ranges."""
    ranges   = getAttribute('ranges', [], **args)
    verbose  = getAttribute('verbose', False, **args)

    if not ranges:
        ranges = getGraphRanges(G, **args)
    result         = G.duplicate()
    result.origMin = G.minpos
    result.origMax = G.maxpos
    # Reset the graph boundaries based on min/max node boundaries
    result.minpos  = G.maxpos
    result.maxpos  = G.minpos
    if verbose:
        sys.stderr.write('Adjusting graph %s (%d-%d) using %d ranges\n' % (G.getName(), G.minpos, G.maxpos, len(ranges)))
    for n in sorted(result.resolvedNodes()):
        modified = False
        for r in ranges:
            if r.overlaps(n):
                n.update(r.shift(n.minpos),r.shift(n.maxpos))
                modified = True
                break

        if not modified:
            # Novel node falls within an intron in an original graph
            loPos = adjustPosition(n.minpos, ranges)
            hiPos = adjustPosition(n.maxpos, ranges)
            if loPos and hiPos:
                n.update(loPos, hiPos)
                modified = True
            else:
                loPos = adjustPosition(n.minpos, ranges)
                hiPos = adjustPosition(n.maxpos, ranges)
                sys.stderr.write('  low position  %s --> %s\n' % (str(n.minpos), str(loPos)))
                sys.stderr.write('  high position %s --> %s\n' % (str(n.maxpos), str(hiPos)))
                modified = True

        if modified:
            result.minpos = min(result.minpos, n.minpos)
            result.maxpos = max(result.maxpos, n.maxpos)
        else:
            raise ValueError('Error shrinking graph edges; unable to adjust node %s\nMake sure splice graph and gene model are both up to date.' % n)

    if verbose:
        sys.stderr.write('Adjusted graph %s (%d-%d) is now (%d-%d)\n' % (G.getName(), G.minpos, G.maxpos, result.minpos, result.maxpos))
    result.annotate()
    return result

def adjustIntron(gap, **args):
    """Adjusts the intron inferred by a gap by shrinking it using
    a constant or a log scale.  The resulting gap is never longer
    than the original."""
    scaleFactor = getAttribute('scaleFactor', MIN_INTRON_SIZE, **args)
    gapSize     = len(gap)
    try:
        return min(logScale(gapSize,scaleFactor), gapSize)
    except ValueError as ve:
        sys.stderr.write('Unable to rescale intron gap %s size %d\n' % (gap, gapSize))
        raise ve

def adjustJunctions(jcts, ranges, **args):
    """Returns a list of junctions adjusted for the given ranges."""
    verbose  = getAttribute('verbose', False, **args)
    if verbose:
        sys.stderr.write('Adjusting %d junctions using %d ranges\n' % (len(jcts), len(ranges)))
    ranges.sort()
    result = []
    for j in jcts:
        loPos = adjustPosition(j.minpos, ranges)
        hiPos = adjustPosition(j.maxpos, ranges)
        if not (hiPos and loPos):
            if verbose:
                sys.stderr.write('Cannot adjust %s; omitting.\n' % j)
            continue
        newJct        = SpliceJunction(j.chromosome, loPos, hiPos, j.anchors, j.sjCode, j.strand)
        newJct.count  = j.count
        newJct.sjCode = j.sjCode
        result.append(newJct)
    return result

def adjustPosition(pos, ranges):
    """Adjusts a position based on the given ranges, if the position
    falls within a range.  If the position is above the top range,
    it adjusts the position based on the top range.  If the position
    falls below the bottom range, the position is unchanged."""
    for r in ranges:
        if r.minpos <= pos <= r.maxpos:
            return r.shift(pos)

    # See if the position falls between ranges
    loRange = closestRangeBelow(ranges, pos)
    hiRange = closestRangeAbove(ranges, pos)
    if hiRange and loRange:
        loDelta = pos - loRange.maxpos
        hiDelta = hiRange.minpos - pos
        # Proportion of range to use is based on its distance from
        # the node; lower delta corresponds to higher proportion
        loProp  = 1.0 - float(loDelta)/(loDelta+hiDelta)
        hiProp  = 1.0 - loProp
        return int(round(loProp*loRange.shift(pos) + hiProp*hiRange.shift(pos)))
    elif pos < ranges[0].minpos:
        return pos
    elif pos > ranges[-1].maxpos:
        return ranges[-1].shift(pos)

def adjustRanges(rangeList, **args):
    """Adjusts ranges by shrinking the gaps between them
    either by a fixed amount or using a log scale."""
    # Infer gaps between ranges
    gaps = getGaps(rangeList)
    pos  = 0
    j    = 0
    for i in range(len(rangeList)):
        r = rangeList[i]
        if i == 0:
            newc = r
        else:
            gap = gaps[j]
            #assert(gap.maxpos < r.minpos)
            intronSize = adjustIntron(gap)
            r.setShift(pos+intronSize)
            #assert(r.delta <= 0)
            j += 1
        pos = r.shiftedMax()

def closestRangeAbove(ranges, pos, **args):
    """Finds the range in the list whose minimum position is closest to
    the given position.  If a range overlaps the position, that is the
    closest.  Returns None if all ranges are below the position.
    Ranges must be sorted."""
    verbose  = getAttribute('verbose', False, **args)
    if ranges[-1].maxpos < pos:
        return None
    best = ranges[-1]
    for r in ranges[-2::-1]:
        if r.minpos <= pos <= r.maxpos:
            return r
        elif pos < r.minpos:
            best = r
        else:
            break
    return best

def closestRangeBelow(ranges, pos, **args):
    """Finds the range in the list whose maximum position is closest to
    the given position.  If a range overlaps the position, that is the
    closest.  Returns None if all ranges are above the position.
    Ranges must be sorted."""
    verbose  = getAttribute('verbose', False, **args)
    if ranges[0].minpos > pos:
        return None
    best = ranges[0]
    for r in ranges[1:]:
        if r.minpos <= pos <= r.maxpos:
            return r
        elif pos > r.maxpos:
            best = r
        else:
            break
    return best

def getGaps(ranges):
    """Returns adjustable range instances for all the gaps between the ranges
    in the given list."""
    gaps = [AdjustableRange(ranges[i-1].maxpos+1, ranges[i].minpos-1) for i in range(1,len(ranges))]
    #assert(len(gaps) == len(ranges)-1)
    return gaps

def getGeneRanges(g, **args):
    """Given a splice graph, returns a list of adjustable ranges
    that represent all the overlapping nodes in the graph."""
    return getRanges(g.exons+g.cds, **args)

def getGraphRanges(G, **args):
    """Given a splice graph, returns a list of adjustable ranges
    that represent all the overlapping nodes in the graph."""
    unresolved = getAttribute('unresolved', False, **args)
    if unresolved:
        return getRanges(G.nodeDict.values(), **args)
    else:
        return getRanges(G.resolvedNodes(), **args)

def getRanges(exonList, **args):
    """Given a list of exons, returns a list of ranges
    that represent all the overlapping regions."""
    result   = []
    exonList.sort()
    for n in exonList:
        container = None
        for r in result:
            if r.overlaps(n):
                r.minpos  = min(r.minpos, n.minpos)
                r.maxpos  = max(r.maxpos, n.maxpos)
                container = r
                break
        if not container:
            result.append(AdjustableRange(n.minpos,n.maxpos))

    result.sort()
    adjustRanges(result, **args)
    return result

def logScale(gapSize, factor=MIN_INTRON_SIZE):
    """Default scaling function for shrinking gaps."""
    return factor + factor*int(math.log(gapSize))
