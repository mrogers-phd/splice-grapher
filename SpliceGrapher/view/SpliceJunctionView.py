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
Module that adds a splice junction view to a matplotlib figure.
"""
from SpliceGrapher.shared.utils      import *
from SpliceGrapher.shared.ShortRead  import KNOWN_JCT, UNKNOWN_JCT, PREDICTED_JCT, UNLABELED_JCT

from matplotlib import patches
from matplotlib.path import Path
from sys        import maxint as MAXINT
import sys

# Color coding:
NOVEL_COLOR = '#00AA88'
LABELS = {KNOWN_JCT:'Known Junction', UNKNOWN_JCT:'Novel Junction', PREDICTED_JCT:'Novel Junction', UNLABELED_JCT:''}
COLORS = {UNLABELED_JCT:'grey', KNOWN_JCT:'grey', UNKNOWN_JCT:NOVEL_COLOR, PREDICTED_JCT:NOVEL_COLOR}
NOVEL_TYPES = [UNKNOWN_JCT, PREDICTED_JCT]

# Colors for showing predicted sites
DONOR_COLOR    = '#D43D1A' # Red = exon end
ACCEPTOR_COLOR = '#3DD41A' # Green = exon start


# Scale for Y-axis (arbitrary)
Y_LIMIT       = 100.0
TRACK_PADDING = 2.8

class SpliceJunctionView(object) :

    def __init__(self, jctList, axis, **args) :
        self.axis      = axis
        self.junctions = jctList
        self.margin    = getAttribute('margin', 16, **args)
        self.yLimit    = getAttribute('yLimit', Y_LIMIT, **args)
        self.verbose   = getAttribute('verbose', False, **args)
        self.overlaps  = {}
        self.patchDict = {}
        self.maxHeight = 1
        self.setOverlaps()

    def drawJunction(self, jct, adjustedY, boxHeight, color, **args) :
        """Draws the actual junction and returns a patch for use in legends."""
        depths  = getAttribute('depths', False, **args)
        xLabels = getAttribute('xLabels', False, **args)

        midX    = 0.5 * (jct.minpos+jct.maxpos)
        minX    = jct.minpos-self.margin
        maxX    = jct.maxpos+self.margin

        halfBox = 0.5 * boxHeight
        lowY    = adjustedY-halfBox
        hiY     = adjustedY+halfBox

        ## Commented out for distribution
        ## assert(lowY >= 0)
        ## if hiY > self.yLimit :
        ##     sys.stderr.write('Junction goes past top of plot: %.3f > %.3f\n' % (hiY, self.yLimit))
        ## assert(hiY <= self.yLimit)

        verts   = [ (jct.minpos, adjustedY), (jct.minpos, hiY), (minX, hiY), (minX, lowY), (jct.minpos, lowY), (jct.minpos, adjustedY), # First block
                    (jct.minpos, adjustedY), (midX, hiY), (midX, hiY), (jct.maxpos, adjustedY), # junction "^" edge
                    (jct.maxpos, adjustedY), (jct.maxpos, hiY), (maxX, hiY), (maxX, lowY), (jct.maxpos, lowY), (jct.maxpos, adjustedY), # Second block
                    ]
        codes   = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY, # First block
                   Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO, # junction edge
                   Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY, # Second block
                   ]

        path  = Path(verts, codes)

        # Unfortunately matplotlib legends won't display the patch, just the color.
        # I may have to write my own extention to the matplotlib Legend class to do that.
        patch = patches.PathPatch(path, facecolor=color, color=color, lw=2)
        self.axis.add_patch(patch)
        bbox = dict(facecolor='lightgrey',edgecolor='lightgrey')

        if xLabels :
            # Show the donor/acceptor site locations
            self.axis.text(jct.minpos, hiY, str(jct.minpos), size='xx-small', ha='right', va='center', bbox=bbox)
            self.axis.text(jct.maxpos, lowY, str(jct.maxpos), size='xx-small', ha='left', va='center', bbox=bbox)

        if depths :
            # Show the junction's read depth at the peak of the "^" intron symbol
            self.axis.text(midX, hiY, str(jct.count), size='xx-small', ha='center', va='top', bbox=bbox)
        return patch


    def junctionsOverlap(self, a,b, margin=0) :
        """
        Returns true if two junctions overlap; false otherwise.  Using a margin permits
        a graph to separate vertically junctions that may appear on top of one another.
        """
        a_min = a.minpos-margin
        a_max = a.maxpos+margin
        b_min = b.minpos-margin
        b_max = b.maxpos+margin
        return a_min <= b_min <= a_max or a_min <= b_max <= a_max or b_min <= a_min <= b_max or b_min <= a_max <= b_max

    def plot(self, **args) :
        """
        Plots all the junctions, keeping track of exon height for overlapping junctions.
        """
        verbose   = getAttribute('verbose', self.verbose, **args)
        #showNovel = getAttribute('showNovel', True, **args)
        dmString  = getAttribute('donors', '', **args)
        amString  = getAttribute('acceptors', '', **args)

        # Ensure boxes are no more than 1/5 of graph
        availHeight = max(5, 1+self.maxHeight)
        boxHeight   = self.yLimit/availHeight
        halfBox     = 0.5 * boxHeight

        # Initialize dictionary for keeping track of used depths
        jctHeight   = dict([(j.id,-1) for j in self.junctions])
        maxRangeSet = set(range(availHeight+1))

        # Iterate over different junction types separately
        for jct in self.junctions :
            ##if (jct.sjCode in NOVEL_TYPES) and not showNovel : continue
            color         = COLORS[jct.sjCode]
            usedPositions = set([jctHeight[j.id] for j in self.overlaps[jct.id]])
            openPositions = maxRangeSet - usedPositions
            if not openPositions :
                if self.verbose : sys.stderr.write('Unable to plot junction %s; ranges = %s; used = %s\n' % (jct, maxRangeSet, usedPositions))
                continue

            height    = min(openPositions)
            adjustedY = halfBox + height*boxHeight
            hiY       = adjustedY + halfBox
            if hiY > self.yLimit :
                if self.verbose : sys.stderr.write('Unable to plot junction %s; y-value = %.1f\n' % (jct, hiY))
                continue

            patch = self.drawJunction(jct, adjustedY, boxHeight, color, **args)

            # Add legend labels only for novel junctions:
            if jct.sjCode in NOVEL_TYPES :
                legendLabel = LABELS[jct.sjCode]
                self.patchDict[legendLabel] = patch

            jctHeight[jct.id] = height

        # Markers used to show predicted sites on junction plots:
        if dmString : self.showMarkers(dmString, DONOR_COLOR, name='donor')
        if amString : self.showMarkers(amString, ACCEPTOR_COLOR, name='acceptor')

        self.axis.set_ylim(0.0, self.yLimit)

    def setOverlaps(self) :
        """
        For each junction, creates a list of other junctions that overlap it.  This
        information is used for establishing the appropriate Y height for plotting a node.
        """
        self.overlaps = {}
        self.minpos   = MAXINT
        self.maxpos   = 0
        maxOverlap    = 0

        for jct1 in self.junctions :
            self.minpos = min(self.minpos, jct1.minpos)
            self.maxpos = max(self.maxpos, jct1.maxpos)
            self.overlaps[jct1.id] = []
            for jct2 in self.junctions :
                if jct2.id == jct1.id : continue
                if self.junctionsOverlap(jct1,jct2) :
                    self.overlaps[jct1.id].append(jct2)
            maxOverlap = max(maxOverlap, len(self.overlaps[jct1.id]))

        # Account for when junction A overlaps B and C, but B and C don't overlap
        maxRangeSet = set(range(maxOverlap))
        if not maxRangeSet : # No overlaps for any junctions
            self.maxHeight = 1
        else :
            tmpHeight   = {}.fromkeys(self.junctions,-1)
            for jct in self.junctions :
                usedPositions  = set([tmpHeight[j] for j in self.overlaps[jct.id]])
                openPositions  = maxRangeSet - usedPositions
                if openPositions :
                    tmpHeight[jct] = min(openPositions)
                else : 
                    tmpHeight[jct] = max(maxRangeSet) + 1
                self.maxHeight = max(self.maxHeight, tmpHeight[jct])

        self.minpos -= self.margin
        self.maxpos += self.margin

    def showMarkers(self, s, color, name='splice site') :
        """Plot vertical lines at the locations given in the
        comma-delimited string."""
        if not s : return

        try :
            positions = [int(x) for x in s.split(',')]
        except Exception :
            raise ValueError('Unrecognized list of %s positions: %s' % (name,s))

        sys.stderr.write('  SpliceJunctionView plotting %s sites from 0 to %d on y-axis\n' % (name, self.yLimit))
        for x in positions :
            sys.stderr.write('  SpliceJunctionView plotting %s %d\n' % (name, x))
            self.axis.plot([x,x], [0,self.yLimit], color=color)
