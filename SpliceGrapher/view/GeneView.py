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
Module that adds a background view of gene boundaries to a figure.
"""
from SpliceGrapher.shared.utils      import *
from SpliceGrapher.formats.GeneModel import *

from pylab import *
import sys

# Fill attributes:
OPAQUENESS    = 0.1
EXON_COLOR    = '#0000FF'
TP_UTR_COLOR  = '#FF0000'
FP_UTR_COLOR  = '#00FF00'

# Area type mappings:
PATCH_NAMES   = {FP_UTR_TYPE:"5' UTR", TP_UTR_TYPE:"3' UTR", CDS_TYPE:"CDS/Exon", EXON_TYPE:"CDS/Exon"}
FILL_COLORS   = {FP_UTR_TYPE:FP_UTR_COLOR, TP_UTR_TYPE:TP_UTR_COLOR, CDS_TYPE:EXON_COLOR, EXON_TYPE:EXON_COLOR}
AREA_VALUE    = {FP_UTR_TYPE:2, TP_UTR_TYPE:2, CDS_TYPE:1, EXON_TYPE:1, None:0}

# Edge attributes
EDGE_COLOR    = '#777777'
LINE_WIDTH    = 1

class GeneView(object) :
    """Class for adding a depiction of a known gene model's exon regions to a plot."""

    def __init__(self, genes, axis, **args) :
        self.axis  = axis
        self.genes = genes if type(genes) == type([]) else [genes]

    def defineRegions(self) :
        """Assign coverage types to regions in the view area."""
        minpos = min([g.minpos for g in self.genes])
        maxpos = max([g.maxpos for g in self.genes])

        # First define individual positions:
        posRange = range(minpos,maxpos+1)
        area     = {}.fromkeys(posRange,None)
        for gene in self.genes :
            features = gene.exons + gene.cds
            for f in features :
                for i in xrange(f.minpos,f.maxpos+1) :
                    if i not in area :
                        raise KeyError('%s position %d is not within %d-%d\n' % (f, i, minpos, maxpos))
                    if AREA_VALUE[f.featureType] > AREA_VALUE[area[i]] :
                        area[i] = f.featureType

        # Next define whole regions as tuples
        regions = {FP_UTR_TYPE:[], TP_UTR_TYPE:[], CDS_TYPE:[], EXON_TYPE:[]}
        value   = 0
        curType = None
        start   = minpos
        for i in posRange :
            if area[i] == curType : continue
            if curType : regions[curType].append((start,i-1))
            curType = area[i]
            start   = i

        # Complete last region
        if curType : regions[curType].append((start,posRange[-1]))

        return regions

    def plot(self, **args) :
        """Plots the model for a gene on the axis and returns a patch dictionary
        that can be used for plot legend information."""

        alpha   = getAttribute('alpha', OPAQUENESS, **args)
        label   = getAttribute('label', 'Gene Model', **args)
        altFill = getAttribute('altFill', {}, **args)

        colors  = dict(FILL_COLORS)
        if altFill : colors.update(altFill)

        regions = self.defineRegions()
        result  = {}
        for rtype in regions :
            rname     = "Gene Model " + PATCH_NAMES[rtype]
            fillColor = colors[rtype]
            for posRange in regions[rtype] :
                result[rname] = self.axis.axvspan(min(posRange), max(posRange),
                                                  ymin=0.0,
                                                  ymax=1.0,
                                                  ec=EDGE_COLOR,
                                                  fc=fillColor,
                                                  lw=LINE_WIDTH,
                                                  alpha=alpha)
        return result
