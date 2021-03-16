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
Module that adds a read depth plot to a matplotlib figure.
"""
from SpliceGrapher.shared.utils import *

# Default threshold for triggering Y-axis log scale
DEFAULT_LOGY_THRESHOLD = 100.0
NORMAL_COLOR    = 'grey'
HIGHLIGHT_COLOR = '#D43D1A'

class ReadDepthView(object) :
    """
    Adds a read depth graph to a plot.
    """
    def __init__(self, depths, axis, **args) :
        if type(depths) != type([]) :
            raise ValueError('Read depths must be given as a list, not %s' % type(depths))
        self.depths   = depths
        self.axis     = axis
        self.xrange   = getAttribute('xrange', [0,len(depths)], **args)
        self.logLimit = getAttribute('loglimit', DEFAULT_LOGY_THRESHOLD, **args)
        if not (type(self.xrange) == type([]) or type(self.xrange) == type(())) :
            raise ValueError('X-range must be given as a tuple or a list, not %s' % type(self.xrange))

    def plot(self, **args) :
        """
        Main method for plotting a read-depth graph.
        """
        hiliteStr = getAttribute('highlight', '', **args)
        verbose   = getAttribute('verbose', False, **args)
        window    = getAttribute('window', 500, **args)
        xLabels   = getAttribute('xLabels', False, **args)
        yLimit    = getAttribute('yLimit', 0.0, **args)
        yTitle    = getAttribute('yTitle', 'Read Depth', **args)

        lastPos = min(len(self.depths), self.xrange[1])
        X       = range(self.xrange[0]-1,lastPos)
        Y       = [self.depths[i] for i in X]
        if not Y : return

        # Clip both ends:
        Y[0]  = 0
        Y[-1] = 0
        maxY  = max(Y)

        # Set log scale on Y axis if max read depth exceeds threshold
        if maxY >= self.logLimit :
            # 0 yields an undefined log value, so use 0.1 instead
            Y = [max(y,0.1) for y in Y]
            self.axis.set_yscale('log', nonposy='clip')
            self.axis.set_ylim(0.3, maxY)
        elif yLimit > 0.0 :
            self.axis.set_ylim(0.0, yLimit)

        #--------------------------
        # Plot the filled graph
        #--------------------------
        self.axis.fill(X,Y, fc=NORMAL_COLOR, ec=NORMAL_COLOR)
        if yTitle : self.axis.set_ylabel(yTitle)

        #--------------------------
        # Plot highlighted regions
        #--------------------------
        # Highlight regions are given by pairs of positions
        # where each pair is hyphenated and pairs are
        # separated by commas:
        #    "123-456" or "123-456,765-987"
        if hiliteStr :
            pairs = hiliteStr.split(',')
            for p in pairs :
                try :
                    duple = tuple([int(k) for k in p.split('-')])
                except ValueError :
                    raise ValueError('Illegal format for highlight region: %s' % hiliteStr)

                if len(duple) != 2 :
                    raise ValueError('Illegal format for highlight region: %s' % hiliteStr)

                Hx     = [i for i in X if min(duple) <= i <= max(duple)]
                Hy     = [self.depths[i] for i in Hx]
                if max(Hy) > 0 :
                    # Highight read depths within region
                    Hy[0]  = 0
                    Hy[-1] = 0
                    if maxY >= self.logLimit :
                        Hy = [max(h,0.1) for h in Hy]
                    self.axis.fill(Hx,Hy, fc=HIGHLIGHT_COLOR, ec=HIGHLIGHT_COLOR)
                else :
                    # Highight region is between clusters; just show boundaries
                    lo = min(Y)
                    self.axis.plot([Hx[0],Hx[0]], [lo,maxY], color=HIGHLIGHT_COLOR)
                    self.axis.plot([Hx[-1],Hx[-1]], [lo,maxY], color=HIGHLIGHT_COLOR)

        #--------------------------
        # Show positions where read depth clusters start and end:
        #--------------------------
        if xLabels :
            bbox = dict(facecolor='lightgrey',edgecolor='lightgrey')
            for i in range(1,len(Y)-1) :
                if Y[i] > 0 and Y[i-1] == 0 :
                    self.axis.text(X[i], maxY, '%d'%X[i], size='xx-small', ha='center', va='center', bbox=bbox)
                elif Y[i] == 0 and Y[i-1] > 0 :
                    self.axis.text(X[i-1], 1, '%d'%X[i-1], size='xx-small', ha='center', va='center', bbox=bbox)
