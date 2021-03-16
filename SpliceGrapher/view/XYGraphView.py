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
Module that adds a simple graph to a plot.
"""
from SpliceGrapher.shared.utils import *
from sys import maxint as MAXINT

BAR_GRAPH     = 'bar'
LINE_GRAPH    = 'line'
DEFAULT_COLOR = 'grey'

class XYGraphView(object) :
    """
    Adds a simple graph to a plot.
    """
    def __init__(self, X, Y, axis, **args) :
        self.X    = X
        self.Y    = Y
        self.axis = axis
        self.minpos = getAttribute('minpos', self.X[0], **args)
        self.maxpos = getAttribute('maxpos', self.X[-1], **args)

    def plot(self, **args) :
        """
        Main method for plotting a read-depth graph.
        """
        plotType   = getAttribute('plottype', BAR_GRAPH, **args)
        color      = getAttribute('color', DEFAULT_COLOR, **args)
        lineStyle  = getAttribute('linestyle', 'dashed', **args)
        lineWidth  = getAttribute('linewidth', 'g-', **args)
        yLabel     = getAttribute('ylabel', 'Y', **args)

        #--------------------------
        # Plot the graph
        #--------------------------
        if plotType == BAR_GRAPH :
            # Don't let matplotlib try to plot a bunch of boxes nobody will see:
            minY      = min(self.Y)
            maxY      = max(self.Y)
            yRange    = maxY-minY
            minHeight = yRange/40.0
            plotX     = []
            plotY     = []
            for i in xrange(len(self.X)) :
                if self.X[i] < self.minpos : continue
                if self.X[i] > self.maxpos : break
                if abs(self.Y[i]) > minHeight :
                    plotX.append(self.X[i])
                    plotY.append(self.Y[i])

            # Leave about 1/4 the distance between X positions as a gap between bars
            delta    = min([self.X[i]-self.X[i-1] for i in xrange(1,len(self.X))])
            barWidth = 2*0.75*delta
            p1 = self.axis.bar(plotX, plotY, barWidth, ec=color, fc=color)
            self.axis.plot([self.minpos,self.maxpos], [0.0,0.0], linestyle='solid', color=color, linewidth=1.0)
        elif plotType == LINE_GRAPH :
            self.axis.plot(self.X,self.Y, linestyle=lineStyle, linewidth=lineWidth, color=color)

        self.axis.set_ylabel(yLabel)
