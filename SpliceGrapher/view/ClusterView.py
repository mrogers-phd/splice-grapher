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
Module for displaying short read clusters within a Matplotlib figure.
"""
from SpliceGrapher.shared.utils import *
from matplotlib                 import patches
from sys import maxsize as MAXINT

HEIGHT_MARGIN = 0.1 # 10%

class ClusterView(object) :
    def __init__(self, clusters, axis, **args) :
        self.axis      = axis
        self.clusters  = sorted(clusters)
        self.xrange    = getAttribute('xrange', (0,MAXINT), **args)

    def plot(self, **args) :
        """
        Plots all the clusters
        """
        labels  = getAttribute('labels', False, **args)
        xLabels = getAttribute('xLabels', False, **args)
        bbox    = dict(facecolor='lightgrey',edgecolor='lightgrey')
        minpos  = min(self.xrange)
        maxpos  = max(self.xrange)

        self.maxHeight = -1
        for cluster in self.clusters :
            if cluster.maxpos < minpos or cluster.minpos > maxpos : continue
            clusterValue   = cluster.avgDepth()
            self.maxHeight = max(self.maxHeight, clusterValue)
            p              = patches.Rectangle((cluster.minpos, 0.0), len(cluster), clusterValue, facecolor='grey', edgecolor='grey')
            self.axis.add_patch(p)

            if xLabels :
                # Show the donor/acceptor site locations
                self.axis.text(cluster.minpos, 0.0, str(cluster.minpos), size='x-small', ha='left', va='center', bbox=bbox)
                self.axis.text(cluster.maxpos, clusterValue, str(cluster.maxpos), size='x-small', ha='right', va='center', bbox=bbox)

            if labels :
                # Show the clusters's id
                middleX = 0.5 * (cluster.minpos+cluster.maxpos)
                middleY = 0.5 * clusterValue
                self.axis.text(middleX, middleY, '%.3f' % clusterValue, color='white', size='small', ha='center', va='center')

        maxY = (1+HEIGHT_MARGIN) * self.maxHeight
        self.axis.set_ylim(0, maxY)
