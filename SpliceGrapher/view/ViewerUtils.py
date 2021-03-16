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
Module of viewer utilities used by several different example scripts.
"""
from SpliceGrapher.shared.config           import *
from SpliceGrapher.shared.utils            import *
from SpliceGrapher.view.SpliceGraphView    import SpliceGraphView
from SpliceGrapher.view.IsoformView        import IsoformView
from SpliceGrapher.view.SpliceJunctionView import SpliceJunctionView
from SpliceGrapher.view.ReadDepthView      import ReadDepthView
from SpliceGrapher.view.ClusterView        import ClusterView
from SpliceGrapher.view.XYGraphView        import XYGraphView
from SpliceGrapher                         import SpliceGraph

from pylab import *
from glob  import glob
import sys, os

ORIGINAL_GRAPH  = 'O'
PREDICTED_GRAPH = 'P'
JUNCTION_GRAPH  = 'J'
DEPTH_GRAPH     = 'R'
XY_GRAPH        = 'X'
ALL_DISPLAYS    = [ORIGINAL_GRAPH, PREDICTED_GRAPH, JUNCTION_GRAPH, DEPTH_GRAPH, XY_GRAPH]
DISPLAY_NAME    = {ORIGINAL_GRAPH:'Original Model', PREDICTED_GRAPH:'Predicted Graph', JUNCTION_GRAPH:'Splice Junctions', DEPTH_GRAPH:'Read Depths', XY_GRAPH:'XY Graph'}

def plotClusters(clusters, displayAxes, minpos, maxpos, **args) :
    """
    Renders a cluster plot within the viewing area.
    """
    prefix    = '' if 'prefix' not in args else args['prefix']+' '
    depthView = ClusterView(clusters, displayAxes, range=(minpos,maxpos), **args)
    depthView.plot(**args)
    displayAxes.set_title('%sRead Clusters' % prefix)

def plotIsoforms(spliceGraph, displayAxes, **args) :
    """
    Renders an isoform plot within the viewing area.
    """
    prefix     = ''    if not 'prefix' in args else args['prefix']+' '
    suffix     = ''    if not 'suffix' in args else ' '+args['suffix']
    adjustment = getAttribute('adjust', 0, **args)
    labels     = getAttribute('labels', False, **args)

    if adjustment : spliceGraph.adjust(adjustment)
    geneName = args['geneName'] if ('geneName' in args and args['geneName']) else spliceGraph.getName()
    formView = IsoformView(spliceGraph, displayAxes, **args)
    formView.plot(**args)
    if 'title' in args :
        displayAxes.set_title('%s%s' % (prefix,args['title']))
    else :
        displayAxes.set_title('%sIsoforms for %s (%s) %s' % (prefix, geneName, spliceGraph.strand, suffix))

    return formView.patchDict, formView.extraPatches

def plotReadDepths(depths, displayAxes, minpos, maxpos, **args) :
    """
    Renders a read depth plot within the viewing area.
    """
    depthView = ReadDepthView(depths, displayAxes, range=(minpos,maxpos), **args)
    depthView.plot(**args)

    prefix = '' if 'prefix' not in args else args['prefix']+' '
    if 'title' in args and args['title'] :
        displayAxes.set_title('%s%s' % (prefix,args['title']))
    else :
        displayAxes.set_title('%sRead depths' % prefix)

def plotSpliceGraph(spliceGraph, displayAxes, **args) :
    """
    Renders a splice graph plot within the viewing area.
    """
    prefix     = ''    if not 'prefix' in args else args['prefix']+' '
    suffix     = ''    if not 'suffix' in args else ' '+args['suffix']
    adjustment = getAttribute('adjust', 0, **args)
    labels     = getAttribute('labels', False, **args)

    if adjustment : spliceGraph.adjust(adjustment)
    geneName  = args['geneName'] if ('geneName' in args and args['geneName']) else spliceGraph.getName()
    graphView = SpliceGraphView(spliceGraph, displayAxes, **args)
    graphView.plot(**args)
    if 'title' in args :
        displayAxes.set_title('%s%s' % (prefix,args['title']))
    else :
        displayAxes.set_title('%sSplice Graph for %s (%s) %s' % (prefix, geneName, spliceGraph.strand, suffix))

    return graphView.patchDict, graphView.extraPatches

def plotSpliceJunctions(junctions, displayAxes, minpos, maxpos, **args) :
    """
    Renders a splice junction plot within the given viewing area.
    """
    mindepth = getAttribute('mindepth', 0, **args) 
    title    = getAttribute('title', '', **args) 
    jctList  = [j for j in junctions if j.maxpos >= minpos and j.minpos <= maxpos and j.count >= mindepth]
    jctView  = SpliceJunctionView(jctList, displayAxes, **args)
    jctView.plot(**args)
    if title :
        displayAxes.set_title(title)
    else :
        prefix   = '' if not 'prefix' in args else args['prefix']+' '
        displayAxes.set_title('%sSplice Junctions with Read Support' % (prefix))

    return jctView.patchDict

def plotXYGraph(X, Y, displayAxes, minpos=0, maxpos=sys.maxint, **args) :
    """
    Renders an X-Y graph within the viewing area.  Assumes X values
    are sorted in ascending order.
    """
    graphView = XYGraphView(X,Y, displayAxes, minpos=minpos, maxpos=maxpos, **args)
    graphView.plot(**args)

    prefix = '' if 'prefix' not in args else args['prefix']+' '
    if 'title' in args :
        displayAxes.set_title('%s%s' % (prefix,args['title']))
    else :
        displayAxes.set_title('%sGraph' % prefix)

def setXticks(minpos, maxpos, minTicks=4, maxTicks=10) :
    """
    Overrides the matplotlib default to use absolute (rather than relative) tick
    mark labels when positions exceed ~10^7
    """
    deltas = [1000000, 500000, 100000, 50000, 10000, 5000, 2000, 1000, 500, 200, 100, 50, 20, 10]
    result = range(minpos,maxpos,10)
    for d in deltas :
        firstTick = d + d*int(int(minpos)/d)
        result    = range(firstTick,maxpos,d)
        if len(result) >= minTicks : break

    # Remove every other tick label
    if len(result) > maxTicks :
        revised = [result[i] for i in range(0,len(result),2)]
        result  = revised

    return result
