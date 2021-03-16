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
from SpliceGrapher.SpliceGraph  import *

import sys

# Map SpliceGraph AS event abbreviations to their associated count functions
COUNT_FUNCTION = {IR_ABBREV:'irCount', ALT5_ABBREV:'alt5Count', ALT3_ABBREV:'alt3Count',
                  ALTB3_ABBREV:'altBCount', ALTB5_ABBREV:'altBCount', ES_ABBREV:'esCount'}

OPTIONAL_FORMS = set([ALTI_ABBREV, ALTT_ABBREV])

def countFunctionName(asType) :
    """Returns the counting function for the given AS type."""
    key = EVENT_ABBREV[asType] if asType in AS_NAMES else asType
    try :
        return COUNT_FUNCTION[key]
    except KeyError :
        raise ValueError('Unrecognized AS type for summary average: %s is not in %s' % (asType, '/'.join(list(EVENT_ABBREVS.keys()))))

class AbstractStatistics(object) :
    """Base class for all statistics classes."""
    def alt3Count(self) : raise NotImplementedError()
    def alt5Count(self) : raise NotImplementedError() 
    def altBCount(self) : raise NotImplementedError() 
    def edgeCount(self) : raise NotImplementedError()
    def esCount(self)   : raise NotImplementedError() 
    def irCount(self)   : raise NotImplementedError() 
    def nodeCount(self) : raise NotImplementedError()

class GraphStatistics(AbstractStatistics) :
    """Encapsulates statistics for a single graph."""

    def __init__(self, graph, **args) :
        annotate = getAttribute('annotate', False, **args)
        self.graph = graph
        if annotate : self.graph.annotate()
        self.nodes      = graph.resolvedNodes()
        self.alt3Events = altSiteEventList(self.nodes, ALT3_ABBREV)
        self.alt5Events = altSiteEventList(self.nodes, ALT5_ABBREV)

    def alt3Count(self) :
        """Returns the number of alternative 3' events in a graph."""
        return sum([(len(e)-1) for e in self.alt3Events])

    def alt5Count(self) :
        """Returns the number of alternative 5' events in a graph."""
        return sum([(len(e)-1) for e in self.alt5Events])

    def altEventCount(self, eventList) :
        """Returns the number of alternative events in a list
        of alt acceptor or donor events."""
        result = 0
        for event in eventList :
            #assert(len(event) > 0)
            result += len(event)-1
        return result

    def altBCount(self) : 
        """Returns the number of simultaneous 3'/5' (alt. both) events in a graph."""
        # Idea is to identify groups of nodes involved in simultaneous events.
        # A distinct complete event is one that contains two or more 5' splice
        # sites and their 3' splice site mates.
        #     A                    B                            C
        #  ||||||||||--------||||||||||||----------------|||||||||||||||
        #  ||||||||||--------||||||||||||||||||||-------------||||||||||
        #  |||||----------------|||||||||----------------|||||||||||||||
        #
        # The example above contains 2 simultaneous 3'/5' events: one between
        # groups A&B and a second between groups B and C.  To count events, we
        # count the number of connections between distinct groups.
        allAB = set([n for n in self.nodes if ALTB3_ABBREV in n.altForms() or ALTB5_ABBREV in n.altForms()])
        if not allAB : return 0

        # Establish overlap groups:
        groups = []
        while allAB :
            n     = allAB.pop()
            group = set([o for o in allAB if overlap(o,n)])
            group.add(n)
            allAB.difference_update(group)
            groups.append(group)

        # Count inter-group connections by identifying nodes in
        # one group that are children of nodes in the other.
        eventCount = 0
        for i in range(len(groups)-1) :
            for j in range(i+1,len(groups)) :
                A     = groups[i]
                B     = groups[j]
                Akids = set([c for n in A for c in n.children])
                Bkids = set([c for n in B for c in n.children])
                if A.intersection(Bkids) or B.intersection(Akids) :
                    eventCount += 1
        return eventCount

    def branchingFactor(self) :
        """Returns the branching factor for a graph."""
        nNodes = self.nodeCount()
        return float(self.edgeCount())/nNodes if nNodes > 0 else 0.0

    def edgeCount(self) :
        """Returns the number of edges in a graph."""
        # Could use either parents or children
        return sum([len(n.parents) for n in self.nodes])

    def esCount(self)   :
        """Returns the number of skipped exons in a graph."""
        result  = 0
        esNodes = [n for n in self.nodes if ES_ABBREV in n.altForms()]
        if esNodes :
            # Unique ES events are associated with distinct edges:
            intronSet = set([])
            for e in edgeSet(self.graph) :
                intronSet.add((e.minpos,e.maxpos))
            for i in intronSet :
                spannedNodes = [n for n in esNodes if min(i) < n.minpos and n.maxpos < max(i)]
                if spannedNodes : result += 1
        return result

    def hasAS(self, useAll=False) :
        """Returns true if the associated graph has AS; false otherwise."""
        if useAll :
            return self.graph.hasAS()
        else :
            return bool(self.graph.altForms() - OPTIONAL_FORMS)

    def irCount(self)   :
        """Returns the number of exons that represent retained introns."""
        allIR = [n for n in self.nodes if IR_ABBREV in n.altForms()]
        if not allIR : return 0
        edges = edgeSet(self.graph)

        # Map nodes to the distinct edges they contain
        edgeMap = {}
        for n in allIR :
            eset = set([e for e in edges if containsEdge(n,e)])
            if eset : edgeMap[n] = tuple(sorted(eset))

        # Now map edge sets to nodes to identify nodes with same edges contained
        nodeMap = {}
        for n in edgeMap :
            key = edgeMap[n]
            nodeMap.setdefault(key,set([]))
            nodeMap[key].add(n)
        return len(nodeMap)

    def maxBranchingFactor(self) :
        """Returns the maximum branching factor for any node in a graph.
        Considers both child edges and parent edges."""
        maxEdges = 0
        for n in self.nodes :
            maxEdges = max(maxEdges, max(len(n.parents), len(n.children)))
        return maxEdges

    def nodeCount(self) :
        """Returns the number of nodes in a graph."""
        return len(self.nodes)

class SummaryStatistics(AbstractStatistics) :
    """Encapsulates a complete summary for a set of graphs."""

    def __init__(self) :
        self.stats  = {}

    def addGraph(self, g, **args) :
        """Adds a single graph to the summary set.  Raises an exception
        if the instance has already loaded the graph."""
        key = g.name.upper() # Case-insensitive keys
        try :
            oldGraph = self.stats[key]
            raise Exception('Summary statistics instance received a duplicate graph: %s' % g.name)
        except KeyError :
            pass

        if len(g.nodeDict) == 0 :
            raise ValueError('Attempted to add empty graph to summary statistics instance: %s' % g.name)

        self.stats[key] = GraphStatistics(g, **args)

    def addGraphFiles(self, fileList, **args) :
        """Adds graphs from a list of files.  Raises an exception if the
        instance has already loaded any of the graphs."""
        verbose = getAttribute('verbose', False, **args)
        if verbose :
            sys.stderr.write('Loading graphs from %d files:\n' % len(fileList))
            initial   = len(self.stats)

        indicator = ProgressIndicator(10000, verbose=verbose)

        for f in fileList :
            indicator.update()
            try :
                g = getFirstGraph(f)
            except KeyError as ke :
                raise KeyError('Invalid graph in file %s: %s' % (f,str(ke)))
            except ValueError :
                continue
            self.addGraph(g, **args)

        indicator.finish()
        if verbose : sys.stderr.write('Summary statistics added %d/%d graphs.\n' % (len(self.stats)-initial, len(fileList)))

    def addGraphs(self, graphs, **args) :
        """Adds a set, dictionary or a list of graphs to the set.  Raises an
        exception if the instance has already loaded any of the graphs."""
        if type(graphs) == type({}) :
            graphList = list(graphs.values())
        elif type(graphs) == type(set([])) :
            graphList = list(graphs)
        elif type(graphs) == type([]) :
            graphList = graphs
        else :
            raise ValueError('Summary statistics instance received a graph bundle that was not a list, set or dictionary: %s.' % type(graphs))

        verbose = getAttribute('verbose', False, **args)
        if verbose : sys.stderr.write('Loading %d graphs\n' % len(graphs))
        indicator = ProgressIndicator(10000, verbose=verbose)
        for g in graphList :
            indicator.update()
            self.addGraph(g, **args)
        indicator.finish()

    def altCount(self, useAll=False) :
        """Returns the total number of graphs stored that have AS.
        Note: ignores alternative initiating or terminating exons
        unless 'useAll' is set."""
        return len([s for s in list(self.stats.values()) if s.hasAS(useAll=useAll)])

    def avg(self, asType) :
        """Returns the average value for the given AS type."""
        return self.computeAverage(self.total(asType))

    def avgEdges(self) :
        """Returns the average number of edges per graph."""
        return self.computeAverage(self.edgeCount())

    def avgBranchingFactor(self) :
        """Returns the average branching factor per graph."""
        total = sum([s.branchingFactor() for s in list(self.stats.values())])
        return self.computeAverage(total)

    def avgMaxBranchingFactor(self) :
        """Returns the average maximum branching factor per graph."""
        total = sum([s.maxBranchingFactor() for s in list(self.stats.values())])
        return self.computeAverage(total)

    def avgNodes(self) :
        """Returns the average number of nodes per graph."""
        return self.computeAverage(self.nodeCount())

    def computeAverage(self, count) :
        """Returns the average value given a sum, based on the number of graphs in this instance."""
        N = len(list(self.stats.values()))
        return float(count)/N if N > 0 else 0.0

    def edgeCount(self) :
        """Returns the total number of edges in all graphs."""
        # Could use either parents or children
        return sum([s.edgeCount() for s in list(self.stats.values())])

    def nodeCount(self) :
        """Returns the number of nodes in a graph."""
        return sum([s.nodeCount() for s in list(self.stats.values())])

    def total(self, asType) :
        """Returns the total for the given AS type."""
        funcName = countFunctionName(asType)
        return sum([getattr(s,funcName)() for s in list(self.stats.values())])
