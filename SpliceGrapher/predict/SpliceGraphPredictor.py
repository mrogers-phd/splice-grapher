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
Augments an existing (reference) splice graph by using short read
evidence to identify novel alternative splicing events.  This script
will only process a single gene at a time, but may be incorporated
into a more general framework for processing splice graphs genome wide.

The predictor uses three stages:
  1. Load an existing splice graph to use as a baseline reference
  2. (Optional) Merge additional splice graphs (e.g., from EST alignments)
  3. Load short-read evidence for identifying novel AS events.  This
     consists of one or more of the following:
     a. ungapped read alignments (single file)
     b. splice junction read data (one or more files)
  4. Build and output a revised splice graph that may predict novel
     exons and introns based on the new evidence.
"""
from SpliceGrapher.shared.utils     import *
from SpliceGrapher.shared.config    import *
from SpliceGrapher.shared.ShortRead import SpliceJunction, PREDICTED_JCT, KNOWN_JCT, depthsToClusters
from SpliceGrapher.shared.GeneModelConverter  import makeSpliceGraph
from SpliceGrapher.formats.sam      import *
from SpliceGrapher.SpliceGraph      import *

import sys, os, re

# Regular expression that accepts a resolvable
# 'language' of acceptor/donor site patterns
RESOLVABLE = re.compile("^d*a(a+|d+a+|d+)?da*$")

# Characters for converting patterns of acceptor/donor
# sites into strings
ACC_CHAR   = 'a'
DON_CHAR   = 'd'

# Special attribute for freezing a node's disposition;
# set to True or False
FROZEN_KEY       = 'frozen'
NOVEL_ATTRIBUTES = [UNRESOLVED_NODE, PREDICTED_NODE]

#----------------------------------------------------------
# Default gap-filling heuristic to use until we develop a
# superior machine-learning method.  For flexibility, this
# function is separate from the prediction class.
def fillGap_default(A, B, graph, **args) :
    """
    Returns True if virtual clusters A and B should be merged;
    false otherwise.  Assumes A and B are ordered by ascending positions.
    """
    gapThreshold = getAttribute('gapThreshold', 0, **args)
    junctions    = getAttribute('junctions', [], **args)
    verbose      = getAttribute('verbose', False, **args)

    # First heuristic: do not merge gaps that are
    # larger than than the gap threshold.
    gap = B.minpos - A.maxpos - 1 # e.g., A=[60-100], B=[101-150] --> gap=0
    if gap > gapThreshold : return False

    # Second heuristic: do not merge if there is a donor
    # donor site in the upstream cluster OR an acceptor
    # site in the downstream cluster.  Relative donor/acceptor
    # positions depend on the strand.
    (X,Y) = (A,B) if graph.strand == '+' else (B,A)
    for j in junctions :
        if X.minpos <= j.donor <= X.maxpos : return False
        if Y.minpos <= j.acceptor <= Y.maxpos : return False

    return True
#----------------------------------------------------------
def addPutativeChildren(node, children) :
    """Adds putative child nodes to an unresolved (frozen) node."""
    newKids = set([c.id for c in children])
    try :
        existing = set(node.attrs[PUTATIVE_CHILDREN].split(','))
    except KeyError :
        existing = set()
    node.attrs[PUTATIVE_CHILDREN] = ','.join(existing|newKids)

def addPutativeParents(node, parents) :
    """Adds putative parent nodes to an unresolved (frozen) node."""
    newFolks = set([p.id for p in parents])
    try :
        existing = set(node.attrs[PUTATIVE_PARENTS].split(','))
    except KeyError :
        existing = set()
    node.attrs[PUTATIVE_PARENTS] = ','.join(existing|newFolks)

def frozenNode(node) :
    """Returns True if the node is Frozen; false otherwise."""
    try :
        return node.attrs[FROZEN_KEY]
    except KeyError :
        return False

def predictGraphs(genes, depths, junctions, dirname, **args) :
    """Makes all predictions for a particular list of gene models.
    Returns the number of graphs predicted."""
    novelOnly      = getAttribute('novelOnly', False, **args)
    verbose        = getAttribute('verbose', False, **args)
    depthThreshold = getAttribute('depthThreshold', 1.0, **args)

    if verbose : sys.stderr.write('  using %d depths and %d junctions\n' % (len(depths), len(junctions)))

    junctions.sort()
    predictedGraphs = set()
    indicator       = ProgressIndicator(10000, verbose=verbose)
    indicator.started = True # forces \n whenever prediction finishes
    for gene in genes :
        indicator.update()
        try :
            origGraph = makeSpliceGraph(gene)
        except ValueError :
            continue

        numNodes  = len(origGraph.nodeDict)
        if numNodes == 0 : continue

        # Look for junctions that overlap the gene
        jctMinpos = origGraph.maxpos
        jctMaxpos = origGraph.minpos
        geneJcts  = []
        for j in junctions :
            if j.maxpos < origGraph.minpos or j.minpos > origGraph.maxpos or j.strand != origGraph.strand : continue
            geneJcts.append(j)
            jctMinpos = min(jctMinpos, j.minpos)
            jctMaxpos = max(jctMaxpos, j.maxpos)

        # Allows prediction algorithm to find novel start/end nodes
        # outside the original gene boundaries:
        (lo,hi) = (origGraph.minpos,origGraph.maxpos)
        # Find appropriate depth range for identifying novel start/end exons
        maxpos = len(depths)-1
        lo = min(origGraph.minpos, min(maxpos, jctMinpos))
        hi = max(origGraph.maxpos, min(maxpos, jctMaxpos))
        while lo > 0 and depths[lo] > 0 : lo -= 1
        while hi < len(depths) and depths[hi] > 0 : hi += 1

        predictor = SpliceGraphPredictor(origGraph)

        predictor.setReadDepths(depths[lo:hi+1], reference=lo, threshold=depthThreshold, useAll=True)
        predictor.addSpliceJunctions(geneJcts, includeOverlaps=True)
        if verbose :
            geneDepths = depths[lo:hi+1]
            sys.stderr.write('    %s: %d-%d\n' % (origGraph.name, origGraph.minpos, origGraph.maxpos))
            sys.stderr.write('        %d depths : min %d / max %d\n' % (len(geneDepths), min(geneDepths), max(geneDepths)))
            sys.stderr.write('        %d junctions:\n' % len(geneJcts))
            for jct in geneJcts :
                sys.stderr.write('      %s\n' % jct)

        filePath = os.path.join(dirname, '%s.gff'%origGraph.name)
        try :
            newGraph = predictor.updateGraph()
            # Always output the graph if it changed
            # and update the prediction count
            if not equivalentGraphs(newGraph, origGraph) :
                newGraph.writeGFF(filePath)
                predictedGraphs.add(newGraph.name)
            elif not novelOnly :
                newGraph.writeGFF(filePath)
        except ValueError :
            pass

    indicator.finish()
    return len(predictedGraphs)

class VirtualCluster(object) :
    """Class that encapsulates a virtual cluster of contiguous
    overlapping short-read clusters and graph nodes."""
    def __init__(self, node=None, cluster=None) :
        self.minpos   = sys.maxint
        self.maxpos   = 0
        self.nodes    = set()
        self.clusters = set()
        if node    : self.addNode(node)
        if cluster : self.addCluster(cluster)

    def addCluster(self, c) :
        """Updates the virtual cluster boundaries based on the given Cluster."""
        self.update(c)
        self.clusters.add(c)

    def addNode(self, n) :
        """Updates the virtual cluster boundaries based on the given SpliceGraphNode."""
        self.update(n)
        self.nodes.add(n)

    def __cmp__(self, o) :
        return self.maxpos-o.maxpos if self.minpos == o.minpos else self.minpos-o.minpos

    def __eq__(self, o) :
        return self.minpos == o.minpos and self.maxpos == o.maxpos

    def maxElement(self) :
        """Returns the element (cluster or node) that has the maximum position."""
        lastNode  = sorted(self.nodes)[-1] if self.nodes else None
        lastClust = sorted(self.clusters)[-1] if self.clusters else None
        #assert(lastNode or lastClust)
        if not lastNode :
            return lastClust
        elif not lastClust :
            return lastNode
        else :
            return lastNode if lastNode.maxpos > lastClust.maxpos else lastClust

    def minElement(self) :
        """Returns the element (cluster or node) that has the minimum position."""
        firstNode  = sorted(self.nodes)[0] if self.nodes else None
        firstClust = sorted(self.clusters)[0] if self.clusters else None
        #assert(firstNode or firstClust)
        if not firstNode :
            return firstClust
        elif not firstClust :
            return firstNode
        else :
            return firstNode if firstNode.minpos < firstClust.minpos else firstClust

    def merge(self, vc) :
        """Merges this cluster with another cluster.  The other cluster remains
        unchanged while this cluster's bounds and node set."""
        self.update(vc)
        self.nodes.update(vc.nodes)

    def overlaps(self, o) :
        """Returns True if two virtual clusters overlap; False otherwise."""
        return o.maxpos > self.minpos and self.maxpos > o.minpos

    def __str__(self) :
        return '%d-%d' % (self.minpos, self.maxpos)

    def update(self, feature) :
        """Updates the virtual cluster boundaries based on the given feature."""
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)

class SimpleJunction(object) :
    """Class that encapsulates only the information necessary to
    identify a junction and its acceptor/donor sites."""
    def __init__(self, acc, don, strand) :
        self.minpos   = min(acc,don)
        self.maxpos   = max(acc,don)
        self.strand   = strand
        self.acceptor = self.maxpos if self.strand == '+' else self.minpos
        self.donor    = self.minpos if self.strand == '+' else self.maxpos

    def __cmp__(self, o) :
        return self.maxpos-o.maxpos if self.minpos == o.minpos else self.minpos-o.minpos

    def __eq__(self, o) :
        return self.minpos == o.minpos and self.maxpos == o.maxpos and self.strand == o.strand

    def __hash__(self) :
        return self.__str__().__hash__()

    def __str__(self) :
        return '%d-%d (%s)' % (self.donor, self.acceptor, self.strand)

class SpliceGraphPredictor(object) :
    nodegen = idFactory(PNODE_PREFIX)

    def __init__(self, referenceGraph, **args) :
        """
        The splice graph builder requires an initial splice graph as a reference,
        This should be a GFF splice graph such as one generated from a gene annotation
        GFF3 file via makeSpliceGraph
        """
        self.talkative       = getAttribute('talkative', False, **args)
        self.verbose         = getAttribute('verbose', False, **args)
        self.clusters        = []
        self.simpleJunctions = set()
        self.novelJunctions  = set()
        self.initializeSpliceGraph(referenceGraph, **args)
        self.verboseMessage('SpliceGraphPredictor for %s instantiated with %d nodes between %d and %d' % \
                (referenceGraph.name, len(self.nodes), self.graph.minpos, self.graph.maxpos))

    def addSpliceGraph(self, auxGraph, **args) :
        """Augments the original splice graph with additional data from another splice graph."""
        mergeEnds = getAttribute('mergeEnds', False, **args)

        # Set attributes prior to merging graphs
        auxNodes  = sorted(auxGraph.nodeDict.values())
        for node in auxNodes :
            node.addAttribute(DISPOSITION_KEY, KNOWN_NODE)

        # Update the global graph
        self.graph = self.graph.union(auxGraph, keepName=True, mergeEnds=mergeEnds)
        self.simpleJunctions.update([SimpleJunction(e.minpos, e.maxpos, self.graph.strand) for e in edgeSet(auxGraph)])

        self.verboseMessage('Added splice graph %s (total nodes = %d, edges = %s)' % (auxGraph, len(self.nodes), len(self.simpleJunctions)))

    def addSpliceJunctions(self, jctList, **args) :
        """
        Imports splice junctions prior to making a prediction.  The caller may add junctions
        from as many files as desired, however only unique junctions will be added.
        """
        useAll = getAttribute('includeOverlaps', False, **args)

        def includeJunction(jct) :
            if useAll :
                return (self.graph.minpos <= j.maxpos and j.minpos <= self.graph.maxpos)
            else :
                return (self.graph.minpos <= j.minpos and j.maxpos <= self.graph.maxpos)

        self.verboseMessage('%s adding junctions from a list of %d records' % (self.graph.name, len(jctList)))
        # Load only novel junctions on the correct strand and within the graph:
        oldLen = len(self.simpleJunctions)
        self.verboseMessage('  (already %d junctions in the graph)' % oldLen)
        for j in jctList :
            #if j.sjCode == KNOWN_JCT or self.graph.strand != j.strand : continue
            if self.graph.strand != j.strand : continue
            if includeJunction(j) :
                newJct = SimpleJunction(j.donor(), j.acceptor(), j.strand)
                self.novelJunctions.add(newJct)
                self.simpleJunctions.add(newJct)
        newLen = len(self.simpleJunctions)
        self.verboseMessage('  %s loaded %d novel junctions' % (self.graph.name, (newLen-oldLen)))

    def constructNodes(self, acc, don) :
        """
        Constructs nodes from all combinations of acceptors
        and donors in the lists.
        """
        result = []
        for a in acc :
            for d in don :
                # Only use donors downstream of acceptor
                if self.graph.strand == '+' and d < a : continue
                if self.graph.strand == '-' and d > a : continue

                newId = self.nodegen.next()
                node  = self.graph.addNode(newId, a, d)
                # If the node already exists (new node id not needed) skip it
                if node.id == newId : result.append(node)
        return result

    def initializeSpliceGraph(self, graph, **args) :
        """Initializes the predictor instance from the given splice graph."""
        self.graph = graph.duplicate()
        self.nodes = sorted(self.graph.resolvedNodes())
        for node in self.nodes :
            node.addAttribute(DISPOSITION_KEY, KNOWN_NODE)
        self.simpleJunctions = set([SimpleJunction(e.minpos, e.maxpos, graph.strand) for e in edgeSet(graph)])
        self.verboseMessage('Initialized splice graph with %d nodes and %d junctions' % (len(self.nodes), len(self.simpleJunctions)))

    def makeSiteString(self, cluster, accSet, donSet) :
        """
        Convert a set of acceptor and donor sites within a virtual cluster
        into a string representation of 'a' and 'd' characters.
        """
        positions = range(cluster.minpos, cluster.maxpos+1)
        symbols   = {}.fromkeys(positions,'')

        # An acceptor and a donor may use the same splice-site location when there are
        # AGGT sequences, based on the inferred adjacent exon (AG|GT).  Thus we put the
        # 'A' immediately before the 'D' by moving it back one position.  (Note that we
        # can do this because it is impossible to have an AG dimer in two consecutive
        # positions).  This might not work for 'AA' dimers because 'AAA' is possible,
        # but this should be extremely(!) rare.
        adjusted = {}
        for a in accSet :
            adjusted[a] = a-1 if a in donSet else a

        for a in accSet :
            if cluster.minpos <= a <= cluster.maxpos :
                pos = adjusted[a]
                symbols[pos] = ACC_CHAR

        for d in donSet :
            if cluster.minpos <= d <= cluster.maxpos :
                symbols[d] = DON_CHAR

        # Build the string
        result = ''.join([symbols[i] for i in positions if symbols[i] in [ACC_CHAR,DON_CHAR]])

        # Correct sequence of 'a' and 'd' is strand-dependent
        return result if self.graph.strand == '+' else result[::-1]

    def setReadDepths(self, depths, **args) :
        """Stores read depth information as clusters."""
        verbose   = getAttribute('verbose', False, **args)
        # minDepth determines at what level clusters should be split
        minDepth  = getAttribute('minDepth', 1, **args)
        # threshold determines the average depth required to accept a cluster
        threshold = getAttribute('threshold', 1, **args)
        reference = getAttribute('reference', 0, **args)
        useAll    = getAttribute('useAll', False, **args)

        allClusters   = depthsToClusters(self.graph.chromosome, depths, minDepth=minDepth, threshold=threshold, reference=reference)
        if useAll :
            self.clusters = allClusters
        else :
            self.clusters = [c for c in allClusters if c.maxpos > self.graph.minpos and c.minpos < self.graph.maxpos]

        self.verboseMessage('%s loaded reads into %d distinct clusters:' % (self.graph.name, len(self.clusters)))

    def talkativeMessage(self, s) :
        """Convenience method for writing messages only in talkative mode."""
        if self.talkative : sys.stdout.write('    %s\n' % s)

    def updateDisposition(self, node) :
        """
        Returns the disposition for a predicted node based on
        whether any of its parent or child edges have been resolved.
        """
        # Known or frozen nodes do not change disposition
        if node.isKnown() or frozenNode(node) : return

        # Default to valid prediction
        node.addAttribute(DISPOSITION_KEY, PREDICTED_NODE)

        # If there are candidate junctions but none of them are linked, node is unresolved
        upstreamCandidates = [j for j in self.simpleJunctions if j.acceptor == node.acceptorEnd()]
        if upstreamCandidates and not node.parents :
            node.addAttribute(DISPOSITION_KEY, UNRESOLVED_NODE)
            return

        downstreamCandidates = [j for j in self.simpleJunctions if j.donor == node.donorEnd()]
        if downstreamCandidates and not node.children :
            node.addAttribute(DISPOSITION_KEY, UNRESOLVED_NODE)

    def updateEdges(self, node, junctions) :
        """Resolves the edges for a node into upstream parents and downstream children."""
        upstream   = [j.donor for j in junctions if j.acceptor == node.acceptorEnd()]
        downstream = [j.acceptor for j in junctions if j.donor == node.donorEnd()]

        if frozenNode(node) :
            # Find all possible parent/child relationships for frozen unresolved nodes:
            parents  = [n for n in self.nodes if n.donorEnd() in upstream]
            children = [n for n in self.nodes if n.acceptorEnd() in downstream]
            addPutativeParents(node,parents)
            addPutativeChildren(node,children)
        else :
            parents  = [n for n in self.nodes if n.donorEnd() in upstream if not frozenNode(n)]
            children = [n for n in self.nodes if n.acceptorEnd() in downstream if not frozenNode(n)]
            for p in parents :
                p.addChild(node)
            for c in children :
                node.addChild(c)

    def updateGraph(self, **args) :
        """
        Updates the splice graph based on the evidence given by short read depths
        and known and putative splice junctions.
        """
        fillGapFunction = getAttribute('gapFunction', fillGap_default, **args)

        savedGraph = self.graph.duplicate()
        savedGraph.annotate()

        # Initialize existing nodes so they are not frozen:
        for n in self.nodes :
            n.addAttribute(FROZEN_KEY,False)

        # Only resolve clusters that are not subsumed by existing nodes:
        unresolvedClusters = []
        for c in self.clusters :
            subsuming = None
            for n in self.nodes :
                if n.minpos <= c.minpos and c.maxpos <= n.maxpos :
                    subsuming = n
                    break
            if not subsuming :
                unresolvedClusters.append(c)

        externalClusters = [c for c in unresolvedClusters if c.maxpos < savedGraph.minpos or c.minpos > savedGraph.maxpos]

        # Merge clusters:
        #   1. convert clusters and nodes into virtual clusters
        #   2. merge clusters that overlap
        #   3. merge clusters when the fill-gap classifier
        #      considers it continuous coverage
        virtual  = [VirtualCluster(cluster=c) for c in unresolvedClusters]
        virtual += [VirtualCluster(node=n) for n in self.nodes]
        virtual.sort()
        distinct = [virtual[0]]
        for v in virtual[1:] :
            if distinct[-1].overlaps(v) :
                distinct[-1].merge(v)
            elif fillGapFunction(distinct[-1], v, self.graph, **args) :
                distinct[-1].merge(v)
            else :
                distinct.append(v)

        # Critical step: create strings from acceptor/donor sites within
        # each cluster and test for membership in regular language.  Create
        # nodes from those patterns the language accepts.
        donors           = set([j.donor for j in self.simpleJunctions])
        acceptors        = set([j.acceptor for j in self.simpleJunctions])

        # Necessary for extending genes by predicting novel exons when a
        # splice junction extends outside gene boundaries.
        for c in externalClusters :
            acc = set([a for a in acceptors if c.minpos <= a <= c.maxpos])
            don = set([d for d in donors if c.minpos <= d <= c.maxpos])
            # No need to add pseudo-sites if an externsl cluster
            # overlaps junctions on both ends
            if acc and don : continue

            if acc :
                pseudoDon = c.maxpos if self.graph.strand == '+' else c.minpos
                donors.add(pseudoDon)

            if don :
                pseudoAcc = c.minpos if self.graph.strand == '+' else c.maxpos
                acceptors.add(pseudoAcc)

        acceptedClusters = []
        novelNodes       = set()
        for c in distinct :
            acc  = set([a for a in acceptors if c.minpos <= a <= c.maxpos])
            don  = set([d for d in donors if c.minpos <= d <= c.maxpos])
            if not (acc or don) : continue
            for n in c.nodes :
                acc.add(n.acceptorEnd())
                don.add(n.donorEnd())

            siteString = self.makeSiteString(c, acc, don)
            if RESOLVABLE.match(siteString) :
                self.talkativeMessage('%s cluster %s site pattern %s good.' % (self.graph.name, c, siteString))
                novelNodes.update(self.constructNodes(acc, don))
            else :
                self.talkativeMessage('%s cluster %s site pattern %s unresolvable.' % (self.graph.name, c, siteString))
                newNodes = self.constructNodes(acc, don)
                for n in newNodes :
                    n.setUnresolved()
                    n.addAttribute(FROZEN_KEY, True)
                # constructNodes only returns novel nodes, so we're safe doing this:
                novelNodes.update(newNodes)

        # Resolve edges for novel nodes
        self.nodes = self.graph.nodeDict.values()
        for node in novelNodes :
            adjacentJunctions = [j for j in self.simpleJunctions if j.acceptor == node.acceptorEnd() or j.donor == node.donorEnd()]
            self.updateEdges(node, adjacentJunctions)
            self.updateDisposition(node)
            self.talkativeMessage('  created new node %s as %s' % (node,node.attrs[DISPOSITION_KEY]))

        # Resolve any unresolved edges
        graphJunctions   = set([SimpleJunction(e.minpos, e.maxpos, self.graph.strand) for e in edgeSet(self.graph)])
        missingJunctions = self.novelJunctions - graphJunctions
        for j in missingJunctions :
            candidates = [n for n in self.graph.nodeDict.values() if n.donorEnd() == j.donor or n.acceptorEnd() == j.acceptor]
            if not candidates : continue
            for node in candidates :
                before    = node.attrs[DISPOSITION_KEY]
                self.updateEdges(node, [j])
                self.updateDisposition(node)

        # Final step: annotate the graph
        self.graph.annotate()

        # Remove temporary internal flags:
        for n in self.graph.nodeDict.values() :
            try :
                del n.attrs[FROZEN_KEY]
            except KeyError :
                continue

        if self.verbose :
            oldSize = len(savedGraph.resolvedNodes())
            newSize = len(self.graph.resolvedNodes())
            oldAS   = len(savedGraph.distinctAltEvents())
            newAS   = len(self.graph.distinctAltEvents())
            self.verboseMessage('Added %d nodes and %d AS events to graph.' % ((newSize-oldSize), (newAS-oldAS)))

        # Write graph to output file
        if 'output' in args :
            self.talkativeMessage('Unresolved nodes in %s:' % self.graph.name)
            for n in self.graph.nodeDict.values() :
                if not n.isUnresolved() : continue

            outStream = args['output']
            if type(outStream) == str :
                self.verboseMessage('Writing new graph to %s' % outStream)
                if self.graph.minpos < savedGraph.minpos or savedGraph.maxpos < self.graph.maxpos :
                    self.talkativeMessage('Updated graph extends original: %d <= %d or %d <= %d' % (self.graph.minpos, savedGraph.minpos, savedGraph.maxpos, self.graph.maxpos))
                outStream = open(outStream, 'w')
            self.graph.writeGFF(outStream)

        return self.graph

    def verboseMessage(self, s) :
        """Convenience method for writing messages only in verbose mode."""
        if self.verbose : sys.stdout.write('%s\n' % s)
