#! /usr/bin/env python
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

ALT_SITE_EVENTS = set([ALT5_ABBREV, ALT3_ABBREV])

class simpleJunction(object) :
    """Class that encapsulates only an edge's start and end
    positions, with hashing for set operations.  These positions
    will be consistent regardless of which nodes they connect."""
    def __init__(self, parent, child) :
        self.parent = parent
        self.child  = child
        self.minpos = self.parent.maxpos if parent.strand == '+' else self.child.maxpos
        self.maxpos = self.child.minpos  if parent.strand == '+' else self.parent.minpos

    def __cmp__(self, o) :
        return self.maxpos-o.maxpos if self.minpos == o.minpos else self.minpos-o.minpos

    def __hash__(self) :
        return self.__str__().__hash__()

    def __str__(self) :
        return '%d,%d' % (self.minpos, self.maxpos)

class uniqueExonFeature(object) :
    """Class that encapsulates a piece of an exon that is
    unique to that exon."""
    def __init__(self, exon, pos) :
        self.minpos = pos
        self.maxpos = pos
        self.exon   = exon

    def descr(self) :
        return 'Feature (%d,%d) --> %s' % (self.minpos, self.maxpos, self.exon)

    def __cmp__(self, o) : return (self.maxpos-o.maxpos) if self.minpos == o.minpos else (self.minpos-o.minpos)
    def __eq__(self, o)  : return (self.minpos,self.maxpos) == (o.minpos,o.maxpos)
    def __hash__(self)   : return self.__str__().__hash__()
    def __len__(self)    : return self.maxpos-self.minpos+1
    def __str__(self)    : return '%d,%d' % (self.minpos,self.maxpos)

    def update(self, pos) :
        self.minpos = min(self.minpos, pos)
        self.maxpos = max(self.maxpos, pos)
        #assert(self.minpos <= self.maxpos)

def getExonFeatures(nodes, **args) :
    """Finds segments that uniquely identify exons in a graph.  Returns
    a dictionary of segments that correspond to exactly one corresponding exon."""
    minOverlap = getAttribute('minOverlap', 1, **args)
    # First find positions that overlap exons
    overlaps = {}
    for n in nodes :
        for i in xrange(n.minpos, n.maxpos+1) :
            overlaps.setdefault(i,[])
            overlaps[i].append(n)

    # Now look for positions that overlap exactly one exon:
    keys    = sorted(overlaps.keys())
    result  = {}
    feature = None
    for k in keys :
        if len(overlaps[k]) != 1 :
            feature = None
        else :
            node = overlaps[k][0]
            if feature and node == feature.exon :
                feature.update(k)
            else :
                feature = uniqueExonFeature(node,k)
                result[feature] = node
    return result

def getGraphDepths(depths, graph) :
    """Returns a list of depths from the given depth list by including
    all depths within the graph.  Depths extend to the nearest end of
    read coverage outside the graph."""
    minpos = graph.minpos
    maxpos = min(len(depths), graph.maxpos)

    # Depths do not reach this graph
    if maxpos < minpos : return [],minpos

    while minpos > 0 and depths[minpos] > 0 : minpos -= 1
    while maxpos < len(depths) and depths[maxpos] > 0 : maxpos += 1
    return depths[minpos:maxpos], minpos

def getIsoformJctSet(isoform, nodeList) :
    """Finds all nodes in the list that belong to the given
    isoform and returns a list of the edges inferred by the nodes."""
    isoformNodes = []
    strand       = None
    for n in nodeList :
        if isoform in n.isoformSet :
            isoformNodes.append(n)
            if not strand : strand = n.strand
    isoformNodes.sort(reverse=(strand=='-'))

    result = set([])
    for i in xrange(1,len(isoformNodes)) :
        result.add(simpleJunction(isoformNodes[i-1], isoformNodes[i]))
    return result

def getSupportedNodeDict(nodes, edgeSet, clusters, **args) :
    """Returns a dictionary of clusters that overlap or contain
    features unique to specific nodes."""
    adjustOverlap = getAttribute('adjustOverlap', False, **args)
    minOverlap    = getAttribute('minOverlap', 1, **args)
    minDepth      = getAttribute('mindepth', 1.0, **args)
    verbose       = getAttribute('verbose', False, **args)

    result   = {}
    features = sorted(getExonFeatures(nodes, **args))
    for f in features :
        # Set overlap size no larger than actual feature if adjustOverlap is set
        overlapSize  = min(len(f), minOverlap) if adjustOverlap else minOverlap
        smallFeature = len(f) < overlapSize

        # We may still allow features smaller than the overlap size to be flagged
        # as supported, provided there are no other nodes/features nearby.  Note
        # that the only way to flag them will be if a cluster contains them completely
        # or the feature contains a cluster and the read length is smaller than the overlap.
        if smallFeature :
            smallFeature = True
            evidence = None
            for n in nodes :
                if n == f.exon : continue
                if f.minpos > n.maxpos and f.minpos - n.maxpos < overlapSize : evidence = n
                if f.maxpos > n.minpos and n.minpos - f.maxpos < overlapSize : evidence = n
                if evidence : break

            if evidence :
                if verbose :
                    sys.stderr.write('Feature %s length %d is smaller than required overlap %d\n' % (f.exon.id, len(f), overlapSize))
                    sys.stderr.write('and too close to node %s\n' % evidence)
                continue
            elif verbose :
                sys.stderr.write('Feature %s length %d is smaller than required overlap %d\n' % (f.exon.id, len(f), overlapSize))
                sys.stderr.write('but no other nodes are close\n')

        for c in clusters :
            # Coverage must exceed the minimum across the unique feature region
            featureDepth = c.avgDepth(minpos=f.minpos,maxpos=f.maxpos)
            if featureDepth < minDepth : continue

            result.setdefault(c,set())
            # Cases split out to make them more readable:
            if f.minpos < c.minpos and c.maxpos < f.maxpos : # Feature contains cluster
                result[c].add(f.exon)
            elif c.minpos < f.minpos and f.maxpos < c.maxpos : # Cluster contains feature
                result[c].add(f.exon)
            elif (f.minpos+overlapSize) <= c.maxpos and c.minpos <= (f.maxpos-overlapSize) : # Cluster overlaps feature by minimum amount
                result[c].add(f.exon)

            if verbose and f.exon in result[c] :
                nodeString  = 'root ' if f.exon.isRoot() else ''
                nodeString += '%s %d-%d (%s)' % (f.exon.id, f.exon.start, f.exon.end, f.exon.isoformString())
                sys.stderr.write('    %s supported by read depth %.1f from %d-%d\n' % (nodeString, featureDepth, f.minpos, f.maxpos))

    return result
    
def makeSupportedGraph(graph, clusters, junctions, **args) :
    """Uses depth and junction information from a SAM file to find those
    splice forms in a graph that are represented by RNA-Seq data.  Returns a graph
    representing only the RNA-Seq supported forms.  If no paths are supported,
    returns None."""
    verbose  = getAttribute('verbose', False, **args)

    if verbose : sys.stderr.write('Finding supported isoforms in %s\n' % graph.getName())

    #-------------------------------------------------
    # Identify isoforms associated with nodes
    nodes             = graph.resolvedNodes()
    graphEdgeSet      = edgeSet(graph)
    supportedNodeDict = getSupportedNodeDict(nodes, graphEdgeSet, clusters, **args)
    if not supportedNodeDict :
        if verbose : sys.stderr.write('  ** no RNA-Seq-supported nodes for %s\n' % graph.getName())
        return

    # Clusters are associated with distinct exon features, but
    # we only choose features associated with unique isoforms:
    isoforms = set()
    for c in supportedNodeDict :
        for n in supportedNodeDict[c] :
            if len(n.isoformSet) == 1 :
                if verbose : sys.stderr.write('  isoform %s supported by node %s\n' % (n.isoformString(), n))
                isoforms.add(n.isoformString())

    #-------------------------------------------------
    # Identify isoforms associated with edges
    # Supported edges are any recapitulated in the RNA-Seq data
    supportedEdgeSet = set([(j.minpos,j.maxpos) for j in junctions])
    graphJct         = [simpleJunction(n,c) for n in nodes for c in n.children]
    validJctSet      = set([e for e in graphJct if (e.minpos,e.maxpos) in supportedEdgeSet])
    knownForms       = graph.isoformDict()

    # Find graph edges associated with unique isoforms
    uniqueIsoformJct = uniqueIsoformEdges(knownForms)

    # Find junctions associated with unique isoforms that are also supported
    if verbose : beforeSize = len(isoforms)
    keysToAdd = [k for k in uniqueIsoformJct.keys() if k in validJctSet]
    isoforms.update([uniqueIsoformJct[e] for e in keysToAdd])
    if verbose :
        afterSize = len(isoforms)
        delta = afterSize-beforeSize
        sys.stderr.write('  added %d isoforms based on junctions:\n' % delta)
        for e in keysToAdd :
            sys.stderr.write('  %s : %s\n' % (e,uniqueIsoformJct[e]))

    if not isoforms :
        if verbose : sys.stderr.write('  no isoforms were uniquely represented in %s\n' % graph.getName())
        return

    sortedIso = sorted(list(isoforms))
    if verbose : sys.stderr.write('  %d isoforms identified for %s:\n  %s\n' % (len(isoforms), graph.getName(), '\n  '.join(sortedIso)))

    # Identify all junctions associated with
    # isoforms to be included in the subgraph
    isoformJct = set([])
    for iso in isoforms :
        isoformJct |= getIsoformJctSet(iso, nodes)
        
    #------------------------
    # Create the new subgraph
    result = SpliceGraph(graph.getName(), graph.chromosome, graph.strand)

    # New graph gets same dimensions as original
    result.minpos = graph.minpos
    result.maxpos = graph.maxpos

    # First add isoform nodes to the new graph:
    referenceNodes = set([])
    for n in nodes :
        # Add a node if it is in one of the isoforms
        if isoforms & n.isoformSet :
            rNode = result.addNode(n.id, n.start, n.end)
            updateAttributes(rNode, n, isoforms)
            referenceNodes.add(n)

    # Next add junctions
    for n in referenceNodes :
        # Add edges only to parents/children that are in the reference set
        node         = result.getNode(n.start,n.end)
        validKids    = set(n.children) & referenceNodes
        for c in validKids :
            child = result.getNode(c.start,c.end)
            edge  = simpleJunction(n,c)
            jct   = (edge.minpos,edge.maxpos)
            if (edge in isoformJct) or (jct in supportedEdgeSet) :
                result.addEdge(node.id,child.id)

        validParents = set(n.parents) & referenceNodes
        for p in validParents :
            parent = result.getNode(p.start,p.end)
            edge   = simpleJunction(p,n)
            jct    = (edge.minpos,edge.maxpos)
            if (edge in isoformJct) or (jct in supportedEdgeSet) :
                result.addEdge(parent.id,node.id)

    return result

def uniqueIsoformEdges(knownForms) :
    """Returns a dictionary of graph edges that participate
    in 1 isoform with the name of the associated form."""
    isoformJcts = dict([(k,getIsoformJctSet(k,knownForms[k])) for k in knownForms])
    jctMap      = {}
    for iso in knownForms :
        # For each junction, establish a set of isoforms it appears in
        for jct in isoformJcts[iso] :
            jctMap.setdefault(jct,set([]))
            jctMap[jct].add(iso)

    return dict([(e,jctMap[e].pop()) for e in jctMap if len(jctMap[e]) == 1])

def updateAttributes(node, referenceNode, isoforms) :
    """Updates a node's attributes using the given reference node
    attributes.  Provides special handling for splice forms and isoforms."""
    for k in referenceNode.attrs.keys() :
        if k == AS_KEY :
            continue # AS annotations may not apply to a subgraph
        elif k == ISO_KEY :
            forms = referenceNode.isoformSet & isoforms
            for f in forms : node.addIsoform(f)
        else :
            node.attrs[k] = referenceNode.attrs[k]
