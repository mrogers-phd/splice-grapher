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
Module for converting GFF gene models into a splice graph.
"""
from SpliceGrapher.shared.config     import *
from SpliceGrapher.SpliceGraph       import *
from SpliceGrapher.shared.utils      import idFactory, getAttribute

def forceAddChild(parent, child):
    """Forces a child node to be added to a parent if the
    child node id is not represented in the child list."""
    duplicate = [c for c in parent.children if c.id == child.id]
    if not duplicate:
        parent.children.append(child)
        forceAddParent(child, parent)

def forceAddParent(child, parent):
    """Forces a parent node to be added to a child if the
    parent node id is not represented in the parent list."""
    duplicate = [p for p in child.parents if p.id == parent.id]
    if not duplicate:
        child.parents.append(parent)

def mergeAdjacentNodes(g, minIntron):
    """Merges nodes when the intron between them is too short."""
    nodes = g.nodeDict.values()
    killList = []
    for n in nodes:
        if n.id in killList:
            continue
        for c in n.children:
            if abs(n.donorEnd()-c.acceptorEnd()) > minIntron:
                continue
            killList.append(c.id)

            # NB: may introduce a duplicate node
            n.update(min(n.minpos,c.minpos), max(n.maxpos,c.maxpos))

            for gc in c.children:
                forceAddChild(n, gc)
                gc.parents = removeNode(gc.parents, c.id)
                #if c in gc.parents:
                gc.removeParent(c)

            for p in c.parents:
                if p.id != n.id:
                    forceAddChild(p, n)
                p.children = removeNode(p.children, c.id)
                #if c in p.children:
                p.removeChild(c)

    for id in killList:
        del g.nodeDict[id]

    result = removeDuplicateNodes(g)

    return result

def removeNode(nodeList, nodeId):
    """Removes the node with the given node ID from the list.
    This is in lieu of n.removeParent() or n.removeChild() that
    rely entirely on start/end positions rather than node ids."""
    return [n for n in nodeList if n.id != nodeId]

def removeDuplicateNodes(graph):
    """Removes duplicate nodes from a graph."""
    # Look for duplicate nodes
    nodes   = graph.nodeDict.values()
    nodeSet = set(nodes)
    if len(nodeSet) == len(nodes):
        return graph

    # Rebuild graph to eliminate duplicate nodes
    result = SpliceGraph(graph.getName(), graph.chromosome, graph.strand)
    idMap  = {}
    for n in graph.nodeDict.values():
        rnode = result.addNode(n.id, n.minpos, n.maxpos)
        # Use addIsoform instead of set operations
        # because other things need to happen:
        for form in n.isoformSet:
            rnode.addIsoform(form)
        idMap[n.id] = rnode.id

    for n in graph.nodeDict.values():
        eid1 = idMap[n.id]
        for c in n.children:
            result.addEdge(eid1, idMap[c.id])
    return result

def geneModelToSpliceGraph(gene, **args):
    """Creates a splice graph from a gene model."""
    useCDS    = getAttribute('useCDS', False, **args)
    minexon   = getAttribute('minexon', 1, **args)
    minintron = getAttribute('minintron', 4, **args)

    def makeKey(e):
        return '%d;%d;%s' % (e.minpos,e.maxpos,e.strand)

    # First pass: create exon nodes
    #   store exon --> exon ID
    graph    = SpliceGraph(gene.id, gene.chromosome, gene.strand)
    graph.attrs.update(gene.attributes)
    graph.attrs.setdefault('Note', gene.featureType)
    exonDict = {}
    exonIds  = idFactory('%s_' % gene.id)

    nodeSet = set()
    if useCDS:
        for key in gene.mrna.keys():
            nodeSet.update(gene.mrna[key].sortedExons())
    else:
        nodeSet.update(gene.exons)

    nodeList = sorted(list(nodeSet)) # Needed for debugging
    badForms = set([])
    for exon in nodeList:
        exon_key = makeKey(exon)
        if exon_key in exonDict:
            continue
        if len(exon) < minexon:
            badForms.update(exon.parents)
            continue

        eid  = exonIds.next()
        node = graph.addNode(eid, exon.minpos, exon.maxpos)
        exonDict[exon_key] = node.id

        for p in exon.parents:
            node.addIsoform(p.id)

    featureList = gene.mrna.values() if useCDS else gene.isoforms.values()
    featureList = [f for f in featureList if f not in badForms]
    if not featureList:
        isoLen = len(gene.isoforms)
        rnaLen = len(gene.mrna)
        if useCDS and rnaLen == 0:
            if isoLen > 0:
                raise ValueError('Gene %s has no mRNA (CDS) records; it has %d isoform records instead' % (gene.name,isoLen))
            else:
                raise ValueError('Gene %s has no mRNA (CDS) records' % gene.name)
        elif isoLen == 0:
            if rnaLen > 0:
                raise ValueError('Gene %s has no isoform records; it has %d mRNA (CDS) records instead' % (gene.name,rnaLen))
            else:
                raise ValueError('Gene %s has no isoform or mRNA records' % gene.name)
        else:
            raise ValueError('All isoforms contain exons shorter than minimum length of %d' % minexon)

    for form in featureList:
        exons = form.sortedExons()
        for i in range(1,len(exons)):
            eid1 = exonDict[makeKey(exons[i-1])]
            eid2 = exonDict[makeKey(exons[i])]
            graph.addEdge(eid1, eid2)

    if graph.minpos < graph.maxpos and len(graph) > 0:
        graph = mergeAdjacentNodes(graph, minintron)

    isoDict = graph.isoformDict()
    # Annotate start/end codons with their associated nodes:
    for formId in isoDict:
        try:
            start = gene.mrna[formId].startCodon()
            end   = gene.mrna[formId].endCodon()
        except KeyError:
            continue
        if not (start or end):
            continue
        for n in isoDict[formId]:
            if start:
                n.addStartCodon(start)
            if end:
                n.addEndCodon(end)

    return graph

def makeSpliceGraph(gene, **args):
    """Converts a gene model to a splice graph  First it attempts to
    use exon records.  If that fails, it tries CDS records."""
    try:
        return geneModelToSpliceGraph(gene, **args)
    except ValueError:
        pass

    try:
        return geneModelToSpliceGraph(gene, useCDS=True, **args)
    except ValueError:
        raise ValueError('Gene model for %s cannot be converted to a splice graph.' % gene.name)
