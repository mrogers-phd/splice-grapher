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
Module that encapsulates a splice graph.
"""
from SpliceGrapher.shared.utils import *
import sys

# Attributes used in standard GFF3 files:
SOURCE_NAME  = 'SpliceGraph'
DEPRECATED_GENE_REC = 'cluster'
GENE_REC       = 'graph'
PARENT_REC     = 'parent'
CHILD_REC      = 'child'
VALID_GENES    = [DEPRECATED_GENE_REC, GENE_REC]
VALID_RECTYPES = [DEPRECATED_GENE_REC, GENE_REC, PARENT_REC, CHILD_REC]

# Graph node attributes for splice graph GFF files
PARENT_ATTR = 'Parent'
ID_ATTR     = 'ID'
KNOWN_ATTRS = [PARENT_ATTR, ID_ATTR]

# Specific node attribute and values for AS forms:
AS_KEY       = 'AltForm'

# Attributes for unresolved nodes:
PUTATIVE_PARENTS  = 'putative_parents'
PUTATIVE_CHILDREN = 'putative_children'

# Attributes for predicting protein sequences:
START_CODON_KEY = 'StartCodon'
END_CODON_KEY   = 'EndCodon'

# Create names and abbreviations for AS forms and map each one to the other
AS_NAMES      = [IR_NAME, ALT5_NAME, ALT3_NAME, ALTB3_NAME, ALTB5_NAME, ES_NAME, ALTI_NAME, ALTT_NAME] \
              = ['Intron Retention',"Alt. 5'","Alt. 3'",'Alt. B3','Alt. B5','Skipped Exon','Alt. Init.','Alt. Term.']

AS_ABBREVS    = [IR_ABBREV, ALT5_ABBREV, ALT3_ABBREV, ALTB3_ABBREV, ALTB5_ABBREV, ES_ABBREV, ALTI_ABBREV, ALTT_ABBREV] \
              = ['IR',"A5","A3",'AB3','AB5','SE','AI','AT']
AS_ABBREV_SET = set(AS_ABBREVS)

EVENT_NAME     = dict([(AS_ABBREVS[i],AS_NAMES[i]) for i in range(len(AS_NAMES))])
EVENT_ABBREV   = dict([(AS_NAMES[i],AS_ABBREVS[i]) for i in range(len(AS_NAMES))])

ALT_35_NAMES   = [ALT3_NAME,ALT5_NAME,ALTB3_NAME,ALTB5_NAME]
ALT_35_ABBREVS = [ALT3_ABBREV,ALT5_ABBREV,ALTB3_ABBREV,ALTB5_ABBREV]
NON_35_NAMES   = set(AS_NAMES) - set(ALT_35_NAMES)
NON_35_ABBREVS = set(AS_ABBREVS) - set(ALT_35_ABBREVS)

# Dispositions for nodes in a graph:
DISPOSITION_KEY = 'disposition'
KNOWN_NODE      = 'known'
PREDICTED_NODE  = 'predicted'
UNRESOLVED_NODE = 'unresolved'

# Edge types for adding edges to new nodes
PARENT_EDGE     = 'parent'
CHILD_EDGE      = 'child'
ALL_EDGE_TYPES  = [PARENT_EDGE, CHILD_EDGE]

# Identifies predicted nodes
PNODE_PREFIX = 'pred_'

# Attributes used for displaying splice graphs:
ACCEPTORS_KEY     = 'acceptors'
DONORS_KEY        = 'donors'
INFORMATION_KEY   = 'as_info'

# Keys for retained intron evidence cross-referencing
IR_TP_EVIDENCE    = 'ir_tp_evidence'
IR_FP_EVIDENCE    = 'ir_fp_evidence'

# Specific node attribute for isoform names
ISO_KEY     = 'Isoforms'

# Specific graph attribute for expanded gene models
ALT_MODEL_KEY = 'Expanded'

VALID_STRANDS = ['+','-']

#========================================================================
# Function definitions
def acceptor(node) :
    """Returns the strand-adjusted position at the start of the acceptor dimer."""
    return node.minpos-2 if node.strand == '+' else node.maxpos

def containsEdge(node, edge) :
    """Returns true if the Node contains the given Edge; false otherwise.  Note
    that the edge must be completely contained within the node (<,>)."""
    return node.minpos < edge.minpos and edge.maxpos < node.maxpos

def containsNode(edge, node) :
    """Returns true if the Edge contains the given Node; false otherwise.  Note
    that the node must be completely contained within the edge (<,>)."""
    return edge.minpos < node.minpos and node.maxpos < edge.maxpos

def childEdges(n) :
    """Returns a set of child edges for a node."""
    return set([Edge(n,c) for c in n.children])

def donor(node) :
    """Returns the strand-adjusted position at the start of the donor dimer."""
    return node.maxpos if node.strand == '+' else node.minpos-2

def edgeSet(G) :
    """Returns a complete set of edges found in splice graph G.
    This includes duplicate edges between distinct nodes as found
    in alternate 3'/5' events."""
    result = set([])
    for n in list(G.nodeDict.values()) :
        result.update(childEdges(n))
    return result

def intronSet(G) :
    """Returns the set of distinct edges found in a splice graph.
    This does not include any duplicate edges."""
    return set([(e.minpos,e.maxpos) for e in edgeSet(G)])

def overlap(a, b) :
    """Returns true if the two nodes overlap and are distinct; false otherwise."""
    return a.id != b.id and a.maxpos > b.minpos and a.minpos < b.maxpos

def overlapsAll(A, nodeList) :
    """Returns true if node A overlaps all nodes in the given node list; false otherwise."""
    for n in nodeList :
        if not overlap(A,n) : return False
    return True

def parentEdges(n) :
    """Returns a set of parent edges for a node."""
    return set([Edge(p,n) for p in n.parents])

#------------------------------------------------------------------------
# AS annotation functions:
def altSiteEventList(nodes, siteType, verbose=False) :
    """Returns a list of alt. 3' or alt. 5' events.  Each element consists
    of a set of nodes involved in the same event.  Note: assumes the graph's
    nodes have already been annotated (see detectAltAcceptor and detectAltDonor)."""
    if siteType not in [ALT3_ABBREV, ALT5_ABBREV] :
        raise ValueError('altNodeList called using non-alternative site type: %s' % siteType)

    def setString(s) :
        return ','.join(['%s'%n.id for n in s])

    def neighborNodes(n) :
        return set(n.parents) if siteType == ALT3_ABBREV else set(n.children)

    # Find nodes annotated with the given event
    eventNodes = [n for n in nodes if siteType in n.altForms()]
    if not eventNodes : return []
    if verbose :
        modifier = 'parents' if siteType == ALT3_ABBREV else 'children'
        sys.stderr.write('\naltSiteEventList found %d %s nodes:\n' % (len(eventNodes), siteType))
        for n in eventNodes :
            sys.stderr.write('  %s\n' % n)
        sys.stderr.write('\n')

    # Sorting nodes ensures that overlapping nodes will
    # not be stored in distinct sets to be merged later
    eventNodes.sort()
    eventList = []
    stored    = set()
    for n in eventNodes :
        eset  = set([n])
        # Grab neighboring nodes (share an edge)
        adj_n = neighborNodes(n)
        if verbose : sys.stderr.write('  Node %s has %s %s\n' % (n.id, modifier, setString(adj_n)))

        # Find other nodes with the same annotation that overlap the current one
        for m in eventNodes :
            if m == n or not overlap(m,n) : continue
            if verbose : sys.stderr.write('  %s overlaps %s\n' % (n.id,m.id))
            # Grab other node's neighbors
            adj_m = neighborNodes(m)

            # If either node overlaps all of the other's parents (alt3)
            # or children (alt5), they are not part of the same event:
            neighborsContained = overlapsAll(m,adj_n) or overlapsAll(n,adj_m)

            # Otherwise, the nodes are part of the same event
            if not neighborsContained :
                eset.add(m)
            elif verbose :
                sys.stderr.write('    %s overlaps %s %s (or vice versa)\n' % (n, m, modifier))

        # If the nodes overlapping the current node also overlap
        # distinct other nodes, merge the events into one set
        if verbose : sys.stderr.write('  %s event set: %s\n' % (n.id, setString(eset)))
        for e in eventList :
            if e & eset :
                if verbose : sys.stderr.write('  Merging %s with %s\n' % (setString(eset), setString(e)))
                e.update(eset)
                stored.update(eset)
                if verbose : sys.stderr.write('  --> %s\n' % setString(e))
                break
            elif verbose :
                sys.stderr.write('    %s does not intersect %s\n' % (setString(eset), setString(e)))

        if not n in stored :
            if verbose : sys.stderr.write('  Storing new event set %s\n' % setString(eset))
            eventList.append(eset)

    if verbose : sys.stderr.write('Final result: %s\n' % ';'.join([setString(e) for e in eventList]))
    return eventList

def detectAltAcceptor(n, nodes, edges) :
    """Looks for evidence that identifies the given node as an alternate acceptor."""
    # Inference rules:
    #  - Establish set of all nodes that overlap the given one and are not graph roots
    #  - Get the set of nodes flanking edges that the node contains (retained introns)
    #  - Get the set of all nodes that overlap parents of the given one (other retained introns)
    #  - Subtract those overlapping parents from the first set
    #  - For those that remain, if any have acceptor sites different from this one, mark it as an alternate acceptor
    if not n.parents : return
    allOverlaps    = set([o for o in nodes if o.parents and overlap(o,n)])
    containedEdges = set([e for e in edges if containsEdge(n,e)])
    flankingNodes  = set([e.child for e in containedEdges if overlapsAll(n,e.child.parents)])
    parentOverlaps = set([o for o in nodes if overlapsAll(o,n.parents)])
    overlapSet = allOverlaps - parentOverlaps - flankingNodes

    for o in overlapSet :
        if o.acceptorEnd() != n.acceptorEnd() :
            n.addAltForm(ALT3_ABBREV)
            return
        
def detectAltBoth(nodes, verbose=False) :
    """Looks for evidence of simultaneous 3'/5' events (alt. both).  Assumes alt. 3' and
    alt. 5' nodes have already been identified."""
    # Using 5' events as a reference... 
    # For each 5' (donor site) event, we look at the nodes to see if any of them share
    # the same child (acceptor site).
    #
    # Examples of alt. both event with 4 nodes:
    #       ||||||||||----------------------|||||||||||||
    #       ||||||||----------------------|||||||||||||
    # Example with 6 nodes:
    #       ||||||||||--------------|||||||||||||
    #       ||||||------------|||||||||||||
    #       |||||||||||||--------------|||||||||||||
    # An alt. both event, and an alt. 5' event:
    #       ||||||||||----------------------|||||||||||||  alt.5' only
    #       ||||||||------------------------|||||||||||
    #       ||||||||||||||----------------------|||||||||||||  alt. both
    #
    # Note that the first two donors have edges to the same acceptor site
    # while the third donor has its own unique acceptor site.  Thus the rule to
    # apply is: any node within an alt 5' group that connects to a unique
    # acceptor that is not shared with any other member of that group becomes
    # an alt. both event.
    #
    # An alt. both event, and an alt. 3' event:
    #       ||||||||||----------------------|||||||||||||  alt.3' only
    #       ||||||||||-------------------|||||||||||
    #       ||||||||||||----------------------|||||||||||||  alt. both
    #
    # Note that the first two acceptors have edges to the same donor site
    # while the third donor has its own unique acceptor site.  Thus the rule to
    # apply is: any node within an alt 3' group that connects to a unique
    # donor that is not shared with any other member of that group becomes
    # an alt. both event.
    alt5Events = altSiteEventList(nodes, ALT5_ABBREV, verbose=verbose)
    for event in alt5Events :
        for n in event :
            nodeAcceptors  = set([c.acceptorEnd() for c in n.children])
            otherAcceptors = set([c.acceptorEnd() for m in event for c in m.children if m != n])
            shared         = nodeAcceptors & otherAcceptors
            if not shared :
                n.removeAltForm(ALT5_ABBREV)
                n.addAltForm(ALTB5_ABBREV)

    alt3Events = altSiteEventList(nodes, ALT3_ABBREV)
    for event in alt3Events :
        for n in event :
            nodeDonors  = set([c.donorEnd() for c in n.parents])
            otherDonors = set([c.donorEnd() for m in event for c in m.parents if m != n])
            shared      = nodeDonors & otherDonors
            if not shared :
                n.removeAltForm(ALT3_ABBREV)
                n.addAltForm(ALTB3_ABBREV) 

def detectAltDonor(n, nodes, edges, verbose=False) :
    """Looks for evidence that identifies the given node as an alternate donor."""
    # Inference rules:
    #  - Establish set of all nodes that overlap the given one
    #  - Get the set of nodes flanking edges that the node contains (retained introns)
    #  - Get the set of all nodes that overlap all children of the given one (other retained introns)
    #  - Subtract those overlapping children from the first set
    #  - For those that remain, if any have acceptor sites different from this one, mark it as an alternate acceptor
    if not n.children : return
    allOverlaps    = set([o for o in nodes if o.children and overlap(o,n)])
    containedEdges = set([e for e in edges if containsEdge(n,e)])
    flankingNodes  = set([e.parent for e in containedEdges if overlapsAll(n,e.parent.children)])
    childOverlaps  = set([o for o in nodes if overlapsAll(o,n.children)])
    overlapSet     = allOverlaps - childOverlaps - flankingNodes

    for o in overlapSet :
        if o.donorEnd() != n.donorEnd() :
            n.addAltForm(ALT5_ABBREV)
            return

def detectRetainedIntron(n, edges) :
    """Looks for evidence in the edge list that identifies the given node as a retained intron."""
    for e in edges :
        if containsEdge(n,e) :
            n.addAltForm(IR_ABBREV)
            return

def detectSkippedExon(n, edges) :
    """Looks for evidence in the edge list that identifies the given node as a skipped exon."""
    annotation = ES_ABBREV
    if not n.parents : annotation = ALTI_ABBREV
    if not n.children : annotation = ALTT_ABBREV
    for e in edges :
        if e.minpos < n.minpos and n.maxpos < e.maxpos :
            n.addAltForm(annotation)
            return

#------------------------------------------------------------------------
# Graph manipulation functions:
def commonAS(A,B) :
    """Returns a graph that contains just the nodes and AS
    events found in both A and B.  Uses the fact that
    A^B = AUB - A\B - B\A"""
    combined      = A.union(B)
    [AminB,BminA] = diffAS(A,B)
    subgraph      = graphMinusAS(combined,AminB)
    return graphMinusAS(subgraph,BminA)

def diffAS(A,B) :
    """Detects alternative splicing differences between two graphs.
    Returns two graphs: one that contains the nodes in A with AS not
    found in B and another that contains nodes in B without events in A."""
    return [graphMinusAS(A,B), graphMinusAS(B,A)]

def equivalentGraphs(A,B) :
    """Returns True if SpliceGraphs A and B are equivalent; False otherwise.
    This is distinct from the SpliceGraph.__eq__ method as it ignores
    attributes that may be changed in the prediction process."""
    Anodes = list(A.nodeDict.values())
    Bnodes = list(B.nodeDict.values())
    if len(Anodes) != len(Bnodes) : return False
    for n in Anodes :
        Aedges = childEdges(n)
        try :
            idx    = Bnodes.index(n)
            o      = Bnodes[idx]
            Bedges = childEdges(o)
            if Aedges != Bedges : return False
        except ValueError :
            return False
    return True

def getFirstGraph(f, **args) :
    """Returns just the first splice graph found in a file."""
    annotate = getAttribute('annotate', False, **args)
    try :
        result = next(SpliceGraphParser(f, **args))
        if annotate : result.annotate()
        return result
    except StopIteration :
        raise ValueError('No graph found in %s' % f)

def graphMinusAS(A,B) :
    """Detects alternative splicing differences between two graphs.
    Returns a graph that contains just the nodes in A that are
    annotated with AS not found in B.  Note that if A is a subgraph
    of B, this will return None."""
    # Get all nodes with AS annotations.  Note that we do
    # not report nodes missing from one graph or the other.
    nodesA = [n for n in A.resolvedNodes() if n.hasAS()]
    result = SpliceGraph(name=A.getName(), chromosome=A.chromosome, strand=A.strand)
    for a in nodesA :
        aAltSet = set(a.altForms())
        if not aAltSet : continue

        b = B.getNode(a.start, a.end)
        if b :
            bAltSet = set(b.altForms())
            aAltSet = aAltSet - bAltSet
            if not aAltSet : continue

        node = result.addNode(a.id, a.minpos, a.maxpos)
        for event in aAltSet :
            node.addAltForm(event)
    return result

def graphSubtract(A,B, **args) :
    """Returns a graph that represents the nodes and edges in A minus those in B."""
    resolvedOnly = getAttribute('resolvedOnly', True, **args)
    name    = '%s\%s' % (A.getName(), B.getName())
    result  = SpliceGraph(name=name, chromosome=A.chromosome, strand=A.strand)
    # Use full graph position range
    result.minpos = A.minpos
    result.maxpos = A.maxpos

    nodesA  = A.resolvedNodes() if resolvedOnly else list(A.nodeDict.values())
    nodesB  = B.resolvedNodes() if resolvedOnly else list(B.nodeDict.values())
    subsetA = []
    for a in nodesA :
        try :
            bidx = nodesB.index(a)
        except ValueError :
            newNode = result.addNode(a.id, a.minpos, a.maxpos)
            newNode.attrs = a.attrs
            subsetA.append(a)

    for a in subsetA :
        for c in a.children :
            if c in subsetA :
                result.addEdge(a.id,c.id)

    return result

def consistencyCoefficients(A,B) :
    """Deprecated.  Use recall(A,B) instead."""
    return recall(A,B)

def jaccardCoefficients(A,B) :
    """Returns the Jaccard coefficients for the similarity between the
    two graphs' vertex and edge sets, respectively.  (Recall that the
    Jaccard coefficient for a set is |A^B|/|AUB|.)  Note that if both
    A and B are the empty set, the Jaccard coefficient is 1."""
    setA  = set(A.resolvedNodes())
    setB  = set(B.resolvedNodes())
    AandB = setA.intersection(setB)
    AorB  = setA.union(setB)
    nodeCoefficient = float(len(AandB)) / len(AorB) if AorB else 1.0

    edgesA = edgeSet(A)
    edgesB = edgeSet(B)
    AandB  = edgesA.intersection(edgesB)
    AorB   = edgesA.union(edgesB)
    edgeCoefficient = float(len(AandB)) / len(AorB) if AorB else 1.0

    return nodeCoefficient, edgeCoefficient

def nodeString(nodeList) :
    """Returns a string representation of all node ids in a list."""
    return ','.join([n.id for n in nodeList])

def recall(A,B) :
    """Returns recall values for the similarity between graphs
    A and B for node and edge sets, respectively.  The recall
    for set A relative to B is |A^B|/|B|.  Note that if A and B
    are both empty, the recall value is 1."""
    nodesB = set(B.resolvedNodes())
    AandB  = set(A.resolvedNodes()) & nodesB
    nodeCoefficient = float(len(AandB)) / len(nodesB) if nodesB else 1.0

    edgesB = edgeSet(B)
    AandB  = edgeSet(A) & edgesB
    edgeCoefficient = float(len(AandB)) / len(edgesB) if edgesB else 1.0

    return nodeCoefficient, edgeCoefficient

def splitFiles(parser, directory=None) :
    """Given a parser on a file that contains multiple models, writes each
    model into its own file prefixed with its graph name.  First removes
    special characters from the name.  For example:
        'AT1G01448|alt 3/5' --> 'AT1G01448_alt_3_5_graph.gff'"""
    import os
    def fixName(s) :
        result = s.replace(' ','_')
        result = result.replace('|','_')
        result = result.replace('/','_')
        result = result.replace('\\','_')
        return result.replace('*','')

    for g in parser.__iter__() :
        outName = '%s_graph.gff' % fixName(g.getName())
        outFile = os.path.join(directry, outName) if directory else outName
        g.writeGFF(outFile)

def updateLeaf(A, B, **args) :
    """Method that updates leaf nodes in graph A with nodes from graph B prior to
    merging the two graphs.  We look for longer nodes with the same acceptor site."""
    uniqueLeaf = getAttribute('uniqueLeaf', False, **args)
    mergeDict = {}
    for leaf in A.getLeaves() :
        longer = [n for n in list(B.nodeDict.values()) if n.acceptorEnd() == leaf.acceptorEnd() and len(n) > len(leaf)]
        if not longer : continue
        if uniqueLeaf and len(longer) > 1 : continue
        maxLen  = max([len(n) for n in longer])
        longest = [n for n in longer if len(n) == maxLen][0]

        # If the new length matches any other nodes,
        # merge them after the loop has finished
        other = A.getNode(longest.minpos, longest.maxpos)
        if other :
            mergeDict[leaf] = other
        else :
            leaf.update(longest.minpos, longest.maxpos)

    # Merge and remove nodes that would have been duplicates
    for node in mergeDict :
        other = mergeDict[node]
        for p in node.parents :
            A.addEdge(p.id, other.id)
            p.removeChild(node)
        del A.nodeDict[node.id]

def updateRoot(A, B, **args) :
    """Method that updates root nodes in graph A using nodes from graph B prior to
    merging the two graphs.  We look for longer nodes that have the same donor site."""
    uniqueRoot = getAttribute('uniqueRoot', False, **args)

    mergeDict = {}
    for root in A.getRoots() :
        longer = [n for n in list(B.nodeDict.values()) if n.donorEnd() == root.donorEnd() and len(n) > len(root)]
        if not longer : continue
        if uniqueRoot and len(longer) > 1 : continue
        maxLen  = max([len(n) for n in longer])
        longest = [n for n in longer if len(n) == maxLen][0]

        # If the new length matches any other nodes,
        # merge them after the loop has finished
        other = A.getNode(longest.minpos, longest.maxpos)
        if other :
            mergeDict[root] = other
        else :
            root.update(longest.minpos, longest.maxpos)

    # Merge and remove nodes that would have been duplicates
    for node in mergeDict :
        other = mergeDict[node]
        for c in node.children :
            A.addEdge(other.id, c.id)
            c.removeParent(node)
        del A.nodeDict[node.id]

#========================================================================
# Class definitions
class Edge(object) :
    """Encapsulates an edge in the graph.  These are not stored with a splice
    graph, but used for splice graph creation and annotation."""

    def __init__(self, parent, child) :
        self.parent = parent
        self.child  = child
        self.pos    = sorted([parent.minpos,parent.maxpos, child.minpos,child.maxpos])
        # minpos/maxpos are related to the edge itself, not the outer positions
        self.minpos = self.pos[1]
        self.maxpos = self.pos[2]

    def __eq__(self, o) :
        # Using minpos/maxpos allows an edge to be compared with an exon.
        return self.minpos == o.minpos and self.maxpos == o.maxpos and self.pos[0] == o.pos[0] and self.pos[3] == o.pos[3]

    def __cmp__(self, o) :
        result = self.pos[0] - o.pos[0]
        if not result : result = self.pos[1] - o.pos[1]
        if not result : result = self.pos[2] - o.pos[2]
        if not result : result = self.pos[3] - o.pos[3]
        return result

    def __getitem__(self, i) :
        return self.pos[i]

    def __hash__(self) :
        return self.__str__().__hash__()

    def __len__(self) :
        return self.maxpos-self.minpos

    def overlaps(self, o) :
        """Returns true if the edge portions overlap; false otherwise."""
        return self.maxpos > o.minpos and self.minpos < o.maxpos

    def sameEdge(self, o) :
        """Returns true if the edge portions are the same; false otherwise."""
        return self.minpos == o.minpos and self.maxpos == o.maxpos

    def __repr__(self) :
        return "%d,%d,%d,%d" % tuple(self.pos)

    def __str__(self) :
        return "%d,%d" % tuple(self.pos[1:3])

class NullNode(object) :
    """Null node object encapsulates the most basic information about a node."""
    def __init__(self, start, end) :
        self.start  = start
        self.end    = end
        self.minpos = min(start,end)
        self.maxpos = max(start,end)

class SpliceGraphNode(object) :
    """This is the node class to use for constructing splice graphs for GFF input/output."""
    def __init__(self, id, start, end, strand, chrom, parents=[], children=[]) :
        self.id         = id
        self.strand     = strand
        self.chromosome = chrom
        self.minpos     = min(start,end)
        self.maxpos     = max(start,end)
        (self.start,self.end) = (self.minpos,self.maxpos) if strand == '+' else (self.maxpos,self.minpos)
        self.parents    = list(parents)
        self.children   = list(children)
        self.attrs      = {}
        self.altFormSet = set()
        self.isoformSet = set()
        # Added for adjustable ranges
        self.origStart  = self.start
        self.origEnd    = self.end

    def acceptorEnd(self) :
        """Returns the position of the exon's 5' (upstream) end."""
        return self.start

    def addAltForm(self, form) :
        """Adds an AS form to the node's list of forms."""
        if len(form.strip()) == 0 : return
        self.altFormSet.add(form)
        for x in self.altFormSet : assert(len(x) > 0)
        self.attrs[AS_KEY] = ','.join(self.altFormSet)

    def addAttribute(self, key, value) :
        """Adds the given key-value pair to a node's attribute list.  If the
        key is already in the attribute dictionary, the new value overwrites the old."""
        self.attrs[key] = value

    def addChild(self, c) :
        """Adds a child node (edge) to a node.  If the child is already known the child
        list is unchanged."""
        if c not in self.children :
            self.children.append(c)
            c.addParent(self)

    def addCodon(self, codon, codonType) :
        """If the codon falls within the node, its start position is added to the set;
        otherwise the list is unchanged."""

        if (type(codon) != tuple) or (len(codon) != 2) :
            raise ValueError('Codons must be provided as a duple of (start,end) positions: received %s' % str(codon))

        pos = min(codon) if self.strand == '+' else max(codon)
        if self.contains(pos) :
            assert(type(pos) == int)
            self.attrs.setdefault(codonType, set())
            self.attrs[codonType].add(pos)

    def addEndCodon(self, codon) :
        """Adds an end codon to the node."""
        self.addCodon(codon, END_CODON_KEY)

    def addFormsFromString(self, s) :
        """Adds AS forms to node list from a string representation."""
        for form in s.split(',') :
            self.addAltForm(form.strip())

    def addIsoformString(self, isoformString) :
        """Adds a set of isoform names to the node's set of isoforms."""
        for iso in isoformString.split(',') :
            self.addIsoform(iso)

    def addIsoform(self, isoform) :
        """Adds an isoform name to the node's set of isoforms."""
        if isoform is None : raise ValueError('Received illegal isoform for %s' % self.id)
        self.isoformSet.add(isoform)
        self.attrs[ISO_KEY] = ','.join(self.isoformSet)

    def addParent(self, p) :
        """Adds a parent node (edge) to a node.  If the parent is already known the parent
        list is unchanged."""
        if p not in self.parents :
            self.parents.append(p)

    def addStartCodon(self, codon) :
        """If the codon position falls within the node it is added to the set;
        otherwise the list is unchanged."""
        self.addCodon(codon, START_CODON_KEY)

    def altForms(self) :
        """Returns a list of AS forms associated with the node."""
        for x in self.altFormSet : assert(len(x) > 0)
        return list(self.altFormSet)

    def altFormString(self) :
        """Returns a string of AS forms associated with the node."""
        try :
            return self.attrs[AS_KEY]
        except KeyError :
            return ''

    def attributeString(self) :
        """Returns a string of all attributes associated with this node, suitable for a
        GFF attributes field.  For example:
                  ID=ei_27;Parent=ei_26;AltForm=Alt. 3',Retained Intron"""
        result = ''
        for k in list(self.attrs.keys()) :
            if k == AS_KEY and not self.altFormSet : continue
            if k == ISO_KEY and not self.isoformSet : continue

            if result : result += ';'
            if k in [START_CODON_KEY, END_CODON_KEY] :
                result += '%s=%s' % (k,self.codonString(k))
            else :
                result += '%s=%s' % (k,self.attrs[k])
        return result

    def branchingFactor(self) :
        return max(len(self.parents), len(self.children))

    def __cmp__(self, other) :
        """Permits sorting based on minimum node position.  Ties are broken by the
        shorter of the two nodes."""
        if self.minpos == other.minpos :
            return self.maxpos - other.maxpos
        else :
            return self.minpos - other.minpos

    def codons(self, codonType) :
        """Returns a list of codon positions within the node, or
        an empty list if there are no codons of the given type."""
        try :
            return sorted(self.attrs[codonType])
        except KeyError :
            return []

    def codonString(self, codonType) :
        """Returns a list of codon positions within the node as a string,
        or None if there are no codons of the given type."""
        try :
            codons = sorted(self.attrs[codonType])
            return ','.join(['%d'%x for x in codons])
        except KeyError :
            pass

    def contains(self, pos) :
        """Returns true if the given position falls within the node; false otherwise."""
        return self.minpos <= pos <= self.maxpos

    def donorEnd(self) :
        """Returns the position of the exon's 3' (downstream) end."""
        return self.end

    def downstreamOf(self, pos) :
        """Returns true if the node is downstream of the position; false otherwise."""
        return (self.minpos > pos) if self.strand == '+' else (self.maxpos < pos)

    def __eq__(self, other) :
        """Node equality is based solely on its start/end positions."""
        return self.minpos == other.minpos and self.maxpos == other.maxpos

    def endCodons(self) :
        """Returns a list of end codon start positions within the node."""
        return self.codons(END_CODON_KEY)

    def endCodonString(self) :
        """Returns a list of end codon start positions within the node as a string."""
        return self.codonString(END_CODON_KEY)

    def gffString(self) :
        """Returns a GFF-formatted string representing the given graph feature."""
        if self.parents :
            # Example: chr1	SpliceGraph	child	37373	37398	.	-	.	ID=gm_20;Parent=gm_19
            result = "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s" % \
                (self.chromosome, SOURCE_NAME, CHILD_REC, self.minpos, self.maxpos, self.strand, self.id, nodeString(self.parents))
        else :
            # Example: chr1	SpliceGraph	parent	37569	37757	.	-	.	ID=gm_19
            result = "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\tID=%s" % \
                (self.chromosome, SOURCE_NAME, PARENT_REC, self.minpos, self.maxpos, self.strand, self.id)
        if self.attrs :
            result += ';%s' % self.attributeString()
        return result

    def hasAS(self) :
        """Returns true if the node has evidence of AS; false otherwise."""
        for x in self.altFormSet : assert(len(x) > 0)
        return len(self.altFormSet & AS_ABBREV_SET) > 0

    def hasDisposition(self, disposition) :
        """Returns true if the node has the given disposition; false otherwise."""
        try :
            return (self.attrs[DISPOSITION_KEY] == disposition)
        except KeyError :
            return False

    def __hash__(self) :
        hashString = '%d,%d' % (self.minpos, self.maxpos)
        return hashString.__hash__()

    def isAltAcceptor(self) :
        """Returns true if the node represents a retained intron; false otherwise."""
        for x in self.altFormSet : assert(len(x) > 0)
        return ALT3_ABBREV in self.altFormSet

    def isAltDonor(self) :
        """Returns true if the node represents a retained intron; false otherwise."""
        for x in self.altFormSet : assert(len(x) > 0)
        return ALT5_ABBREV in self.altFormSet

    def isKnown(self) :
        """Returns true if the node is known from the gene model; false otherwise."""
        return self.hasDisposition(KNOWN_NODE)

    def isLeaf(self) :
        """Returns true if the node is a leaf node; false otherwise."""
        return len(self.children) == 0

    def isPredicted(self) :
        """Returns true if the node is predicted; false otherwise."""
        return self.hasDisposition(PREDICTED_NODE)

    def isRetainedIntron(self) :
        """Returns true if the node represents a retained intron; false otherwise."""
        for x in self.altFormSet : assert(len(x) > 0)
        return IR_ABBREV in self.altFormSet

    def isRoot(self) :
        """Returns true if the node is a root node; false otherwise."""
        return len(self.parents) == 0

    def isSkippedExon(self) :
        """Returns true if the node represents a retained intron; false otherwise."""
        for x in self.altFormSet : assert(len(x) > 0)
        return ES_ABBREV in self.altFormSet

    def isUnresolved(self) :
        """Returns true if the node is unresolved; false otherwise."""
        return self.hasDisposition(UNRESOLVED_NODE)

    def isoformList(self) :
        """Returns a list of isoforms associated with the node."""
        return list(self.isoformSet)

    def isoformString(self) :
        """Returns a string representation of isoforms associated with the node."""
        try :
            return self.attrs[ISO_KEY]
        except KeyError :
            return None

    def __len__(self) :
        """Returns the length of the genomic region represented by the node."""
        return self.maxpos-self.minpos+1

    def __lt__(self, other) :
        """Permits sorting based on minimum node position.  Ties are broken by the
        shorter of the two nodes."""
        if self.minpos == other.minpos :
            return self.maxpos < other.maxpos
        else :
            return self.minpos < other.minpos

    def putativeChildren(self) :
        """Returns a set of putative children for an unresolved node,
        or an empty set if there are none or the node is not unresolved."""
        try :
            return self.attrs[PUTATIVE_CHILDREN].split(',')
        except KeyError :
            return set()

    def putativeParents(self) :
        """Returns a set of putative parents for an unresolved node,
        or None if there are none or the node is not unresolved."""
        try :
            return self.attrs[PUTATIVE_PARENTS].split(',')
        except KeyError :
            return set()

    def removeAltForm(self, form) :
        """Adds an AS form to the node's list of forms."""
        self.altFormSet.discard(form)
        for x in self.altFormSet : assert(len(x) > 0)
        self.attrs[AS_KEY] = ','.join(self.altFormSet)

    def removeChild(self, child) :
        """Removes the child from this node's list."""
        self.children.remove(child)

    def removeParent(self, parent) :
        """Removes the parent from this node's list."""
        self.parents.remove(parent)

    def __repr__(self) :
        if self.parents :
            return '%s %d-%d children=%s parents=%s' % (self.id, self.start, self.end, nodeString(self.children), nodeString(self.parents))
        else :
            return 'Root %s %d-%d children=%s' % (self.id, self.start, self.end, nodeString(self.children))

    def setUnresolved(self, preserve=[], acceptors=[], donors=[]) :
        """
        Makes a node unresolved by setting its disposition attribute and
        by removing edges between the node and parents and children in the graph.
        Caller may preserve some display edges by listing nodes in a preserve list.
        """
        self.addAttribute(DISPOSITION_KEY, UNRESOLVED_NODE)
        if acceptors : self.addAttribute(ACCEPTORS_KEY, listString(acceptors))
        if donors    : self.addAttribute(DONORS_KEY, listString(donors))

        for c in self.children :
            if c.id in preserve : continue
            c.removeParent(self)
            self.children.remove(c)

        for p in self.parents :
            if p.id in preserve : continue
            p.removeChild(self)
            self.parents.remove(p)

    def startCodons(self) :
        """Returns a list of start codon start positions within the node."""
        return self.codons(START_CODON_KEY)

    def startCodonString(self) :
        """Returns a list of start codon start positions within the node as a string."""
        return self.codonString(START_CODON_KEY)

    def __str__(self) :
        if self.parents :
            return '%s %d-%d (parents:%d/children:%d)' % (self.id, self.start, self.end, len(self.parents), len(self.children))
        else :
            return 'Root %s %d-%d (children:%d)' % (self.id, self.start, self.end, len(self.children))

    def update(self, minpos, maxpos) :
        """Mutates the node based on the revised min/max positions."""
        assert(minpos <= maxpos)
        self.minpos = minpos
        self.maxpos = maxpos
        self.start  = minpos if self.strand == '+' else maxpos
        self.end    = maxpos if self.strand == '+' else minpos

    def upstreamOf(self, pos) :
        """Returns true if the node is upstream of the position; false otherwise."""
        return (self.maxpos < pos) if self.strand == '+' else (self.minpos > pos)

class SpliceGraph(object) :
    """Main class for storing splice graph data."""
    def __init__(self, name, chromosome, strand) :
        """
        Instantiate a splice graph.

        :Parameters:
           'name'       - gene name to associate with the graph
           'chromosome' - chromosome/scaffold name to associate with the graph (GFF format requirement)
           'strand'     - strand associated with the graph
        """
        self.chromosome = chromosome
        self.strand     = strand
        self.nodeDict   = {}
        self.attrs      = {}
        self.minpos     = MAXINT
        self.maxpos     = 0
        # NB: guarantees attributes exist for the graph
        self.setName(name)

    def addNode(self, newId, start, end) :
        """Add a node to the splice graph if it doesn't already exist."""
        try :
            # Look for another node with the same start/end positions
            tmpNode  = NullNode(start,end)
            allNodes = list(self.nodeDict.values())
            idx      = allNodes.index(tmpNode)
            node     = allNodes[idx]
        except ValueError :
            self.nodeDict[newId] = SpliceGraphNode(newId, start, end, self.strand, self.chromosome)
            node              = self.nodeDict[newId]
            self.minpos       = min(self.minpos, node.minpos)
            self.maxpos       = max(self.maxpos, node.maxpos)
     
        return node

    def addCodons(self, codonList, codonType) :
        """Adds codons to every node in the graph."""
        for codon in codonList :
            for node in list(self.nodeDict.values()) :
                node.addCodon(codon, codonType)

    def addEdge(self, pid, cid) :
        """Add an edge to the splice graph if it doesn't already exist."""
        try :
            parent = self.nodeDict[pid]
        except KeyError :
            raise Exception("Error adding edge from node %s: node not found in graph" % pid)

        try :
            child = self.nodeDict[cid]
        except KeyError :
            raise Exception("Error adding edge to node %s: node not found in graph" % cid)

        parent.addChild(child)

    def addEndCodons(self, codonList) :
        """Adds end codons to all nodes in the graph."""
        self.addCodons(codonList, END_CODON_KEY)

    def addStartCodons(self, codonList) :
        """Adds start codons to all nodes in the graph."""
        self.addCodons(codonList, START_CODON_KEY)

    def adjust(self, adjustment) :
        """
        Shifts the positions of all nodes in the graph by the given amount.
        This is useful for example for converting a graph based on positions
        [0,n-1] for an environment that uses positions [1,n].
        """
        for n in list(self.nodeDict.values()) :
            n.start  += adjustment
            n.end    += adjustment
            n.minpos += adjustment
            n.maxpos += adjustment

        self.minpos += adjustment
        self.maxpos += adjustment

    def altForms(self) :
        """Returns a set of all AS forms found in the graph."""
        result = set()
        for n in list(self.nodeDict.values()) :
            result.update(n.altFormSet)
        for x in result : assert(len(x) > 0)
        return result

    def annotate(self, verbose=False) :
        """Updates all AS annotations for the graph."""
        edges = edgeSet(self)
        nodes = self.resolvedNodes()
        for n in nodes :
            # First reset the node AS annotations
            n.altFormSet    = set()
            n.attrs[AS_KEY] = ''
            for x in n.altFormSet : assert(len(x) > 0)

            # Assign annotations to the node
            detectSkippedExon(n, edges)
            detectRetainedIntron(n, edges)
            if n.children : detectAltDonor(n, nodes, edges, verbose=verbose)
            if n.parents  : detectAltAcceptor(n, nodes, edges)
            
        # Removed (11/21/2011) after conversation with Asa:
        ## detectAltBoth(nodes)
            
    def attributeString(self) :
        """Returns a string of all graph attributes, for GFF attributes fields."""
        return ';'.join(['%s=%s'%(k,self.attrs[k]) for k in sorted(self.attrs.keys())])

    def branchingStats(self) :
        """Returns the min, max and average branching factor for the nodes in the graph."""
        branches = [n.branchingFactor() for n in self.resolvedNodes()]
        # NB: len(branches) has to be at least 1 since there must be at least 1 node to have a graph
        try :
            avg = float(sum(branches))/len(branches)
        except ZeroDivisionError :
            raise ValueError('Attempted to compute branching factor on an empty graph')

        return min(branches), max(branches), avg

    def __cmp__(self, other) :
        """Permits sorting graphs based on minimum position.  Ties are broken by the
        shorter of the two graphs."""
        if self.minpos == other.minpos :
            return self.maxpos - other.maxpos
        else :
            return self.minpos - other.minpos

    def deleteNode(self, n) :
        """Removes a node from a graph, along with all edges attached to it.
        Returns the node that was deleted."""
        # Caller must catch KeyError if this fails
        node = self.nodeDict[n.id]
        for p in node.parents :
            p.removeChild(node)
        for c in node.children :
            c.removeParent(node)
        del self.nodeDict[n.id]
        return node

    def distinctAltEvents(self) :
        """Returns a list of distinct AS events found in the graph.
        Each event is reported as a tuple: (node-start,node-end,event-type,count).
        Alt. 5'/Alt. 3' events are the only ones with a count greater than 1 and
        are reported using the shortest node involved."""
        result    = []
        alt5Nodes = [n for n in self.resolvedNodes() if ALT5_ABBREV in n.altFormSet]
        while alt5Nodes :
            n      = alt5Nodes.pop()
            others = [o for o in alt5Nodes if o.acceptorEnd() == n.acceptorEnd()]
            for o in others :
                if len(o) < len(n) : n = o
                alt5Nodes.remove(o)
            result.append((n.start, n.end, ALT5_ABBREV, len(others)+1))

        alt3Nodes = [n for n in self.resolvedNodes() if ALT3_ABBREV in n.altFormSet]
        while alt3Nodes :
            n      = alt3Nodes.pop()
            others = [o for o in alt3Nodes if o.donorEnd() == n.donorEnd()]
            for o in others :
                if len(o) < len(n) : n = o
                alt3Nodes.remove(o)
            result.append((n.start, n.end, ALT3_ABBREV, len(others)+1))

        for n in self.resolvedNodes() :
            for asType in NON_35_ABBREVS :
                if asType in n.altFormSet :
                    result.append((n.start, n.end, asType,1))

        return result

    def downstreamOf(self, a, b) :
        """Returns true if a is downstream of b; false otherwise."""
        return (a>b) if self.strand == '+' else (b>a)

    def duplicate(self) :
        """Returns a duplicate of the current graph.  This performs
        the same function as copy.deepcopy(), but without the risk
        of stack overflow for large graphs."""
        result = SpliceGraph(self.getName(), self.chromosome, self.strand)
        result.minpos = self.minpos
        result.maxpos = self.maxpos
        result.attrs.update(self.attrs)

        # First pass: create nodes in graph
        for n in list(self.nodeDict.values()) :
            newNode = result.addNode(n.id, n.start, n.end)
            newNode.attrs.update(n.attrs)
            newNode.isoformSet     = set(n.isoformSet)
            newNode.attrs[ISO_KEY] = ','.join(newNode.isoformSet)
            
        # Second pass: create edges in graph
        for n in list(self.nodeDict.values()) :
            for c in n.children :
                result.addEdge(n.id, c.id)

        return result

    def __eq__(self, other) :
        """Returns true if two splice graphs are identical; false otherwise."""
        nodes      = list(self.nodeDict.values())
        otherNodes = list(other.nodeDict.values())
        if len(nodes) != len(otherNodes) : return False
        for n in nodes :
            edges = childEdges(n)
            try :
                idx        = otherNodes.index(n)
                o          = otherNodes[idx]
                otherEdges = childEdges(o)
                if edges != otherEdges : return False
                if n.attrs != o.attrs : return False
            except ValueError :
                return False
        return True

    def expandedModel(self, name, oldRange, newRange) :
        """Marks a graph when it expands a gene model beyond its original bounds.
        oldRange and newRange must be (minpos,maxpos) pairs."""
        if len(oldRange) != 2 or type(oldRange[0]) != int or type(oldRange[1]) != int :
            raise ValueError('old range must be 2 int values')

        if len(newRange) != 2 or type(newRange[0]) != int or type(newRange[1]) != int :
            raise ValueError('new range must be 2 int values')

        self.attrs[ALT_MODEL_KEY] = "%s from (%d,%d) to (%d,%d)" % \
                (name, oldRange[0], oldRange[1], newRange[0], newRange[1])

    def getAcceptors(self) :
        """Returns a list of distinct acceptor sites in the graph."""
        result = set([acceptor(n) for n in list(self.nodeDict.values()) if not n.isRoot()])
        return list(result)

    def getDonors(self) :
        """Returns a list of distinct donor sites in the graph."""
        result = set([donor(n) for n in list(self.nodeDict.values()) if not n.isLeaf()])
        return list(result)

    def getNode(self, start, end) :
        """Returns the node represented by the given positions, if it exists."""
        # NB: May want to add option to locate only resolved nodes
        try :
            # Look for a node with the same start/end positions
            tmpNode  = NullNode(start,end)
            allNodes = list(self.nodeDict.values())
            idx      = allNodes.index(tmpNode)
            return allNodes[idx]
        except ValueError :
            return None

    def getLeaves(self) :
        """Returns a list of nodes that have no children."""
        return [n for n in list(self.nodeDict.values()) if not n.children]

    def getName(self) :
        """Central method for retrieving the name.  This is a temporary
        hack until a better method is implemented."""
        return self.attrs[ID_ATTR]

    def getRoots(self) :
        """Returns a list of nodes that have no parents."""
        return [n for n in list(self.nodeDict.values()) if not n.parents]

    def gffString(self, node=None) :
        """Returns a GFF-formatted string representing the given graph feature."""
        result = "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s" % \
                (self.chromosome, SOURCE_NAME, GENE_REC, self.minpos, self.maxpos, self.strand, self.attributeString())
        return result

    def hasAS(self) :
        """Returns true if the graph has evidence of AS; false otherwise."""
        for n in self.resolvedNodes() :
            if n.hasAS() : return True
        return False

    def isEmpty(self) :
        """Returns true if the graph is empty (has no nodes or edges); false otherwise."""
        return (len(self.nodeDict) == 0)

    def isoformDict(self) :
        """Returns a dictionary where each isoform name is
        keyed to a sorted list of nodes in the isoform."""
        result = {}
        nodes  = self.resolvedNodes()
        nodes.sort(reverse=(self.strand=='-'))
        for n in nodes :
            for iso in n.isoformSet :
                result.setdefault(iso,[])
                result[iso].append(n)
        return result

    def isoformGraph(self, isoNames, **args) :
        """Returns a splice graph that represents just the given isoforms.
        Isoforms may be given as a comma-separated string, a list or a set.
        Setting 'fuzzyMatch' to True allows you to specify part of an isoform
        name in the list; it will match all isoforms that contain that part."""
        fuzzyMatch = getAttribute('fuzzyMatch', False, **args)

        isoNames = asList(isoNames)
        result   = None
        for name in isoNames :
            if fuzzyMatch :
                nodes = []
                query = name.upper()
                for n in self.resolvedNodes() :
                    for x in n.isoformSet :
                        if x.upper().find(query) >= 0 :
                            nodes.append(n)
            else :
                nodes = [n for n in self.resolvedNodes() if name in n.isoformSet]

            if not nodes :
                raise ValueError('%s not found in graph %s (isoforms %s)\n' % (name, self.getName(), ','.join(self.isoforms())))

            graph = SpliceGraph(name, self.chromosome, self.strand)
            prev  = None
            for n in nodes :
                curr                = graph.addNode(n.id, n.minpos, n.maxpos)
                curr.attrs          = dict(n.attrs)
                curr.isoformSet     = set([name])
                curr.attrs[ISO_KEY] = ','.join(curr.isoformSet)
                if prev : graph.addEdge(prev.id,curr.id)
                prev = curr
            result = graph if result is None else result.union(graph)
        if result : result.setName('|'.join(isoNames))
        return result

    def isoforms(self) :
        """Returns a list of all isoform strings associated with nodes in the graph."""
        isoSet = set()
        for n in self.resolvedNodes() :
            isoSet = isoSet|n.isoformSet
        return sorted(list(isoSet))

    def __len__(self) :
        return self.maxpos-self.minpos+1 if self.maxpos > self.minpos else 0

    def mergeNodes(self, a, b) :
        """Merges two nodes by using the lowest min position and the highest
        max position of the two nodes.  The new node will have the same identifier
        as the first original node.  All parents and all children are combined and the
        original nodes will be deleted from the graph.  Returns the original nodes."""
        node1    = self.nodeDict[a]
        node2    = self.nodeDict[b]
        minpos   = min(node1.minpos, node2.minpos)
        maxpos   = max(node1.maxpos, node2.maxpos)

        # Caller must catch KeyError if this fails
        self.deleteNode(node1)
        self.deleteNode(node2)

        # Either create a new node or grab an existing one
        newNode = self.addNode(node1.id, minpos, maxpos)

        # Only add parents/children that are still in the graph
        for c in node1.children+node2.children :
            if c.id in self.nodeDict :
                newNode.addChild(c)
        for p in node1.parents+node2.parents :
            if p.id in self.nodeDict :
                p.addChild(newNode)

        self.nodeDict[newNode.id] = newNode
        return (node1,node2)

    def nameString(self) :
        """Returns the graph name.  Graphs that have been merged will have more than
        one name given as a hyphenated list."""
        names = set(self.getName().split(','))
        if len(names) > 1 :
            return 'Merged ' + '+'.join(names)
        else :
            return self.getName()

    def __ne__(self, other) :
        """Returns true if two splice graphs are different; false otherwise."""
        return not self.__eq__(other)

    def __str__(self) :
        return '%s (%s) %d-%d (%d nodes)' % (self.getName(), self.strand, self.minpos, self.maxpos, len(self.nodeDict))

    def resolvedNodes(self) :
        """Returns a list of all resolved nodes: those actually used in the graph."""
        return [x for x in list(self.nodeDict.values()) if not x.isUnresolved()]

    def setName(self, name) :
        """Provides a central method for storing the name in both places.
        This is a temporary hack until a better method is implemented."""
        self.name           = name
        self.attrs[ID_ATTR] = name

    def union(self, other, **args) :
        """Merges two graphs into one.  Adds nodes that are distinct and merges nodes
        that have the same start/end positions.  Merged nodes will share AS
        information.  New attributes may be added, but conflicts will be resolved
        in favor of the existing graph."""
        keepName  = getAttribute('keepName', False, **args)
        mergeEnds = getAttribute('mergeEnds', False, **args)
        idgen     = getAttribute('idGenerator', idFactory('%s_U'%self.getName()), **args)

        # establish correct strand
        strand = self.strand
        if self.strand in VALID_STRANDS :
            if other.strand in VALID_STRANDS and self.strand != other.strand :
                raise AttributeError('Cannot merge graphs with conflicting strands.\n')
        elif other.strand in VALID_STRANDS :
            strand = other.strand

        # update graph name
        if keepName :
            newName = self.getName()
        else :
            names = set(self.getName().split(','))
            names.add(other.getName())
            newName = ','.join(names)

        if mergeEnds :
            updateRoot(other, self)
            updateLeaf(other, self)

        result   = SpliceGraph(newName, self.chromosome, self.strand)
        allNodes = self.resolvedNodes() + other.resolvedNodes()
        for node in allNodes :
            newNode = result.addNode(next(idgen), node.minpos, node.maxpos)
            for c in node.children :
                cNode = result.addNode(next(idgen), c.minpos, c.maxpos)
                newNode.addChild(cNode)
            for p in node.parents :
                pNode = result.addNode(next(idgen), p.minpos, p.maxpos)
                pNode.addChild(newNode)
            for iso in node.isoformSet :
                newNode.addIsoform(iso)
            for k in node.attrs :
                if k == AS_KEY :
                    for f in node.altForms() :
                        newNode.addAltForm(f)
                elif k != ISO_KEY :
                    newNode.attrs[k] = node.attrs[k]
        return result

    def unresolvedNodes(self) :
        """Returns a list of all unresolved nodes in the graph."""
        return [x for x in list(self.nodeDict.values()) if x.isUnresolved()]

    def upstreamOf(self, a, b) :
        """Returns true if a is upstream of b; false otherwise."""
        return (a<b) if self.strand == '+' else (b<a)

    def validate(self, halt=False) :
        """Returns None if the splicegraph is valid; otherwise returns a reason it is invalid."""
        reason     = None
        allNodes   = list(self.nodeDict.values())
        ## allNodes   = self.resolvedNodes()
        nodeSet    = set(allNodes)
        nodeList   = list(nodeSet)
        allNodeIds = [x.id for x in allNodes]
        uniqueIds  = [x.id for x in nodeSet]
        roots      = self.getRoots()
        leaves     = self.getLeaves()
        if len(roots) == 0 or len(leaves) == 0 :
            reason = 'Graph is missing roots or leaves (%d roots, %d leaves)' % (len(roots), len(leaves))

        for n in allNodes :
            if reason : break
            # Detect duplicate nodes: only one will appear in set
            if not n.id in uniqueIds :
                other  = nodeList[nodeList.index(n)]
                reason = 'Duplicate node %s (%d-%d) matches %s (%d-%d)' % (n.id, n.minpos, n.maxpos, other.id, other.minpos, other.maxpos)

            # Detect invalid child/parent ids
            for o in n.children :
                if o.id not in allNodeIds :
                    reason = 'Node %s: child %s not in graph' % (n.id, o.id)
                elif o.id not in uniqueIds :
                    twin   = nodeList[nodeList.index(o)]
                    reason = 'Node %s: child %s (%d-%d) is not unique (%s %d-%d)' % \
                            (n.id, o.id, o.minpos, o.maxpos, twin.id, twin.minpos, twin.maxpos)

            for o in n.parents :
                if o.id not in allNodeIds :
                    reason = 'Node %s parent %s not in graph' % (n.id, o.id)
                elif o.id not in uniqueIds :
                    twin   = nodeList[nodeList.index(o)]
                    reason = 'Node %s parent %s (%d-%d) is not unique (%s %d-%d)' % \
                            (n.id, o.id, o.minpos, o.maxpos, twin.id, twin.minpos, twin.maxpos)

        if reason and halt :
            import sys
            sys.stderr.write('\n** Error in graph:\n')
            sys.stderr.write('  %s\n' % '\n  '.join(['%s'%x for x in allNodes]))
            raise ValueError('Illegal graph:\n%s' % reason)
        elif reason :
            return reason

    def writeGFF(self, fileRef, haltOnError=False) :
        """Writes a splice graph to a file in GFF format.  The file may be given
        either as an output stream or as a file path.
        """
        reason = self.validate()
        if reason :
            if haltOnError :
                raise ValueError('Cannot write invalid splice graph %s to file:\n"%s"\n' % (self.getName(),reason))
            else :
                sys.stderr.write('** Warning: splice graph for %s is invalid:\n"%s"\n' % (self.getName(),reason))

        roots     = set(self.getRoots())
        other     = set(self.nodeDict.values()) - roots
        newFile   = (type(fileRef) == str)
        outStream = open(fileRef,'w') if newFile else fileRef

        # First write the graph, then parents, then remaining nodes
        outStream.write('%s\n' % self.gffString())
        for p in roots  :
            outStream.write('%s\n' % p.gffString())
        for o in other  :
            outStream.write('%s\n' % o.gffString())

        if newFile :
            outStream.close()
        return True

class SpliceGraphParser(object) :
    """Class that parses a GFF file filled with splice graphs and provides an
    iterator over each graph in the file."""

    def __init__(self, fileRef, **args) :
        """Parses a GFF file filled with splice graphs and returns each graph in an iterator.
        The file may be given either as an input stream or as a file path."""
        self.verbose = getAttribute('verbose', False, **args)

        if type(fileRef) == type('') :
            self.instream = ezopen(fileRef)
        else :
            self.instream = fileRef

        if self.instream is None :
            raise ValueError('No input file stream given.')

        self.graphDict = {}
        self.loadFromFile()

    def __iter__(self) :
        """Iterator implementation."""
        return self

    def __next__(self) :
        """Iterator implementation."""
        try :
            key = list(self.graphDict.keys())[self.graphId]
            self.graphId += 1
            return self.graphDict[key]
        except Exception :
            raise StopIteration

    def __len__(self) :
        """Returns the number of nodes in the graph."""
        return len(self.graphDict)

    def loadFromFile(self) :
        """Loads all graphs stored in a GFF file."""
        lineNo = 0
        graph  = None
        # aliases keep track of nodes with different ids but identical start/end positions
        alias  = {}
        edges  = set([])
        indicator = ProgressIndicator(100000, verbose=self.verbose)
        for line in self.instream :
            indicator.update()
            lineNo += 1
            if line.startswith('#') : continue
            s       = line.strip()
            parts   = s.split('\t')

            try :
                recType = parts[2]
                start   = int(parts[3])
                end     = int(parts[4])
            except IndexError :
                raise ValueError('Illegal record in splice graph file at line %d:\n\t%s' % (lineNo,s))

            if recType not in VALID_RECTYPES :
                raise ValueError('Illegal record type in splice graph file at line %d:\n\t%s' % (lineNo,s))

            try :
                # Convert 'ID=ABC;Parent=X,Y,Z' into {'ID':'ABC', 'Parent':'X,Y,Z'}
                attrs = dict([tuple(p.split('=')) for p in parts[-1].split(';') if p])
            except Exception as eee :
                raise ValueError("Illegal attribute field '%s' at line %d in GFF file." % (parts[-1],lineNo))

            try :
                id = attrs['ID']
            except KeyError :
                raise ValueError("GFF attribute field '%s' has no ID at line %d" % (parts[-1],lineNo))

            if recType.lower() in VALID_GENES :
                # Add edges to previous graph
                if graph is not None :
                    for e in edges :
                        graph.addEdge(alias[e[0]], alias[e[1]])
                # Start new graph
                graph        = SpliceGraph(name=id, chromosome=parts[0], strand=parts[6])
                graph.minpos = min(start,end)
                graph.maxpos = max(start,end)
                for k in attrs :
                    if k not in KNOWN_ATTRS :
                        graph.attrs[k] = attrs[k]
                self.graphDict[id] = graph
                edges              = set([])
                alias              = {}
            elif graph is None :
                raise ValueError("Graph feature found before graph header at line %d" % lineNo)
            else :
                node      = graph.addNode(id, start, end)
                alias[id] = node.id
                for k in attrs :
                    if k == AS_KEY :
                        node.addFormsFromString(attrs[k])
                    elif k == ISO_KEY and attrs[k] :
                        node.addIsoformString(attrs[k])
                    elif k in [START_CODON_KEY, END_CODON_KEY] and attrs[k] :
                        node.attrs[k] = set([int(x) for x in attrs[k].split(',')])
                    elif k not in KNOWN_ATTRS :
                        node.addAttribute(k, attrs[k])

                if PARENT_ATTR in attrs :
                    parents = attrs[PARENT_ATTR].split(',')
                    for p in parents :
                        edges.add((p,id))

        if graph is not None :
            for e in edges :
                graph.addEdge(alias[e[0]], alias[e[1]])

        indicator.finish()
        # Initialize iterator counter
        self.graphId = 0
