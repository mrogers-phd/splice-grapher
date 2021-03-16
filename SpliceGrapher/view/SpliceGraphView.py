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
Module that adds a splice graph depiction to a matplotlib figure.
"""
from SpliceGrapher.shared.utils                 import *
from SpliceGrapher.SpliceGraph                  import *
from SpliceGrapher.shared.SpliceGraphPath       import *
from SpliceGrapher.predict.SpliceGraphPredictor import DISPOSITION_KEY, UNRESOLVED_NODE, ACCEPTORS_KEY, DONORS_KEY, INFORMATION_KEY

# paths and patches for bezier splines
from pylab           import *
from matplotlib.path import Path
import matplotlib.patches as patches

import sys, getopt, os

# AS types and relationship to node edge/body highlighting
EDGE_SET  = set([ES_ABBREV, IR_ABBREV])
BODY_SET  = set([ALT5_ABBREV, ALT3_ABBREV])
LABEL_SET = BODY_SET.union(EDGE_SET)

# Colors for each AS type
LABEL_COLOR  = 'red'
NORMAL_COLOR = 'grey'
ES_COLOR     = '#00EE00'
ALT3_COLOR   = 'orange'
ALT5_COLOR   = '#DD77FF'
ALTB3_COLOR  = '#00DD00'
ALTB5_COLOR  = '#00DD00'

# Colors for codon annotations
START_CODON_COLOR = '#55FF55'
END_CODON_COLOR   = '#FF5555'

FILL_COLORS  = {ES_ABBREV:'grey', ALT5_ABBREV:ALT5_COLOR, ALT3_ABBREV:ALT3_COLOR, ALTB3_ABBREV:ALTB3_COLOR, ALTB5_ABBREV:ALTB5_COLOR, IR_ABBREV:'grey'}
EDGE_COLORS  = {ES_ABBREV:ES_COLOR,  ALT5_ABBREV:ALT5_COLOR, ALT3_ABBREV:ALT3_COLOR, ALTB3_ABBREV:ALTB3_COLOR, ALTB3_ABBREV:ALTB5_COLOR, IR_ABBREV:'blue'}

GENE_LABEL_FILL = '#D0D0D0'

# Drawing properties for normal/unresolved nodes
NORMAL_STYLE     = {'fc':NORMAL_COLOR, 'ec':NORMAL_COLOR, 'ls':'solid',  'lw':2, 'weight':'bold',  'style':'normal'}
UNRESOLVED_STYLE = {'fc':'#EEEEEE',     'ec':'#777777',   'ls':'dotted', 'lw':1, 'weight':'light', 'style':'italic'}

# Scale for Y-axis (arbitrary)
Y_LIMIT     = 100.0

# Proportion of limit to use as empty border
BORDER      = 0.08

# Proportion of x-range to enforce separation of root/leaf nodes
END_MARGIN_PCT = 0.05

# Bezier spline path description
CURVE_CODES = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,]

# Nodes 'overlap' if they are within this many positions/nt of each other:
NODE_MARGIN = 4

# Maximum number of nodes that still allow us to make a 'pretty' graph:
MAX_PRETTY_NODES = 90

class SpliceGraphView(object) :

    def __init__(self, graph, axis, **args) :
        self.graph             = graph
        self.axis              = axis
        self.yLimit            = getAttribute('yLimit', Y_LIMIT, **args)
        self.xLimits           = getAttribute('xLimits', (self.graph.minpos,self.graph.maxpos), **args)
        self.includeUnresolved = getAttribute('unresolved', False, **args)
        self.verbose           = getAttribute('verbose', False, **args)
        self.xMin              = min(self.xLimits)
        self.xMax              = max(self.xLimits)
        self.border            = BORDER*self.yLimit
        self.nodes             = self.graph.resolvedNodes()
        self.unresolvedNodes   = self.graph.unresolvedNodes() if self.includeUnresolved else []
        self.maxHeight         = 0
        self.patchDict         = {}
        self.extraPatches      = {}
        self.level             = {}
        self.edgeWidth         = getAttribute('minwidth', 1, **args)

        self.setLevels()

    def addCurve(self, vertices, color=NORMAL_COLOR, lw=1) :
        """
        Add a Bezier curve using the given 4-vertex specification.
        """
        path   = Path(vertices, CURVE_CODES)
        patch  = patches.PathPatch(path, fill=False, lw=lw, ec=color)
        self.axis.add_patch(patch)

    def arrowWidth(self) :
        """
        Finds a 'nice' arrow width for most purposes with a maximum
        size of 8% of the graph height.
        """
        return min(0.8*(self.yLimit-self.border), 0.6*self.trackHeight())

    def getExonSkippingPartners(self, node) :
        """
        Find all children whose edges skip an exon in the graph.
        """
        # Find child edges that contain skipped nodes:
        result = {}
        for c in node.children :
            # Get start/end positions of edge to child
            (edgeStart,edgeEnd) = (node.maxpos, c.minpos) if node.strand == '+' else (c.maxpos, node.minpos)

            # Establish min/max positions covered by child edge
            for s in self.skippedExons :
                if edgeStart <= s.minpos and s.maxpos <= edgeEnd :
                    (edgeMin,edgeMax) = result.setdefault(c, (s.minpos, s.maxpos))
                    result[c]         = (min(edgeMin,s.minpos), max(edgeMax,s.maxpos))
        return result

    def isSkippedExon(self, node) :
        """Returns true if the node represents a skipped exon; false otherwise."""
        return AS_KEY in node.attrs and ES_ABBREV in node.attrs[AS_KEY]

    def nodeStyles(self, node) :
        """Determines a node's appropriate edge and body color."""
        # Special handling for unresolvable nodes:
        if node.isUnresolved() :
            return dict(UNRESOLVED_STYLE)

        result = dict(NORMAL_STYLE)
        if node.altFormSet :
            # Display only recognized AS types
            bodyType = node.altFormSet.intersection(BODY_SET)
            edgeType = node.altFormSet.intersection(EDGE_SET)
            if bodyType : result['fc'] = FILL_COLORS[bodyType.pop()]
            result['ec'] = EDGE_COLORS[edgeType.pop()] if edgeType else result['fc']
        return result

    def nodesOverlap(self, a,b, margin=0) :
        """
        Returns true if two nodes overlap; false otherwise.  Using a margin permits
        a graph to separate vertically nodes that may appear on top of one another.
        """
        if (a.isRoot() and b.isRoot()) or (a.isLeaf() and b.isLeaf()) :
            margin = (self.xMax-self.xMin+1)*END_MARGIN_PCT
        else :
            margin = NODE_MARGIN
        a_min = a.minpos-margin
        a_max = a.maxpos+margin
        b_min = b.minpos-margin
        b_max = b.maxpos+margin
        return a_min <= b_min <= a_max or a_min <= b_max <= a_max or b_min <= a_min <= b_max or b_min <= a_max <= b_max

    def plot(self, **args) :
        """Main method for plotting a graph."""
        self.xLabels  = getAttribute('xLabels', False, **args)
        self.textSize = getAttribute('textSize', 'x-small', **args)

        self.plotNodes(**args)
        self.plotEdges(**args)
        self.plotGeneLabels(**args)
        self.axis.set_ylim(0, self.yLimit)

    def plotEdges(self, **args) :
        """Plots Bezier edges between exons in the graph."""
        highlightAS = getAttribute('highlightAS', True, **args)
        skipColor   = ES_COLOR if highlightAS else NORMAL_COLOR

        # 2-point spline between two endpoints:
        for n in self.nodes :
            if n.isUnresolved() and not self.includeUnresolved : continue

            # Invoke special handling for skipped exons as needed
            skipPartners = set([])
            if self.skippedExons :
                skipPartners = self.getExonSkippingPartners(n)

            startY = self.nodeY[n]
            if startY < 0 :
                if self.subsumed[n] : continue
                sys.stderr.write('Bad Y position set for node %s: %d\n' % (n.id,startY))

            for c in n.children :
                if c.isUnresolved() and not self.includeUnresolved : continue
                endY = self.nodeY[c]
                if endY < 0 :
                    if self.subsumed[c] : continue
                    sys.stderr.write('Bad Y position set for node %s: %d\n' % (c.id,endY))

                skipEdge  = (c in skipPartners)
                avoidExon = False
                if c in skipPartners :
                    region  = skipPartners[c]
                    skipped = [e for e in self.skippedExons if self.nodeY[e] == startY and region[0] <= e.minpos and e.maxpos <= region[1]]
                    if skipped : avoidExon = True

                if skipEdge and startY == endY and avoidExon :
                    midY   = startY + max(10,0.5*self.trackHeight())
                    #
                    # Use intervening exon to determine spline control points, so that
                    # lines (hopefully) don't get hidden behind intervening exons.
                    #   1st control point: between upstream exon and intervening exon
                    #   2nd control point: between intervening exon and downstream exon
                    # Best-looking results are when we use consistent splines on both ends,
                    # so we use minimum distance between a flanking exon and the intervening exon
                    #
                    if self.graph.strand == '+' :
                        mindelta = min(abs(region[0]-n.end), abs(c.start-region[1]))
                    else :
                        mindelta = -min(abs(region[1]-n.end), abs(c.start-region[0]))
                    ctrl1       = n.end + 0.5*mindelta
                    ctrl2       = c.start - 0.5*mindelta
                    offset1     = n.end + mindelta
                    offset2     = c.start - mindelta
                    firstCurve  = [(n.end, startY), (ctrl1, startY), (ctrl1, midY), (offset1, midY),]
                    secondCurve = [(offset2, midY), (ctrl2, midY), (ctrl2, endY), (c.start, endY),]

                    self.addCurve(firstCurve, color=skipColor, lw=self.edgeWidth)
                    # connect upstream "S" to downstream "S"
                    self.axis.plot([offset1,offset2], [midY,midY], lw=self.edgeWidth, linestyle='-', color=skipColor)
                    self.addCurve(secondCurve, color=skipColor, lw=self.edgeWidth)
                else :
                    midpt = 0.5 * (c.start + n.end)
                    verts = [(n.end, startY), (midpt, startY), (midpt, endY), (c.start, endY),]
                    edgeColor = skipColor if skipEdge else NORMAL_COLOR
                    self.addCurve(verts, color=edgeColor, lw=self.edgeWidth)

    def plotGeneLabels(self, **args) :
        """Adds gene labels to the bottom of the graph, if available."""
        genes  = getAttribute('genes', [], **args)
        text_y = (0.5*BORDER)*self.yLimit
        for gene in genes :
            # Ensure label does not extend beyond plot boundaries:
            rectMinpos = max(gene.minpos, self.xMin)
            rectMaxpos = min(gene.maxpos, self.xMax)
            if rectMaxpos < gene.minpos or rectMinpos > gene.maxpos : continue

            rectLength = rectMaxpos-rectMinpos+1
            midpt      = (gene.minpos+gene.maxpos)/2
            if rectMinpos < midpt < rectMaxpos and rectLength > 1 :
                p     = patches.Rectangle((rectMinpos, 0.0), rectLength, 2*text_y, fill=True, fc=GENE_LABEL_FILL, ec=GENE_LABEL_FILL, lw=1.0)
                self.axis.add_patch(p)
                self.axis.text(midpt, text_y, gene.id, weight='bold', size='x-small', color='black', ha='center', va='center')

    def plotNodes(self, **args) :
        """Plots all the nodes in the list."""
        labels      = getAttribute('labels', False, **args)
        labelColor  = getAttribute('labelColor', LABEL_COLOR, **args)
        highlightAS = getAttribute('highlightAS', True, **args)
        showCodons  = getAttribute('showCodons', False, **args)
        urmargin    = getAttribute('urmargin', 0, **args)
        graphWidth  = self.xMax - self.xMin + 1
        arrWidth    = self.arrowWidth()
        maxHeadLen  = 0.01*graphWidth
        middleY     = 0.5 * self.yLimit
        offset      = -0.5 * (self.maxHeight % 2)
        urCounter   = 0

        keys       = sorted(self.level.keys())
        self.subsumed = {}.fromkeys(keys,False)
        resolvedNodes = self.graph.resolvedNodes()
        for node in keys :
            # omit unresolved nodes subsumed by other nodes
            if node.isUnresolved() :
                if not self.includeUnresolved : continue
                self.subsumed[node] = [o for o in resolvedNodes if (o.minpos-urmargin) <= node.minpos and (o.maxpos+urmargin) >= node.maxpos]
                if self.subsumed[node] : continue

            self.nodeY[node] = middleY + (self.level[node]+offset)*self.trackHeight()
            arrowLen         = len(node)-1 if node.strand == '+' else -len(node)+1
            styles           = self.nodeStyles(node) if highlightAS else NORMAL_STYLE

            if abs(arrowLen) > 0 :
                # Arrow head should never be longer than exon; at most half its length
                headLength = min(abs(arrowLen/2), maxHeadLen)

                #assert(self.nodeY[node] >= 0)
                patch = self.axis.arrow(node.start, self.nodeY[node], arrowLen, 0.0,
                                    fc=styles['fc'], ec=styles['ec'], ls=styles['ls'], lw=styles['lw'],
                                    width=arrWidth,
                                    head_width=arrWidth,
                                    head_length=headLength,
                                    shape='full',
                                    length_includes_head=True)
            elif len(node) > 0 :
                loY   = self.nodeY[node]
                hiY   = loY + arrWidth
                patch = patches.Rectangle((node.start, loY), len(node), hiY, fill=True,
                        fc=styles['fc'], ec=styles['ec'], ls=styles['ls'], lw=styles['lw'])
            else :
                raise ValueError('Node %s has non-positive length' % n)

            if showCodons :
                loY = self.nodeY[node] - arrWidth/2
                hiY = self.nodeY[node] + arrWidth/2
                for c in node.startCodons() :
                    self.axis.plot([c,c],[loY,hiY], color=START_CODON_COLOR, lw=styles['lw'])
                for c in node.endCodons() :
                    self.axis.plot([c,c],[loY,hiY], color=END_CODON_COLOR, lw=styles['lw'])

            # Create legend keys:
            if node.isUnresolved() :
                urCounter += 1
                self.extraPatches['Unresolved'] = patch
            else :
                asAbbrevs = list(LABEL_SET.intersection(node.altForms()))
                asTypes   = [EVENT_NAME[a] for a in asAbbrevs]
                key       = None
                if asTypes :
                    asTypes.sort()
                    key = '/'.join(asTypes)
                    self.patchDict[key] = patch

            # Show exon labels (always on for unresolved nodes)
            if labels or node.isUnresolved() :
                midpt     = 0.5*(node.start+node.end)
                textColor = styles['ec'] if node.isUnresolved() else labelColor
                nodeLabel = 'U%d' % urCounter if node.isUnresolved() else node.id
                self.axis.text(midpt, self.nodeY[node], nodeLabel,
                               weight=styles['weight'], style=styles['style'],
                               size=self.textSize, color=textColor, ha='center', va='center')

            # Show X positions
            if self.xLabels :
                bbox = dict(facecolor='lightgrey', edgecolor='lightgrey')
                self.axis.text(node.start, self.nodeY[node]+arrWidth/2, '%d'%node.origStart,
                               size='xx-small', color='#000000', ha='right', va='center', bbox=bbox)
                self.axis.text(node.end, self.nodeY[node]-arrWidth/2, '%d'%node.origEnd,
                               size='xx-small', color='#000000', ha='left', va='center', bbox=bbox)

    def setLevels(self) :
        """
        Sets the display Y level for each node in the graph.
        """
        if len(self.nodes) < MAX_PRETTY_NODES :
            self.paths = getAllPaths(self.graph)
            self.paths.sort(reverse=True)
        else :
            # Do the best we can
            longest    = getLongestPath(self.graph)
            self.paths = [longest]
            others     = [n for n in self.nodes if n not in longest.nodes]
            if others :
                others.sort(reverse=(self.graph.strand=='-'))
                fakePath   = SpliceGraphPath(others[0])
                for o in others[1:] :
                    fakePath.append(o)
                self.paths.append(fakePath)

        # Convenience method for computing total parent-child deltas
        def totalDelta() :
            result = 0
            for n in self.level :
                if n.isLeaf() : continue
                result += sum([abs(self.level[n]-self.level[c]) for c in n.children])
            return result

        # prepare for plotting
        self.nodeY        = {}.fromkeys(self.graph.nodeDict.values(), -1)
        self.skippedExons = set([])

        # Initialize nodes in main path, defined as
        # the one with the most exons
        distance = {}
        overlaps = {}
        allNodes = set(self.nodes + self.unresolvedNodes)
        for n in allNodes :
            overlaps[n] = set([])

        mainPath = self.paths[0]
        for n in mainPath.nodes :
            distance[n] = 0
            if self.isSkippedExon(n) : self.skippedExons.add(n)

        endMargin = (self.xMax-self.xMin+1)*END_MARGIN_PCT
        allEdges  = edgeSet(self.graph)
        # First organize nodes bottom-to-top
        for p in self.paths[1:] :
            # all edges on same level from previous paths
            levelEdges  = [e for e in allEdges if e.parent in distance and e.child in distance and distance[e.parent] == distance[e.child]]
            for n in p.nodes :
                if n in distance : continue
                overlaps[n] = set([o for o in distance if self.nodesOverlap(n,o)])
                for o in overlaps[n] :
                    overlaps.setdefault(o,set([]))
                    overlaps[o].add(n)
                if self.isSkippedExon(n) : self.skippedExons.add(n)
                if overlaps[n] :
                    distance[n] = max([(1+distance[o]) for o in overlaps[n]])
                else :
                    edges       = [e for e in levelEdges if containsNode(e,n)]
                    distance[n] = max([(1+distance[e.parent]) for e in edges]) if edges else 0

        self.maxHeight = max(distance.values())
        if self.unresolvedNodes :
            self.maxHeight += 1
            for n in self.unresolvedNodes :
                distance[n] = self.maxHeight

        # Convert bottom-to-top to center-out
        self.level = {}
        for n in distance :
            k = distance[n]+1
            # f(0)=0, f(1)=1, f(2)=-1, f(3)=2, ...
            self.level[n] = ((-1)**k)*(k/2)

        ## # VALIDATION/DEBUGGING:
        ## setNodes = set(self.level.keys())
        ## unset    = allNodes-setNodes
        ## assert(len(unset) == 0)

        # Swap nodes to minimize parent-child distances: a simple, greedy heuristic
        levelEdges = [e for e in allEdges if self.level[e.parent] == self.level[e.child]]
        total      = totalDelta()
        for n in self.level :
            if n.isUnresolved() : continue
            for c in n.children :
                # If the child overlaps nothing on parent's level, simply move it
                delta = abs(self.level[n]-self.level[c])
                if delta == 0 : continue

                # try swapping with other nodes that this child overlaps
                # if a swap improves the total parent-child differences,
                # keep it; otherwise revert back
                overlapping = overlaps[c] if c in overlaps else []
                for o in overlapping :
                    if o in mainPath.nodes or o.isUnresolved() : continue # leave unresolved and main path alone
                    if abs(self.level[n]-self.level[o]) >= delta : continue
                    # swap levels
                    self.level[o],self.level[c] = self.level[c],self.level[o]
                    # update total
                    revised    = totalDelta()
                    # a swap may cause a retained intron node to be on top of one of
                    # its overlapping nodes; for now we just undo the swap
                    violations  = [x.id for x in overlaps[c] if self.level[x] == self.level[c]]
                    violations += [x.id for x in overlaps[o] if self.level[x] == self.level[o]]
                    if violations or revised >= total :
                        self.level[o],self.level[c] = self.level[c],self.level[o]
                    else :
                        levelEdges = [e for e in allEdges if self.level[e.parent] == self.level[e.child]]
                        total      = revised
                        delta      = abs(self.level[n]-self.level[c])

    def trackHeight(self) :
        return (self.yLimit-self.border)/(2+self.maxHeight)
