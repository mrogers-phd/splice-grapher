"""
Module that adds a graph's annotated isoforms to a matplotlib figure
using a separate graph for each isoform.
"""
from SpliceGrapher.shared.utils                 import *
from SpliceGrapher.shared.subgraph              import *
from SpliceGrapher.view.SpliceGraphView         import *
from SpliceGrapher.SpliceGraph                  import *
from SpliceGrapher.predict.SpliceGraphPredictor import DISPOSITION_KEY, UNRESOLVED_NODE, ACCEPTORS_KEY, DONORS_KEY, INFORMATION_KEY

from pylab import *
import sys

HIGHLIGHT_OUTLINE    = '#000000'
HIGHLIGHT_COLOR      = '#FFD700'
ISOLABEL_COLOR       = 'black'
DEFAULT_ISOLABEL_POS = 'top'
MINIMUM_WEIGHT       = 0.001

START_CODON_COLOR    = '#55FF55'
END_CODON_COLOR      = '#FF5555'

class IsoformView(SpliceGraphView) :

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
        self.nodes             = graph.resolvedNodes()
        self.patchDict         = {}
        self.extraPatches      = {}
        self.level             = {}
        self.edgeWidth         = getAttribute('minwidth', 1, **args)
        self.isoformMap        = self.graph.isoformDict()
        self.maxHeight         = len(self.isoformMap)

        if self.graph.maxpos < self.xMin or self.graph.minpos > self.xMax :
            raise ValueError('Splice graph %d-%d does not fall within specified x-range %d-%d' % \
                    (self.graph.minpos, self.graph.maxpos, self.xMin, self.xMax))

        if not self.isoformMap :
            sys.stderr.write('IsoformView warning: graph %s has no annotated isoforms\n' % self.graph.getName())
            sys.exit(1)

        # A caller can turn highlighting on either by supplying a
        # set of highlight features or by setting highlight to True
        self.features          = getAttribute('features', set(), **args)
        self.highlight         = getAttribute('highlight', bool(self.features), **args)

        if self.highlight and not self.features :
            self.features  = set(self.isoformMap.keys())
            self.features |= set([n for n in self.nodes if len(n.isoformSet) == 1])

    def nodeStyles(self, node, level) :
        """Determines a node's appropriate edge and body color."""
        # Special handling for unresolvable nodes:
        if node.isUnresolved() :
            return dict(UNRESOLVED_STYLE)

        result = dict(NORMAL_STYLE)
        if node in self.features :
            result['fc'] = HIGHLIGHT_COLOR
            result['ec'] = HIGHLIGHT_OUTLINE
        elif node.altFormSet :
            # Display only recognized AS types
            bodyType = node.altFormSet.intersection(BODY_SET)
            edgeType = node.altFormSet.intersection(EDGE_SET)
            if bodyType : result['fc'] = FILL_COLORS[bodyType.pop()]
            result['ec'] = EDGE_COLORS[edgeType.pop()] if edgeType else result['fc']
        return result

    def plot(self, **args) :
        """Main method for plotting a graph."""
        self.xLabels  = getAttribute('xLabels', False, **args)
        self.textSize = getAttribute('textSize', 'x-small', **args)

        self.plotNodes(**args)
        self.plotGeneLabels(**args)
        self.axis.set_ylim(0, self.yLimit)

    def plotNodes(self, **args) :
        """Plots all the nodes in the list."""
        isoLabelPos = getAttribute('isoformLabelPosition', DEFAULT_ISOLABEL_POS, **args)
        isoLabels   = getAttribute('isoformLabels', False, **args)
        isoWeights  = getAttribute('isoformWeights', None, **args)
        labelColor  = getAttribute('labelColor', LABEL_COLOR, **args)
        labels      = getAttribute('labels', False, **args)
        highlightAS = getAttribute('highlightAS', True, **args)
        senseStrand = getAttribute('senseStrand', True, **args)
        showCodons  = getAttribute('showCodons', False, **args)
        sortByName  = getAttribute('sortByName', False, **args)
        urmargin    = getAttribute('urmargin', 0, **args)

        # Identify all edges that skip an exon (if there are any)
        allEdges      = edgeSet(self.graph)
        skippedExons  = [n for n in self.graph.nodeDict.values() if ES_ABBREV in n.altFormSet]
        skippingEdges = set([])
        for skipped in skippedExons :
            skippingEdges.update([e for e in allEdges if containsNode(e,skipped)])

        keys       = sorted(self.level.keys())
        isoLengths = dict([(n,acceptor(self.isoformMap[n][0])) for n in self.isoformMap])
        isoKeys    = isoLengths.keys()
        if isoWeights :
            # Account for any isoforms missing from the dictionary
            for k in isoKeys :
                isoWeights.setdefault(k,0.0)
            isoKeys = [x for x in isoKeys if isoWeights[x] >= MINIMUM_WEIGHT]
            isoKeys.sort(cmp=lambda x,y : int(10000*isoWeights[x])-int(10000*isoWeights[y]))
            self.maxHeight = len(isoKeys)
        elif sortByName :
            isoKeys.sort()
        else :
            isoKeys.sort(cmp=lambda x,y : isoLengths[x]-isoLengths[y])

        graphWidth  = self.xMax - self.xMin + 1
        arrWidth    = self.arrowWidth()
        maxHeadLen  = 0.01 * graphWidth
        offset      = -0.5 * (self.maxHeight % 2)
        urCounter   = 0

        level = 0
        for formName in isoKeys :
            level += 1
            nodes = self.isoformMap[formName]

            for i in range(len(nodes)) :
                node       = nodes[i]
                Ypos       = offset + level*self.trackHeight()
                #assert(Ypos >= 0)
                Ymiddle    = Ypos + arrWidth/2
                arrowLen   = len(node)-1 if node.strand=='+' else -len(node)+1
                headLength = min(abs(arrowLen/2), maxHeadLen)
                styles     = self.nodeStyles(node, level) if highlightAS else NORMAL_STYLE

                # Show isoform label
                if i == 0 and isoLabels :
                    textColor = styles['ec'] if node.isUnresolved() else ISOLABEL_COLOR
                    nodeLabel = 'U%d' % urCounter if node.isUnresolved() else node.id
                    if isoLabelPos == DEFAULT_ISOLABEL_POS :
                        textStart = node.start
                        justify   = 'left' if senseStrand else 'right'
                        textY     = Ypos+0.8*arrWidth
                        valign    = 'top'
                    else :
                        textStart = node.start-headLength if node.strand=='+' else node.start+headLength
                        justify   = 'right' if senseStrand else 'left'
                        textY     = Ypos
                        valign    = 'center'

                    labelText = '%s (%.1f%%)'%(formName,100.0*isoWeights[formName]) if isoWeights else formName
                    self.axis.text(textStart, textY, labelText,
                                   weight=styles['weight'], style=styles['style'],
                                   size=self.textSize, color=textColor, ha=justify, va=valign)
                # draw edge:
                if i > 0 :
                    edge      = Edge(nodes[i-1],nodes[i])
                    edgeColor = NORMAL_COLOR
                    styles['lw'] = min(styles['lw'], self.edgeWidth)
                    if edge in self.features :
                        edgeColor = HIGHLIGHT_COLOR
                        #styles['lw'] += 1 # increase edge width
                    elif edge in skippingEdges :
                        edgeColor = ES_COLOR
                    self.axis.plot([donor(nodes[i-1]), acceptor(nodes[i])], [Ypos,Ypos], color=edgeColor, lw=styles['lw'])

                if abs(arrowLen) > 0 :
                    patch = self.axis.arrow(node.start, Ypos, arrowLen, 0.0,
                                            fc=styles['fc'], ec=styles['ec'], ls=styles['ls'], lw=styles['lw'],
                                            width=arrWidth,
                                            head_width=arrWidth,
                                            head_length=headLength,
                                            shape='full',
                                            length_includes_head=True)
                elif len(node) > 0 :
                    #loY   = self.nodeY[node]
                    loY   = Ypos
                    hiY   = loY + arrWidth
                    patch = patches.Rectangle((node.start, loY), len(node), hiY, fill=True,
                                              fc=styles['fc'], ec=styles['ec'], ls=styles['ls'], lw=styles['lw'])
                if showCodons :
                    loY = Ypos - arrWidth/2
                    hiY = Ypos + arrWidth/2
                    for c in node.startCodons() :
                        self.axis.plot([c,c],[loY,hiY], color=START_CODON_COLOR, lw=styles['lw'])
                    for c in node.endCodons() :
                        self.axis.plot([c,c],[loY,hiY], color=END_CODON_COLOR, lw=styles['lw'])

                asAbbrevs = list(LABEL_SET.intersection(node.altForms()))
                asTypes   = [EVENT_NAME[a] for a in asAbbrevs]
                key       = None
                if node in self.features :
                    key = 'Unique to Isoform'
                    self.patchDict[key] = patch
                elif asTypes and highlightAS :
                    asTypes.sort()
                    key = '/'.join(asTypes)
                    self.patchDict[key] = patch

                # Show exon labels (always on for unresolved nodes)
                if labels :
                    midpt     = 0.5*(node.start+node.end)
                    textColor = styles['ec'] if node.isUnresolved() else labelColor
                    nodeLabel = 'U%d' % urCounter if node.isUnresolved() else node.id
                    self.axis.text(midpt, Ypos+arrWidth/2, nodeLabel,
                                   weight=styles['weight'], style=styles['style'],
                                   size=self.textSize, color=textColor, ha='center', va='center')

                # Show X positions
                if self.xLabels :
                    bbox = dict(facecolor='lightgrey', edgecolor='lightgrey')
                    self.axis.text(node.start, Ypos+arrWidth/2, '%d'%node.origStart, size='xx-small',
                            color='#000000', ha='right', va='center', bbox=bbox)
                    self.axis.text(node.end, Ypos-arrWidth/2, '%d'%node.origEnd, size='xx-small',
                            color='#000000', ha='left', va='center', bbox=bbox)

    def trackHeight(self) :
        availHeight = self.yLimit-self.border
        return availHeight/(self.maxHeight+1)
