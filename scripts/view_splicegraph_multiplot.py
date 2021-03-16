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
"""
Viewer capable of showing one or more forms of evidence used to
predict splice graphs along with the predicted graph itself:
  1. Splice graphs
  2. Splice junctions
  3. Short read depths
"""
from SpliceGrapher.shared.config             import *
from SpliceGrapher.shared.utils              import *
from SpliceGrapher.shared.ShortRead          import *
from SpliceGrapher.shared.GeneModelConverter import *
from SpliceGrapher.view.ViewerUtils          import *
from SpliceGrapher.view.GeneView             import GeneView
from SpliceGrapher.formats.loader            import *
from SpliceGrapher.formats.sam               import *
from SpliceGrapher                           import SpliceGraph
from SpliceGrapher.formats                   import gtf, wig, bed, xydata

from optparse import OptionParser, OptionGroup

import sys, os

# Arbitrary Y limit for all but read depths:
Y_LIMIT = 100.0

def splitString(s, maxwidth=20) :
    """Inserts newline characters in a string to maintain the given width."""
    parts  = s.split()
    result = parts[0]
    length = len(result)
    for i in range(1,len(parts)) :
        newlength = length + len(parts[i]) + 1
        if newlength <= maxwidth :
            result += ' ' + parts[i]
            length = newlength
        else :
            result += '\n' + parts[i]
            length  = len(parts[i])
    return result

def refinePatchDict(d, maxwidth=20) :
    """Replaces long keys with multi-line keys in a legend patch dictionary."""
    result = {}
    for k in d.keys() :
        newkey = splitString(k, maxwidth=maxwidth)
        result[newkey] = d[k]
    return result

#==========================================================================
# Incompatibilities in some matplotlib versions may yield runtime warnings
import warnings
warnings.filterwarnings('ignore')
#==========================================================================

## DISPLAY_HEIGHT  = {ORIGINAL_GRAPH:0.20, PREDICTED_GRAPH:0.22, JUNCTION_GRAPH:0.10, DEPTH_GRAPH:0.15}
DISPLAY_HEIGHT  = {ORIGINAL_GRAPH:0.18, PREDICTED_GRAPH:0.18, JUNCTION_GRAPH:0.08, DEPTH_GRAPH:0.12, XY_GRAPH:0.12}
MAX_HEIGHT      = sum(DISPLAY_HEIGHT.values())
X_PAD_FRACTION  = 0.01

#==========================================================================
# Main program
#==========================================================================
USAGE = """%prog [options] gene-name

Interface for viewing alternative splicing information for a gene."""

#==========================================================================
# Initialize command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-v', dest='verbose',      default=False,       action='store_true', help='use verbose output [default: %default]')

fileGroup  = OptionGroup(parser, 'File Options')
fileGroup.add_option('-m', dest='model',        default=SG_GENE_MODEL, help='GFF gene model reference [default: %default]')
fileGroup.add_option('-o', dest='output',       default=None,          help='Output file (extension determines format) [default: screen]')
fileGroup.add_option('-d', dest='depth_file',   default=None,          help='SAM alignment file (gapped and spliced reads) [default: %default]')
fileGroup.add_option('-s', dest='splice_graph', default=None,          help='GFF splice graph file [default: %default]')
fileGroup.add_option('-G', dest='orig_graph',   default=None,          help='Baseline graph file [default: none]')
fileGroup.add_option('-X', dest='xydata',       default=None,          help='File of X,Y value pairs for X-Y plot [default: %default]')
parser.add_option_group(fileGroup)

displayGroup = OptionGroup(parser, 'Display Options')
displayGroup.add_option('-c', dest='jctcover', default=False, help='Display read coverage on junctions [default: %default]', action='store_true')
displayGroup.add_option('-x', dest='xLabels',  default=False, help='Add genomic position labels to plots [default: %default]', action='store_true')
displayGroup.add_option('-E', dest='edge',     default=2,     help='Minimum edge thickness [default: %default]', type='int')
displayGroup.add_option('-J', dest='minjct',   default=2,     help='Minimum coverage for displaying junctions [default: %default]', type='int')
displayGroup.add_option('-L', dest='legend',   default=False, help='Add legend to graph [default: %default]', action='store_true')
displayGroup.add_option('-H', dest='height',   default=8.5,   help='Display window height (inches) [default: %default]', type='float')
displayGroup.add_option('-W', dest='width',    default=11,    help='Display window width (inches) [default: %default]', type='float')
displayGroup.add_option('-F', dest='fontsize', default=12,    help='Font size for plot titles [default: %default]', type='int')
displayGroup.add_option('-S', dest='shrink',   default=False, help='Shrink introns relative to exons [default: %default]', action='store_true')
displayGroup.add_option('-D', dest='display',  default='OPJR',
        help="Display order string for plots (O=original model, P=predicted graph, J=splice junctions, R=read depth, X=XY graph).  Examples: 'OPR', 'RJP' [default: '%default']")
parser.add_option_group(displayGroup)

# Deprecated to simplify interface: 
## fileGroup.add_option('-b', dest='bed_file',     default=None,          help='BED predicted junctions file [default: %default]')
## fileGroup.add_option('-j', dest='junctions',    default=None,          help='SAM spliced alignment file (spliced reads only) [default: %default]')
## fileGroup.add_option('-t', dest='gtf_file',     default=None,          help='GTF file [default: %default]')
## fileGroup.add_option('-w', dest='wig_file',     default=None,          help='WIG depth file [default: %default]')
## displayGroup.add_option('-l', dest='labels',   default=False, help='Include exon labels [default: %default]', action='store_true')
## displayGroup.add_option('-A', dest='adjust',   default=False, help='Adjust splice graph to match gene boundaries [default: %default]', action='store_true')
## displayGroup.add_option('-U', dest='urmargin', default=0,     help='Margin for subsuming unresolved node into known exons [default: %default]', type='int')
## displayGroup.add_option('-C', dest='clusters', default=False, help='Show read clusters instead of read depths [default: %default]', action='store_true')
## displayGroup.add_option('-T', dest='threshold',default=1,     help='Minimum threshold for clusters (-C option) [default: %default]', type='int')
## displayGroup.add_option('--titles', dest='titles',   default=None,  help='Alternate title for each display [default: %default]')

#==========================================================================
# Process command-line options:
opts, args   = parser.parse_args(sys.argv[1:])
opts.display = opts.display.upper()

# Removed after version 0.1.0 to simplify interface:
opts.wig_file  = False
opts.bed_file  = False
opts.gtf_file  = None
opts.urmargin  = sys.maxint
opts.junctions = None
opts.adjust    = False
opts.clusters  = False
opts.threshold = 1
opts.titles    = None
opts.labels    = False
#

if not opts.output :
    try :
        from pylab import *
    except Exception :
        sys.stderr.write('\n** Error initiating matplotlib/pylab modules\n')
        sys.stderr.write('Check your $DISPLAY or try using the -o option.\n\n')
        sys.exit(1)

if len(args) != 1 :
    parser.print_help()
    if len(args) > 1 :
        sys.stderr.write('Too many gene names on command line (expected 1):\n    %s\n' % '\n    '.join(args))
    sys.exit(1)

if not opts.model :
    parser.print_help()
    sys.stderr.write('** No GFF gene model specified.  Use the -m option or set SG_GENE_MODEL in your SpliceGrapher configuration.\n')
    sys.exit(1)

if opts.clusters and not opts.depth_file :
    parser.print_help()
    if opts.wig_file :
        sys.stderr.write('Clusters not implemented for .wig files yet.\n')
    else :
        sys.stderr.write('You must specify a SAM read depth file if you wish to plot clusters.\n')
    sys.exit(1)

titleStrings = opts.titles.split(',') if opts.titles else []

if opts.titles and len(titleStrings) < len(opts.display) :
    parser.print_help()
    sys.stderr.write('** You must enter as many titles (%d) as there are plots (%d)\n' % (len(titleStrings), len(opts.display)))
    sys.exit(1)

if (opts.xydata == '') ^ (XY_GRAPH in opts.display) :
    parser.print_help()
    sys.stderr.write('You specified an X-Y data file without including %s in the display.\n' % XY_GRAPH)
    sys.exit(1)

validateFile(opts.model)
if opts.depth_file   : validateFile(opts.depth_file)
if opts.orig_graph   : validateFile(opts.orig_graph)
if opts.junctions    : validateFile(opts.junctions)
if opts.splice_graph : validateFile(opts.splice_graph)
if opts.gtf_file     : validateFile(opts.gtf_file)
if opts.bed_file     : validateFile(opts.bed_file)
if opts.wig_file     : validateFile(opts.wig_file)

# Required parameter:
geneName = args[0]

writeStartupMessage()

#==========================================================================
# Determine what plots to include:
DISPLAY = {}
DISPLAY[ORIGINAL_GRAPH]  = True
DISPLAY[PREDICTED_GRAPH] = (opts.splice_graph or opts.gtf_file)
DISPLAY[JUNCTION_GRAPH]  = (opts.depth_file or opts.junctions or opts.bed_file)
DISPLAY[DEPTH_GRAPH]     = (opts.depth_file or opts.wig_file)
DISPLAY[XY_GRAPH]        = (opts.xydata is not None)
displayList              = [d for d in opts.display if DISPLAY[d]]

if not displayList :
    parser.print_help()
    sys.stderr.write('\nNo data specified via -s, -r, -w or -t option.  Nothing to do.\n')
    sys.exit(1)
elif opts.verbose :
    sys.stderr.write('Building the following displays:\n')
    sys.stderr.write('  %s\n' % ', '.join([DISPLAY_NAME[d] for d in displayList if DISPLAY[d]]))

if not (opts.gtf_file or opts.splice_graph) and DISPLAY[PREDICTED_GRAPH] :
    raise ValueError('No graph provided for predicted graph plot.')
#==========================================================================

#==========================================================================
# Load gene model information -- absolutely required!
geneModel = loadGeneModels(opts.model, verbose=opts.verbose)
gene      = geneModel.getGeneByName(geneName)
if not gene :
    raise Exception('Unable to locate gene model for gene %s\n' % geneName)

# Establish graph X boundaries (plus small padding) based on gene
minpos = gene.minpos
maxpos = gene.maxpos
# Find splice graph boundaries
predictedGraph = None
if opts.gtf_file :
    predictedGraph = gtf.getFirstGraph(opts.gtf_file)
elif opts.splice_graph :
    predictedGraph = SpliceGraph.getFirstGraph(opts.splice_graph)

if opts.orig_graph :
    origGraph = SpliceGraph.getFirstGraph(opts.orig_graph)
else :
    origGraph = geneModelToSpliceGraph(gene)
    origGraph.annotate()

if predictedGraph :
    if opts.adjust :
        adjustment = gene.minpos - predictedGraph.minpos
        predictedGraph.adjust(adjustment)
    else :
        minpos = min(minpos, predictedGraph.minpos)
        maxpos = max(maxpos, predictedGraph.maxpos)
padding = X_PAD_FRACTION*(maxpos-minpos)

depths    = []
junctions = []
clusters  = []
if opts.depth_file :
    if opts.junctions :
        depthDict = getSamDepths(opts.depth_file, maxpos=maxpos, verbose=opts.verbose)
    else :
        depthDict,jctDict = getSamReadData(opts.depth_file, maxpos=maxpos, minjct=opts.minjct, verbose=opts.verbose)
    try :
        depths = depthDict[gene.chromosome]
    except KeyError :
        sys.stderr.write('** Warning: no depths found for %s in %s\n' % (gene.chromosome, opts.depth_file))

if opts.clusters :
    clusters = depthsToClusters(gene.chromosome, depths, threshold=opts.threshold, verbose=opts.verbose)

if opts.junctions :
    jctDict = getSamJunctions(opts.junctions, maxpos=maxpos, minjct=opts.minjct, verbose=opts.verbose)

if opts.junctions or opts.depth_file :
    try :
        junctions = jctDict[gene.chromosome]
    except KeyError :
        jctFile = opts.junctions if opts.junctions else opts.depth_file
        sys.stderr.write('** Warning: no junctions found for %s in %s' % (gene.chromosome, jctFile))

if opts.xydata :
    (Xvalues,Yvalues) = xydata.getXYData(opts.xydata, minpos=minpos, maxpos=maxpos, verbose=opts.verbose)

neighborGenes = geneModel.getGenesInRange(gene.chromosome, minpos, maxpos, strand=gene.strand)

#==========================================================================
# Draw the graphs:
# Initialize matplotlib view area:
rcParams['figure.figsize'] = opts.width, opts.height
rcParams['font.size']      = opts.fontsize
rcParams['font.weight']    = 'bold'

# Establish reference height for subplots
totalHeight = sum([DISPLAY_HEIGHT[d] for d in displayList])
topLine     = 0.94 if opts.legend else 0.92
axLeft      = 0.05
axWidth     = 0.90 

# Establish reference height for subplots
patchDict   = {}
for i in range(len(displayList)) :
    displayType = displayList[i]
    titleString = titleStrings[i] if titleStrings else None
    height      = MAX_HEIGHT * DISPLAY_HEIGHT[displayType]/totalHeight
    curAxes     = axes([axLeft, topLine-height, axWidth, height])
    topLine     = topLine - height - 0.06

    # Place gene view in background of all graphs
    gv = GeneView(neighborGenes, curAxes)
    gv.plot()

    #------------------------------------------------------------
    # Plot appropriate view
    tmpPatches = {}
    if displayType in [ORIGINAL_GRAPH, PREDICTED_GRAPH] :
        if displayType == ORIGINAL_GRAPH :
            graph = origGraph
            title = titleString if titleString else 'Gene Models for %s (%s)' % (geneName, gene.strand)
        else :
            graph = predictedGraph
            title = titleString if titleString else 'Predicted Splice Graph for %s (%s)' % (geneName, gene.strand)

        if opts.verbose : sys.stderr.write('Rendering splice graph\n')
        tmpPatches, extraPatches = plotSpliceGraph(graph, curAxes,
                                     labels=opts.labels,
                                     xLabels=opts.xLabels,
                                     unresolved=False,
                                     urmargin=opts.urmargin,
                                     geneName=geneName,
                                     title=title,
                                     minwidth=opts.edge,
                                     gtf=(opts.gtf_file is not None))

    elif displayType == JUNCTION_GRAPH :
        if opts.verbose : sys.stderr.write('Rendering splice junctions\n')
        tmpPatches = plotSpliceJunctions(junctions, curAxes, minpos, maxpos,
                                         depths=opts.jctcover,
                                         mindepth=opts.minjct,
                                         title=titleString,
                                         xLabels=opts.xLabels)

    elif displayType == DEPTH_GRAPH :
        if opts.clusters :
            if opts.verbose : sys.stderr.write('Rendering read clusters\n')
            plotClusters(clusters, curAxes, minpos, maxpos,
                         labels=opts.labels,
                         title=titleString,
                         xLabels=opts.xLabels)
        else :
            if opts.verbose : sys.stderr.write('Rendering read depths\n')
            plotReadDepths(depths, curAxes, minpos, maxpos,
                           labels=opts.labels,
                           title=titleString,
                           xLabels=opts.xLabels)

    elif displayType == XY_GRAPH :
        if opts.verbose : sys.stderr.write('Rendering X-Y graph\n')
        title = titleString if titleString else os.path.basename(opts.xydata)
        plotXYGraph(Xvalues, Yvalues, curAxes, minpos, maxpos, title=titleString)

    #------------------------------------------------------------
    # Tasks common to all graphs:
    if tmpPatches :
        patchDict.update(tmpPatches)
    
    # curAxes.grid(True)
    if gene.strand == '+' :
        curAxes.set_xlim(minpos-padding, maxpos+padding)
    else :
        curAxes.set_xlim(maxpos+padding, minpos-padding)

    # Display x positions only on last (bottom) plot
    if displayType == displayList[-1] :
        xvalues = setXticks(int(minpos-padding), int(round(maxpos+padding)))
        curAxes.set_xticks(xvalues)
        curAxes.set_xticklabels(['%d'%x for x in xvalues])
    else :
        curAxes.set_xticklabels([])

    # Display Y positions only on depth graph or XY graphs
    if displayType not in [DEPTH_GRAPH, XY_GRAPH] :
        curAxes.set_yticklabels([])
        curAxes.set_yticks([])
#==========================================================================

if opts.legend and patchDict :
    # Values found via trial and error:
    rcParams['legend.borderaxespad'] = 0.01
    rcParams['legend.handlelength']  = 0.02
    rcParams['legend.handletextpad'] = 0.01
    rcParams['legend.labelspacing']  = 0.008
    rcParams['legend.fontsize']      = max(8,0.75*opts.fontsize)
    patchDict = refinePatchDict(patchDict, maxwidth=24)
    numCols   = max(3,len(patchDict))
    keys      = sorted(patchDict.keys())
    patches   = [patchDict[k] for k in keys]
    placement = 'lower center'
    figlegend(patches, keys, placement, ncol=numCols)

if opts.output :
    if opts.verbose : sys.stderr.write('Writing graph output to %s\n' % opts.output)
    ext = opts.output.split('.')[-1]
    dpi = 100 if ext.lower() == 'png' else 400
    savefig(opts.output, dpi=dpi, format=ext.lower())
else :
    show()
