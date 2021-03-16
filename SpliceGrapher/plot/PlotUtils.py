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
from SpliceGrapher.formats.loader            import *
from SpliceGrapher.shared.GeneModelConverter import *
from SpliceGrapher.formats.xydata            import *
from SpliceGrapher.formats.sam               import *
from SpliceGrapher.SpliceGraph               import *
from SpliceGrapher.shared.utils              import *
from SpliceGrapher.shared.adjust             import *
from SpliceGrapher.shared.ShortRead          import depthsToClusters
from SpliceGrapher.SpliceGraph               import getFirstGraph
from SpliceGrapher.view.ViewerUtils          import *
from SpliceGrapher.view.ReadDepthView        import DEFAULT_LOGY_THRESHOLD
from SpliceGrapher.view.GeneView             import GeneView
from SpliceGrapher.plot.PlotterConfig        import *

from sys import maxint as MAXINT
import matplotlib, re, sys, math

# Formatting directives for inserting values into title strings
DIRECTIVE_MATCH = re.compile('%\w+\S')
GENE_DIRECTIVE   = '%gene'
PATH_DIRECTIVE   = '%path'
STRAND_DIRECTIVE = '%strand'

# determined via trial and error:
AXIS_LEFT      = 0.08
AXIS_WIDTH     = 0.84 
DISPLAY_HEIGHT = 0.99
X_PAD_FRACTION = 0.01
Y_PAD_FRACTION = 0.0029

# rcParam values found through trial and error:
PLOT_AREA_SPECS   = {'font.weight':'bold', 'legend.borderaxespad': 0.01, 'legend.handlelength': 0.02, 'legend.handletextpad': 0.01, 'legend.labelspacing': 0.008}

# Legend values
LEGEND_PADDING    = 0.04
LEGEND_PLACEMENT  = 'lower center'
LEGEND_FONT_RATIO = 0.8
LEGEND_MAX_FONT_SIZE = 8.0

# Order for loading plot source data
DATA_LOAD_ORDER  = [SPLICEGRAPH_PLOT, ISOFORM_PLOT, GENE_MODEL_PLOT, READ_DEPTH_PLOT, JUNCTION_PLOT, XY_PLOT]
SHRINKABLE_PLOTS = [SPLICEGRAPH_PLOT, ISOFORM_PLOT, GENE_MODEL_PLOT, READ_DEPTH_PLOT, JUNCTION_PLOT]

class GenePlotSpecs(object) :
    """Simple data structure that contains information common to a
    set of plots for a single gene."""
    def __init__(self, gene, bgndGenes, data, strand, chromosome, minpos, maxpos, **args) :
        self.gene      = gene
        self.bgndGenes = bgndGenes
        self.data      = data
        self.strand    = strand
        self.chrom     = chromosome.lower()
        self.minpos    = minpos
        self.maxpos    = maxpos
        self.origMin   = getAttribute('origMin', minpos, **args)
        self.origMax   = getAttribute('origMax', maxpos, **args)
        self.fontSize  = getAttribute('fontSize', 12, **args)

    def setFontSize(self, fontSize) :
        self.fontSize = fontSize

def formatTitleString(p, **args) :
    """Formats a title string based on formatting directives."""
    geneName   = getAttribute('geneName', None, **args)
    graph      = getAttribute('graph', None, **args)
    result     = p.title_string
    try :
        envList    = DIRECTIVE_MATCH.findall(result)
    except TypeError :
        raise ValueError('Bad title string %s' % str(result))

    sourcePath = resolveSourcePath(p, **args)
    base       = os.path.basename(sourcePath)
    for env in envList :
        name = env.strip()
        if name == GENE_DIRECTIVE :
            if geneName :
                gname = geneName
            elif graph :
                gname = graph.getName()
            else :
                gname,ext = os.path.splitext(base)
            result = result.replace(env,gname)
        elif name == STRAND_DIRECTIVE :
            if not graph : raise ValueError('No graph provided for strand information')
            result = result.replace(env,graph.strand)
        elif name == PATH_DIRECTIVE :
            result   = result.replace(env,base)
    return result

def generatePlot(p, plotSpecs, plotAxes, **args) :
    """Generates the appropriate plot based on the plot type.
    Returns patch dictionaries for a legend for splice graph plots."""
    verbose     = getAttribute('verbose', False, **args)
    shrink      = getAttribute('shrink', False, **args)
    xTickLabels = getAttribute('xTickLabels', False, **args)

    if p.hide : return {},{}

    data        = getSourceData(p, plotSpecs, **args)
    patches     = {}
    extra       = {}
    padding     = X_PAD_FRACTION * (plotSpecs.maxpos-plotSpecs.minpos)

    # Place gene view in background
    if p.background : GeneView(plotSpecs.gene, plotAxes).plot()

    graph = data if p.plot_type == GENE_MODEL_PLOT else None
    title = formatTitleString(p, graph=graph, **args)

    # Only show gene labels at bottom if labels are turned on:
    geneList = plotSpecs.bgndGenes if p.labels else []

    if p.plot_type in [GENE_MODEL_PLOT, SPLICEGRAPH_PLOT] :
        patches,extra = plotSpliceGraph(data, plotAxes,
                                        geneName=data.name,
                                        genes=geneList,
                                        title=title,
                                        xLimits=[plotSpecs.gene.minpos, plotSpecs.gene.maxpos],
                                        minwidth=p.edge_weight,
                                        labels=p.labels,
                                        showCodons=p.codons,
                                        xLabels=p.x_labels)
    elif p.plot_type == ISOFORM_PLOT :
        # Nodes to highlight in isoform plot:
        weightDict    = loadIsoformWeights(data.name,p.iso_weights) if p.iso_weights else None
        features      = [n for n in data.nodeDict.values() if n.id in p.features] if p.features else []
        patches,extra = plotIsoforms(data, plotAxes,
                                     geneName=data.name,
                                     genes=geneList,
                                     xLimits=[plotSpecs.gene.minpos, plotSpecs.gene.maxpos],
                                     title=title,
                                     minwidth=p.edge_weight,
                                     labels=p.labels,
                                     isoformLabels=p.iso_labels,
                                     isoformWeights=weightDict,
                                     features=features,
                                     highlight=bool(p.highlight),
                                     sortByName=(not weightDict),
                                     showCodons=p.codons,
                                     xLabels=p.x_labels)
    elif p.plot_type == JUNCTION_PLOT :
        message = ''
        if not data : message = 'No spliced reads for %s' % title
        junctions = [jct for jct in data if jct.strand == plotSpecs.gene.strand]
        if not (junctions or message) : message = 'All junctions on wrong strand for %s' % title
        junctions = [jct for jct in junctions if plotSpecs.minpos <= jct.minpos and jct.maxpos <= plotSpecs.maxpos]
        if not (junctions or message) : message = 'All junctions outside gene bounds for %s' % title
        if verbose and message : sys.stderr.write('%s\n' % message)

        jpatch = plotSpliceJunctions(junctions, plotAxes, plotSpecs.minpos, plotSpecs.maxpos,
                            depths=p.labels,
                            mindepth=p.min_coverage,
                            xLabels=p.x_labels,
                            showNovel=p.novel_jct,
                            acceptors=p.acceptors,
                            donors=p.donors,
                            title=title)
        if p.novel_jct : patches = jpatch
    elif p.plot_type == READ_DEPTH_PLOT :
        if p.clusters :
            clusters = depthsToClusters(plotSpecs.gene.chromosome, data, **args)
            plotClusters(clusters, plotAxes, plotSpecs.minpos, plotSpecs.maxpos,
                           labels=p.labels,
                           xLabels=p.x_labels,
                           yLimit=p.y_limit,
                           loglimit=p.log_threshold,
                           title=title)
        else :
            plotReadDepths(data, plotAxes, plotSpecs.minpos, plotSpecs.maxpos,
                           labels=p.labels,
                           xLabels=p.x_labels,
                           yLimit=p.y_limit,
                           loglimit=p.log_threshold,
                           highlight=p.highlight,
                           title=title)
    elif p.plot_type == XY_PLOT :
        from SpliceGrapher.view.XYGraphView import BAR_GRAPH, LINE_GRAPH
        style = LINE_GRAPH if p.line else BAR_GRAPH
        plotXYGraph(data[0], data[1], plotAxes, plotSpecs.minpos, plotSpecs.maxpos, plottype=style, title=title)
    else :
        raise ValueError('Unrecognized plot type for %s: %s\n' % (p.name, p.plot_type))
    
    if plotSpecs.strand == '+' :
        plotAxes.set_xlim(plotSpecs.minpos-padding, plotSpecs.maxpos+padding)
    else :
        plotAxes.set_xlim(plotSpecs.maxpos+padding, plotSpecs.minpos-padding)

    if xTickLabels :
        if shrink and p.plot_type in SHRINKABLE_PLOTS :
            xvalues = [int(plotSpecs.minpos), int(round(plotSpecs.maxpos))]
            xlabels = ['%d'%x for x in [plotSpecs.origMin,plotSpecs.origMax]]
        else :
            xvalues = setXticks(int(plotSpecs.minpos-padding), int(round(plotSpecs.maxpos+padding)))
            xlabels = ['%d'%x for x in xvalues]
        plotAxes.set_xticks(xvalues)
        plotAxes.set_xticklabels(xlabels)
    else :
        plotAxes.set_xticklabels([])

    # Display Y positions only on depth graph or XY graphs
    if p.plot_type not in [READ_DEPTH_PLOT, XY_PLOT] :
        plotAxes.set_yticklabels([])
        plotAxes.set_yticks([])

    return patches,extra

def getAvailableHeight(plots, fontSize, **args) :
    """Determines the total height available for actual plotting areas."""
    legend    = getAttribute('legend', False, **args)
    legendPad = getAttribute('legendPadding', LEGEND_PADDING, **args)
    padding = (1 + len(plots)) * getTitlePadding(fontSize, **args)
    result  = DISPLAY_HEIGHT - padding
    if legend : result -= legendPad
    return result

def getCommonFiles(plots, **args) :
    """Returns a set of files common to all directories given for a set of plots."""
    verbose = getAttribute('verbose', False, **args)
    result  = set([])
    for i in xrange(len(plots)) :
        if verbose : sys.stderr.write('Looking for files in %s\n' % plots[i].source_file)
        if i == 0 :
            result = set(getFileList(plots[i], **args))
        else :
            result &= set(getFileList(plots[i], **args))
    return result

def getFileList(p, **args) :
    """Returns a list of all files relevant to the given plot."""
    from glob import glob
    verbose = getAttribute('verbose', False, **args)
    if p.file_format not in MULTI_FILE_FORMATS :
        raise ValueError('Illegal format for parsing a directory: %s' % p.file_format)
    result = []
    for t in VALID_FILE_FORMAT[p.plot_type] :
        for ext in VALID_FILE_EXTS[t] :
            target = os.path.join(p.source_file, '*%s'%ext)
            if verbose : sys.stderr.write('  %s\n' % target)
            files  = [os.path.basename(f) for f in glob(target)]
            result += files
    return result

def getSourceData(p, geneSpecs, **args) :
    """Returns the appropriate form of data for the given
    plot type from the source data dictionary."""
    sourcePath = resolveSourcePath(p, **args)
    if p.plot_type == GENE_MODEL_PLOT :
        result = geneSpecs.data[sourcePath][SPLICEGRAPH_PLOT]
    elif p.plot_type in [JUNCTION_PLOT, READ_DEPTH_PLOT] :
        try :
            result = geneSpecs.data[sourcePath][p.plot_type][geneSpecs.chrom]
        except KeyError :
            return []
    else :
        result = geneSpecs.data[sourcePath][p.plot_type]
    return result

def getTitlePadding(fontSize, **args) :
    """Returns the amount of space to leave for title text."""
    padFraction = getAttribute('padFraction', Y_PAD_FRACTION, **args)
    return padFraction * fontSize

def initializePlots(geneSpecs, width, height, plots, **args) :
    """Takes the relative plot sizes for all plots and normalizes
    them so they sum to 1.  Sets min and max plot boundaries if
    they were not in the configuration."""
    # Initialize matplotlib view area:
    rcParams.update(PLOT_AREA_SPECS)
    rcParams['figure.figsize']  = width, height
    rcParams['font.size']       = geneSpecs.fontSize
    rcParams['legend.fontsize'] = max(LEGEND_MAX_FONT_SIZE, LEGEND_FONT_RATIO*geneSpecs.fontSize)

    # Normalize plot relative sizes
    visible      = [p for p in plots if not p.hide]
    totalHeights = sum([p.relative_size for p in visible])
    scaleFactor  = DISPLAY_HEIGHT/totalHeights
    for p in visible :
        p.relative_size *= scaleFactor

def loadData(plot, dataDict, **args) :
    """Loads the given data file based on the plot type."""
    gene       = getAttribute('gene', None, **args)
    maxSAM     = getAttribute('maxSAM', MAXINT, **args)
    minjct     = getAttribute('minjct', 2, **args)
    storedDict = getAttribute('stored', None, **args)
    verbose    = getAttribute('verbose', False, **args)

    if plot.plot_type == GENE_MODEL_PLOT and not gene :
        raise ValueError('No gene provided for %s %s' % (plot.name, plot.plot_type))

    # First resolve source file path
    sourcePath = resolveSourcePath(plot, **args)

    # Data may already have been loaded
    dataDict.setdefault(sourcePath,{})
    if plot.file_format != GENE_MODEL_FORMAT :
        if plot.plot_type in dataDict[sourcePath] :
            return sourcePath
        elif storedDict :
            try :
                dataDict[sourcePath][plot.plot_type] = storedDict[sourcePath][plot.plot_type]
                return sourcePath
            except Exception :
                pass

    if plot.file_format not in VALID_FILE_FORMAT[plot.plot_type] :
        raise ValueError('Unrecognized file format %s for %s plot' % (plot.file_format, plot.plot_type))

    # Load the data
    if verbose : sys.stderr.write('Loading data for %s file %s\n' % (plot.plot_type, sourcePath))
    if plot.file_format == GENE_MODEL_FORMAT :
        try :
            graph = geneModelToSpliceGraph(gene)
        except ValueError :
            graph = geneModelToSpliceGraph(gene, useCDS=True)
        graph.annotate()
        dataDict[sourcePath][SPLICEGRAPH_PLOT] = graph
        dataDict[sourcePath][GENE_MODEL_PLOT]  = gene
    elif plot.file_format == SPLICEGRAPH_FORMAT :
        graph = getFirstGraph(sourcePath)
        if plot.annotate : graph.annotate()
        dataDict[sourcePath][plot.plot_type] = graph
    elif plot.file_format == CSV_FORMAT :
        dataDict[sourcePath][plot.plot_type] = getXYData(sourcePath)
    elif plot.file_format == SAM_FORMAT :
        if not gene : raise ValueError('Must provide a gene to load gene-specific SAM data.')
        depthDict,jctDict = getSamReadData(sourcePath, maxpos=maxSAM, minjct=minjct, chromosomes=gene.chromosome, verbose=verbose)
        dataDict[sourcePath][READ_DEPTH_PLOT] = depthDict
        dataDict[sourcePath][JUNCTION_PLOT]   = jctDict

    return sourcePath

def loadGeneData(plots, plotConfig, maxSAM=0, **args) :
    """Goes through each plot type and loads data and plot
    information associated with the each type."""
    geneName = getAttribute('geneName', None, **args)
    models   = getAttribute('models', None, **args)
    verbose  = getAttribute('verbose', False, **args)

    ranges   = []
    dataDict = {}
    bgndGene = None
    gene     = None
    strand   = None
    chrom    = None
    minpos   = MAXINT
    maxpos   = 0
    adjMin   = minpos
    adjMax   = maxpos
    for plotType in DATA_LOAD_ORDER :
        for p in plots :
            if p.plot_type != plotType : continue

            if p.plot_type == GENE_MODEL_PLOT :
                model = models[p.source_file]
                name  = geneName if geneName else p.gene_name
                gene  = model.getGeneByName(name)

            # Limit for SAM files
            samLimit   = max(maxSAM, maxpos)
            sourcePath = loadData(p, dataDict, gene=gene, maxSAM=samLimit, **args)

            if p.plot_type in MINMAX_PLOT_TYPES :
                # Note: MINMAX_PLOT_TYPES data types must have minpos/maxpos
                #       attributes as well as strand and chromosome
                obj    = dataDict[sourcePath][p.plot_type]
                maxpos = max(maxpos, obj.maxpos)
                minpos = min(minpos, obj.minpos)
                if not strand : strand = obj.strand
                if not chrom  : chrom  = obj.chromosome

            if plotConfig.shrink_introns and p.plot_type in SHRINKABLE_PLOTS :
                if not ranges :
                    if SPLICEGRAPH_PLOT in dataDict[sourcePath] :
                        ranges = getGraphRanges(dataDict[sourcePath][SPLICEGRAPH_PLOT], scaleFactor=plotConfig.shrink_factor)
                    elif ISOFORM_PLOT in dataDict[sourcePath] :
                        ranges = getGraphRanges(dataDict[sourcePath][ISOFORM_PLOT], scaleFactor=plotConfig.shrink_factor)
                    else :
                        raise ValueError('No gene model, splice graph or isoform provided')

                if p.plot_type in [SPLICEGRAPH_PLOT, ISOFORM_PLOT] :
                    dataDict[sourcePath][p.plot_type] = adjustGraph(dataDict[sourcePath][p.plot_type], ranges=ranges)
                elif p.plot_type == GENE_MODEL_PLOT :
                    dataDict[sourcePath][GENE_MODEL_PLOT] = adjustGene(dataDict[sourcePath][GENE_MODEL_PLOT], ranges=ranges)
                    dataDict[sourcePath][SPLICEGRAPH_PLOT] = adjustGraph(dataDict[sourcePath][SPLICEGRAPH_PLOT], ranges=ranges)
                elif p.plot_type == READ_DEPTH_PLOT and chrom in dataDict[sourcePath][READ_DEPTH_PLOT] :
                    dataDict[sourcePath][READ_DEPTH_PLOT][chrom] = adjustDepths(dataDict[sourcePath][READ_DEPTH_PLOT][chrom], ranges=ranges)
                elif p.plot_type == JUNCTION_PLOT and chrom in dataDict[sourcePath][JUNCTION_PLOT] :
                    dataDict[sourcePath][JUNCTION_PLOT][chrom] = adjustJunctions(dataDict[sourcePath][JUNCTION_PLOT][chrom], ranges=ranges)

                if p.plot_type in [SPLICEGRAPH_PLOT, ISOFORM_PLOT, GENE_MODEL_PLOT] :
                    adjMax = max(adjMax, dataDict[sourcePath][p.plot_type].maxpos)
                    adjMin = min(adjMin, dataDict[sourcePath][p.plot_type].minpos)
                
            if p.plot_type == GENE_MODEL_PLOT :
                bgndGene = dataDict[sourcePath][GENE_MODEL_PLOT]

    if plotConfig.shrink_introns :
        return GenePlotSpecs(bgndGene, [bgndGene], dataDict, strand, chrom, adjMin, adjMax, origMin=minpos, origMax=maxpos)
    else :
        return GenePlotSpecs(bgndGene, [bgndGene], dataDict, strand, chrom, minpos, maxpos)

def loadIsoformWeights(geneName, path) :
    """Loads isoform weights from a PSGinfer output file
    and stores just the values for the given gene.  Expected
    format:
        gene-id<TAB>isoform-id<TAB>frequency
    """
    result = {}
    for line in ezopen(path) :
        if not line.startswith(geneName) : continue
        parts = line.strip().split('\t')
        result[parts[1]] = float(parts[2])
    return result if result else None

def getGeneModels(plots, **args) :
    """Returns a dictionary of GeneModel instances for the
    source files for all gene model plots."""
    verbose   = getAttribute('verbose', False, **args)
    genePlots = [p for p in plots if p.plot_type == GENE_MODEL_PLOT]
    result    = {}
    for p in genePlots :
        if p.source_file not in result :
            result[p.source_file] = loadGeneModels(p.source_file, verbose=verbose)
        elif verbose :
            sys.stderr.write('Gene models already loaded from %s\n' % p.source_file)
    return result

def plotLegend(patches, **args) :
    """Adds a legend to the current plot."""
    keys      = getAttribute('keys', None, **args)
    if not keys : keys = sorted(patches.keys())
    patchList = [patches[k] for k in keys]
    numCols   = int(math.ceil(math.sqrt(len(patches))))
    figlegend(patchList, keys, LEGEND_PLACEMENT, ncol=numCols)

def resolveSourcePath(p, **args) :
    """Resolves a source directory to a given file name if necessary."""
    fileName = getAttribute('fileName', None, **args)
    result   = p.source_file
    if os.path.isdir(result) :
        if p.file_format in MULTI_FILE_FORMATS :
            if fileName :
                result = os.path.join(result, fileName)
            else :
                raise ValueError('No file name specified for multi-plot %s : %s\n' % (p.plot_type, p.source_file))
        else :
            raise ValueError('Directory used when source file required for %s : %s\n' % (p.plot_type, p.source_file))
    return result
