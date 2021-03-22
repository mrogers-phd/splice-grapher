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
Plotter that uses a config file to define
a sequence of plots to appear on one page.
"""
from SpliceGrapher.shared.utils       import *
from SpliceGrapher.plot.PlotterConfig import *
from optparse import OptionParser
import os,sys

#==========================================================================
# Incompatibilities in some matplotlib versions may yield runtime warnings
import warnings
warnings.filterwarnings('ignore')
#==========================================================================

def helpFormat(strings, offset="    ", count=6):
    """Formats a list of strings for extended help."""
    result = ""
    for i in range(0,len(strings),count):
        if result:
            result += '\n'
        result += offset + ', '.join(strings[i:i+count])
    return result

EXTENDED_HELP="""

---------------------- Extended Help ---------------------- 
Configuration file starts with main heading as follows
(values below are examples only):

    [SinglePlotConfig]
    legend           = True
    output_file      = plot.pdf
    shrink_introns   = True

Valid main section tags are:
%s

This is followed by one or more plot sections for the same gene,
each with a unique name of your choosing:

    [GeneModelGraph]
    plot_type     = gene
    file_format   = gene_model
    source_file   = MX2_model.gff
    gene_name     = MX2
    relative_size = 1.0
    title_string  = Gene Model for MX2

    [SpliceGrapherPrediction]
    relative_size = 1.5
    plot_type     = splice_graph
    source_file   = MX2_pred.gff
    title_string  = SpliceGrapher MX2 Prediction

Valid plot section tags:
%s

Valid plot_type names:
%s

Valid file_format names:
%s

See the SpliceGrapher tutorial for more information on these features.
""" % (helpFormat(VALID_SINGLE_TAGS), helpFormat(VALID_PLOT_TAGS), helpFormat(VALID_PLOT_TYPES), helpFormat(VALID_FILE_FORMATS))

DEFAULT_FONT   = 12
DEFAULT_HEIGHT = 11.0
DEFAULT_WIDTH  = 8.5

VALID_FORMATS = ['emf', 'eps', 'pdf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz']

def setValue(argv, key, cfgValue, default):
    """Simple method that performs the logic for assigning a value to a
    parameter.  If there is no configuration value, returns the default.
    If there is a configuration value, the result depends on whether the
    user entered something on the command line."""
    if not cfgValue:
        return default
    return default if key in argv else cfgValue

USAGE = """%prog config-file [options]

Produces a plot for a single gene following guidelines specified in a
configuration file.  If the -o option is used, the output format is
determined by the file extension.  Valid file extensions are:
  """ + ', '.join(VALID_FORMATS)

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-F', dest='fontsize',default=DEFAULT_FONT,   help='Font size (points) [default: %default-point]', type='float')
parser.add_option('-H', dest='height',  default=DEFAULT_HEIGHT, help='Plotting area width (inches) [default: %default]', type='float')
parser.add_option('-W', dest='width',   default=DEFAULT_WIDTH,  help='Plotting area width (inches) [default: %default]', type='float')
parser.add_option('-o', dest='output',  default=None,           help='Output file [default: stdout]')
parser.add_option('-v', dest='verbose', default=False,          help='Verbose mode [default: %default]', action='store_true')
parser.add_option('--extended', dest='extended', default=False,  help='Output extended help text and quit [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if opts.extended:
    parser.print_help()
    sys.stderr.write(EXTENDED_HELP)
    sys.exit(1)

if len(args) != 1:
    parser.print_help()
    sys.exit(1)

cfgFile     = args[0]
validateFile(cfgFile)
config      = SinglePlotConfig(cfgFile)
displayList = config.getPlotList()
if not displayList:
    sys.stderr.write('No plots specified in config file; nothing to do!\n')
    sys.exit(1)
plotConfig = config.getConfiguration()
for p in displayList:
    validateFile(p.source_file)

#=============================================================
# Save matplotlib includes until after parsing arguments
import matplotlib

outputFile = opts.output if opts.output else plotConfig.output_file
if outputFile:
    foo,ext = os.path.splitext(outputFile)
    ext     = ext.replace('.','')
    if ext.lower() not in VALID_FORMATS:
        raise ValueError('Unrecognized output format: %s' % ext.upper())
    matplotlib.use('agg')

from SpliceGrapher.plot.PlotUtils import *
#=============================================================

writeStartupMessage()

geneModels = getGeneModels(displayList, verbose=opts.verbose)
geneSpecs  = loadGeneData(displayList, plotConfig, models=geneModels, verbose=opts.verbose)

if geneSpecs.strand is None:
    raise Exception('No strand found for plots')
if geneSpecs.chrom is None:
    raise Exception('No chromosome found for plots')
if geneSpecs.minpos == MAXINT or geneSpecs.maxpos == 0:
    raise Exception('Unable to establish min/max plot boundaries.')

# Set values based on command-line > configuration > defaults relation
width    = setValue(sys.argv, '-W', plotConfig.width, opts.width)
height   = setValue(sys.argv, '-H', plotConfig.height, opts.height)
fontsize = setValue(sys.argv, '-F', plotConfig.fontsize, opts.fontsize)
geneSpecs.setFontSize(fontsize)

if opts.verbose:
    sys.stderr.write('Display area %.1f x %.1f; %.1fpt font size\n' % (width, height, fontsize))

initializePlots(geneSpecs, width, height, displayList)
titlePadding = getTitlePadding(fontsize)
topLine      = DISPLAY_HEIGHT - titlePadding
availHeight  = getAvailableHeight(displayList, fontsize, legend=plotConfig.legend)

patchDict = {}
for p in displayList:
    if p.hide:
        continue
    height = availHeight * p.relative_size
    curAxes = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
    topLine = topLine - height - titlePadding
    patches,ignore = generatePlot(p, geneSpecs, curAxes, shrink=plotConfig.shrink_introns, xTickLabels=(p==displayList[-1]), verbose=True)
    patchDict.update(patches)

if plotConfig.legend and patchDict:
    plotLegend(patchDict)

if outputFile:
    if opts.verbose:
        sys.stderr.write('Writing graph output to %s\n' % outputFile)
    savefig(outputFile, dpi=400)
else:
    show()
