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
Example script that demonstrates the splice graph prediction module.
"""
from SpliceGrapher.shared.config                import *
from SpliceGrapher.predict.SpliceGraphPredictor import *
from SpliceGrapher.shared.GeneModelConverter    import *
from SpliceGrapher.formats.sam                  import *
from SpliceGrapher.formats.loader               import *
from SpliceGrapher                              import SpliceGraph

from optparse import OptionParser
import sys, os

USAGE = """%prog gene-or-file [options]

Predicts an output splice graph given a baseline graph and additional input such as short read
data and putative splice junction data.  The baseline may be given either as a gene name or a
file (if a file with the given name is not found, it is assumed to be a gene name).

If a gene name is given, a gene model must also be provided using SG_GENE_MODEL or the -m option."""

parser = OptionParser(usage=USAGE)
parser.add_option('-d', dest='alignments',default=None, help='SAM file containing RNA-Seq alignment data [default: %default]')
parser.add_option('-j', dest='junctions', default=None, help='Alternate SAM file with spliced alignments (deprecated)')
parser.add_option('-J', dest='jmindepth', default=1,    help='Minimum depth required for junction evidence [default: %default]', type='int')
parser.add_option('-m', dest='model',     default=SG_GENE_MODEL, help='File for output splice graph (GFF) [default: %default]')
parser.add_option('-M', dest='minanchor', default=0,    help='Minimum anchor size required for junction evidence [default: %default]', type='int')
parser.add_option('-o', dest='output',    default=sys.stdout, help='File for output splice graph (GFF) [default: stdout]')
parser.add_option('-s', dest='graphs',    default=None, help='Comma-separated list of splice graph files to augment baseline graph [default: %default]')
parser.add_option('-T', dest='threshold', default=1,    help='Minimum depth threshold for identifying clusters [default: %default]', type='int')
parser.add_option('-v', dest='verbose',   default=False, action='store_true', help='use verbose output [default: %default]')
parser.add_option('--CDS', dest='useCDS', default=False, action='store_true', help='use mRNA/CDS records to infer isoforms [default: %default]')
opts, args   = parser.parse_args(sys.argv[1:])

if len(args) != 1 :
    parser.print_help()
    sys.exit(1)

if opts.junctions : sys.stderr.write('\n** Note: -j option will be deprecated in later versions.  Use a single SAM file instead.\n\n')

useInputFile = os.path.exists(args[0])

if useInputFile :
    graphFile = args[0]
    validateFile(graphFile)
else :
    geneName = args[0]
    if opts.verbose : sys.stderr.write('No file found, assuming %s is a gene name.\n' % geneName)
    if not opts.model :
        parser.print_help()
        sys.stderr.write('You must specify a gene model if you enter a gene name.\n')
        sys.exit(1)
    validateFile(opts.model)

if opts.alignments    : validateFile(opts.alignments)
if opts.junctions : validateFile(opts.junctions)
if opts.graphs :
    for s in opts.graphs.split(',') :
        validateFile(s)
writeStartupMessage()

if opts.verbose :
    sys.stderr.write('Predicting splice graph with the following inputs:\n')
    if useInputFile :
        sys.stderr.write('  Reference baseline:  %s\n' % graphFile)
    else :
        sys.stderr.write('  Gene name:           %s\n' % geneName)
    if opts.graphs    : sys.stderr.write('  Additional graphs:   %s\n' % opts.graphs)
    if opts.alignments    : sys.stderr.write('  Short reads:         %s\n' % opts.alignments)
    if opts.junctions : sys.stderr.write('  Splice junctions:    %s\n' % opts.junctions)

if useInputFile :
    origGraph = SpliceGraph.getFirstGraph(graphFile)
else :
    geneModel = loadGeneModels(opts.model, verbose=opts.verbose)
    gene      = geneModel.getGeneByName(geneName)
    if not gene : raise ValueError('Gene %s not found in %s' % (geneName, opts.model))
    origGraph = geneModelToSpliceGraph(gene,useCDS=opts.useCDS)

predictor = SpliceGraphPredictor(origGraph, verbose=opts.verbose)

if opts.graphs :
    for f in opts.graphs.split(',') :
        graph = SpliceGraph.getFirstGraph(f)
        if graph :
            predictor.addSpliceGraph(graph, mergeEnds=True)
        else :
            raise Exception('No graphs found in %s' % f)

depths    = []
junctions = []
if opts.alignments :
    depthDict, jctDict = getSamReadData(opts.alignments, maxpos=origGraph.maxpos, minjct=opts.jmindepth, minanchor=opts.minanchor, verbose=opts.verbose)
    try :
        depths = depthDict[origGraph.chromosome]
    except KeyError :
        sys.stderr.write('** Warning: no alignment information found for %s in %s' % (origGraph.chromosome, opts.alignments))

    predictor.setReadDepths(depths, minDepth=opts.threshold, verbose=opts.verbose)

    if jctDict :
        try :
            junctions = jctDict[origGraph.chromosome]
        except KeyError :
            sys.stderr.write('** Warning: SAM file %s does not contain any spliced alignments for chromosome %s\n' % (opts.alignments, origGraph.chromosome))
            junctions = []

if opts.junctions :
    jctDict = getSamJunctions(opts.junctions, maxpos=origGraph.maxpos, minanchor=opts.minanchor, minjct=opts.jmindepth, verbose=opts.verbose)
    try :
        junctions = jctDict[origGraph.chromosome]
    except KeyError :
        raise ValueError('** Warning: SAM file %s does not contain any spliced alignments for chromosome %s' % (opts.junctions, origGraph.chromosome))

if junctions :
    predictor.addSpliceJunctions(junctions, verbose=opts.verbose)

predictor.updateGraph(output=opts.output)
