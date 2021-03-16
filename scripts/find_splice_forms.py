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
from SpliceGrapher.shared.utils              import *
from SpliceGrapher.shared.subgraph           import *
from SpliceGrapher.shared.ShortRead          import *
from SpliceGrapher.SpliceGraph               import *
from SpliceGrapher.shared.GeneModelConverter import *
from SpliceGrapher.formats.loader            import *
from SpliceGrapher.formats.sam               import *
from optparse                                import OptionParser
import os,sys

USAGE = """%prog SAM-file [options]

Uses existing gene models to predict splice forms represented in a set of RNA-Seq
alignments.  For each gene, it produces a splice graph to represent the splice forms
based on spliced alignments and read depths found in the SAM file.  If there is no
unique evidence to identify any splice form, no graph will be produced."""

parser = OptionParser(usage=USAGE)
parser.add_option('-d', dest='outdir',    default=None,  help='Output directory (overrides -o) [default: %default]')
parser.add_option('-D', dest='avgdepth',  default=1,     help='Minimum average depth for accepting a cluster [default: %default]', type='float')
parser.add_option('-j', dest='minjct',    default=2,     help='Minimum junction threshold [default: %default]', type='int')
parser.add_option('-m', dest='model',     default=SG_GENE_MODEL, help='Gene model file [default: %default]')
parser.add_option('-o', dest='output',    default=None,  help='Output file [default: stdout]')
parser.add_option('-O', dest='overlap',   default=1,     help='Minimum number of bases that a read cluster must overlap a feature [default: %default]', type='int')
parser.add_option('-S', dest='list',      default=None,  help='List of splice graphs to use as a reference (overrides -m) [default: %default]')
parser.add_option('-t', dest='threshold', default=1,     help='Minimum read coverage threshold for predicting exons [default: %default]', type='float')
parser.add_option('-v', dest='verbose',   default=False, help='Verbose mode [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 1 :
    parser.print_help()
    sys.exit(1)

if not (opts.model or opts.list) :
    parser.print_help()
    sys.stderr.write('** No GFF gene model specified.  Use the -m option or set SG_GENE_MODEL in your SpliceGrapher configuration.\n')
    sys.exit(1)

samFile = args[0]
validateFile(samFile)
if opts.model : validateFile(opts.model)
if opts.list  : validateFile(opts.list)

# Grab SAM alignments
# Note that gene selection is driven by data in SAM file
if opts.verbose : sys.stderr.write('Loading SAM file\n')
depthDict,jctDict = getSamReadData(samFile, minjct=opts.minjct, verbose=opts.verbose)
maxposDict = {}
chromSet   = set(depthDict.keys())
for chrom in chromSet :
    maxposDict[chrom] = len(depthDict[chrom])
    # In case there are no spliced alignments for one of the chromosomes:
    jctDict.setdefault(chrom,[])

# Load gene models
graphDict  = {}
indicator  = ProgressIndicator(10000, verbose=opts.verbose)
foundChrom = set()
if opts.list :
    if opts.verbose : sys.stderr.write('Loading reference splicegraphs\n')
    for line in ezopen(opts.list) :
        indicator.update()
        graph = getFirstGraph(line.strip())
        chrom = graph.chromosome
        graphDict.setdefault(chrom,{})
        graphDict[chrom][graph.getName()] = graph
        maxposDict.setdefault(chrom,0)
        maxposDict[chrom] = max(maxposDict[chrom], graph.maxpos)
        if chrom in chromSet :
            foundChrom.add(chrom)
else :
    geneModel  = loadGeneModels(opts.model, verbose=opts.verbose)
    # Grab chromosome-specific graphs from gene model
    useCDS     = False
    if opts.verbose : sys.stderr.write('Loading gene model splicegraphs\n')
    for chrom in chromSet :
        try :
            geneRecs = geneModel.getGeneRecords(chrom, geneFilter=gene_type_filter)
        except KeyError :
            sys.stderr.write('Warning: chromosome %s not found in gene model.\n' % chrom)
            continue

        foundChrom.add(chrom)
        for gene in geneRecs :
            if gene.minpos > maxposDict[chrom] : continue
            indicator.update()
            try :
                graph  = geneModelToSpliceGraph(gene, useCDS=useCDS)
            except KeyError :
                try :
                    useCDS = not useCDS
                    graph  = geneModelToSpliceGraph(gene, useCDS=useCDS)
                except KeyError :
                    continue
            except ValueError as ve :
                sys.stderr.write('Cannot build graph for gene %s: %s\n' % (gene.name, str(ve)))
                continue

            graphDict.setdefault(chrom,{})
            graphDict[chrom][graph.getName()] = graph
            maxposDict.setdefault(chrom,0)
            maxposDict[chrom] = max(maxposDict[chrom], graph.maxpos)
indicator.finish()

if not foundChrom :
    sys.stderr.write('** Error: gene models/splice graphs do not contain any of the chromosomes identified in the SAM file.\n')
    sys.exit(1)

missing = chromSet - foundChrom
if missing :
    sys.stderr.write('** Warning: not all chromosomes were found in gene models/splice graphs:\n')
    sys.stderr.write('     missing: %s\n' % ','.join(missing))
    sys.stderr.write('     found:   %s\n' % ','.join(foundChrom))

outStream = sys.stdout if not opts.output else open(opts.output,'w')
for c in graphDict :
    if opts.verbose : sys.stderr.write('Generating splice form graphs in chromosome %s\n' % c)

    for g in graphDict[c] :
        graph    = graphDict[c][g]
        chrom    = graph.chromosome.lower()
        depths,minpos = getGraphDepths(depthDict[chrom], graph)
        clusters = depthsToClusters(chrom, depths, minDepth=opts.avgdepth, threshold=opts.threshold, reference=minpos)

        junctions = []
        for jct in jctDict[chrom] :
            if jct.strand != graph.strand : continue
            if jct.maxpos < graph.minpos : continue
            if jct.minpos > graph.maxpos : break
            junctions.append(jct)

        if opts.verbose :
            sys.stderr.write('Identifying forms using %d clusters / %d junctions.\n' % (len(clusters), len(junctions)))
            sys.stderr.write('  clusters:\n')
            for ccc in clusters :
                sys.stderr.write('    %d-%d; avg depth=%.2f\n' % (ccc.minpos, ccc.maxpos, ccc.avgDepth()))
            sys.stderr.write('  junctions:\n')
            for jjj in junctions :
                sys.stderr.write('    %s; coverage=%d\n' % (jjj, jjj.count))

        # Core of the program: construct a subgraph that represents only
        # transcripts supported by clusters and junctions from the RNA-Seq data.
        newGraph = makeSupportedGraph(graph, clusters, junctions, verbose=opts.verbose, minOverlap=opts.overlap, mindepth=opts.avgdepth)

        if newGraph is not None : 
            if opts.outdir :
                if not os.path.isdir(opts.outdir) :
                    if opts.verbose : sys.stderr.write('Directory %s not found; creating it\n' % opts.outdir)
                    os.makedirs(opts.outdir)
                outFile = os.path.join(opts.outdir, '%s.gff'%newGraph.getName())
                if opts.verbose : sys.stderr.write('  creating %s\n' % outFile)
                outStream = open(outFile,'w')
            newGraph.annotate()
            newGraph.writeGFF(outStream)
            if opts.outdir :
                outStream.close()
        elif opts.verbose :
            sys.stderr.write('  insufficient splice form evidence for %s\n' % graph.getName())

if opts.verbose : sys.stderr.write('Done.\n')
