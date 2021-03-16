#! /usr/bin/env python
# 
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
Script for converting a GFF gene model into a splice graph.
"""
from SpliceGrapher.shared.config     import *
from SpliceGrapher.shared.utils      import *
from SpliceGrapher.shared.GeneModelConverter import *
from SpliceGrapher.formats.loader    import *
from SpliceGrapher.formats.GeneModel import GENE_TYPE, defaultGeneFilter, gene_type_filter

from optparse    import OptionParser

import sys, os

USAGE = """%prog [options]

Converts GFF gene models into splice graphs.  Produces a splice graph
for every gene in its own separate file."""

parser = OptionParser(usage=USAGE)
parser.add_option('-A',    dest='alltypes',  default=False, help='Include all gene types [default: protein-coding only]', action='store_true')
parser.add_option('-a',    dest='annotate',  default=False, help='annotate any AS found in the model [default: %default]', action='store_true')
parser.add_option('-c',    dest='chrom',     default=None,  help='generate graphs only for genes in the given chromosome [default: %default]')
parser.add_option('--CDS', dest='cds',       default=False, help='use CDS/UTR records for exon boundaries [default: exon records]', action='store_true')
parser.add_option('-d',    dest='dir',       default=None,  help='top-level directory for gene file(s) [default: %default]')
parser.add_option('-g',    dest='genes',     default=None,  help='generate graphs for genes given either in a comma-separated list or a file [default: %default]')
parser.add_option('-e',    dest='minexon',   default=1,     help='minimum allowable exon size [default: %default]', type='int')
parser.add_option('-i',    dest='minintron', default=1,     help='minimum allowable intron size [default: %default]', type='int')
parser.add_option('-m',    dest='model',     default=SG_GENE_MODEL, help='optional GFF3/GTF gene model file [default: %default]')
parser.add_option('-n',    dest='name',      default=False, help='use gene name instead of gene id for output filename [default: %default]', action='store_true')
parser.add_option('-o',    dest='output',    default=None,  help='file for output splice graph [default: %default]')
parser.add_option('-S',    dest='subdirs',   default=False, help='create chromosome subdirectories (only with -d option) [default: %default]', action='store_true')
parser.add_option('-v',    dest='verbose',   default=False, help='use verbose output [default: %default]', action='store_true')

opts, args = parser.parse_args(sys.argv[1:])

if len(args) > 0 :
    parser.print_help()
    if args : sys.stderr.write('\nReceived %d unexpected parameters:\n  %s\n' % (len(args), '\n  '.join(args)))
    sys.exit(1)

if not opts.model :
    parser.print_help()
    sys.stderr.write('** No GFF gene model specified.  Use the -m option or set SG_GENE_MODEL in your SpliceGrapher configuration.\n')
    sys.exit(1)

validateFile(opts.model)
writeStartupMessage()
gene_filter = defaultGeneFilter if opts.alltypes else gene_type_filter

# Load gene model
geneModel = loadGeneModels(opts.model, verbose=opts.verbose, alltypes=True)
geneList  = []
if opts.genes :
    geneString = opts.genes
    if os.path.isfile(opts.genes) :
        allGenes   = [s.strip() for s in ezopen(opts.genes)]
        geneString = ','.join(allGenes)

    for gid in geneString.split(',') :
        gene = geneModel.getGeneByName(gid)
        if not gene :
            sys.stderr.write('** Warning: gene %s not found\n' % gid)
        else :
            geneList.append(gene)
elif opts.chrom :
    geneList = geneModel.getGeneRecords(opts.chrom, geneFilter=gene_filter, verbose=opts.verbose)
else :
    geneList = geneModel.getAllGenes(geneFilter=gene_filter, verbose=opts.verbose)

if not geneList :
    sys.stderr.write('No genes found; exiting.\n')
    sys.exit(1)

geneList.sort()
if opts.verbose : sys.stderr.write('Generating graphs for %d gene models\n' % len(geneList))

outStream = sys.stdout
if opts.output :
    outFile = os.path.join(opts.dir, opts.output) if opts.dir else opts.output
    outStream = open(outFile,'w')

indicator = ProgressIndicator(10000, description='genes', verbose=opts.verbose)
for gene in geneList :
    indicator.update()

    try :
        graph = makeSpliceGraph(gene, minexon=opts.minexon, verbose=opts.verbose)
    except ValueError :
        sys.stderr.write('Warning: failed to convert gene %s to splice graph\n' % gene.name)
        #sys.stderr.write('  there are %d isoforms and %d mrna\n' % (len(gene.isoforms), len(gene.mrna)))
        continue

    if graph.minpos > graph.maxpos or len(graph) == 0 :
        sys.stderr.write('Warning: failed to convert gene %s to splice graph\n' % gene.name)
        continue

    # Annotate any AS events found in the gene model
    if opts.annotate : graph.annotate()

    if opts.dir and not opts.output :
        outDir = opts.dir
        if opts.subdirs :
            base = os.path.basename(opts.dir)
            if base != gene.chromosome :
                outDir = os.path.join(opts.dir, gene.chromosome)
                if not os.path.exists(outDir) :
                    os.makedirs(outDir)
                    if opts.verbose :
                        indicator.reset()
                        sys.stderr.write('created directory %s\n' % outDir)
        fileName  = gene.name if opts.name else gene.id
        outFile   = os.path.join(outDir, '%s.gff' % fileName.upper())
        outStream = open(outFile, 'w')

    graph.writeGFF(outStream)

indicator.finish()
