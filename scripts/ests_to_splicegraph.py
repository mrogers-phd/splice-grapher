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
Script for converting EST/cDNA alignments into splice graphs.
"""
from SpliceGrapher.shared.config     import *
from SpliceGrapher.shared.utils      import *
from SpliceGrapher.formats.psl       import *
from SpliceGrapher.formats.loader    import *
from SpliceGrapher.predict.SpliceSiteValidator import *
from SpliceGrapher.SpliceGraph       import SpliceGraph, updateRoot, updateLeaf

from optparse import OptionParser
import sys, os

def pslCmp(a,b) :
    return (a.Tend-b.Tend) if a.Tstart == b.Tstart else (a.Tstart-b.Tstart)

USAGE = """%prog [options] psl-file

Converts EST alignments from a PSL alignment file into splice graphs by using ESTs
that fall within a gene on the same strand.  Target sequences referenced in the PSL
file must match names in the gene model file for graphs to be produced.  Only ESTs
that match a gene's strand will be included for that gene, so reverse-complement
matches should be resolved prior to running this script."""

parser = OptionParser(usage=USAGE)
parser.add_option('-a', dest='annotate',default=False,   help='Annotate the output graphs [default: %default]', action='store_true')
parser.add_option('-A', dest='adjust',  default=False,   help='Adjust graphs to gene boundaries (requires gene model) [default: %default]', action='store_true')
parser.add_option('-d', dest='dir',     default='.',     help='Output directory for splice graph files (overrides -o option) [default: %default]')
parser.add_option('-f', dest='fasta',   default=None,    help='FASTA reference to validate splice sites [default: %default]')
parser.add_option('-g', dest='gene',    default=None,    help='Gene name [default: infer from PSL]')
parser.add_option('-i', dest='intron',  default=4,       help='Minimum intron length allowed [default: %default]', type='int')
parser.add_option('-m', dest='model',   default=SG_GENE_MODEL, help='Gene model file [default: %default]')
parser.add_option('-s', dest='singles', default=False,   help='Allow single-exon ESTs to be included [default: %default]', action='store_true')
parser.add_option('-S', dest='suffix',  default='_ests', help='File name suffix for ouput files [default: %default]')
parser.add_option('-V', dest='valid',   default=False,   help='Only output files for valid graphs [default: %default]', action='store_true')
parser.add_option('-v', dest='verbose', default=False,   help='Use verbose output [default: %default]', action='store_true')

opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 1 :
    parser.print_help()
    sys.exit(1)

if not opts.model :
    parser.print_help()
    sys.stderr.write('** You must specify a reference gene model (-m).\n')
    sys.exit(1)

pslFile = args[0]
validateFile(pslFile)
validateFile(opts.model)
if opts.fasta : validateFile(opts.fasta)

writeStartupMessage()

geneModel = loadGeneModels(opts.model, verbose=opts.verbose)
if not geneModel :
    raise Exception('Unable to load gene model from %s\n' % opts.model)

validator = None
if opts.fasta :
    validator = SpliceSiteValidator(opts.fasta, verbose=opts.verbose)

# Store PSL records by target (usually chromosome) and strand:
pslRecords = loadPSLRecords(pslFile)
pslDict    = {}
for r in pslRecords :
    key = r.Tname.lower()
    pslDict.setdefault(key,{})
    pslDict[key].setdefault(r.strand,[])
    pslDict[key][r.strand].append(r)

# Sort each list by start/end positions:
for k1 in pslDict :
    for k2 in pslDict[k1] :
        pslDict[k1][k2].sort(cmp=pslCmp)

counts = {'-':0,'+':0}
for chrom in pslDict :
    geneList = geneModel.getGeneRecords(chrom)
    for gene in geneList :
        if gene.strand not in pslDict[chrom] :
            continue

        # Find PSL records on the same chromosome and strand that overlap the given gene boundaries.
        # estsToSpliceGraph() will only use the exons that overlap the gene.
        ## pslRecs = [r for r in pslDict[chrom][gene.strand] if r.Tstart >= gene.minpos and r.Tend <= gene.maxpos]
        pslRecs = [r for r in pslDict[chrom][gene.strand] if r.Tend > gene.minpos and r.Tstart < gene.maxpos]

        if not pslRecs :
            if opts.verbose : sys.stderr.write('** No PSL records for gene %s; skipping.\n' % gene.id)
            continue

        if opts.verbose :
            idString = ', '.join([r.Qname for r in pslRecs])
            sys.stderr.write('Creating graph for gene %s (%s) using %d PSL records: %s.\n' % (gene.id,gene.strand,len(pslRecs),idString))

        # additional options:
        #     includeSingles = single-exon ESTs
        if opts.verbose and len(gene.id) > 20 :
            sys.stderr.write('** Warning, possible bad gene name: %s\n' % gene.id)
        graph = estsToSpliceGraph(gene.id, pslRecs, chrom,
                                  geneBounds=[gene.minpos,gene.maxpos],
                                  min_intron=opts.intron,
                                  strand=gene.strand,
                                  validator=validator,
                                  verbose=opts.verbose,
                                  includeSingles=opts.singles)

        if opts.annotate :
            graph.annotate()

        try :
            outPath   = os.path.join(opts.dir, '%s%s.gff' % (gene.id,opts.suffix))
            graph.writeGFF(outPath, haltOnError=True)
            counts[gene.strand] += 1
        except ValueError, ve :
            if len(graph.nodeDict) == 0 :
                if opts.singles :
                    sys.stderr.write('Unable to create graph for %s: created 0 nodes from %d PSL records.\n' % (gene.id, len(pslRecs)))
                else :
                    sys.stderr.write('Unable to create graph for %s: %d records may only be single-exon ESTs.\n' % (gene.id, len(pslRecs)))
            else :
                raise ValueError(ve)

total = sum(counts.values())
sys.stdout.write('Created a total of %d graphs (%d + strand/%d - strand) in %s\n' % (total, counts['+'], counts['-'], os.path.abspath(opts.dir)))
