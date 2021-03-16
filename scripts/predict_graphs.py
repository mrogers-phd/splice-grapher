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
from SpliceGrapher.shared.config     import *
from SpliceGrapher.shared.utils      import *
from SpliceGrapher.shared.ShortRead  import *
from SpliceGrapher.formats.loader    import *
from SpliceGrapher.formats.sam       import *
from SpliceGrapher.predict.SpliceGraphPredictor import *

from optparse import OptionParser
import os,sys

def depthPredictions(chromosomes, dDict, jDict, model, opts) :
    """Run predictions based on alignment information from a SpliceGrapher depths file."""
    nGenes = nPred  = 0
    sys.stderr.write('Generating predictions:\n')
    chromList = sorted(chromosomes)
    for c in chromList :
        sys.stderr.write('  starting chromosome %s' % c)
        geneList  = model.getGeneRecords(c)
        nGenes   += len(geneList)
        geneList.sort()
        depths = dDict[c] if c in dDict else []
        jcts   = jDict[c] if c in jDict else []
        outDir = os.path.join(opts.outdir, c)
        if not os.path.isdir(outDir) : os.makedirs(outDir)
        nPred += predictGraphs(geneList, depths, jcts, outDir, depthThreshold=opts.avgdepth, novelOnly=opts.novel, minAvgDepth=opts.mindepth)
    return nPred, nGenes

def getCommonChromosomes(geneSet, depthSet, depthsFile, opts) :
    result = geneSet & depthSet
    if result :
        sys.stderr.write('Found %d chromosomes in common between %s and %s:\n' % \
                (len(result), opts.model, depthsFile))
        sys.stderr.write('  %s\n' % ', '.join(sorted(result)))
    else :
        sys.stderr.write('Found no chromosomes in common between %s and %s:\n' % (opts.model, depthsFile))
        sys.stderr.write('  Gene models: %s\n' % ', '.join(sorted(geneSet)))
        sys.stderr.write('  Depths file: %s\n' % ', '.join(sorted(depthSet)))
        sys.exit(1)
    return result

USAGE = """%prog SAM-file [options]

Predicts splice graphs for all genes based on evidence in the given SAM file.
Predicted graphs will be put into separate files under the given output directory,
creating a subdirectory for each chromosome."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-a', dest='anchor',   default=8,   help='Minimum anchor required for junctions [default: %default]', type='int')
parser.add_option('-d', dest='outdir',   default='.', help='Output file [default: current directory]')
parser.add_option('-j', dest='minjct',   default=2,   help='Minimum coverage required for splice junctions [default: %default]', type='int')
parser.add_option('-D', dest='avgdepth', default=1,   help='Minimum average read depth required for clusters [default: %default]', type='float')
parser.add_option('-m', dest='model',    default=SG_GENE_MODEL, help='Gene models file [default: %default]')
parser.add_option('-N', dest='novel',    default=False, help='Only write splice graphs that were modified [default: %default]', action='store_true')
parser.add_option('-v', dest='verbose',  default=False, help='Verbose mode [default: %default]', action='store_true')
parser.add_option('--mindepth', dest='mindepth',   default=None,
        help='Minimum gene-wise read-depth for making predictions [default: %default]', type='int')
opts, args = parser.parse_args(sys.argv[1:])

MIN_ARGS = 1
if len(args) != MIN_ARGS :
    parser.print_help()
    if args : sys.stderr.write('\nExpected %d parameters; received %d:\n  %s\n' % (MIN_ARGS, len(args), '\n  '.join(args)))
    sys.exit(1)

if not opts.model :
    parser.print_help()
    if args : sys.stderr.write('\nYou must supply gene models (-m option) to make predictions\n')
    sys.exit(1)

if opts.mindepth is None :
    opts.mindepth = 0

dataFile = args[0]
validateFile(dataFile)
validateFile(opts.model)

# Load reference information from gene models
model        = loadGeneModels(opts.model, verbose=opts.verbose)
geneChromSet = set(model.getChromosomes())

if isDepthsFile(dataFile) :
    if opts.verbose : sys.stderr.write('Loading depth information from %s\n' % dataFile)
    depthDict,jctDict = readDepths(dataFile, minjct=opts.minjct, minanchor=opts.anchor, verbose=opts.verbose)
    depthChromSet     = set(depthDict.keys())
    commonChromSet    = getCommonChromosomes(geneChromSet, depthChromSet, dataFile, opts)
    predCount, totalGenes = depthPredictions(commonChromSet, depthDict, jctDict, model, opts)
else :
    # SAM-based predictions:
    # Load header information about the chromosomes stored in the SAM
    # file and get the first (seed) alignment.
    if opts.verbose : sys.stderr.write('Initializing SAM input from %s\n' % dataFile)
    samStream      = ezopen(dataFile)
    cset, seedLine = getSamHeaderInfo(samStream)
    depthChromSet  = set([c.lower() for c in cset])

    # Find chromosomes common to both gene models and alignment records:
    commonChromSet = getCommonChromosomes(geneChromSet, depthChromSet, dataFile, opts)

    # Track chromosomes that already have been processed
    usedChrom  = {}.fromkeys(depthChromSet,False)
    totalGenes = predCount  = 0
    chrBuffer  = []
    sys.stderr.write('Generating predictions:\n')

    # We will have a seed line (first alignment) for each new chromosome
    while seedLine :
        c = None
        # Only want data relevant to gene models
        while c not in commonChromSet :
            chrBuffer, seedLine = getNextSamChromosome(samStream, seedLine, verbose=opts.verbose)
            if len(chrBuffer) <= 1 : break
            c = chrBuffer[0].split('\t')[2].lower()
            if c not in commonChromSet :
                sys.stderr.write('  ** skipping unrecognized chromosome %s in %s **\n' % (c, dataFile))

        if len(chrBuffer) <= 1 : break

        try :
            if usedChrom[c] : raise ValueError('Unsorted file: chromosome %s found twice\n' % c)
            usedChrom[c] = True
        except KeyError :
            raise ValueError('Chromosome %s not found in %s\n' % dataFile)

        sys.stderr.write('  starting chromosome %s\n' % c)
        sys.stderr.write('    loading alignment records from %s\n' % dataFile)
        depthDict, jctDict = getSamReadData(chrBuffer, minjct=opts.minjct, minanchor=opts.anchor)

        sys.stderr.write('    predicting splice graphs ')
        geneList    = model.getGeneRecords(c)
        totalGenes += len(geneList)
        geneList.sort()

        depths = depthDict[c] if c in depthDict else []
        jcts   = jctDict[c] if c in jctDict else []
        outDir = os.path.join(opts.outdir, c)
        if not os.path.isdir(outDir) : os.makedirs(outDir)
        predCount += predictGraphs(geneList, depths, jcts, outDir, depthThreshold=opts.avgdepth, novelOnly=opts.novel, minAvgDepth=opts.mindepth)
        remaining = [c for c in commonChromSet if not usedChrom[c]]
        if not remaining : break

sys.stderr.write('Finished: splice graphs were modified for %s / %s genes.\n' % \
        (commaFormat(predCount), commaFormat(totalGenes)))
