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
Python script for classifying putative splice sites within a FASTA file.
"""
from SpliceGrapher.shared.config            import *
from SpliceGrapher.predict.SiteClassifier   import *
from SpliceGrapher.shared.utils             import *
from SpliceGrapher.shared.ShortRead         import *
from SpliceGrapher.shared.streams           import *
from SpliceGrapher.shared.dna               import *
from SpliceGrapher.formats.FastaLoader      import FastaLoader
from SpliceGrapher.formats.loader           import *
from SpliceGrapher.predict.SpliceSite       import *
from SpliceGrapher.predict.ClassifierConfig import *

from optparse import OptionParser
import sys, os.path

try :
    from PyML import *
    from PyML.containers import SequenceData
    from PyML.containers import Labels
    from PyML.classifiers import SVM
    from PyML.classifiers.svm import loadSVM
except Exception :
    sys.stderr.write('\n** Unable to import PyML modules required for this script.\n')
    sys.exit(1)

DEFAULT_DIMERS = ['gt', 'gc', 'ag']
DIMER_STRING   = ','.join(DEFAULT_DIMERS)

#------------------------
# Methods
#------------------------
def logMsg(s) :
    """
    Convenience method for outputting log messages.
    """
    sys.stderr.write('%s\n' % s)

def classifyGenes(classifiers, geneModel, seqDict, **args) :
    """
    Main method for classifying all locations within known genes.
    """
    posRange  = getAttribute('posRange', None, **args)
    clusters  = getAttribute('clusters', None, **args)
    outStream = getAttribute('outStream', sys.stdout, **args)
    outputAll = getAttribute('outputAll', False, **args)
    verbose   = getAttribute('verbose', False, **args)

    if len(seqDict) != 1 : raise ValueError('Must have exactly 1 sequence: found %d' % len(seqDict))
    keys  = seqDict.keys()
    chrom = keys[0]

    # Track performance on known sites
    truePos  = {}.fromkeys(classifiers.keys(), 0)
    falseNeg = {}.fromkeys(classifiers.keys(), 0)

    chromosome = geneModel.getChromosome(chrom)
    if not chromosome :
        raise Exception('Chromosome %s not found in gene model.' % chrom)

    # Absolute minimum window size required by classifiers
    # is the maximum intron/exon size used by any of them.
    minWindow = 0
    for svm in classifiers.values() :
        minWindow = max(minWindow, max(svm.config.exonSize(), svm.config.intronSize()))
    minWindow += 1

    genes = geneModel.getGeneRecords(chrom)
    genes.sort()
    seq   = seqDict[chrom]
    if posRange :
        positionRange = posRange
    else :
        positionRange = [0, len(seq)]

    for i in range(len(genes)) :
        # Note that we may include a gene that starts before the given range or that
        # ends after the range.  By allowing overlaps we ensure that all genes will
        # be covered, even if some predictions might be duplicated.
        gene = genes[i]
        if gene.maxpos < positionRange[0] or gene.minpos > positionRange[1] : continue

        #----------------------------------------------------------
        # First check intergenic regions that contain clusters:
        if clusters :
            lower = positionRange[0] if i == 0 else genes[i-1].maxpos
            lower = max(lower, minWindow)
            upper = min(gene.minpos, positionRange[1]-minWindow)
            if lower < upper :
                hasReads = False
                for c in clusters :
                    if lower <= c.maxpos and c.minpos <= upper :
                        hasReads = True
                        break
                if hasReads :
                    classifySequence(chrom, '-', classifiers, outStream, minpos=lower, maxpos=upper, knownDict=None, verbose=verbose)
                    classifySequence(chrom, '+', classifiers, outStream, minpos=lower, maxpos=upper, knownDict=None, verbose=verbose)

        #----------------------------------------------------------
        # Now check sites within the gene.
        # Storing known sites allows us to track classifier
        # performance for sites we already know about.
        knownDict                = {}
        knownDict[ACCEPTOR_SITE] = gene.acceptorList()
        knownDict[DONOR_SITE]    = gene.donorList()

        if verbose : logMsg('  Classifying sites for gene %s (%d-%d)' % (gene.id, gene.minpos, gene.maxpos))
        (tp,fn) = classifySequence(chrom, gene.strand, classifiers, outStream,
                                   minpos=gene.minpos, maxpos=gene.maxpos, knownDict=knownDict, verbose=verbose, outputAll=outputAll)

        for t in classifiers.keys() :
            truePos[t]  += tp[t]
            falseNeg[t] += fn[t]

    if verbose :
        for t in classifiers.keys() :
            total = truePos[t] + falseNeg[t]
            if total > 0 :
                logMsg('Performance on %s sites: %d TP, %d FN, ACC = %.5f' % (t, truePos[t], falseNeg[t], float(truePos[t])/(truePos[t]+falseNeg[t])))
            else :
                logMsg('Performance on %s sites: %d TP, %d FN, ACC = 0.0' % (t, truePos[t], falseNeg[t]))

def classifySequence(chrom, strand, classifiers, outStream, **args) :
    """
    Classifies a single DNA sequence and writes any positive classifications to the given output stream.
    (Note: DNA sequence is obtained through the 'getSequence' method.)
    """
    minpos     = getAttribute('minpos', 0, **args)
    maxpos     = getAttribute('maxpos', sys.maxint, **args)
    knownDict  = getAttribute('knownDict', None, **args)
    clusters   = getAttribute('clusters', None, **args)
    verbose    = getAttribute('verbose', False, **args)
    outputAll  = getAttribute('outputAll', False, **args)

    # Track known site performance
    truePos  = {}.fromkeys(classifiers.keys(), 0)
    falseNeg = {}.fromkeys(classifiers.keys(), 0)

    # Turn off PyML messages:
    hideStdout()

    k = 0 # k, for'kluster' index
    for i in range(minpos, maxpos) :
        if clusters :
            # Use clusters to direct classification:
            while i > clusters[k].maxpos :
                k += 1
                if k >= len(clusters) : break
            if not clusters[k].contains(i) : continue

        for cName in classifiers.keys() :
            svm       = classifiers[cName]
            siteDimer = svm.config.dimer()
            siteType  = ACCEPTOR_SITE if svm.config.acceptor() else DONOR_SITE

            (exonSeq, intronSeq, dimer) = getSpliceSiteParts(chrom, i, siteType, strand, getSequence,
                                                             exonWindow=svm.config.exonSize(),
                                                             intronWindow=svm.config.intronSize())
            if dimer.lower() == siteDimer.lower() :
                newSeq = intronSeq + exonSeq if svm.config.acceptor() else exonSeq + intronSeq

                if svm.config.mismatchProfile() :
                    ssData = positionalKmerData(svm.config, [newSeq])
                else :
                    ssData = SequenceData([newSeq], mink=svm.config.mink(), maxk=svm.config.maxk(), maxShift=svm.config.maxShift())

                if svm.config.normalize()       : ssData.attachKernel('cosine')
            
                (ssClass,score) = svm.classify(ssData,0)

                if knownDict and (i in knownDict[siteType]) :
                    if ssClass == 1 :
                        truePos[cName] += 1
                    else :
                        falseNeg[cName] += 1
                    continue
            
                # Output the site and its score
                if ssClass == 1 or outputAll :
                    siteCode = ACCEPTOR_CODE if svm.config.acceptor() else DONOR_CODE
                    outStream.write('%s\t%s\t%d\t%s\t%.6f\t%s\n' % (chrom, strand, i, siteDimer, score, siteCode))

    showStdout()
    return (truePos,falseNeg)

def getSequence(chrom, pos1, pos2, strand) :
    return seqDict.subsequence(chrom, pos1, pos2, reverse=(strand=='-'))

USAGE="""%prog chromosome [options]

Applies one or more SVMs to the chromosome from the given FASTA
reference sequence file to classify novel splice sites."""

parser = OptionParser(usage=USAGE)
parser.add_option('-A', dest='all',       default=False, help='Output all results, not just positive scores [default: %default]', action='store_true')
parser.add_option('-c', dest='cfg_files', default=None,  help='list of dimer configuration files to load [default=%default]')
parser.add_option('-f', dest='fasta',     default=SG_FASTA_REF,  help='FASTA reference file reference [default=%default]')
parser.add_option('-m', dest='model',     default=SG_GENE_MODEL, help='GFF genome reference [default=%default]')
parser.add_option('-o', dest='output',    default=None,  help='output file [default=stdout]')
parser.add_option('-r', dest='range',     default=None,  help='prediction position range (e.g., "10000,20000") [default=whole chromosome]')
parser.add_option('-v', dest='verbose',   default=False, help='verbose mode [default: %default]', action='store_true')

#=======================================================
# Main program
opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 1 :
    parser.print_help()
    sys.exit(1)

if not opts.cfg_files :
    parser.print_help()
    sys.stderr.write('You must enter at least one classifier configuration file to use.\n')
    sys.exit(1)

errStrings = []
if not opts.model : errStrings.append('** No GFF gene model specified.  Use the -m option or set SG_GENE_MODEL in your SpliceGrapher configuration.\n')
if not opts.fasta : errStrings.append('** No FASTA reference specified.  Use the -f option or set SG_FASTA_REF in your SpliceGrapher configuration.\n')
if errStrings :
    parser.print_help()
    sys.stderr.write('\n%s\n' % '\n'.join(errStrings))
    sys.exit(1)

chrom = args[0].lower()
#-------------------------------------------------------
# Required parameters:
validateFile(opts.model)
validateFile(opts.fasta)
writeStartupMessage()

#-------------------------------------------------------
# Optional parameters:
if opts.model : validateFile(opts.model)
## if opts.depths  : validateFile(opts.depths)

posRange = [0,sys.maxint]
if opts.range :
    parts    = [int(x) for x in opts.range.split(',')]
    posRange = (min(parts), max(parts))

#-------------------------------------------------------
# Load the splice site classifiers
classifiers = {}
for ctype in opts.cfg_files.split(',') :
    if opts.verbose : logMsg("\nLoading '%s' classifier:" % ctype)
    classifiers[ctype] = SiteClassifier(ctype, verbose=opts.verbose)

#-------------------------------------------------------
# Load the FASTA reference to be classified:
if opts.verbose :
    logMsg('SiteClassifier running in verbose mode')
    logMsg('Loading FASTA reference from %s' % opts.fasta)

seqDict = FastaLoader(opts.fasta, sequenceID=chrom, verbose=opts.verbose)

# Load the cluster file for guiding chromosome-wide classification (deprecated)
clusters = None

# Load reference gene model
geneModel = loadGeneModels(opts.model, verbose=opts.verbose)

outStream = file(opts.output, 'w', 0) if opts.output else sys.stdout

# Perform the classification:
classifyGenes(classifiers, geneModel, seqDict,
              posRange=posRange,
              clusters=clusters,
              outStream=outStream,
              outputAll=opts.all,
              verbose=opts.verbose)
