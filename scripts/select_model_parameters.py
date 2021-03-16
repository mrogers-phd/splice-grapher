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
Script that tries a variety of parameters to identify an optimal set
for classifying a particular splice-site dimer.
"""
from SpliceGrapher.shared.config            import *
from SpliceGrapher.formats.fasta            import *
from SpliceGrapher.shared.dna               import *
from SpliceGrapher.shared.utils             import *
from SpliceGrapher.shared.streams           import *
from SpliceGrapher.predict.SpliceSite       import truncateSequences
from SpliceGrapher.predict.ClassifierConfig import ClassifierConfig

from optparse import OptionParser
import sys, os

try :
    from PyML import *
    from PyML.containers  import SequenceData, Labels
    from PyML.classifiers import SVM, modelSelection
except Exception :
    sys.stderr.write('\n** Unable to import PyML modules required for this script.\n')
    sys.exit(1)

C_DEFAULTS = [10.0**x for x in range(-2,3)]
C_DEF_STR  = ','.join(['%.4f'%x for x in C_DEFAULTS])

# Sizes larger than these shouldn't be necessary and
# may cause training to fail to converge
MIN_EXON_SIZE    = 5
MAX_EXON_SIZE    = 50
MIN_INTRON_SIZE  = 4
MAX_INTRON_SIZE  = 50
MIN_MINK         = 1
MAX_MAXK         = 10
MIN_SHIFT        = 0
MAX_SHIFT        = 10

# Hard-coded no-shift zones within intron/exon regions
NO_SHIFT_INTRON  = 10
NO_SHIFT_EXON    = 3

# Keep user combinations somewhat sane.
MAX_COMBINATIONS = 100

def getShiftProfile(exonsize, intronsize, acceptor) :
    """
    Shift profile given by locations where no shifts are permitted.
    """
    if acceptor :
        noShiftStart = intronsize - NO_SHIFT_INTRON
        noShiftEnd   = intronsize + NO_SHIFT_EXON
    else :
        noShiftStart = exonsize - NO_SHIFT_EXON
        noShiftEnd   = exonsize + NO_SHIFT_INTRON
    return (noShiftStart, noShiftEnd)

def makeConfig(dimer, exonSize, intronSize, mink, maxk, maxShift, C, opts, ROC) :
    """
    Creates a new ClassifierConfig instance using the given parameters.
    """
    result = ClassifierConfig()
    result.setAcceptor(opts.acceptor)
    result.setDimer(dimer)
    result.setC(C)
    result.setExonSize(exonSize)
    result.setIntronSize(intronSize)
    result.setMink(mink)
    result.setMaxk(maxk)
    result.setMaxShift(maxShift)
    if opts.profile :
        (noShiftStart, noShiftEnd) = getShiftProfile(exonSize, intronSize, opts.acceptor)
        result.setNoShiftStart(noShiftStart)
        result.setNoShiftEnd(noShiftEnd)
    result.setMismatchProfile(opts.profile)
    result.setNormalize(not opts.nonnorm)
    result.setROC(ROC)
    siteType = 'acc' if opts.acceptor else 'don'
    result.setFastaPath('%s_%s.fa'%(dimer,siteType))
    result.setSVMPath('%s_%s.svm'%(dimer,siteType))
    result.addComment(timeString('Generated automatically'))
    return result

def makePositionalKmerData(trainingSet, exonsize, intronsize, **args) :
    """
    Returns a positional kmer kernel based on the given command-line options.
    """
    acceptor      = getAttribute('acceptor', False, **args)
    headerHandler = getAttribute('headerHandler', process_labeled_fasta_header, **args)
    maxk          = getAttribute('maxk', 1, **args)
    maxShift      = getAttribute('maxShift', 0, **args)
    mink          = getAttribute('mink', 1, **args)

    # Shift profile given by locations where no shifts are permitted
    (noShiftStart, noShiftEnd) = getShiftProfile(exonsize, intronsize, acceptor)

    # Mismatch profile constructed in positionalKmer method:
    return positionalKmer(trainingSet, exonsize, intronsize, noShiftStart=noShiftStart, noShiftEnd=noShiftEnd, **args)

def positionalKmer(trainingSet, exonSize, intronSize, **args) :
    """
    Creates a positional kmer with no mismatches allowed for the 13 positions flanking a splice site,
    in which 10 positions are on the intron side and 3 on the exon side.  If the intron is smaller than
    10 nt, the entire region allows no mismatches.  A similar rule applies to exons shorter than 3nt
    (though that hardly seems like a useful size).
    """
    acceptor      = getAttribute('acceptor', False, **args)
    headerHandler = getAttribute('headerHandler', process_labeled_fasta_header, **args)
    maxk          = getAttribute('maxk', 0, **args)
    mismatches    = getAttribute('mismatches', 1, **args)
    mink          = getAttribute('mink', 0, **args)

    # If no mismatches are allowed the mismatch profile is all 0's
    if mismatches == 0 :
        mismatchProfile = [0 for i in xrange(mink, maxk+1)]
    elif acceptor :
        zeroStart        = max(0, intronSize-10)
        mismatchProfile  = [mismatches for i in xrange(zeroStart)] + [0 for i in xrange(zeroStart,intronSize)]
        zeroEnd          = min(3, exonSize)
        mismatchProfile += [0 for i in xrange(zeroEnd)] + [mismatches for i in xrange(zeroEnd,exonSize)]
    else : # donor
        zeroStart        = max(0, exonSize-3)
        mismatchProfile  = [mismatches for i in xrange(zeroStart)] + [0 for i in xrange(zeroStart,exonSize)]
        zeroEnd          = min(10, intronSize)
        mismatchProfile += [0 for i in xrange(zeroEnd)] + [mismatches for i in xrange(zeroEnd,intronSize)]

    return SequenceData(trainingSet, mismatchProfile=mismatchProfile, **args)

USAGE="""%prog FASTA-train dimer [options]

Builds a weighted-degree kernel SVM for labeled FASTA sequences and runs cross-validation
on a variety of parameter settings.  Saves the best SVM, configuration and data for later use.

Note: some parameter settings may cause an SVM to fail to converge.  Often this may be corrected
by adding training data or by decreasing intron or exon size."""

#--------------------------------------------
# Command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-a', dest='acceptor',   default=False,      help='Treat sites as acceptors [default = donor sites]', action='store_true')
parser.add_option('-e', dest='exonsize',   default='8,12,16',  help='List of exon sizes to try [default = %default]')
parser.add_option('-C', dest='Clist',      default=C_DEF_STR,  help='List of regularization constants for SVM [default = %default]')
parser.add_option('-i', dest='intronsize', default='15,20,25', help='List of intron sizes to try [default = %default]')
parser.add_option('-l', dest='logfile',    default=None,       help='Optional logfile for tracking performance [default = dimer]')
parser.add_option('-k', dest='kfolds',     default=5,          help='Number of folds to run cross-validation [default = %default]', type='int')
parser.add_option('-m', dest='mink',       default='1',        help='Minimum k value for kernel [default = %default]')
parser.add_option('-M', dest='maxk',       default='1',        help='Maximum k value for kernel [default = %default]')
parser.add_option('-N', dest='nonnorm',    default=False,      help='Turn normalization OFF [default = %default]', action='store_true')
parser.add_option('-p', dest='profile',    default=False,      help='Use mismatch profile [default = %default]', action='store_true')
parser.add_option('-P', dest='prefix',     default=None,       help='Optional prefix for output files [default = dimer]')
parser.add_option('-S', dest='shift',      default='0',        help='List of maximum shift values [default = %default]')
parser.add_option('-v', dest='verbose',    default=False,      help='Verbose mode [default = %default]', action='store_true')

opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 2 :
    parser.print_help()
    sys.exit(1)

trainingFile  = args[0]
validateFile(trainingFile)
writeStartupMessage()

dimer = args[1]
if not opts.prefix :
    opts.prefix = '%s_acc'%dimer if opts.acceptor else '%s_don'%dimer

logStream = None
if opts.logfile :
    logStream = open(opts.logfile, 'w', 0)

#-----------------------------------------------------------
# Convert string lists to numbers
# Outer loop: intron/exon sizes
exonList   = [int(x) for x in opts.exonsize.split(',')]
intronList = [int(x) for x in opts.intronsize.split(',')]
minkList   = [int(x) for x in opts.mink.split(',')]
maxkList   = [int(x) for x in opts.maxk.split(',')]
shiftList  = [int(x) for x in opts.shift.split(',')]

#-----------------------------------------------------------
# Validate lists
if min(exonList)   < MIN_EXON_SIZE   : raise ValueError('Exon lengths must be at least %d' % MIN_EXON_SIZE)
if max(exonList)   > MAX_EXON_SIZE   : raise ValueError('Exon lengths cannot exceed %d' % MAX_EXON_SIZE)

if min(intronList) < MIN_INTRON_SIZE : raise ValueError('Intron lengths must be at least %dnt' % MIN_INTRON_SIZE)
if max(intronList) > MAX_INTRON_SIZE : raise ValueError('Intron lengths cannot exceed %dnt' % MAX_INTRON_SIZE)

if min(minkList)   < MIN_MINK        : raise ValueError('mink must be at least %d' % MIN_MINK)
if max(maxkList)   > MAX_MAXK        : raise ValueError('maxk cannot exceed %d' % MAX_MAXK)
if max(minkList)   > max(maxkList)   : raise ValueError('Maximum mink value cannot exceed the maximum maxk value')

if min(shiftList)  < MIN_SHIFT       : raise ValueError('Shift values must be at least %d' % MIN_SHIFT)
if max(shiftList)  > MAX_SHIFT       : raise ValueError('Shift values cannot exceed %d' % MAX_SHIFT)

#-----------------------------------------------------------
# Make sure user doesn't get lost in combinatorial explosion
combinations = len(exonList) * len(intronList) * len(minkList) * len(maxkList) * len(shiftList)
if combinations > MAX_COMBINATIONS :
    raise ValueError('Too many parameter combinations: %d.  Reduce values to get at most %d combinations.' % (combinations, MAX_COMBINATIONS))

#-----------------------------------------------------------
# Generator for all parameter combinations:
paramGenerator = ((e,i,lo,hi,s) for e in exonList for i in intronList for lo in minkList for hi in maxkList for s in shiftList if hi >= lo)

# C values handled by PyML modelSelector
Cvals = C_DEFAULTS
if opts.Clist :
    Cvals = [float(x) for x in opts.Clist.split(',')]

if opts.verbose :
    logMessage(timeString('Building weighted-degree kernel with the following parameters:\n'), logstream=logStream)
    logMessage('  dimer:              %s\n' % dimer, logstream=logStream)
    logMessage('  acceptor site:      %s\n' % ['False', 'True'][opts.acceptor], logstream=logStream)
    logMessage('  C:                  %s\n' % listString(Cvals), logstream=logStream)
    logMessage('  exon size range:    %s\n' % listString(exonList), logstream=logStream)
    logMessage('  intron size range:  %s\n' % listString(intronList), logstream=logStream)
    logMessage('  minimum k-mer size: %s\n' % listString(minkList), logstream=logStream)
    logMessage('  maximum k-mer size: %s\n' % listString(maxkList), logstream=logStream)
    logMessage('  maximum shift size: %s\n' % listString(shiftList), logstream=logStream)
    logMessage('  CV folds:           %d\n' % opts.kfolds, logstream=logStream)
    if opts.nonnorm   : logMessage('  Normalization turned OFF\n', logstream=logStream)
    if opts.profile   : logMessage('  using profiles for mismatches and shifts\n', logstream=logStream)
    logMessage('Comparing %d different parameter combinations.\n' % combinations, logstream=logStream)
    logMessage('Creating files %(out)s.svm, %(out)s.fa and %(out)s.cfg\n' % {'out':opts.prefix}, logstream=logStream)

#-----------------------------------------------------------
# Start main processing loops
bestSVM    = None
bestROC    = 0.0
bestConfig = None
bestData   = None

hideStdout()
for (exonSize,intronSize,mink,maxk,maxShift) in paramGenerator :
    if opts.verbose :
        logMessage(timeString("exon size=%d, intron size=%d, mink=%d, maxk=%d, shift=%d\n" % (exonSize, intronSize, mink, maxk, maxShift)), logstream=logStream)

    #--------------------------------------------------
    # extract training examples of the correct size
    tempTrain = "%s_tmp_e%d_i%d.fa" % (dimer, exonSize, intronSize)
    truncateSequences(trainingFile, exonSize, intronSize, tempTrain, acceptor=opts.acceptor)

    #--------------------------------------------------
    # Load sequences & create kernel
    if opts.profile :
        data = makePositionalKmerData(tempTrain, exonSize, intronSize,
                                      mink=mink, maxk=maxk, maxShift=maxShift, headerHandler=process_labeled_fasta_header)
    else :
        data = SequenceData(tempTrain, mink=mink, maxk=maxk, maxShift=maxShift, headerHandler=process_labeled_fasta_header)

    if not opts.nonnorm :
        data.attachKernel('cosine')

    #--------------------------------------------------
    # Use PyML framework to find optimal model:
    param    = modelSelection.Param(svm.SVM(), 'C', Cvals)
    selector = modelSelection.ModelSelector(param, numFolds=opts.kfolds, foldsToPerform=opts.kfolds, measure='roc')
    selector.train(data)

    if selector.log.maxSuccessRate > bestROC :
        if opts.verbose :
            logMessage(timeString("  ** New best ROC: %.4f (C=%.4f)\n" % (selector.log.maxSuccessRate, selector.classifier.C)), logstream=logStream)
        bestROC    = selector.log.maxSuccessRate
        bestSVM    = selector.classifier
        bestData   = data
        bestConfig = makeConfig(dimer, exonSize, intronSize, mink, maxk, maxShift, selector.classifier.C, opts, bestROC)
        bestSVM.save('%s.svm' % opts.prefix)
        os.rename(tempTrain, '%s.fa' % opts.prefix)
        bestConfig.save('%s.cfg' % opts.prefix)
    elif opts.verbose :
        logMessage(timeString("  ROC=%.4f at C=%.4f\n" % (selector.log.maxSuccessRate, selector.classifier.C)), logstream=logStream)

showStdout()
msgValues = (bestROC, bestConfig.exonSize(), bestConfig.intronSize(), bestConfig.mink(), bestConfig.maxk(), bestConfig.maxShift())
logMessage(timeString("Overall best SVM with ROC %.4f (exon=%d, intron=%d, mink=%d, maxk=%d, shift=%d)\n" % msgValues), logstream=logStream)
