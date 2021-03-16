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
Script that demonstrates the entire process of generating splice-junction
sequences prior to running alignments.
"""
from SpliceGrapher.shared.config import *
from SpliceGrapher.shared.utils  import *
from SpliceGrapher.formats.fasta import *
from optparse                    import OptionParser
import os,sys,subprocess

LEGAL_DIMERS    = ['%s%s'%(a,b) for a in 'acgt' for b in 'acgt']
DON_STRINGS     = 'gt'
ACC_STRINGS     = 'ag'
TRAINING_MIN    = 10
LOGFILE_FORMAT  = '%s_training.log'
SELECT_COMMAND  = 'select_model_parameters.py %s %s -l %s -P %s %s'
ROC_COMMAND     = 'generate_roc.py %s -o %s'
REPORT_FILE     = 'splice_site_report.txt'

def fileExists(path, opts, alert=True) :
    """Convenience function that loosely infers whether a file has been created
    and writes a message if desired."""
    if opts.overwrite : return False
    exists = os.path.exists(path) and os.path.getsize(path) > 0
    if exists and alert :
        sys.stderr.write(' + File %s already exists.\n' % path)
    return exists

def fileLen(path):
    """Simple function to get an exact file length."""
    with open(path) as f:
        for i, l in enumerate(f): pass
    return i+1

def generateTrainingData(dimer, opts, subsetSize, acceptor=False, logstream=None) :
    """Generates data necessary for training an SVM for one splice-site dimer."""
    accFlag = '-a' if acceptor else ''
    suffix  = 'acc' if acceptor else 'don'
    lcDimer = dimer.lower()
    posFile = '%s_%s_training.fa' % (lcDimer, suffix)

    # 2 subsets, 2 lines per sequence:
    expectedSize = 4*subsetSize 

    if not fileExists(posFile, opts, alert=False) or fileLen(posFile) != expectedSize :
        rptOption = '' if os.path.exists(REPORT_FILE) else '-r %s' % REPORT_FILE
        command   = 'generate_splice_site_data.py %s -f %s -m %s -d %s -n %d -o %s %s' % \
                (accFlag, opts.fasta, opts.model, dimer, subsetSize, posFile, rptOption)
        runCommand(command, logstream=logstream, debug=opts.debug)

        negFile = '%s_%s_neg.fa' % (lcDimer, suffix)
        command = 'generate_splice_site_data.py %s -f %s -m %s -d %s -n %d -o %s -N' % \
                (accFlag, opts.fasta, opts.model, dimer, subsetSize, negFile)
        runCommand(command, logstream=logstream, debug=opts.debug)

        command = 'cat %s >> %s' % (negFile, posFile)
        runCommand(command, logstream=logstream, debug=opts.debug)
    else :
        sys.stderr.write(' + File %s already exists.\n' % posFile)

    return posFile

def selectModelParameters(dimer, trainingFile, opts, acceptor=False, logstream=None) :
    accFlag       = '-a'  if acceptor else ''
    suffix        = 'acc' if acceptor else 'don'
    logFile       = LOGFILE_FORMAT % dimer
    fullDimerName = '%s_%s' % (dimer, suffix)
    dimerFile     = fullDimerName + '.cfg'
    if not fileExists(dimerFile, opts) :
        command = SELECT_COMMAND % (trainingFile, dimer, logFile, fullDimerName, accFlag)
        runCommand(command, logstream=logstream, debug=opts.debug, stderr=nullStream, stdout=nullStream)
    return dimerFile

USAGE = """%prog [options]

Script for running the full splice junction prediction pipeline.
This will perform the following steps:
    1. Generate known and recombined junction sequences
    2. Create a database of predicted splice sites
       a. Generate training data for each splice site dimer
       b. Find optimal parameters for each splice-site dimer SVM
       c. Apply SVMs to dimers in the reference genome to create the database
       d. Use predicted sites to generate predicted junction sequences (--gendb option)"""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-a', dest='acceptors',  default=ACC_STRINGS,   help='Acceptor site dimers to predict [default: %default]')
parser.add_option('-D', dest='debug',      default=False,         help="Show but don't run commands [default: %default]", action='store_true')
parser.add_option('-d', dest='donors',     default=DON_STRINGS,   help='Donor site dimers to predict [default: %default]')
parser.add_option('-f', dest='fasta',      default=SG_FASTA_REF,  help='Fasta reference file [default: %default]')
parser.add_option('-m', dest='model',      default=SG_GENE_MODEL, help='Gene model annotations (GFF3) [default: %default]')
parser.add_option('-l', dest='logfile',    default=None,          help='Optional log file [default: %default]')
parser.add_option('-o', dest='overlap',    default=10,            help='Required overlap on either side of a junction [default: %default]', type='int')
parser.add_option('-O', dest='overwrite',  default=False,         help='Over-write any existing files [default: %default]', action='store_true')
parser.add_option('-r', dest='readlen',    default=36,            help='Read length [default: %default]', type='int')
parser.add_option('-t', dest='training',   default=1000,          help='Training data set size [default: %default]', type='int')
parser.add_option('-v', dest='verbose',    default=False,         help='Verbose mode [default: %default]', action='store_true')
parser.add_option('--gendb', dest='gendb', default=False,         help="Use predicted splice sites to generate splice junction sequence database (can take a long time) [default: %default]", action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

#-------------------------------------------------------------------------------------------------
# Parse the command line and make sure everything looks OK
errStrings = []
if not opts.model : errStrings.append('** No GFF gene model specified.  Use the -m option or set SG_GENE_MODEL in your SpliceGrapher configuration.')
if not opts.fasta : errStrings.append('** No FASTA reference specified.  Use the -f option or set SG_FASTA_REF in your SpliceGrapher configuration.')
if errStrings :
    parser.print_help()
    sys.stderr.write('\n%s\n' % '\n'.join(errStrings))
    sys.exit(1)

validateFile(opts.model)
validateFile(opts.fasta)
writeStartupMessage()

if opts.overlap < 0 :
    raise ValueError('Overlap size (%d) must be non-negative' % opts.overlap)

if opts.readlen <= 0 :
    raise ValueError('Read length (%d) must be greater than zero' % opts.readlen)

if opts.training < TRAINING_MIN :
    raise ValueError('Training set size (%d) must be at least %d' % (opts.training, TRAINING_MIN))

if opts.overlap >= opts.readlen :
    raise ValueError('Overlap size (%d) must be less than read length (%d)' % (opts.overlap, opts.readlen))

DONOR_DIMERS = opts.donors.split(',')
for dimer in DONOR_DIMERS :
    if len(dimer) != 2 : raise ValueError("Illegal donor dimer '%s': dimers must have exactly 2 chars" % dimer)
    if dimer.lower() not in LEGAL_DIMERS :
        sys.stderr.write("** Warning: donor dimer '%s' contains an unexpected character\n" % dimer)

ACCEPTOR_DIMERS = opts.acceptors.split(',')
for dimer in ACCEPTOR_DIMERS :
    if len(dimer) != 2 : raise ValueError("Illegal acceptor dimer '%s': dimers must have exactly 2 chars" % dimer)
    if dimer.lower() not in LEGAL_DIMERS :
        sys.stderr.write("** Warning: acceptor dimer '%s' contains an unexpected character\n" % dimer)

# Establish output streams for different purposes:
nullStream = open(os.devnull, 'a')
logStream  = open(opts.logfile, 'w') if opts.logfile else nullStream
commandLog = logStream if opts.verbose else nullStream

# Determine pos/neg training subset, sequence window sizes:
subsetSize = int(0.5 * opts.training)
window     = opts.readlen-opts.overlap

#-------------------------------------------------------------------------------------------------
# Start the process:
#    1. Generate known and recombined junction sequences
outFile = 'known_jct_%d.fa' % opts.readlen
if not fileExists(outFile, opts) :
    logMessage('\nGenerating known, recombined junction sequences\n', logstream=logStream)
    command    = 'generate_known_junctions.py -f %s -m %s -w %d -o %s' % (opts.fasta, opts.model, window, outFile)
    commandLog = logStream if opts.verbose else nullStream
    runCommand(command, logstream=commandLog, debug=opts.debug)

#-------------------------------------------------------------------------------------------------
#    2. Generate predicted junction sequences
#       a. Generate training data for each splice site dimer
#       b. Find optimal parameters for an SVM based on these data
logMessage('\nCreating SVM models for splice sites\n', logstream=logStream)

if opts.verbose : logMessage('   Donor dimers (%s)\n' % ','.join(DONOR_DIMERS), logstream=logStream)
dimerNames = []
for dimer in DONOR_DIMERS :
    if opts.verbose : logMessage('     generating %s dimer training data (%d examples)\n' % (dimer, opts.training), logstream=logStream)
    trainingFile  = generateTrainingData(dimer, opts, subsetSize, acceptor=False, logstream=commandLog)
    fullDimerName = selectModelParameters(dimer, trainingFile, opts)
    dimerNames.append(fullDimerName)

if opts.verbose : logMessage('   Acceptor dimers (%s)\n' % ','.join(ACCEPTOR_DIMERS), logstream=logStream)
for dimer in ACCEPTOR_DIMERS :
    if opts.verbose : logMessage('     generating %s dimer training data (%d examples)\n' % (dimer, opts.training), logstream=logStream)
    trainingFile  = generateTrainingData(dimer, opts, subsetSize, acceptor=True, logstream=commandLog)
    fullDimerName = selectModelParameters(dimer, trainingFile, opts, acceptor=True)
    dimerNames.append(fullDimerName)

#-------------------------------------------------------------------------------------------------
# c. Apply SVMs to all sequences in the main FASTA file
#   i. split out sequences from FASTA file
logMessage('\nUsing SVM classifiers to identify predicted sites.\n', logstream=logStream)

fileDict = {}
if not opts.debug :
    if opts.verbose : logMessage('  Creating separate FASTA file for each chromosome\n', logstream=logStream)
    for rec in fasta_itr(opts.fasta) :
        chrom    = process_fasta_header(rec.header)
        filePath = '%s.fa' % chrom.lower()
        fileDict[chrom] = filePath
        if fileExists(filePath, opts) : continue

        if opts.verbose : logMessage('    creating %s.fa\n' % chrom, logstream=logStream)
        fileStream = open(filePath, 'w')
        fileStream.write(str(rec))

#   ii. run classification on each file
if opts.verbose : logMessage('\n  Applying classifiers to each chromosome\n')
dimerString = ','.join(dimerNames)
chromList   = sorted(fileDict.keys())
predictions = {}
for chrom in chromList :
    chromFile = '%s_sites.dat' % chrom
    predictions[chrom] = chromFile
    if fileExists(chromFile, opts) : continue
    if opts.verbose : logMessage('    classifying sites in %s' % chrom, logstream=logStream)
    command = 'classify_sites.py %s -o %s -c %s -f %s -m %s -v' % (chrom, chromFile, dimerString, fileDict[chrom], opts.model)
    runCommand(command, logstream=logStream, debug=opts.debug, stderr=nullStream, stdout=nullStream)

#  d. Optional: use predicted sites to generate splice junction sequence database
if opts.gendb :
    logMessage('\nUsing predicted sites to construct a splice-junction sequence database.\n', logstream=logStream)
    for chrom in chromList :
        predictedFasta = '%s_pred_%d.fa' % (chrom, opts.readlen)
        if fileExists(predictedFasta, opts) : continue
        if opts.verbose : logMessage('  generating sequences for %s' % chrom, logstream=logStream)
        logStream = open('%s_pred.out' % chrom, 'w')
        errStream = open('%s_pred.err' % chrom, 'w')
        command   = 'generate_predicted_junctions.py %s -f %s -m %s -o %s -W %d -d %s -a %s' % (predictions[chrom], fileDict[chrom], opts.model, predictedFasta, window, opts.donors, opts.acceptors)
        runCommand(command, logstream=logStream, debug=opts.debug, stderr=nullStream, stdout=nullStream)

logMessage('\nFinished.\n', logstream=logStream)
