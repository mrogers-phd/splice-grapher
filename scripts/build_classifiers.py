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
Script that demonstrates the entire process of creating splice
site models given a reference genome and its corresponding
gene models.
"""
from SpliceGrapher.shared.config import *
from SpliceGrapher.shared.utils  import *
from SpliceGrapher.formats.fasta import *
from optparse                    import OptionParser
import os,sys,subprocess

C_DEFAULTS      = [10.0**x for x in range(-2,3)]
C_DEF_STR       = ','.join(['%.2f'%x for x in C_DEFAULTS])
EXON_DEFAULT    = '8,12,16'
INTRON_DEFAULT  = '15,20,25'
MAXK_DEFAULT    = '1'
SHIFT_DEFAULT   = '0'

LEGAL_DIMERS    = ['%s%s'%(a,b) for a in 'acgt' for b in 'acgt']
DON_STRINGS     = 'gt'
ACC_STRINGS     = 'ag'
TRAINING_MIN    = 10
LOGFILE_FORMAT  = '%s_training.log'
SELECT_COMMAND  = 'select_model_parameters.py %s %s -l %s -P %s %s %s'
REPORT_FILE     = 'splice_site_report.txt'
ZIP_COMMAND     = 'zip classifiers.zip ??_???.cfg ??_???.svm ??_???.fa'
CLEAN_COMMAND   = 'rm ??_tmp_e*_i*.fa'

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
    accFlag = ' -a' if acceptor else ''
    suffix  = 'acc' if acceptor else 'don'
    lcDimer = dimer.lower()
    posFile = '%s_%s_training.fa' % (lcDimer, suffix)

    # 2 subsets, 2 lines per sequence:
    expectedSize = 4*subsetSize 
    fasta_part   = '' if opts.fasta == SG_FASTA_REF else ' -f %s' % opts.fasta
    model_part   = '' if opts.model == SG_GENE_MODEL else ' -m %s' % opts.model

    if not fileExists(posFile, opts, alert=False) or fileLen(posFile) != expectedSize :
        rptOption = '' if os.path.exists(REPORT_FILE) else '-r %s' % REPORT_FILE
        command   = 'generate_splice_site_data.py%s%s%s -d %s -n %d -o %s %s' % \
                (accFlag, fasta_part, model_part, dimer, subsetSize, posFile, rptOption)
        runCommand(command, logstream=logstream, debug=opts.commands)

        negFile = '%s_%s_neg.fa' % (lcDimer, suffix)
        command = 'generate_splice_site_data.py%s%s%s -d %s -n %d -o %s -N' % \
                (accFlag, fasta_part, model_part, dimer, subsetSize, negFile)
        runCommand(command, logstream=logstream, debug=opts.commands)

        command = 'cat %s >> %s' % (negFile, posFile)
        runCommand(command, logstream=logstream, debug=opts.commands)

        # Use command instead of os.remove() so user can see what's happening
        command = 'rm %s' % negFile
        runCommand(command, logstream=logstream, debug=opts.commands)
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
        wrappedOptions = ''
        if opts.Clist != C_DEF_STR : wrappedOptions += ' -C %s' % opts.Clist
        if opts.exonsize != EXON_DEFAULT : wrappedOptions += ' -e %s' % opts.exonsize
        if opts.intronsize != INTRON_DEFAULT : wrappedOptions += ' -i %s' % opts.intronsize
        if opts.maxk != MAXK_DEFAULT : wrappedOptions += ' -M %s' % opts.maxk
        if opts.shift != SHIFT_DEFAULT : wrappedOptions += ' -S %s' % opts.shift
        if opts.profile : wrappedOptions += ' -p'
        command = SELECT_COMMAND % (trainingFile, dimer, logFile, fullDimerName, accFlag, wrappedOptions)
        runCommand(command, logstream=logstream, debug=opts.commands)
    return dimerFile

USAGE = """%prog [options]

Creates an organism's splice-site classifiers by performing the following steps:
   1. Generates training data for each splice site dimer
   2. Selects optimal parameters for each splice-site-dimer classifier
   3. Creates a .zip file that contains the resulting classifiers"""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-a', dest='acceptors',  default=ACC_STRINGS,   help='Acceptor site dimers to predict [default: %default]')
parser.add_option('-d', dest='donors',     default=DON_STRINGS,   help='Donor site dimers to predict [default: %default]')
parser.add_option('-f', dest='fasta',      default=SG_FASTA_REF,  help='FASTA genomic reference [default: %default]')
parser.add_option('-m', dest='model',      default=SG_GENE_MODEL, help='Gene model annotations (GFF3/GTF) [default: %default]')
parser.add_option('-n', dest='examples',   default=2000,          help='Number of examples in training set [default: %default]', type='int')
parser.add_option('-l', dest='logfile',    default=None,          help='Optional log file [default: %default]')
parser.add_option('-v', dest='verbose',    default=False,         help='Verbose mode [default: %default]', action='store_true')
parser.add_option('--commands', dest='commands', default=False,   help="Show but don't run commands [default: %default]", action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

# Omitted to simplify interface
#parser.add_option('-C', dest='Clist',      default=C_DEF_STR,     help='List of regularization constants for SVM [default = %default]')
#parser.add_option('-e', dest='exonsize',   default=EXON_DEFAULT,  help='List of exon sizes to try [default = %default]')
#parser.add_option('-i', dest='intronsize', default=INTRON_DEFAULT, help='List of intron sizes to try [default = %default]')
#parser.add_option('-M', dest='maxk',       default=MAXK_DEFAULT,  help='Maximum k value for kernel [default = %default]')
#parser.add_option('-O', dest='overwrite',  default=False,         help='Over-write any existing files [default: %default]', action='store_true')
#parser.add_option('-p', dest='profile',    default=False,         help='Use mismatch profile [default = %default]', action='store_true')
#parser.add_option('-S', dest='shift',      default=SHIFT_DEFAULT, help='List of maximum shift values [default = %default]')
#

# If options are commented back in above, comment them out here:
opts.Clist      = C_DEF_STR
opts.exonsize   = EXON_DEFAULT
opts.intronsize = INTRON_DEFAULT
opts.maxk       = MAXK_DEFAULT
#opts.examples   = 2000
opts.overwrite  = True
opts.profile    = False
opts.shift      = SHIFT_DEFAULT
#

#-------------------------------------------------------------------------------------------------
# Parse the command line and make sure everything looks OK
errStrings = []
if not opts.model : errStrings.append('** No gene models specified.  Use the -m option or set SG_GENE_MODEL in your SpliceGrapher configuration.')
if not opts.fasta : errStrings.append('** No FASTA reference specified.  Use the -f option or set SG_FASTA_REF in your SpliceGrapher configuration.')
if errStrings :
    parser.print_help()
    sys.stderr.write('\n%s\n' % '\n'.join(errStrings))
    sys.exit(1)

validateFile(opts.model)
validateFile(opts.fasta)
writeStartupMessage()

if opts.examples < TRAINING_MIN :
    raise ValueError('Training set size (%d) must be at least %d' % (opts.examples, TRAINING_MIN))

#if opts.overlap >= opts.readlen :
#    raise ValueError('Overlap size (%d) must be less than read length (%d)' % (opts.overlap, opts.readlen))

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

# Determine pos/neg training subset, sequence window sizes:
subsetSize = int(0.5 * opts.examples)

#-------------------------------------------------------------------------------------------------
# 1. Generate training data for each splice site dimer
# 2. Find optimal parameters for an SVM based on these data
logMessage('\nCreating SVM models for splice sites\n', logstream=logStream)

if opts.verbose : logMessage('   Donor dimers (%s)\n' % ','.join(DONOR_DIMERS), logstream=logStream)
dimerNames = []
for dimer in DONOR_DIMERS :
    if opts.verbose : logMessage('     generating %s dimer training data (%d examples)\n' % (dimer, opts.examples), logstream=logStream)
    trainingFile  = generateTrainingData(dimer, opts, subsetSize, acceptor=False, logstream=logStream)
    fullDimerName = selectModelParameters(dimer, trainingFile, opts)
    dimerNames.append(fullDimerName)

if opts.verbose : logMessage('   Acceptor dimers (%s)\n' % ','.join(ACCEPTOR_DIMERS), logstream=logStream)
for dimer in ACCEPTOR_DIMERS :
    if opts.verbose : logMessage('     generating %s dimer training data (%d examples)\n' % (dimer, opts.examples), logstream=logStream)
    trainingFile  = generateTrainingData(dimer, opts, subsetSize, acceptor=True, logstream=logStream)
    fullDimerName = selectModelParameters(dimer, trainingFile, opts, acceptor=True)
    dimerNames.append(fullDimerName)

# 3. Store classifiers in a zip file for later use
if opts.verbose : logMessage('Create .zip file\n')
runCommand(ZIP_COMMAND, logstream=logStream, debug=opts.commands, stderr=nullStream, stdout=nullStream)
if opts.verbose : logMessage('Removing temporary files\n')
runCommand(CLEAN_COMMAND, logstream=logStream, debug=opts.commands, stderr=nullStream, stdout=nullStream)

logMessage('\nFinished.\n', logstream=logStream)

