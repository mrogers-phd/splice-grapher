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
Script that runs the PSGInfer software (Dewey et al., 2010)
using SpliceGrapher predictions as a reference.
"""
from SpliceGrapher.shared.config import *
from SpliceGrapher.shared.utils  import *
from SpliceGrapher.formats.fasta import *

from optparse import OptionParser
from glob import glob
import os,sys,subprocess

# Determined via simulations
DEFAULT_THRESH   = 0.01

MAKE_FASTA_CMD   = 'cat %s/*.fa > %s'
# Only genes with unresolved/novel paths:
PUTATIVE_GTF_CMD = 'generate_putative_sequences.py %s -f %s -M %s -m %s -UA'
UPDATE_CMD       = 'psginfer_update_graphs.py %s %s %s -vd %s -t %s'

# PSGInfer 1.1.3 software:
PREPARE_LOG      = 'psg_prepare_reference.log'
INFER_LOG        = 'psg_infer_frequencies.log'
PREPARE_CMD      = 'psg_prepare_reference.py -g %s %s < %s'
INFER_CMD        = 'psg_infer_frequencies.py %s %s %s %s'

PSGINFER_DIR     = 'psginfer_data'

FASTA_DIR        = 'FASTA_REF'
OUT_SUFFIX       = 'psg_reference'
RESULTS_SUFFIX   = 'psg_results'

PUTATIVE_BED     = 'putative_forms.bed'
PUTATIVE_MAP     = 'putative_forms.map'
PUTATIVE_GTF     = 'putative_forms.gtf'
REFERENCE_DIR    = 'putative_forms_psg_reference'
RESULTS_DIR      = 'putative_forms_psg_results'
PRED_DIR         = 'psgpred'
ISOFORM_RESULTS  = 'isoform.results.txt'

USAGE = """%prog graph-list fastq1 fastq2 [options]

Example:
    %prog initial_predictions.lis reads_1.fq reads_2.fq

Runs the PSGInfer software (Dewey et al., 2010) using SpliceGrapher
predictions as a reference."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-d', dest='outdir',  default='.',          help='Output directory [default: %default]')
parser.add_option('-D', dest='debug',   default=False,        help='Debug mode (just show commands) [default: %default]', action='store_true')
parser.add_option('-f', dest='fasta',   default=SG_FASTA_REF, help='Genome reference FASTA [default: %default]')
parser.add_option('-l', dest='logfile', default=None,         help='Optional log file [default: %default]')
parser.add_option('-t', dest='thresh',  default=DEFAULT_THRESH, help='Probability threshold [default: %default]', type='float')
parser.add_option('-v', dest='verbose', default=False,        help='Verbose mode [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

#-------------------------------------------------------------------------------------------------
# Parse the command line and make sure everything looks OK
MIN_ARGS = 3
if len(args) != MIN_ARGS :
    parser.print_help()
    if args : sys.stderr.write('\nExpected %d parameters; received %d:\n  %s\n' % (MIN_ARGS, len(args), '\n  '.join(args)))
    sys.exit(1)

# Check required files
samFile = args[0]
for f in args : validateFile(f)

writeStartupMessage()

logStream = open(opts.logfile,'w') if opts.logfile else None
# Use absolute paths 
graphList = args[0]
fq1File   = os.path.abspath(args[1])
fq2File   = os.path.abspath(args[2])

# Since we're calling another package, set up log files:
psgPrepLog  = open(PREPARE_LOG,'w')
psgInferLog = open(INFER_LOG,'w')

# PSGInfer does its work in the same directory as the separate FASTA files:
workDir = os.path.abspath(FASTA_DIR)
if not os.path.exists(workDir) :
    os.makedirs(workDir)

#-------------------------------------------------------------------------------------------------
# 1. split FASTA reference into FASTA directory:
for rec in fasta_itr(opts.fasta) :
    prefix  = rec.header.split()[0]
    outFile = os.path.join(workDir, '%s.fa'%prefix)
    if os.path.exists(outFile) : continue
    if opts.verbose : sys.stderr.write('  creating %s\n' % outFile)
    outStream = open(outFile,'w')
    outStream.write(str(rec))
    outStream.close()

#-------------------------------------------------------------------------------------------------
# 2. generate_putative_sequences.py initial_predictions.lis -f Genome_reference.fa --gtf -M put.gtf -m put.map -o put.fa 
cmd = PUTATIVE_GTF_CMD % (graphList, opts.fasta, PUTATIVE_GTF, PUTATIVE_MAP)
runCommand(cmd, logstream=logStream, debug=opts.debug, stdout=logStream, stderr=logStream)
gtfFile = os.path.abspath(PUTATIVE_GTF)
if not (opts.debug or os.path.exists(gtfFile)) : raise Exception('%s was not created' % gtfFile)

#-------------------------------------------------------------------------------------------------
#  3. psg_prepare_reference.py -g FASTA_REF psginfer_data < putative_forms.gtf'
cmd = PREPARE_CMD % (FASTA_DIR, PSGINFER_DIR, gtfFile)
runCommand(cmd, logstream=logStream, debug=opts.debug, stdout=psgPrepLog, stderr=psgPrepLog)
if not (opts.debug or os.path.isdir(PSGINFER_DIR)) :
    raise Exception('Error running psg_prepare_reference.py; check %s for details' % PREPARE_LOG)

#-------------------------------------------------------------------------------------------------
#  4. psg_infer_frequencies.py reference_name sample_name upstream_read_file downstream_read_file
cmd = INFER_CMD % (PSGINFER_DIR, PSGINFER_DIR, fq1File, fq2File)
runCommand(cmd, logstream=logStream, debug=opts.debug, stdout=psgInferLog, stderr=psgInferLog)
isoformFile = os.path.join(PSGINFER_DIR, ISOFORM_RESULTS)
if not (opts.debug or os.path.exists(isoformFile)) :
    raise Exception('Error running psg_infer_frequencies.py; check %s for details' % INFER_LOG)

#-------------------------------------------------------------------------------------------------
#  5. psginfer_update_graphs.py initial_predictions.lis putative_forms.map FASTA_REF/putative_forms_psg_results/isoform.results.txt -vd psgpred
cmd = UPDATE_CMD % (graphList, PUTATIVE_MAP, isoformFile, opts.outdir, str(opts.thresh))
runCommand(cmd, logstream=logStream, debug=opts.debug, stdout=None, stderr=None)
if not (opts.debug or os.path.isdir(opts.outdir)) : raise Exception('%s was not created' % opts.outdir)

logMessage('\nFinished.\n', logstream=logStream)
