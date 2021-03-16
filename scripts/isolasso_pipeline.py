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
Script that runs the IsoLasso software (Li et al., 2011)
using SpliceGrapher predictions as a reference.
"""
from SpliceGrapher.shared.config import *
from SpliceGrapher.shared.utils  import *
from SpliceGrapher.formats.fasta import *

from optparse import OptionParser
from glob import glob
import os,sys,subprocess

# Tricky bit for users:
REQUIRED_PROGS   = ['gtfToGenePred', 'genePredToBed']

UCSC_PROG_MESSAGE = """
This pipeline requires the following programs:
    %s

They should be available from the UCSC website:

    http://hgdownload.cse.ucsc.edu/admin/exe

Alternately, go to the UCSC Genome Bioinformatics website (genome.ucsc.edu)
and click on the 'Downloads' link.  Scroll down to the 'Utilities Downloads'
section and look for the 'other utilities' link.  Select your platform and
look for these tools there.

""" % '\n    '.join(REQUIRED_PROGS)

# Main steps in the pipeline:
# -A -- all genes but no unresolved nodes
PUTATIVE_GTF_CMD        = 'generate_putative_sequences.py %s -A -f %s -M %s -m %s'
# -AU -- include unresolved
PUTATIVE_UNRESOLVED_CMD = 'generate_putative_sequences.py %s -AU -f %s -M %s -m %s'

GTF_CONVERT      = 'gtfToGenePred %s stdout | genePredToBed | sort -k1,1 -k2,2n > %s'
LASSO_CMD        = 'runlasso.py -x %s --forceref %s -o %s %s'
UPDATE_CMD       = 'isolasso_update_graphs.py %s %s %s -t %.2f -d %s'

def attributeDict(line) :
    """Parses a GTF line and returns the attributes as a dictionary."""
    tokens     = line.strip().split('\t')
    attrToken  = tokens[-1]
    # e.g., ['gene_id "Inst1"', 'transcript_id "AT1G01100_p14"', 'FPKM  "0"', 'frac "1.000000"', 'conf_lo "0.0"', 'conf_hi "2.0"', 'cov "0.1"', '']
    attributes = attrToken.split(';')
    result     = {}
    for pair in attributes :
        # e.g., 'transcript_id "AT1G01100_p14"'
        pair = pair.replace('"','')
        if not pair : continue
        key,val = pair.split()
        result[key] = val
    return result

def checkIsoLassoValues(isoFile) :
    """Does a quick validation of an IsoLasso output GTF file."""
    fpkmVals = set()
    for line in ezopen(isoFile) :
        attrs = attributeDict(line)
        fpkmVals.add(float(attrs['FPKM']))
    if not fpkmVals :
        raise Exception('No FPKM values found in %s' % isoFile)

    return 'range of FPKM values was %.3f to %.3f\n' % (min(fpkmVals), max(fpkmVals))

USAGE = """%prog graph-files SAM-file [options]

Where:
    graph-files    is either a file of paths to splice graph files
                   or a directory of splice graphs with the structure:
                   top-level/chromosome-dir/graph-file
    SAM-file       is a SAM file to be used with IsoLasso

Example:
    %prog my_predictions filtered.sam

Runs the IsoLasso software (Li et al., 2011) using SpliceGrapher predictions as a reference."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-C', dest='cem',       default=False,  help='Use CEM instead of LASSO [default: %default]', action='store_true')
parser.add_option('-d', dest='outdir',    default='isolasso', help='Output directory [default: %default]')
parser.add_option('-p', dest='params',    default='',     help='Additional parameters for IsoLasso [default: %default]')
parser.add_option('-f', dest='fasta',     default=SG_FASTA_REF, help='FASTA genome reference [default: %default]')
parser.add_option('-t', dest='threshold', default=1.0,    help='Minimum FPKM threshold [default: %default]', type='float')
parser.add_option('-U', dest='unresolved',default=False,  help='Include unresolved transcripts [default: %default]', action='store_true')
parser.add_option('-v', dest='verbose',   default=False,  help='Verbose mode [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

#-------------------------------------------------------------------------------------------------
# Parse the command line and make sure everything looks OK
MIN_ARGS = 2
if len(args) != MIN_ARGS :
    parser.print_help()
    if args : sys.stderr.write('\nExpected %d parameters; received %d:\n  %s\n' % (MIN_ARGS, len(args), '\n  '.join(args)))
    sys.exit(1)

# Make sure the user has downloaded the required UCSC programs:
for p in REQUIRED_PROGS :
    fullPath = findFile(p, os.environ['PATH'])
    if not fullPath :
        sys.stderr.write("\nUnable to find the program '%s' in your path.\n" % p)
        sys.stderr.write(UCSC_PROG_MESSAGE)
        sys.exit(1)

# Check required files
for f in args : validateFile(f)

if not opts.fasta :
    parser.print_help()
    sys.stderr.write('\nYou must specify a genome reference sequence (-f option).\n')
    sys.exit(1)

validateFile(opts.fasta)
fileDB     = args[0]
samFile    = args[1]

if opts.cem :
    opts.params += '--useem'

PUTATIVE_CMD = PUTATIVE_UNRESOLVED_CMD if opts.unresolved else PUTATIVE_GTF_CMD

base       = os.path.basename(fileDB)
filePrefix,ext = os.path.splitext(base)
graphList  = makeGraphListFile(fileDB) if os.path.isdir(fileDB) else fileDB

# Use directory name to make (somewhat) unique, meaningful file names
putativeBed = '%s_forms.bed' % filePrefix
putativeMap = '%s_forms.map' % filePrefix
putativeGTF = '%s_forms.gtf' % filePrefix

writeStartupMessage()

#-------------------------------------------------------------------------------------------------
# 1. generate_putative_sequences.py initial_predictions.lis -f Genome_reference.fa -A --gtf -M put.gtf -m put.map -o put.fa 
cmd = PUTATIVE_CMD % (graphList, opts.fasta, putativeGTF, putativeMap)
runCommand(cmd, stdout=None, stderr=None)
if not os.path.exists(putativeGTF) : raise Exception('%s was not created' % putativeGTF)

#-------------------------------------------------------------------------------------------------
#  2. gtf2bed putative_forms.gtf

# Tricky bit for users; need to integrate this more seamlessly:
cmd   = GTF_CONVERT % (putativeGTF, putativeBed)
##cmd = GTF_TO_BED % putativeGTF
runCommand(cmd, stdout=None, stderr=None)
if not os.path.exists(putativeBed) : raise Exception('%s was not created' % putativeBed)

#-------------------------------------------------------------------------------------------------
#  3. runlasso.py -x putative_forms.bed --forceref -o isolasso filtered.sam'
cmd = LASSO_CMD % (putativeBed, opts.params, opts.outdir, samFile)
isoLassoGTF = '%s.pred.gtf' % opts.outdir
runCommand(cmd, stdout=None, stderr=None)
if not os.path.exists(isoLassoGTF) : raise Exception('%s was not created' % isoLassoGTF)

#-------------------------------------------------------------------------------------------------
#  4. isolasso_update_graphs.py initial_predictions.lis putative_transcripts.map isolasso.pred.gtf -v -t 400.0 -d isolasso_pred
cmd = UPDATE_CMD % (graphList, putativeMap, isoLassoGTF, opts.threshold, opts.outdir)
runCommand(cmd, stdout=None, stderr=None)
if not os.path.isdir(opts.outdir) :
    # (if we got this far, the file must exist)
    message = checkIsoLassoValues(isoLassoGTF)
    raise Exception('%s was not created\n%s\n' % (opts.outdir, message))

logMessage('\nFinished.\n')
