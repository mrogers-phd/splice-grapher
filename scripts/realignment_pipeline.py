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
Script that handles the process of predicting and updating
splice graphs by realigning reads to a database of putative transcripts.
"""
from SpliceGrapher.shared.config import *
from SpliceGrapher.shared.utils  import *
from SpliceGrapher.formats.fasta import *

from optparse import OptionParser
import os,sys,subprocess,multiprocessing

GENE_MODEL_ERR = '** No gene models specified.  Use the -m option or set SG_GENE_MODEL in your SpliceGrapher configuration.'
FASTA_FILE_ERR = '** No FASTA reference specified.  Use the -f option or set SG_FASTA_REF in your SpliceGrapher configuration.'

PUTATIVE_FASTA    = 'putative.fa'
PREDICTION_SCRIPT = 'predict_graphs.py'
TRANSCRIPT_SCRIPT = 'generate_putative_sequences.py'
TRANSCRIPT_MAP    = 'putative_transcripts.map'
BWA_SAM           = 'filtered_bwa.sam'

# File extensions used to validate BWA alignments 
BWA_EXTS  = ['amb', 'ann', 'bwt', 'pac', 'sa']

# Novel transcripts denoted [gene-id]_u[#]
NOVEL_TRANS_ID    = '_u'

def fileExists(path) :
    """Convenience function that loosely infers whether a file has been created."""
    return os.path.exists(path) and os.path.getsize(path) > 0

def writeNovelAlignments(inputSam, outputSam, initiate=False) :
    """Finds novel alignments in the input SAM file and appends them
    to the output file."""
    outStream = open(outputSam,'w') if initiate else open(outputSam,'a')
    for line in ezopen(inputSam) :
        if line.startswith('@') : continue
        parts = line.split('\t')
        if parts[2].find(NOVEL_TRANS_ID) > 0 :
            outStream.write(line)

def getBWAAlnFileName(fqFile) :
    """Given a path to a FASTQ file, returns the expected BWA alignment file name."""
    if not fqFile : return ''
    prefix,ext = os.path.splitext(fqFile)
    return fqFile.replace(ext, '.aln')

def getReadLength(fqFile) :
    """Naive heuristic that determines fastq read length for a file based on the first read in the file."""
    result = 0
    ctr    = 0
    for line in ezopen(fqFile) :
        ctr += 1
        if ctr == 2 :
            result = len(line.strip())
            break
    return result

USAGE = """%prog original-dir [options]

Where:
    original-dir = directory containing original predictions
                   created with predict_graphs.py

Single-end read example:
    %prog original_predictions -1 myreads.fq 

Paired-end example:
    %prog original_predictions -1 myreads_1.fq -2 myreads_2.fq

Runs the realignment procedure to resolve exons in a set of splice graphs.
First it generates putative transcripts, then uses BWA to align reads
to the putative transcripts, and places updated graphs in an output
directory (local directory by default)."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-1', dest='first',   default=None,
        help='FASTQ file containing first half of mate pairs for paired-end reads, or the file containing single-end reads [default: %default]')
parser.add_option('-2', dest='second',  default=None,
        help='FASTQ file containing second half of mate pairs (paired-end reads only) [default: %default]')
parser.add_option('-d', dest='outdir',  default='.',           help='Output directory [default: %default]')
parser.add_option('-f', dest='fasta',   default=SG_FASTA_REF,  help='Genome reference file [default: %default]')
parser.add_option('-m', dest='model',   default=SG_GENE_MODEL, help='Gene model annotations (GFF3/GTF) [default: %default]')
parser.add_option('-v', dest='verbose', default=False,         help='Verbose mode [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

#-------------------------------------------------------------------------------------------------
# Parse the command line and make sure everything looks OK
MIN_ARGS = 1
if len(args) != MIN_ARGS :
    parser.print_help()
    if args : sys.stderr.write('\nExpected %d parameters; received %d:\n  %s\n' % (MIN_ARGS, len(args), '\n  '.join(args)))
    sys.exit(1)

errStrings = []
if not opts.model : errStrings.append(GENE_MODEL_ERR)
if not opts.fasta : errStrings.append(FASTA_FiLE_ERR)
if errStrings :
    parser.print_help()
    sys.stderr.write('\n%s\n' % '\n'.join(errStrings))
    sys.exit(1)

if not opts.first :
    parser.print_help()
    sys.stderr.write('\nYou must specify at least one FASTQ file (-1).\n')
    sys.exit(1)

bwaPath = findFile('bwa', os.environ['PATH'])
if not bwaPath :
    parser.print_help()
    sys.stderr.write('\nBWA executable not found in your PATH.\n')
    sys.exit(1)

# Check required files
origDir = args[0]
for f in [opts.model, opts.fasta] : validateFile(f)

sourceDir = os.path.abspath(origDir)
destDir   = os.path.abspath(opts.outdir)
if sourceDir == destDir :
    parser.print_help()
    sys.stderr.write('\nYou must specify different input and output splice graph directories:\n')
    sys.stderr.write('  input:  %s\n' % sourceDir)
    sys.stderr.write('  output: %s\n' % destDir)
    sys.exit(1)

firstFiles  = opts.first.split(',')
secondFiles = opts.second.split(',') if opts.second else []
if secondFiles and (len(firstFiles) != len(secondFiles)) :
    raise ValueError('You must provide the same numbers of paired-end mate files')
for f in firstFiles + secondFiles : validateFile(f)

writeStartupMessage()

#-------------------------------------------------------------------------------------------------
# 1. generate putative transcripts (only for updated/predicted genes)
graphFile = makeGraphListFile(sourceDir)
if not fileExists(graphFile) : raise Exception('No graphs created in %s' % sourceDir)

logMessage('\nGenerating putative transcript database\n')
cmd = '%s "%s" -f "%s" -U -o %s -m %s' % (TRANSCRIPT_SCRIPT, graphFile, opts.fasta, PUTATIVE_FASTA, TRANSCRIPT_MAP)
runCommand(cmd)
if not fileExists(PUTATIVE_FASTA) : raise Exception('No putative transcripts created in %s' % PUTATIVE_FASTA)

#-------------------------------------------------------------------------------------------------
# 2. run BWA on putative transcripts
#    a. create BWA database
#    b. run alignments (twice for PE data)
dbExists = True
for ext in BWA_EXTS :
    path = '%s.%s' % (PUTATIVE_FASTA, ext)
    if not os.path.exists(path) :
        dbExists = False
        break

if not dbExists :
    logMessage('\nCreating BWA indexes:\n')
    cmd = '%s index %s' % (bwaPath, PUTATIVE_FASTA)
    runCommand(cmd)

for ext in BWA_EXTS :
    path = '%s.%s' % (PUTATIVE_FASTA, ext)
    if not os.path.exists(path) : raise Exception("BWA indexing failed; missing %s" % path)

allCPU = multiprocessing.cpu_count()
useCPU = allCPU/2 if allCPU > 3 else 1
logMessage('\nRunning BWA alignments on FASTQ files using %d/%d CPUs\n' % (useCPU, allCPU))

for fastqFile in firstFiles + secondFiles :
    readLength = getReadLength(fastqFile)
    alnFile    = getBWAAlnFileName(fastqFile)
    options    = "-f %s -t %d -M 8 -R 2 -o 0 -e 0 -d 0 -i %d" % (alnFile, useCPU, readLength)
    cmd        = '%s aln %s "%s" "%s"' % (bwaPath, options, PUTATIVE_FASTA, fastqFile)
    runCommand(cmd)
    if not fileExists(alnFile) : raise Exception('BWA failed for %s: %s not found' % (fastqFile, alnFile))

if secondFiles :
    bwaLog = open('bwa_pe.log','w')
    for i in xrange(len(firstFiles)) :
        file1    = firstFiles[i]
        file2    = secondFiles[i]
        alnFile1 = getBWAAlnFileName(file1)
        alnFile2 = getBWAAlnFileName(file2)
        prefix   = 'pair_%d' % (i+1)
        rawSam   = '%s.raw' % prefix
        cmd      = '%s sampe -n 2 -N 2 -P -f "%s" "%s" "%s" "%s" "%s" "%s"' % (bwaPath, rawSam, PUTATIVE_FASTA, alnFile1, alnFile2, file1, file2)
        runCommand(cmd, stdout=bwaLog, stderr=bwaLog)
        if not fileExists(rawSam) : raise Exception('BWA paired-end alignments failed: %s not created.' % rawSam)

        # Filter out incorrect mate-pair alignments
        goodPairs = '%s.pairs' % prefix
        cmd       = 'get_good_pairs.py %s %s' % (rawSam, goodPairs)
        runCommand(cmd)
        if not os.path.exists(goodPairs) : raise Exception('get_good_pairs.py failed: %s not created.' % goodPairs)

        # Grab only unresolved/putative transcript alignments
        writeNovelAlignments(goodPairs, BWA_SAM, initiate=(i==0))

else : # single-end reads
    bwaLog = open('bwa_se.log','w')
    for i in xrange(len(firstFiles)) :
        file1    = firstFiles[i]
        alnFile1 = getBWAAlnFileName(file1)
        prefix   = 'single_%d' % (i+1)
        rawSam   = '%s.raw' % prefix
        cmd      = 'bwa samse -n 2 -f "%s" "%s" "%s" "%s"' % (rawSam, PUTATIVE_FASTA, alnFile1, file1)
        runCommand(cmd, stdout=bwaLog, stderr=bwaLog)
        if not fileExists(rawSam) : raise Exception('BWA single-end alignments failed: %s not created.' % rawSam)
        writeNovelAlignments(rawSam, BWA_SAM, initiate=(i==0))

#-------------------------------------------------------------------------------------------------
# 4. run update script
# Construct list of graphs in initial prediction dir
updateLog = open('fix_unresolved.log','w')
logMessage('\nUpdating original predictions.\n')

# Coverage-based update procedure
cmd = 'fix_unresolved.py "%s" "%s" "%s" -f "%s" -c -d "%s" -v' % (graphFile, TRANSCRIPT_MAP, BWA_SAM, PUTATIVE_FASTA, destDir)
runCommand(cmd, stdout=updateLog, stderr=updateLog)

logMessage('\nFinished.\n')
