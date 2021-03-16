#!/usr/bin/python
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
from SpliceGrapher.shared.config import *
from SpliceGrapher.shared.utils  import *
from SpliceGrapher.formats.sam   import *
from SpliceGrapher.shared.ShortRead import *
from optparse                    import OptionParser
import os,sys

USAGE = """%prog SAM-file [options]

Converts a SAM/BAM file into a SpliceGrapher depth file."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-o', dest='output',  default=None,  help='Output file [default: same file with .depths extension]')
parser.add_option('-v', dest='verbose', default=False, help='Verbose mode [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 1 :
    parser.print_help()
    sys.exit(1)

samFile = args[0]
validateFile(samFile)

if not opts.output :
    pfx,ignore  = os.path.splitext(samFile)
    opts.output = '%s.depths' % pfx

if os.path.isfile(opts.output) :
    sys.stderr.write('%s will be over-written after SAM data has been loaded\n' % opts.output)

outStream = open(opts.output, 'w')
oldChrom  = None
chrLines  = []
used      = set()

sys.stderr.write(timeString('Parsing records\n'))
invalidRecs = 0
totalRecs   = 0
validChrom  = set()
for line in samIterator(samFile) :
    # Ignore headers
    if line.startswith('@') :
        if line.startswith('@SQ') :
            parts = line.strip().split('\t')
            # @SQ	SN:20	LN:64444167
            [ignore,chrom] = parts[1].split(':')
            validChrom.add(chrom)
        continue
    totalRecs += 1

    s     = line.strip()
    parts = s.split('\t')
    chrom = parts[2]
    if chrom not in validChrom :
        if opts.verbose : sys.stderr.write('  skipping record with invalid chromosome %s\n' % chrom)
        invalidRecs += 1
        continue

    if chrom != oldChrom :
        # write current set of records
        if chrLines :
            if opts.verbose : sys.stderr.write('  converting SAM records to depths\n')
            depthDict,jctDict = getSamReadData(chrLines)
            if opts.verbose : sys.stderr.write('  writing to %s\n' % opts.output)
            writeDepths(outStream, depthDict=depthDict, jctDict=jctDict)
            used.add(oldChrom)

        # initialize new chromosome
        chrLines = []
        if chrom in used :
            raise ValueError("%s appears to be unsorted (chromosome '%s' found twice)" % (samFile, chrom))
        elif opts.verbose :
            sys.stderr.write('chromosome %s:\n' % chrom)
            sys.stderr.write('  reading records\n')

    chrLines.append(s)
    oldChrom = chrom

if chrLines :
    if opts.verbose : sys.stderr.write('  converting SAM records to depths\n')
    depthDict,jctDict = getSamReadData(chrLines)
    if opts.verbose : sys.stderr.write('  writing to %s\n' % opts.output)
    writeDepths(outStream, depthDict=depthDict, jctDict=jctDict)

if invalidRecs > 0 :
    sys.stderr.write('Omitted %s invalid records out of %s\n' % (commaFormat(invalidRecs), commaFormat(totalRecs)))
if opts.verbose : sys.stderr.write('Finished.\n')
