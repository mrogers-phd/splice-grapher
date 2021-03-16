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
from SpliceGrapher.shared.utils import *
from optparse                   import OptionParser
import os,sys,gzip

USAGE = """%prog SAM-files [options]

Collates a list of SAM files into distinct files for each chromosome."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-d', dest='outdir',  default='.',   help='Output directory [default: %default]')
parser.add_option('-v', dest='verbose', default=False, help='Verbose mode [default: %default]', action='store_true')
parser.add_option('-z', dest='gzip',    default=False, help='Use gzip compression on output files [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if len(args) < 1 :
    parser.print_help()
    sys.exit(1)

for f in args : validateFile(f)
headerSet = set()
headers   = []
streams   = {}
ctr       = 0
totalRecs = 0
for f in args :
    ctr += 1
    if opts.verbose : sys.stderr.write('Parsing file %d/%d: %s\n' % (ctr,len(args),f))

    indicator = ProgressIndicator(1000000, verbose=opts.verbose)
    for line in ezopen(f) :
        indicator.update()
        totalRecs += 1
        parts      = line.strip().split('\t')
        if line.startswith('@') :
            headerSet.add(line)
            headers = sorted(headerSet)
            continue

        chrom = parts[2].capitalize()
        if chrom not in streams :
            if not os.path.isdir(opts.outdir) :
                if opts.verbose : sys.stderr.write('Creating directory %s\n' % opts.outdir)
                os.makedirs(opts.outdir)
            outPath = os.path.join(opts.outdir, '%s.sam'%chrom)
            if opts.gzip : outPath += '.gz'
            streams[chrom] = gzip.open(outPath,'w') if opts.gzip else open(outPath,'w')
            for h in headers :
                streams[chrom].write(h)
        streams[chrom].write(line)
    indicator.finish()

sys.stderr.write('Finished processing %d files and %s SAM records.\n' % (ctr, commaFormat(totalRecs)))
