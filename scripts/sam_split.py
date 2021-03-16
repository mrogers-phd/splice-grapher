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

USAGE = """%prog SAM-file [options]

Splits a SAM file into individual files for each target
sequence (usually chromosomes).  Output files will have
the names [pre][chromosome][suf].sam, where pre and suf
are optional prefix and suffix strings."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-p', dest='prefix',  default='',    help='Optional prefix string for output files [default: none]')
parser.add_option('-s', dest='suffix',  default='',    help='Optional suffix string for output files [default: none]')
parser.add_option('-v', dest='verbose', default=False, help='Verbose mode [default: %default]', action='store_true')
parser.add_option('-z', dest='gzip',    default=False, help='Use gzip compression on output [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 1 :
    parser.print_help()
    sys.exit(1)

validateFile(args[0])
headers   = []
streams   = {}
indicator = ProgressIndicator(1000000, verbose=opts.verbose)
for line in ezopen(args[0]) :
    indicator.update()
    if line.startswith('@') :
        headers.append(line)
    else :
        parts = line.split('\t')
        key   = parts[2]
        try :
            streams[key].write(line)
        except KeyError :
            outPath = '%s%s%s.sam' % (opts.prefix, key, opts.suffix)
            if opts.gzip : outPath += '.gz'
            indicator.finish()
            sys.stderr.write('Creating %s\n' % outPath)
            indicator.reset()
            streams[key] = gzip.open(outPath,'w') if opts.gzip else open(outPath,'w')
            streams[key].writelines(headers)
            streams[key].write(line)

indicator.finish()
