#! /usr/bin/env python
# Copyright (C) 2015 by Colorado State University
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
from SpliceGrapher.shared.utils               import *
from SpliceGrapher.statistics.GraphStatistics import *
from SpliceGrapher.SpliceGraph                import *
from optparse                                 import OptionParser
import os,sys

USAGE = """%prog file-list [options]

Shows statistics for each of the splice graphs given by the file list."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-A', dest='all',     default=False, help='Show all genes [default: only AS genes]', action='store_true')
parser.add_option('-C', dest='csv',     default=False, help='Use CSV spreadsheet output format [default: tab-delimited]', action='store_true')
parser.add_option('-o', dest='output',  default=None,  help='Output file [default: stdout]')
parser.add_option('-v', dest='verbose', default=False, help='Verbose mode [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 1 :
    parser.print_help()
    sys.exit(1)

writeStartupMessage()

if os.path.isdir(args[0]) :
    fileList  = [os.path.join(args[0],f) for f in os.listdir(args[0])]
else :
    fileList  = [s.strip() for s in ezopen(args[0])]

fileList.sort()
graphList = [getFirstGraph(f) for f in fileList]

outStream = sys.stdout if not opts.output else open(opts.output,'w')
if opts.csv :
    outStream.write('Gene,Alt5,Alt3,ES,IR,Total\n')
    lineFmt = '%s,%d,%d,%d,%d,%d\n'
else :
    outStream.write('Gene\tAlt5\tAlt3\tES\tIR\tTotal\n')
    lineFmt = '%s\t%d\t%d\t%d\t%d\t%d\n'

for g in graphList :
    stats   = GraphStatistics(g)
    total   = stats.alt5Count()+stats.alt3Count()+stats.esCount()+stats.irCount()
    if not opts.all and (total == 0) : continue
    outStream.write(lineFmt % (g.name, stats.alt5Count(), stats.alt3Count(), stats.esCount(), stats.irCount(), total))
outStream.close()
