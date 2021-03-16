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
from SpliceGrapher.shared.utils               import *
from SpliceGrapher.statistics.GraphStatistics import *
from SpliceGrapher.SpliceGraph                import *
from optparse                                 import OptionParser
import os,sys

STAT_TYPES = [IR_ABBREV, ES_ABBREV, ALT5_ABBREV, ALT3_ABBREV]

def addLatexTableTitles(outStream, titles, cols=1) :
    outStream.write('  \\multicolumn{%d}{|c|}{%s} ' % (cols,titles[0]))
    for t in titles[1:-1] :
        outStream.write('\n  & \\multicolumn{%d}{|c|}{%s}' % (cols,t))
    outStream.write('\n  & %s' % titles[-1])
    outStream.write('  \\\\ \n')

def endLatexTableLine(outStream) :
    outStream.write('\\\\ \n')

def startLatexTable(outStream, title, asTypes) :
    colString = '|rr'*asTypes + '|r|'
    outStream.write('\\begin{tabular}{%s}\n' % colString)
    outStream.write('  \\multicolumn{%d}{|c|}{%s}' % (2*asTypes+1, title))
    endLatexTableLine(outStream)

def addLatexValues(outStream, count, percent) :
    outStream.write('%s   &    (%.1f\\%%)  &  ' % (commaFormat(count),percent))

def finishLatexTable(outStream) :
    outStream.write('\\end{tabular}\n')

USAGE = """%prog file-list [options]

Computes statistics for all the splice graphs given by the file list.
The file list may either be a list of splicegraph GFF files or a directory
that contains splicegraph GFF files."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-a', dest='annotate',default=False, help='Annotate graphs before generating statistics [default: %default]', action='store_true')
parser.add_option('-o', dest='output',  default=None,  help='Output file [default: stdout]')
parser.add_option('-v', dest='verbose', default=False, help='Verbose mode [default: %default]', action='store_true')
parser.add_option('-L', dest='latex',   default=False, help='Write table for import into a Latex document [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 1 :
    parser.print_help()
    sys.exit(1)

writeStartupMessage()

if os.path.isdir(args[0]) :
    fileList  = [os.path.join(args[0],f) for f in os.listdir(args[0])]
else :
    fileList  = [s.strip() for s in ezopen(args[0])]


# 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
#           Alternative Splicing Statistics for 28,412 Graphs (4,029 with AS)
# Intron Retention    Skipped Exon       Alt. 5'            Alt. 3'             Total
#    1,987 (35.9%)        550 ( 9.9%)      1,092 (19.8%)      1,899 (34.4%)       5,528

summary   = SummaryStatistics()
summary.addGraphFiles(fileList, verbose=opts.verbose, annotate=opts.annotate)
asTotal   = sum([summary.total(t) for t in STAT_TYPES])

outStream = sys.stdout if not opts.output else open(opts.output,'w')
titleString = 'Alternative Splicing Statistics for %s Graphs (%s with AS)' % (commaFormat(len(fileList)), commaFormat(summary.altCount()))
if opts.latex :
    startLatexTable(outStream, titleString, len(STAT_TYPES))
    titles = ["Intron Retention","Skipped Exon","Alt. 5'","Alt. 3'","Total"]
    addLatexTableTitles(outStream, titles, cols=2)
else :
    outStream.write('          %s\n' % titleString)
    outStream.write("Intron Retention    Skipped Exon       Alt. 5'            Alt. 3'             Total\n")

for asType in STAT_TYPES :
    typeTotal = summary.total(asType)
    typePct   = float(100*typeTotal)/asTotal if asTotal > 0 else 0.0
    if opts.latex :
        addLatexValues(outStream, typeTotal, typePct)
    else :
        outStream.write('%7s (%4.1f%%)    ' % (commaFormat(typeTotal), typePct))

if opts.latex :
    outStream.write('%8s  ' % commaFormat(asTotal))
    endLatexTableLine(outStream)
    finishLatexTable(outStream)
else :
    outStream.write('%8s\n' % commaFormat(asTotal))
