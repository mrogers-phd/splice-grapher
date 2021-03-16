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
from SpliceGrapher.formats.gtf  import *
from optparse                   import OptionParser
from sys import maxint as MAXINT
import os,sys,gzip

USAGE = """%prog GTF-file [options]

Converts an ENSEMBL GTF gene model annotation file into its GFF3 equivalent.
Note that it only accepts records with the protein_coding tag by default."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-A', dest='alltypes', default=False, help='Accept all source types (overrides -E and -S) [default: %default]', action='store_true')
parser.add_option('-o', dest='output',   default=None,  help='Output file [default: stdout]')
parser.add_option('-E', dest='ensembl',  default=False, help='Accept all ENSEMBL source types (see --show-types) [default: %default]', action='store_true')
parser.add_option('-S', dest='sources',  default=PROTEIN_CODING,  help='Comma-separated list of GTF source types to accept [default: %default]')
parser.add_option('-v', dest='verbose',  default=False, help='Verbose mode [default: %default]', action='store_true')
parser.add_option('-z', dest='gzip',     default=False, help='Use gzip compression on output [default: %default]', action='store_true')
parser.add_option('--show-types', dest='showtypes', default=False, help='Outputs ENSEMBLE source types and exits. [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if opts.showtypes :
    print("Known ENSEMBL source types:")
    for t in ALL_ENSEMBL_SOURCES :
        print("  ", t)
    sys.exit(0)

if len(args) != 1 :
    parser.print_help()
    sys.exit(1)

outStream = sys.stdout
if opts.output :
    outStream = gzip.open(opts.output,'w') if opts.gzip else open(opts.output,'w')

KNOWN_SOURCES = opts.sources.split(',') if not opts.ensembl else ALL_ENSEMBL_SOURCES

# First load all records from the input GTF file and
# store them as GTF_Line records:
indicator = ProgressIndicator(1000000, verbose=opts.verbose)
chromDict = {}
skipped   = {}
lineNo    = 0
for line in ezopen(args[0]) :
    indicator.update()
    lineNo += 1
    if line.startswith('#') : continue
    rec = GTF_Line(line)

    # If filtering is on, ignore records associated with unrecognized gene types
    if (not opts.alltypes) and (rec.source() not in KNOWN_SOURCES) :
        try :
            skipped[rec.source()] += 1
        except KeyError :
            skipped[rec.source()] = 1
            if opts.verbose : sys.stderr.write('skipping %s records\n' % rec.source())
        continue 

    chrom = rec.seqname()
    chromDict.setdefault(chrom, GTF_Chromosome(chrom))

    # Ignore any records not associated with a gene (regardless of type)
    key = rec.gene_name()
    if not key : continue

    chromDict[chrom].append(rec)
indicator.finish()

# Sort records by chromosome before writing them to GFF3 file
keys = sorted(chromDict.keys())
if opts.verbose :
    print "Stored information for %d chromosomes:" % len(keys)
    if skipped :
        print "  skipped the following record types:"
        for k in sorted(skipped.keys()) :
            print "   %s (%s records)" % (k, commaFormat(skipped[k]))

for k in keys :
    # Search chromosome for invalid genes
    invalid = False
    for g in chromDict[k].keys() :
        if not chromDict[k][g].valid :
            invalid = True
            break

    if opts.verbose : 
        if invalid :
            print "  chromosome %s: %12d genes (%d invalid **)" % (k, len(chromDict[k]), invalid)
        else :
            print "  chromosome %s: %12d genes" % (k, len(chromDict[k]))

    # Write out entire chromosome
    chromDict[k].writeGFF3(outStream)
