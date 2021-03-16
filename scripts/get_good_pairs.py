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
from SpliceGrapher.formats.sam  import *
from optparse                   import OptionParser
import os,sys,math,functools,re

NONSENSE_CHROM = ''
NONSENSE_FLAG  = -1

def getStrand(s) :
    """Uses the SAM flag format to determine the strand."""
    return '+' if int(s)&16 else '-'

def validCigar(seq, cigar) :
    ungap = '%dM' % len(seq)
    if cigar in [ungap, EXACT_CIGAR] : return True
    try :
        # Convert S,H,I and D tokens into M or N:
        converted = convertCigarString(cigarTokens(cigar), len(seq))
        return True
    except ValueError :
        return False

USAGE = """%prog input-SAM-file output-SAM-file [options]

Identifies valid (concordant) paired-end read alignments in a SAM file
and places them in the given output file."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-i', dest='indicator', default='/',   help="Read mate indicator (e.g., '/' for /1, /2 style) [default: '%default']")
parser.add_option('-v', dest='verbose',   default=False, help='Verbose mode [default: %default]', action='store_true')
parser.add_option('-I', dest='pairids',   default=False, help='Look for pair ids [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 2 :
    parser.print_help()
    sys.exit(1)

samInput  = args[0]
samOutput = args[1]
validateFile(samInput)

reads     = {}
sys.stderr.write('First phase: identifying good pairs\n')
indicator = ProgressIndicator(10000000, verbose=opts.verbose)
badMatch  = 0
badCigar  = {}
headers   = []
for line in ezopen(samInput) :
    indicator.update()
    if line.startswith('@') :
        headers.append(line)
        continue
    parts = line.strip().split('\t')
    # MapSplice example:
    # 0                1   2       3       4   5   6   7       8   9            10           11              12
    # HWUSI-...976#0/1 83  Chr3    8011911 121 88M =   8011884 -27 AGC...TAG    ___..._bb    NM:i:0IH:i:1    HI:i:1

    # For a read pair to be a good match, we require:
    #  1. same chromosome
    try :
        chrom = parts[2].lower()
        seq   = parts[9]
    except IndexError :
        raise ValueError('** Illegal SAM record at line %d in %s\noffending line:\n%s' % (indicator.ctr, samInput, line))

    #  2. good cigar string 
    cigar = parts[5]
    if not validCigar(seq, cigar) :
        badCigar[cigar] = badCigar.setdefault(cigar,0) + 1
        continue

    readId = parts[0]
    if opts.pairids :
        pos = readId.rfind(opts.indicator)
        if pos < 0 : raise ValueError('Invalid read ID at line %d: %s\n' % (indicator.count(), readId))
        key = readId[:pos]
        reads.setdefault(key,{})
        ctr = readId[pos+1:]
    else :
        pos = len(readId)
        key = readId
        reads.setdefault(key,{})
        ctr = '2' if reads[key] else '1'

    flags = int(parts[1])
    try :
        ignore = reads[key][ctr]
        reads[key][ctr] = (NONSENSE_FLAG,NONSENSE_CHROM)
    except KeyError :
        reads[key][ctr] = (flags,chrom)
        if len(reads[key]) > 2 : raise ValueError('Too many identifiers for same prefix at line %d: %s\n' % (indicator.count(), readId))

indicator.finish()
badCigarReads = sum(badCigar.values())
badCigarTypes = len(badCigar)
acceptedPairs = 2*len(reads)
sys.stderr.write('  first pass: accepted %s/%s records; %s records with %s distinct bad CIGAR strings\n' % \
        (commaFormat(len(reads)), commaFormat(indicator.count()), commaFormat(badCigarReads), commaFormat(badCigarTypes)))

# We also require one part on the + strand and the other on the - strand
# Next, find reads with only 2 parts, each aligned exactly once, with opposite strands
sys.stderr.write('Second phase: removing multireads and mismatched pairs\n')
goodReads = {}
indicator.reset()
incomplete = 0
multiread  = 0
badStrands = 0
for key in reads :
    indicator.update()
    # Identify reads with missing mates
    if len(reads[key]) < 2 :
        incomplete += 1
        continue

    # Identify duplicate entries
    (a,b) = reads[key].values()
    (aflag,achr) = a
    (bflag,bchr) = b
    if NONSENSE_FLAG in [aflag,bflag] :
        multiread += 1
        continue

    # Alternate method for identifying chromosome mismatch
    if achr != bchr :
        #if opts.verbose : sys.stderr.write('%s chromosome mismatch: %s != %s\n' % (key, achr, bchr))
        badMatch += 1
        continue 

    # Reject matching strands (16='-', 0='+')
    if aflag&16 == bflag&16 :
        #if opts.verbose : sys.stderr.write('%s strands match: %d == %d\n' % (key, aflag, bflag))
        badStrands += 1
        continue

    goodReads[key] = True
indicator.finish()
sys.stderr.write('  second pass: %s good alignments; %s missing mates; %s multireads; %s chromosome mismatches; %s same-strand errors\n' % \
        (commaFormat(len(goodReads)), commaFormat(incomplete), commaFormat(multiread), commaFormat(badMatch), commaFormat(badStrands)))

if not goodReads :
    sys.stderr.write('No good read pairs found; exiting.\n')
    sys.exit(0)

sys.stderr.write('Third phase: writing good pairs to %s\n' % samOutput)
outStream = open(samOutput,'w')
if headers : outStream.writelines(headers)
indicator.reset()
for line in ezopen(samInput) :
    indicator.update()
    if line.startswith('@') : continue
    parts  = line.strip().split('\t')
    readId = parts[0]
    if opts.pairids :
        pos = readId.rfind(opts.indicator)
        key = readId[:pos]
    else :
        key = readId

    try :
        if goodReads[key] : outStream.write(line)
    except KeyError :
        pass

indicator.finish()
sys.stderr.write('Finished.\n')
