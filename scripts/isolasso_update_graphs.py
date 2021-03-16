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
from SpliceGrapher.formats.FastaLoader  import *
from SpliceGrapher.SpliceGraph  import *
from optparse                   import OptionParser
import os,sys

# IsoLasso GTF example (last field only):
#   gene_id "Inst1"; transcript_id "AT1G01100_p14"; FPKM  "0"; frac "1.000000"; conf_lo "0.0"; conf_hi "2.0"; cov "0.1";
# Required IsoLasso attributes:
TRANS_KEY = 'transcript_id'
FPKM_KEY  = 'FPKM'

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

USAGE = """%prog graph-list transcript-map IsoLasso-GTF [options]

Updates SpliceGrapher predictions using transcripts confirmed by IsoLasso."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-d', dest='outdir',    default='.',  help='Top-level output directory [default: %default]')
parser.add_option('-t', dest='threshold', default=1.0,  help='FPKM minimum threshold [default: %default]', type='float')
parser.add_option('-v', dest='verbose',   default=False,help='Verbose mode [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 3 :
    parser.print_help()
    sys.exit(1)

if opts.threshold < 0.0 :
    parser.print_help()
    sys.stderr.write('\nThreshold must be a non-negative value.\n')
    sys.exit(1)

for f in args : validateFile(f)
graphList   = args[0]
transMap    = args[1]
isoLassoGTF = args[2]

if opts.verbose : sys.stderr.write('Loading IsoLasso results from %s\n' % isoLassoGTF)
transDict = {}
for line in ezopen(isoLassoGTF) :
    attrs = attributeDict(line)
    transId = attrs[TRANS_KEY]
    fpkm    = float(attrs[FPKM_KEY])
    transDict[transId] = fpkm

formSet         = set([t for t in transDict if transDict[t] > opts.threshold])
formVals        = [transDict[t] for t in formSet]
minFPKM         = min(formVals)
maxFPKM         = max(formVals)
totalCandidates = len(formSet)
if opts.verbose : sys.stderr.write('  loaded %s forms (FPKM range %.2f to %.0f)\n' % (commaFormat(totalCandidates), minFPKM, maxFPKM))

nodePattern = {}
allIds      = set()
indicator   = ProgressIndicator(100000, verbose=opts.verbose)
if opts.verbose : sys.stderr.write('Loading exon pattern information from %s\n' % transMap)
for line in ezopen(transMap) :
    indicator.update()
    parts = line.strip().split('\t')
    key   = parts[0]
    allIds.add(key)
    if key not in formSet : continue
    nodePattern[key] = parts[1]
indicator.finish()
missing = formSet - set(nodePattern.keys())

if nodePattern :
    sys.stderr.write('  found exon patterns for %s/%s putative transcripts in %s\n' % \
            (commaFormat(len(nodePattern)), commaFormat(totalCandidates), transMap))
else :
    sys.stderr.write('No putative transcripts found in %s; exiting\n' % transMap)
    sys.stderr.write('formSet: %s\n' % str(formSet))
    sys.stderr.write('all known transcripts: %s\n' % str(allIds))
    sys.exit(0)

# Now update the graphs
allUpdatedNodes = set()
transcripts     = sorted(nodePattern.keys())
sys.stderr.write('Adding isoforms\n')
indicator.reset()
resolvedSet = set()
newTranscripts = 0
updatedGraphs  = 0
for line in ezopen(graphList) :
    # Load the original prediction
    graph  = getFirstGraph(line.strip())
    geneId = graph.name
    indicator.update()

    geneForms        = set([s for s in transcripts if s.startswith(graph.name)])
    validTranscripts = geneForms & formSet
    if validTranscripts : updatedGraphs += 1

    # Ensure that parent-child relationships exist for entire
    # transcript and resolve any unresolved nodes.
    for transId in validTranscripts :
        newTranscripts += 1
        resolvedSet.add(transId)
        nodes = [graph.nodeDict[n] for n in nodePattern[transId].split(',')] 
        prev  = nodes[0]
        if prev.isUnresolved() :
            resolvedSet.add(prev)
            prev.attrs[DISPOSITION_KEY] = PREDICTED_NODE

        for curr in nodes[1:] :
            if curr.isUnresolved() :
                resolvedSet.add(curr)
                curr.attrs[DISPOSITION_KEY] = PREDICTED_NODE
            prev.addChild(curr)
            curr.addIsoform(transId)
            prev = curr

    outDir  = os.path.join(opts.outdir,graph.chromosome)
    if not os.path.exists(outDir) : os.makedirs(outDir)
    outFile = os.path.join(outDir, '%s.gff'%graph.name.upper())
    graph.annotate()
    graph.writeGFF(outFile)

if opts.verbose :
    sys.stderr.write('Added %s novel transcripts to %s graphs\n' % (commaFormat(newTranscripts), commaFormat(updatedGraphs)))
    sys.stderr.write('Resolved %s previously unresolved nodes\n' % commaFormat(len(resolvedSet)))
sys.stderr.write('Finished.\n')
