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

USAGE = """%prog graph-list transcript-map psginfer-isoform-table [options]

Updates SpliceGrapher predictions using transcripts confirmed by psgInfer."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-d', dest='outdir',    default='.',  help='Top-level output directory [default: %default]')
parser.add_option('-t', dest='threshold', default=0.1,  help='Threshold for accepting isoforms (between 0 and 1)[default: %default]', type='float')
parser.add_option('-v', dest='verbose',   default=False,help='Verbose mode [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 3 :
    parser.print_help()
    sys.exit(1)

if opts.threshold < 0.0 or opts.threshold > 1.0 :
    parser.print_help()
    sys.stderr.write('\nThreshold must be between 0 and 1.\n')
    sys.exit(1)

for f in args : validateFile(f)
graphList = makeGraphListFile(args[0]) if os.path.isdir(args[0]) else args[0]
transMap  = args[1]
psgTable  = args[2]

# PSGInfer isoform table:
# gene_id     transcript_id   prob
# AT1G01560   AT1G01560.1     1.0
# AT1G06960   AT1G06960.1     0.36
# AT1G06960   AT1G06960_p1    0.64
if opts.verbose : sys.stderr.write('Loading psgInfer results from %s\n' % psgTable)
putativeForms = set()
ctr     = 0
for line in ezopen(psgTable) :
    if line.startswith('gene_id') : continue
    ctr += 1
    parts  = line.strip().split('\t')
    try :
        prob   = float(parts[2])
    except IndexError :
        sys.stderr.write('line %d: bad record in psgInfer table:\n%s' % (ctr, line))
        sys.exit(1)

    if prob < opts.threshold : continue

    # Only store putative/unresolved transcripts:
    formId = parts[1]
    if formId.find('_p') >= 0 or formId.find('_u') >= 0 :
        putativeForms.add(formId)
totalCandidates = len(putativeForms)
if opts.verbose : sys.stderr.write('  loaded %s forms\n' % commaFormat(totalCandidates))

# Map transcript ids such as 'AT1G06960_p1' to lists of exon ids
nodePattern = {}
indicator   = ProgressIndicator(100000, verbose=opts.verbose)
if opts.verbose : sys.stderr.write('Loading transcript node information from %s\n' % transMap)
for line in ezopen(transMap) :
    indicator.update()
    parts = line.strip().split('\t')
    key   = parts[0]
    if key not in putativeForms : continue
    nodePattern[key] = parts[1]
indicator.finish()
missing = putativeForms - set(nodePattern.keys())

if nodePattern :
    sys.stderr.write('  found patterns for %s/%s putative transcripts in %s\n' % \
            (commaFormat(len(nodePattern)), commaFormat(totalCandidates), transMap))
    if missing :
        sys.stderr.write('Examples of missing ids:\n')
        missingList = list(missing)
        for x in missingList[:5] :
            sys.stderr.write('     %s\n' % x)
else :
    sys.stderr.write('No putative transcript patterns found in %s; exiting\n' % transMap)
    sys.exit(0) # Not really an error, but no reason to continue

# Now update the graphs
allUpdatedNodes = set()
transcripts     = sorted(nodePattern.keys())
resolvedSet     = set()
newTranscripts  = 0
updatedGraphs   = 0
sys.stderr.write('Adding isoforms\n')
indicator.reset()
for line in ezopen(graphList) :
    # Load the original prediction
    graph  = getFirstGraph(line.strip())
    geneId = graph.name
    indicator.update()

    # Note that for species like human, the only transcripts that match the
    # gene name will be those that SpliceGrapher created.  In Arabidopsis,
    # both known and predicted transcripts will match.
    geneForms        = set([s for s in transcripts if s.startswith(graph.name)])

    # Reduce set to putative/predicted forms:
    validTranscripts = geneForms & putativeForms

    graphResolved = set()
    if validTranscripts :
        # Ensure that parent-child relationships exist for entire
        # transcript and resolve any unresolved nodes.
        for transId in validTranscripts :
            newTranscripts += 1
            nodes = [graph.nodeDict[n] for n in nodePattern[transId].split(',')] 
            prev  = nodes[0]
            if prev.isUnresolved() :
                graphResolved.add(prev)
                resolvedSet.add(prev)
                old = prev.attrs[DISPOSITION_KEY]
                prev.attrs[DISPOSITION_KEY] = PREDICTED_NODE
                new = prev.attrs[DISPOSITION_KEY]

            prev.addIsoform(transId)
            for curr in nodes[1:] :
                if curr.isUnresolved() :
                    graphResolved.add(curr)
                    resolvedSet.add(curr)
                    curr.attrs[DISPOSITION_KEY] = PREDICTED_NODE
                curr.addIsoform(transId)
                prev.addChild(curr)
                prev = curr
        updatedGraphs += 1

    outDir  = os.path.join(opts.outdir,graph.chromosome)
    if not os.path.exists(outDir) : os.makedirs(outDir)
    outFile = os.path.join(outDir, '%s.gff'%graph.name.upper())
    graph.annotate()
    graph.writeGFF(outFile)

if opts.verbose :
    sys.stderr.write('Added %s novel transcripts to %s graphs\n' % (commaFormat(newTranscripts), commaFormat(updatedGraphs)))
    sys.stderr.write('Resolved %s previously unresolved nodes\n' % commaFormat(len(resolvedSet)))
