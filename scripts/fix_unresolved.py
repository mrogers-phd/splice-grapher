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

PUTATIVE_PARENTS  = 'putative_parents'
PUTATIVE_CHILDREN = 'putative_children'

class PseudoNode(object) :
    def __init__(self, node, minpos, maxpos) :
        self.node   = node
        self.id     = node.id
        self.absmin = node.minpos
        self.absmax = node.maxpos
        self.minpos = minpos
        self.maxpos = maxpos

    def __eq__(self,o) :
        return self.id == o.id

    def __len__(self) :
        return self.maxpos-self.minpos+1

    def longString(self) :
        return '%s = %s' % (self, self.node)

    def __str__(self) :
        return '%s %d-%d' % (self.id, self.minpos, self.maxpos)

def delAttributes(node, keyList) :
    """Deletes the list of attributes from the given node."""
    for key in keyList :
        try :
            del node.attrs[key]
        except KeyError :
            continue

def getPseudoNodes(nodes) :
    """Returns a list of pseudo-nodes that consist of start/end positions
    relative to the start of a transcript.  Assumes the nodes are sorted."""
    total  = -1
    result = []
    for n in nodes :
        minpos = total+1
        maxpos = total+len(n)
        result.append(PseudoNode(n,minpos,maxpos))
        total  = maxpos
    return result

def getNodeList(unode, key, nodeDict) :
    """Returns a list of putative parent/child nodes for the given unresolved node."""
    result = []
    try :
        ids    = set(unode.attrs[key].split(','))
        keys   = set(nodeDict.keys())
        common = keys & ids
        for nid in common :
            result.append(nodeDict[nid])
    except KeyError : # key not in unode.attrs
        pass
    return result

def loadPutativeSamRecords(path, verbose=False) :
    """Loads records from a SAM file that are recognized as matches to
    a putative transcript (have the form <trans-id> <node-list>).
    Returns a dictionary that maps transcript ids to strings containing node lists."""
    indicator = ProgressIndicator(1000000, verbose=verbose)
    if verbose : sys.stderr.write('Loading SAM records from %s\n' % path)
    result = {}
    for line in ezopen(path) :
        indicator.update()
        if line.startswith('@') : continue

        # 0                                     1   2               3 
        # HWUSI-EAS1760_0012:8:120:5891:20838#0	99	AT5G14550_35	1162	...
        transId = line.split('\t')[2]
        result.setdefault(transId,[])
        result[transId].append(line.strip())
    indicator.finish()
    return result

def seqString(s, indent=0, width=60) :
    """Returns a DNA sequence string formatted to appear on multiple lines."""
    pfx    = ' '*indent if indent else ''
    result = ''
    for i in range(0,len(s),width) :
        if result : result += pfx
        result += s[i:i+width+1] + '\n'
    return result.strip()

USAGE = """%prog graph-list transcript-map SAM-file [options]

Updates SpliceGrapher predictions using ungapped read alignments
to putative transcripts.  Works with either single-end or paired-end reads."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-c', dest='coverage',    default=False, help='Use all coverage to accept nodes [default: %default]', action='store_true')
parser.add_option('-d', dest='outdir',      default='.',   help='Top-level output directory [default: %default]')
parser.add_option('-f', dest='fasta',       default=None,  help='Optional FASTA file for validating sequences [default: %default]')
parser.add_option('-O', dest='min_overlap', default=10,    help='Minimum overlap required to resolve a node [default: %default]', type='int')
parser.add_option('-R', dest='report',      default=None,  help='Report update statistics to file [default: %default]')
parser.add_option('-S', dest='sensitive',   default=False, help='Accepting nodes with overlapping reads [default: %default]', action='store_true')
parser.add_option('-v', dest='verbose',     default=False, help='Verbose mode [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 3 :
    parser.print_help()
    sys.exit(1)

if opts.sensitive and opts.coverage :
    parser.print_help()
    sys.stderr.write('\nYou cannot choose both sensitive (-S) and coverage (-c) options.\n')
    sys.exit(1)

if opts.sensitive : sys.stderr.write('Using sensitive resolution: accepting all nodes with overlapping reads\n')

for f in args : validateFile(f)
graphList = args[0]
transMap  = args[1]
samFile   = args[2]

# Load all relevant SAM records
samRecords = loadPutativeSamRecords(samFile, verbose=opts.verbose)
if not samRecords :
    sys.stderr.write('  no SAM records found; exiting\n')
    sys.exit(0)
sys.stderr.write('  found SAM records for %s transcripts in %s\n' % (commaFormat(len(samRecords)), samFile))

fastaDB = FastaLoader(opts.fasta) if opts.fasta else None

nodePattern = {}
indicator   = ProgressIndicator(100000, verbose=opts.verbose)
if opts.verbose : sys.stderr.write('Loading putative transcript information from %s\n' % transMap)
for line in ezopen(transMap) :
    indicator.update()
    if line.startswith('chr') : continue
    parts = line.strip().split('\t')
    key   = parts[0]
    if not key in samRecords : continue
    nodePattern[key] = parts[1]
indicator.finish()

if nodePattern :
    sys.stderr.write('  found patterns for %s putative transcripts in %s\n' % (commaFormat(len(nodePattern)), transMap))
else :
    sys.stderr.write('No putative transcript patterns found in %s; exiting\n' % transMap)
    sys.exit(0) # Not really an error, but no reason to continue

# Now update the graphs
allUnresolved   = set()
allUpdatedNodes = set()
transcripts     = sorted(nodePattern.keys())
sys.stderr.write('Resolving nodes\n')
indicator.reset()
for line in ezopen(graphList) :
    # Load the original prediction
    graph  = getFirstGraph(line.strip())
    geneId = graph.name
    indicator.update()

    allUnresolved.update(graph.unresolvedNodes())
    geneForms       = [s for s in transcripts if s.startswith(graph.name)]
    updatedNodeSet  = set()
    hasUpdatedNodes = {}.fromkeys(geneForms,False)
    for transId in geneForms :
        # Get transcript nodes and sort them
        transcriptNodes = [graph.nodeDict[n] for n in nodePattern[transId].split(',')]
        transcriptNodes.sort(reverse=(graph.strand=='-'))

        # Convert real nodes into pseudo-nodes with positions
        # relative to the start of the transcript
        transcriptPseudoNodes = getPseudoNodes(transcriptNodes)

        transcriptLength      = sum([len(n) for n in transcriptPseudoNodes])
        #if fastaDB : assert(transcriptLength == len(fastaDB[transId]))

        # Form read pairs from SAM records for current transcript
        readPair = {}
        for s in samRecords[transId] :
            r   = SAMRecord(s)
            key = r.read()
            readPair.setdefault(key,[])
            readPair[key].append(r)

        resolvedNodeSet = set()
        if opts.coverage :
            coverage = set()
            for pid in readPair :
                pairMin = sys.maxint
                pairMax = 0
                for r in readPair[pid] :
                    last    = r.pos() + len(r.query()) - 1
                    pairMin = min(pairMin, r.pos())
                    pairMax = max(pairMax, last)
                coverage.update(xrange(pairMin,pairMax+1))

            # Find unresolved nodes that overlap the range
            for n in transcriptPseudoNodes :
                if n.node.isUnresolved() :
                    # No overlaps used for first/last node:
                    loOffset = 0 if n == transcriptPseudoNodes[0] else opts.min_overlap
                    hiOffset = 0 if n == transcriptPseudoNodes[-1] else opts.min_overlap+1
                    missingCoverage = set(xrange(n.minpos-loOffset, n.maxpos+hiOffset)) - coverage
                    if missingCoverage : continue
                    resolvedNodeSet.add(n.node)

        else :
            # Find unresolved nodes that fall within any read pairs for the transcript
            for pid in readPair :
                # Find all nodes associated with a read pair.  Note that read
                # positions are relative to the start of a transcript.
                # Track coverage across full pair:
                pairMin = sys.maxint
                pairMax = 0
                for r in readPair[pid] :
                    last    = r.pos() + len(r.query()) - 1
                    pairMin = min(pairMin, r.pos())
                    pairMax = max(pairMax, last)

                # Find unresolved nodes that overlap the range
                for n in transcriptPseudoNodes :
                    if n.maxpos < pairMin or n.minpos > pairMax : continue
                    if not n.node.isUnresolved() : continue
                    overlapSize      = min(len(n), min(n.maxpos-pairMin+1, pairMax-n.minpos+1))
                    required_overlap = min(opts.min_overlap, len(n)) if opts.sensitive else len(n)
                    if overlapSize < required_overlap : continue
                    resolvedNodeSet.add(n.node)

        if not resolvedNodeSet : continue
        hasUpdatedNodes[transId] = True

        # Only allow other nodes from the same transcript to become
        # parents/children of the resolved node.  Add other nodes
        # from the transcripts to a set of candidates for that node.
        transcriptDict = {}
        for n in transcriptNodes :
            transcriptDict[n.id] = n

        for node in resolvedNodeSet :
            for p in getNodeList(node, PUTATIVE_PARENTS, transcriptDict) :
                p.addChild(node)
            for c in getNodeList(node, PUTATIVE_CHILDREN, transcriptDict) :
                node.addChild(c)
            # Add to resolved node set for changing the
            # node disposition after all updates are done
            updatedNodeSet.add(node)

    if updatedNodeSet and opts.verbose : sys.stderr.write('  resolved %d nodes in %s\n' % (len(updatedNodeSet), graph.name))

    # Change the newly resolved nodes to predicted nodes
    # and remove putative information
    allUpdatedNodes |= updatedNodeSet
    for n in updatedNodeSet :
        n.addAttribute(DISPOSITION_KEY, PREDICTED_NODE)
        delAttributes(n, [PUTATIVE_PARENTS, PUTATIVE_CHILDREN])

    outDir  = os.path.join(opts.outdir,graph.chromosome)
    if not os.path.exists(outDir) : os.makedirs(outDir)
    outFile = os.path.join(outDir, '%s.gff'%graph.name.upper())
    graph.annotate()
    graph.writeGFF(outFile)

resolvedCount  = len(allUpdatedNodes)
origUnresolved = len(allUnresolved)
pct = resolvedCount*100.0/origUnresolved if origUnresolved else 0.0
sys.stderr.write('Resolved %s/%s previously unresolved nodes (%.1f%%)\n' % (commaFormat(resolvedCount), commaFormat(origUnresolved), pct))

if opts.report :
    rptStream = open(opts.report,'w')
    rptStream.write('%s\n' % ' '.join(args))
    rptStream.write('Resolved %s/%s previously unresolved nodes (%.1f%%)\n' % (commaFormat(resolvedCount), commaFormat(origUnresolved), pct))
