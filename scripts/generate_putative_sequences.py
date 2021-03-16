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
from SpliceGrapher.shared.config          import *
from SpliceGrapher.shared.utils           import *
from SpliceGrapher.formats.fasta          import FastaRecord
from SpliceGrapher.formats.FastaLoader    import *
from SpliceGrapher.formats.GeneModel      import *
from SpliceGrapher.SpliceGraph            import *
from SpliceGrapher.shared.SpliceGraphPath import *

from optparse import OptionParser
import os,sys

PUTATIVE_CHILDREN = 'putative_children'
PUTATIVE_PARENTS  = 'putative_parents'
MARKER_KEY        = 'unresolved_marker'

def addGraph(G, geneModel) :
    """Creates a gene from the given splice graph and adds the gene to the gene model instance."""
    G.attrs.setdefault('ID',G.name)
    G.attrs.setdefault('Name',G.name)
    gene = Gene(G.name, None, start(G), end(G), G.chromosome, G.strand, attr=G.attrs)

    geneModel.addChromosome(1, G.maxpos, G.chromosome)
    geneModel.addGene(gene)

    isoDict = G.isoformDict()
    if not isoDict : raise ValueError('No isoforms found in %s' % G.name)

    keys    = sorted(isoDict.keys())
    for k in keys :
        nodes = isoDict[k]
        form  = Isoform(k, start(nodes[0]), end(nodes[-1]), G.chromosome, G.strand)
        if not nodes : raise ValueError('No nodes in isoform %s' % k)
        for node in nodes :
            exon = Exon(start(node), end(node), G.chromosome, G.strand)
            gene.addExon(form, exon)

def isMarked(node) :
    """Returns True if a node is marked; False otherwise."""
    try :
        return node.attrs[MARKER_KEY]
    except KeyError :
        return False

def hasMarkedNodes(nodeList) :
    """Returns True if a path has any marked nodes; False otherwise."""
    for n in nodeList :
        if isMarked(n) : return True
    return False

def hasPredictedNodes(nodeList) :
    """Returns True if a graph has any predicted nodes; False otherwise."""
    for n in nodeList :
        if n.isPredicted() : return True
    return False

def hasUnresolvedNodes(nodeList) :
    """Returns True if a graph has any unresolved nodes; False otherwise."""
    for n in nodeList :
        if n.isUnresolved() : return True
    return False

def putativeChildren(node) :
    """Returns a list of the putative child nodes for an unresolved node."""
    try :
        if node.attrs[PUTATIVE_CHILDREN] :
            return node.attrs[PUTATIVE_CHILDREN].split(',')
    except KeyError :
        pass
    return []

def putativeParents(node) :
    """Returns a list of the putative parent nodes for an unresolved node."""
    try :
        if node.attrs[PUTATIVE_PARENTS] :
            return node.attrs[PUTATIVE_PARENTS].split(',')
    except KeyError :
        pass
    return []

def nodeString(nodeList) :
    """Provides a consistent way to represent nodes in a list."""
    return ','.join([n.id for n in nodeList])

def end(x) :
    return x.maxpos if x.strand == '+' else x.minpos

def start(x) :
    return x.minpos if x.strand == '+' else x.maxpos


USAGE = """%prog graph-files [options]

Where:
    graph-files  is a file containing a list of splice graphs
                 or a top-level directory with the structure:
                 top-level-dir/chromosome-dir/gff-files

Generates putative sequences for all possible paths through each
of the splice graphs provided."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-A', dest='allgenes',   default=False,        help='Include all genes in list [default: %default]', action='store_true')
parser.add_option('-f', dest='fasta',      default=SG_FASTA_REF, help='FASTA reference sequence [default: %default]')
parser.add_option('-l', dest='seqlimit',   default=1000,         help='Exclude graphs with more paths than this [default: %default]', type='int')
parser.add_option('-o', dest='output',     default=None,         help='Output FASTA file [default: one per chromosome]')
parser.add_option('-U', dest='unresolved', default=False,        help='Include unresolved-node transcripts [default: %default]', action='store_true')
parser.add_option('-m', dest='mapfile',    default=None,         help='Create a file that maps transcript ids to exon ids [default: %default]')
parser.add_option('-M', dest='savemodel',  default=None,         help='Convert the putative graphs to a gene model [default: %default]')
parser.add_option('-v', dest='verbose',    default=False,        help='Verbose mode [default: %default]', action='store_true')
parser.add_option('--gff3', dest='gffmodel',default=False,help='Save model using GFF3 format [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

MIN_ARGS = 1
if len(args) != MIN_ARGS :
    parser.print_help()
    if args : sys.stderr.write('\nExpected %d parameters; received %d:\n  %s\n' % (MIN_ARGS, len(args), '\n  '.join(args)))
    sys.exit(1)

if not opts.fasta :
    parser.print_help()
    sys.stderr.write('\nYou must provide a FASTA reference (-f)\n')
    sys.exit(1)

if opts.gffmodel and not opts.savemodel :
    parser.print_help()
    sys.stderr.write('\nYou must provide an output file (-M) to save as GFF3 gene models (--gff)\n')
    sys.exit(1)

fileDB    = makeGraphListFile(args[0]) if os.path.isdir(args[0]) else args[0]
fileList  = [s.strip() for s in ezopen(fileDB)]

# Loading graphs should be quick; don't make the user
# wait if there are errors
if opts.verbose : sys.stderr.write('Loading graphs:\n')
indicator  = ProgressIndicator(10000, verbose=opts.verbose)
graphs     = {}
loaded     = 0
unresolved = 0
for f in fileList :
    indicator.update()
    graph         = getFirstGraph(f)
    hasUnresolved = hasUnresolvedNodes(graph.nodeDict.values())
    if not (opts.allgenes or hasUnresolved) : continue

    if hasUnresolved :
        # update parent/child relationships for
        # unresolved nodes and "pretend" they are predicted
        unresolved += 1
        nodes = graph.nodeDict.values()
        for n in nodes :
            if not n.isUnresolved() : continue

            parents = putativeParents(n)
            for pid in parents :
                pnode = graph.nodeDict[pid]
                pnode.addChild(n)

            children = putativeChildren(n)
            for cid in children :
                cnode = graph.nodeDict[cid]
                n.addChild(cnode)

            n.addAttribute(DISPOSITION_KEY, PREDICTED_NODE)
            n.addAttribute(MARKER_KEY, True)

    graphs.setdefault(graph.chromosome,{})
    graphs[graph.chromosome][graph.name] = graph
    loaded += 1

indicator.finish()
sys.stderr.write("Loaded %s graphs from of %s files (%s with unresolved nodes)\n" % (commaFormat(loaded), commaFormat(len(fileList)), commaFormat(unresolved)))
if loaded == 0 :
    sys.stderr.write('Nothing to do; exiting\n')
    sys.exit(0)

# Next load FASTA sequences
db = FastaLoader(opts.fasta, verbose=opts.verbose)

if opts.verbose : sys.stderr.write('Generating sequences:\n')
if opts.output : outStream = open(opts.output,'w')
if opts.savemodel : model = GeneModel(None)

# Take graphs in order by chromosome and gene name
mapStream  = open(opts.mapfile,'w') if opts.mapfile else None
chromList  = sorted(graphs.keys())
fastaList  = []
seqCount   = 0
graphCount = 0
tooLong    = 0
for c in chromList :
    graphList   = sorted(graphs[c].keys())
    graphCount +=  len(graphList)

    if opts.verbose : sys.stderr.write('  %s graphs in %s ' % (commaFormat(len(graphList)), c))
    indicator = ProgressIndicator(10000, verbose=opts.verbose)
    for gid in graphList :
        indicator.update()
        graph   = graphs[c][gid]
        reverse = (graph.strand == '-')
        nodes   = graph.nodeDict.values()
        chrom   = graph.chromosome

        # Grab the DNA sequence for every distinct exon
        seqDict = {}
        for n in nodes :
            # Include dimers before and after each exon to make it easy to verify
            (a,d) = (n.maxpos+1, n.minpos-3) if reverse else (n.minpos-3, n.maxpos+1)
            seq   = db.subsequence(chrom, a, d, reverse=reverse)
            seqDict[n.id] = seq[2:-2]

        # Create a DNA sequence for every path through the graph
        paths = []
        try :
            paths = getAllPaths(graph, limit=opts.seqlimit)
            paths.sort()
        except ValueError, ve :
            tooLong += 1
            continue

        isoDict    = graph.isoformDict()
        isoStrings = {}
        # Map node strings back to isoform name
        for k in isoDict.keys() :
            nodes = isoDict[k]
            key   = nodeString(nodes)
            isoStrings[key] = k
            
        isoCtr = 0
        for path in paths :
            marked = hasMarkedNodes(path.nodes)
            if marked and not opts.unresolved : continue

            header  = None
            nodeStr = nodeString(path.nodes)
            try :
                header = isoStrings[nodeStr]
            except KeyError : # novel paths
                # Do not create novel paths through graphs that have
                # no predictions; just output their known isoforms.
                if not (marked or hasPredictedNodes(path.nodes)) : continue
                isoCtr += 1
                tChar   = 'u' if marked else 'p'
                # produces strings like: AT1G01020_p3 or AT1G01060_u1
                header  = '%s_%s%d' % (graph.name, tChar, isoCtr)

            if opts.mapfile and header : mapStream.write('%s\t%s\n' % (header, nodeStr))

            seq  = ''.join([seqDict[n.id] for n in path.nodes])
            if opts.output : outStream.write('>%s\n%s\n' % (header, seq))
            seqCount += 1
            if opts.savemodel :
                for n in path.nodes :
                    n.addIsoform(header)

        if opts.savemodel :
            addGraph(graph, model)

    indicator.finish()

if opts.output : outStream.close()

if opts.savemodel :
    if opts.gffmodel :
        model.writeGFF(opts.savemodel)
    else :
        model.writeGTF(opts.savemodel)

if fastaList :
    sys.stderr.write('Created %s sequences from %s graphs in the following files:\n' % (commaFormat(seqCount), commaFormat(graphCount)))
    sys.stderr.write('  %s\n' % '\n  '.join(fastaList))
else :
    sys.stderr.write('Created %s sequences from %s graphs\n' % (commaFormat(seqCount), commaFormat(graphCount)))

if tooLong > 0 :
    sys.stderr.write('There were %s genes with more than %s paths.\n' % (commaFormat(tooLong), commaFormat(opts.seqlimit)))
sys.stderr.write('\nFinished.\n')
