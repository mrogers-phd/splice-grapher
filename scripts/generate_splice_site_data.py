#! /usr/bin/env python
# 
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
"""
Program for generating FASTA examples for a given splice site type.
"""
from SpliceGrapher.shared.GeneModelConverter import *
from SpliceGrapher.shared.config       import *
from SpliceGrapher.shared.utils        import *
from SpliceGrapher.predict.SpliceSite  import *
from SpliceGrapher.formats.fasta       import FastaRecord
from SpliceGrapher.formats.FastaLoader import FastaLoader
from SpliceGrapher.formats.loader      import *
from SpliceGrapher.SpliceGraph         import *

from optparse import OptionParser, OptionGroup
from glob     import glob
from sys      import maxint as MAXINT
import os, random, sys, warnings

def trainingFilter(g) :
    """Filter for accepting/rejecting training examples."""
    for k in TRAIN_ACCEPT.keys() :
        if k not in g.attributes : return False
        target = TRAIN_REJECT[k].lower()
        value  = g.attributes[k].lower()
        if value.find(target) < 0 : return False

    for k in TRAIN_REJECT.keys() :
        if k not in g.attributes : continue
        target = TRAIN_REJECT[k].lower()
        value  = g.attributes[k].lower()
        if value.find(target) >= 0 : return False

    return gene_type_filter(g)

def getSequence(chrom, pos1, pos2, strand) :
    # To be used with SpliceSite.getSpliceSiteParts()
    return fastaDB.subsequence(chrom, pos1, pos2, reverse=(strand=='-'))

def randomWriteStringBuffer(stringBuffer, outStream, limit=MAXINT, verbose=False) :
    """Writes a random subset of the given string list to the output stream provided.
    Returns the number of strings written."""
    if verbose : sys.stderr.write('Shuffling %d FASTA records in buffer\n' % (len(stringBuffer)))

    permutedBuffer = random.sample(stringBuffer, limit) if limit < len(stringBuffer) else stringBuffer

    if verbose :
        sys.stderr.write('Writing buffer\n')
        indLimit  = max(100, min(100000, len(permutedBuffer)))
        indicator = ProgressIndicator(indLimit)

    for s in permutedBuffer :
        if verbose : indicator.update()
        outStream.write(s)
    if verbose : indicator.finish()

    return len(permutedBuffer)

def writeReport(options, jctDict) :
    """Writes a report detailing statistics on the dimers encountered
    for donor and acceptor sites along with statistics on junction pairs."""
    donorDict    = {}
    acceptorDict = {}
    total        = 0
    for d in jctDict :
        donorDict[d] = 0
        for a in jctDict[d] :
            acceptorDict[a] = acceptorDict.setdefault(a,0) + jctDict[d][a]
            donorDict[d]   += jctDict[d][a]
            total          += jctDict[d][a]

    # First part of report: dimer frequency
    outStream = open(opts.report, 'w')
    outStream.write('Breakdown of donor sites for %d introns:\n' % total)
    pairs = list(donorDict.items())
    pairs.sort(cmp=lambda a,b : b[1]-a[1])
    for p in pairs :
        pct = (100.0*p[1])/total
        outStream.write(' %s: %7d (%5.2f%%)\n' % (p[0],p[1],pct))

    outStream.write('Breakdown of acceptor sites for %d introns:\n' % total)
    pairs = list(acceptorDict.items())
    pairs.sort(cmp=lambda a,b : b[1]-a[1])
    for p in pairs :
        pct = (100.0*p[1])/total
        outStream.write(' %s: %7d (%5.2f%%)\n' % (p[0],p[1],pct))

    # Second part of report: splice junction frequency
    totalJct = 0
    for d in jctDict.keys() :
        for a in jctDict[d].keys() :
            totalJct += jctDict[d][a]

    jctRecs = []
    for d in jctDict.keys() :
        for a in jctDict[d].keys() :
            pct = (100.0*jctDict[d][a])/totalJct
            jctRecs.append((d, a, jctDict[d][a], pct))

    # Sort in decreasing order of frequency
    jctRecs.sort(cmp=lambda a,b : b[2]-a[2])
    outStream.write('Breakdown of splice junctions:\n')
    for rec in jctRecs :
        outStream.write(" %s-%s: %7d (%7.4f%%)\n" % rec)

# How far to extend sequences around splice sites:
WINDOW_SIZE = 50

USAGE="""%prog [options]"""

DESCR="""Generates positive or negative examples of splice sites based on
an organism's GFF gene model annotations or a set of splice graphs, and its
FASTA reference sequences.  Default splice site type is 'donor' and the
default dimer is 'GT'."""

parser = OptionParser(usage=USAGE, description=DESCR)
required = OptionGroup(parser, 'Required', "Note: use a gene model (-m), a list of splice graphs (-S), or both, but at least one is required.  A FASTA reference is always required.")
required.add_option('-f', dest='fasta',        default=SG_FASTA_REF,  help='(Required) FASTA reference file [default: %default]')
required.add_option('-m', dest='model',        default=SG_GENE_MODEL, help='GFF gene model file [default: %default]')
required.add_option('-S', dest='splicegraphs', default=None,          help='File containing a list of splice graph file paths [default: %default]')
parser.add_option_group(required)

optional = OptionGroup(parser, 'Other')
optional.add_option('-a', dest='acceptor',     default=False,         help='Look for dimer in acceptor sites [default: %default]', action='store_true')
optional.add_option('-D', dest='add_dimer',    default=False,         help='Include dimer in sequence output [default: %default]', action='store_true')
optional.add_option('-d', dest='dimers',       default='GT',          help='Dimers to look for [default: %default]')
optional.add_option('-i', dest='minintron',    default=4,             help='Minimum allowed intron length [default: all]', type='int')
optional.add_option('-n', dest='limit',        default=0,             help='Generate given number of examples [default: all]', type='int')
optional.add_option('-W', dest='window',       default=WINDOW_SIZE,   help='Set window size on either side of splice site [default: %default]', type='int')
optional.add_option('-o', dest='outfile',      default=None,          help='Output file [default: %default]')
optional.add_option('-r', dest='report',       default=None ,         help='Produce splice site frequency report [default: %default]')
optional.add_option('-N', dest='negative',     default=False,         help='Generate negative examples [default: %default]', action='store_true')
optional.add_option('-v', dest='verbose',      default=False,         help='Verbose mode [default: %default]', action='store_true')
parser.add_option_group(optional)

#-------------------------------------------------------
# Main program
#-------------------------------------------------------
opts, args = parser.parse_args(sys.argv[1:])

errStrings = []
if not (opts.model or opts.splicegraphs) :
    errStrings.append('** No GFF gene model or splice graphs specified.  Use the -m or -S option, or set SG_GENE_MODEL in your SpliceGrapher configuration.')
if not opts.fasta :
    errStrings.append('** No FASTA reference specified.  Use the -f option or set SG_FASTA_REF in your SpliceGrapher configuration.')
if errStrings :
    parser.print_help()
    sys.stderr.write('\n%s\n' % '\n'.join(errStrings))
    sys.exit(1)

validateFile(opts.fasta)
if opts.model : validateFile(opts.model)
if opts.splicegraphs : validateFile(opts.splicegraphs)

writeStartupMessage()

if opts.verbose and TRAIN_ACCEPT or TRAIN_REJECT :
    sys.stderr.write('Using accept/reject filters from config.py to select genes.\n')

#-------------------------------------------------------
# Interpret splice site options:
siteType  = 'acceptor' if opts.acceptor else 'donor'
dimerList = [s.upper() for s in opts.dimers.split(',')]
#-------------------------------------------------------

seqLimit  = opts.limit if opts.limit > 0 else MAXINT
# Negative examples are selected at random and can be written
# on the fly; positive examples are buffered first, then written
# in random order
useBuffer = opts.limit and not opts.negative

if opts.negative :
    classLabel = 0
    classType  = 'negative'
else :
    classLabel = 1
    classType  = 'positive'

#-------------------------------------------------------
# Load reference data:
fastaDB = FastaLoader(opts.fasta, verbose=opts.verbose)

graphDict = {}
badChrom  = set([])
if opts.model :
    geneModel = loadGeneModels(opts.model, verbose=opts.verbose)
    genes     = geneModel.getAllGenes(trainingFilter)
    for g in genes :
        if g.chromosome.lower() not in fastaDB.keys() :
            if opts.verbose and g.chromosome.lower() not in badChrom :
                sys.stderr.write('Ignoring gene %s: chromosome %s not in FASTA reference\n' % (g.id, g.chromosome))
                badChrom.add(g.chromosome.lower())
            continue

        try :
            graph = geneModelToSpliceGraph(g, useCDS=(len(g.isoforms)==0), minintron=opts.minintron)
        except ValueError :
            if opts.verbose :
                sys.stderr.write('Warning: unable to create splice graph for %s: no isoforms found\n' % g.id)
            continue

        key   = graph.getName().upper()
        # Initialize before adding any auxiliary graphs
        graphDict[key] = [graph]

    if opts.verbose : sys.stderr.write('Retained %d genes for training\n' % len(graphDict))

if opts.splicegraphs :
    sys.stderr.write("Loading splice graphs listed in %s\n" % opts.splicegraphs)
    indicator = ProgressIndicator(10000, verbose=opts.verbose)
    ctr = 0
    for f in ezopen(opts.splicegraphs) :
        indicator.update()
        try :
            graph = getFirstGraph(f.strip())
        except ValueError as ve :
            sys.stderr.write('Warning: %s\n' % ve)
            continue
        if graph.isEmpty() : continue

        if graph.chromosome.lower() not in fastaDB.keys() :
            if opts.verbose and graph.chromosome.lower() not in badChrom :
                sys.stderr.write('Ignoring gene %s: chromosome %s not in FASTA reference\n' % (graph.getName(), graph.chromosome))
                badChrom.add(graph.chromosome.lower())
            continue

        key = graph.getName().upper()
        if opts.model :
            try :
                graphDict[key].append(graph)
            except KeyError :
                sys.stderr.write('** Warning: skipping auxiliary graph %s: does not match any gene model entry\n' % graph.getName())
                continue
        else :
            graphDict.setdefault(key,[])
            graphDict[key].append(graph)
        ctr += 1

    indicator.finish()
    sys.stderr.write("Loaded %d auxiliary splice graphs\n" % ctr)

if not graphDict : raise Exception('No graphs loaded; cannot continue.')
#-------------------------------------------------------

outStream = open(opts.outfile,'w') if opts.outfile else sys.stdout

if opts.verbose : sys.stderr.write('Generating %s examples of %s sequences for %s dimers:\n' % (classType, siteType, opts.dimers))

# Iterator is just the graph list for positive examples,
# or a random iterator over negative examples
graphIter = RandomListIterator(graphDict.keys()) if opts.negative else graphDict.keys()

seqCtr    = 0
geneCount = 0
jctDict   = {}
recBuffer = []
geneSet   = set([])

indicator = ProgressIndicator(100000, verbose=opts.verbose)
for gid in graphIter :
    geneCount += 1
    geneSet.add(gid)

    knownAcceptors = set([])
    knownDonors    = set([])
    edges          = set([])
    for g in graphDict[gid] :
        knownAcceptors.update(g.getAcceptors())
        knownDonors.update(g.getDonors())
        edges.update(edgeSet(g))

    # Use first graph in list as reference
    graph = graphDict[gid][0]

    if opts.report :
        for x in edges :
            (estr, istr, dondimer) = getSpliceSiteParts(graph.chromosome, donor(x.parent), 'gt', graph.strand, getSequence, exonWindow=2, intronWindow=2)
            (estr, istr, accdimer) = getSpliceSiteParts(graph.chromosome, acceptor(x.child), 'acc', graph.strand, getSequence, exonWindow=2, intronWindow=2)
            if not (dondimer and accdimer) :
                seq = fastaDB.sequence(graph.chromosome)
                raise ValueError("Gene %s bad donor '%s' or acceptor '%s' found in edge '%s' on %s (%s %d nt)" % (graph.getName(), dondimer, accdimer, x, graph.chromosome, graph.strand, len(seq)))
            jctDict.setdefault(dondimer,{})
            jctDict[dondimer].setdefault(accdimer,0)
            jctDict[dondimer][accdimer] += 1

    if opts.negative :
        # Randomly permute sites from a list of unknown splice sites
        allSites     = set(range(graph.minpos,graph.maxpos+1))
        unknownSites = allSites-knownAcceptors if opts.acceptor else allSites-knownDonors
        siteList     = random.sample(unknownSites,len(unknownSites))
    else :
        siteList = knownAcceptors if opts.acceptor else knownDonors

    try :
        # adjust for window length + dimer length at upper end
        siteLimit = fastaDB.sequenceLength(graph.chromosome) - opts.window - 2
        siteList  = [x for x in siteList if x < siteLimit]
    except KeyError : # unknown chromosome
        continue

    for pos in siteList :
        try :
            (exonSeq, intronSeq, dimer) = getSpliceSiteParts(graph.chromosome, pos, siteType, graph.strand, getSequence,
                                                             exonWindow=opts.window,
                                                             intronWindow=opts.window)
        except KeyError :
            warnings.warn('No sequence information for chromosome %s' % graph.chromosome)
            continue

        #--------------------------------------
        # Dimer match --> write out the example
        if dimer.upper() in dimerList :
            header = '%s;%s;%s;%s;%d label=%d' % (graph.chromosome, graph.getName(), graph.strand, siteType, pos, classLabel)
            (firstPart, secondPart) = (intronSeq, exonSeq) if opts.acceptor else (exonSeq, intronSeq)
            seq = firstPart + dimer + secondPart if opts.add_dimer else firstPart + secondPart

            fastaRec = FastaRecord(header,seq)
            if useBuffer :
                # NB: strings use less memory than full records
                recBuffer.append(str(fastaRec))
            else :
                outStream.write(str(fastaRec))
            seqCtr += 1
            indicator.update()

            # Take only first example for negative data
            if opts.negative : break
        #-----------------------------------
    if opts.negative and seqCtr >= seqLimit : break
indicator.finish()

if useBuffer and recBuffer :
    seqCtr = randomWriteStringBuffer(recBuffer, outStream, seqLimit, verbose=opts.verbose)
elif useBuffer :
    raise Exception('** No sequences found')

if opts.report :
    writeReport(opts, jctDict)

if opts.verbose :
    if opts.outfile :
        sys.stderr.write('Wrote %d %s sequences from %d genes with %s dimer to %s\n' % (seqCtr, siteType, len(geneSet), opts.dimers, opts.outfile))
    else :
        sys.stderr.write('Wrote %d %s sequences from %d genes with %s dimer\n' % (seqCtr, siteType, len(geneSet), opts.dimers))
