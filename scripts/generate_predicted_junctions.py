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
Python script that generates splice-junction sequences for short-read
alignments based on predicted splice sites.  Output includes combinations
of known splice sites with predicted sites.
"""
from SpliceGrapher.shared.config          import *
from SpliceGrapher.shared.utils           import *
from SpliceGrapher.formats.loader         import *
from SpliceGrapher.formats.FastaLoader    import FastaLoader
from SpliceGrapher.formats.fasta          import FastaRecord
from SpliceGrapher.predict.PredictedSites import *
from SpliceGrapher.predict.SpliceSite     import *

from optparse import OptionParser
import gzip, os, sys

DEFAULT_ACCEPTORS = 'ag'
DEFAULT_DONORS    = 'gt,gc'
DEFAULT_WINDOW    = 30

HEADER_FORMAT     = "%s;%d;%s;%d;%d;%d;%d;%s;%s"

def getSequence(chrom, pos1, pos2, strand) :
    # NOTE: relies on global 'seqDict' FastaLoader instance
    return seqDict.subsequence(chrom, pos1, pos2, reverse=(strand=='-'))

USAGE="""Usage: %prog predictions-file [options]

Python script that uses splice site predictions to generate splice-junction
sequences for short-read alignments.  Outputs only novel splice junctions
where at least one of the two splice sites is not previously documented.

Where: 'predictions-file' contains predicted splice sites using the format:
chromosome	strand	position	dimer	score	site-code"""

parser = OptionParser(usage=USAGE)
parser.add_option('-a', dest='acceptors', default=DEFAULT_ACCEPTORS, help='Comma-separated list of acceptor dimers [default: %default]')
parser.add_option('-d', dest='donors',    default=DEFAULT_DONORS,    help='Comma-separated list of donor dimers [default: %default]')
parser.add_option('-f', dest='fasta',     default=SG_FASTA_REF,      help='Fasta reference file [default: %default]')
parser.add_option('-g', dest='gene',      default=None,              help='Generate sequences only for the given gene [default: %default]')
parser.add_option('-m', dest='model',     default=SG_GENE_MODEL,     help='Gene model GFF3 file [default: %default]')
parser.add_option('-o', dest='output',    default=None,              help='FASTA output file [default: stdout]')
parser.add_option('-W', dest='window',    default=DEFAULT_WINDOW,    help='Sequence window either side of splice site [default: %default]', type='int')
parser.add_option('-v', dest='verbose',   default=False,             help='Verbose mode [default: %default]', action='store_true')

opts, args = parser.parse_args(sys.argv[1:])

if len(args) < 1 :
    parser.print_help()
    sys.exit(1)

validateFile(opts.model)
validateFile(opts.fasta)
writeStartupMessage()

#---------------------------------------------------------------
# Required args
predictionsFile = args[0]

#---------------------------------------------------------------
# Optional args
# Establish list of dimers to process
acceptorList = [s.upper() for s in opts.acceptors.split(',')]
donorList    = [s.upper() for s in opts.donors.split(',')]

#---------------------------------------------------------------
# Load the gene model
geneModel = loadGeneModels(opts.model, verbose=opts.verbose)

#---------------------------------------------------------------
# Load predicted sites
if opts.verbose : sys.stderr.write('Loading predicted sites from %s\n' % predictionsFile)
predictions = ChromosomeSiteMap(predictionsFile, verbose=opts.verbose)
chrom       = predictions.chromosome

#---------------------------------------------------------------
# Find the specific gene
specificGene = None
if opts.gene :
    specificGene = geneModel.getGeneByName(opts.gene)
    if not specificGene :
        raise Exception("Gene '%s' was not found" % opts.gene)
    geneList = [specificGene]
else :
    geneIds  = sorted(geneModel.getGenes(chrom))
    geneList = [geneModel.getGene(chrom,id) for id in geneIds]

#---------------------------------------------------------------
# Load FASTA sequences
if opts.verbose : sys.stderr.write('Loading FASTA sequences from %s\n' % opts.fasta)
seqDict = FastaLoader(opts.fasta, verbose=opts.verbose)

#---------------------------------------------------------------
# Start output stream
outStream = sys.stdout
if opts.output :
    ext = opts.output[opts.output.rfind('.'):]
    if '.gzip'.find(ext.lower()) >= 0 :
        outStream = gzip.open(opts.output, 'w')
    else :
        outStream = file(opts.output, 'w')

#---------------------------------------------------------------
# Generate sequences
seqID = 0
for gene in geneList :
    reversed = (gene.strand == '-')

    # Get predicted donor and acceptor sites
    predictedAcceptors = set(predictions.getSites(ACCEPTOR_SITE, gene.strand, minpos=gene.minpos, maxpos=gene.maxpos))
    predictedDonors    = set(predictions.getSites(DONOR_SITE, gene.strand, minpos=gene.minpos, maxpos=gene.maxpos))

    # Get previously-known donors and acceptors
    knownAcceptors = set(gene.acceptorList())
    knownDonors    = set(gene.donorList())

    # Ensure predicted sites and known sites have nothing in common
    accIntersection = predictedAcceptors.intersection(knownAcceptors)
    donIntersection = predictedDonors.intersection(knownDonors)

    # Distinguish between novel sites and known sites
    novelDonors    = predictedDonors - knownDonors
    novelAcceptors = predictedAcceptors - knownAcceptors

    # Add novel sites to known sites
    acceptorSet = knownAcceptors.union(predictedAcceptors)
    donorSet    = knownDonors.union(predictedDonors)

    # Sets are unordered: change to sorted lists
    donors    = list(donorSet)
    acceptors = list(acceptorSet)
    donors.sort(reverse=reversed)
    acceptors.sort(reverse=reversed)
    
    knownJunctions = gene.getJunctions()
    geneSequences  = 0
    for d in donors :
        (exon1, intron1, dimer1) = getSpliceSiteParts(chrom, d, DONOR_SITE, gene.strand, getSequence,
                                                      exonWindow=opts.window,
                                                      intronWindow=opts.window)
        # Ensure donor site is legal
        if dimer1 not in donorList :
            if d not in knownDonors : # Should never happen, given nature of prediction software:
                sys.stderr.write('Illegal PREDICTED donor site at %s %d %s: %s|%s|%s\n' % (chrom, d, gene.strand, exon1, dimer1, intron1))
            continue

        for a in acceptors :
            # These pairs are already taken care of in initial FASTA reference
            if (d in knownDonors) and (a in knownAcceptors) : continue

            # Ensure acceptor is downstream of donor
            if reversed :
                if a >= d : continue
            else :
                if a <= d : continue

            (exon2, intron2, dimer2) = getSpliceSiteParts(chrom, a, ACCEPTOR_SITE, gene.strand, getSequence,
                                                          exonWindow=opts.window,
                                                          intronWindow=opts.window)

            # Ensure acceptor site is legal
            if dimer2 not in acceptorList :
                if a not in knownAcceptors :
                    sys.stderr.write('Illegal PREDICTED acceptor site at %s %d %s: %s|%s|%s\n' % (chrom, a, gene.strand, intron2, dimer2, exon2))
                continue

            # Header should always be donor 5'/donor 3'/acceptor 5'/acceptor 3' order
            seqID += 1
            if gene.strand == '+' :
                ## header = HEADER_FORMAT % (chrom.capitalize(), seqID, gene.id, d-opts.window+2, d+1, a, a+opts.window-1, gene.strand, PREDICTED_JCT)
                header = HEADER_FORMAT % (chrom.capitalize(), seqID, gene.id, d-opts.window+1, d, a+2, a+opts.window+1, gene.strand, PREDICTED_JCT)
            else : 
                header = HEADER_FORMAT % (chrom.capitalize(), seqID, gene.id, d+opts.window, d+2, a, a-opts.window+1, gene.strand, PREDICTED_JCT)

            seq = exon1 + exon2
            outStream.write(str(FastaRecord(header,seq)))
            geneSequences += 1
            # end for a in acceptors
        # end for d in donors

    if opts.verbose :
        sys.stderr.write('Gene %s: %d known donors/%d known acceptors, %d predicted donors/%d predicted acceptors, %d sequences written\n' % \
            (gene.id, len(knownDonors), len(knownAcceptors), len(predictedDonors), len(predictedAcceptors), geneSequences))
