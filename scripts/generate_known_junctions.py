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
Program for generating known splice junction sequences.
"""
from SpliceGrapher.shared.config       import *
from SpliceGrapher.shared.utils        import *
from SpliceGrapher.predict.SpliceSite  import *
from SpliceGrapher.formats.loader      import *
from SpliceGrapher                     import SpliceGraph
from SpliceGrapher.formats.FastaLoader import FastaLoader
from SpliceGrapher.formats.fasta       import FastaRecord

from glob     import glob
from optparse import OptionParser
import os, sys, warnings

# Window around splice sites:
WINDOW_SIZE = 100

def makeHeader(chrom, name, strand, ctr, don, acc, code) :
    """Creates a FASTA header for a splice junction sequence."""
    don3prime = don   if strand == '+' else don+2
    acc5prime = acc+2 if strand == '+' else acc
    don5prime = don-opts.window+1 if strand == '+' else don+opts.window+1
    acc3prime = acc+opts.window+2 if strand == '+' else acc-opts.window+1
    return '%s;%d;%s;%d;%d;%d;%d;%s;%s' % (chrom, ctr, name, don5prime, don3prime, acc5prime, acc3prime, strand, code)

class InferredIntron(object) :
    """
    Simple class that stores information about an inferred intron.
    """
    def __init__(self, chr, strand, donPos, accPos) :
        self.donPos = donPos
        self.accPos = accPos
        self.chr    = chr
        self.strand = strand

    def acceptor(self) :
        return self.accPos

    def donor(self) :
        return self.donPos

    def __eq__(self, other) :
        return self.chr == other.chr and self.strand == other.strand and self.donPos == other.donPos and self.accPos == other.accPos

    def __len__(self) :
        return abs(self.accPos-self.donPos)+1

    def __hash__(self) :
        return self.__str__().__hash__()

    def __str__(self) :
        return "%s,%s,%d,%d" % (self.chr, self.strand, self.donPos, self.accPos)

def getSequence(chr, pos1, pos2, strand) :
    # Used with SpliceSite.getSpliceSiteParts()
    return loader.subsequence(chr, pos1, pos2, reverse=(strand=='-'))

def getGeneIntrons(gene, min_intron=4) :
    """Uses a gene's isoforms to infer all the introns in a gene."""
    chr    = gene.chromosome
    strand = gene.strand
    result = set([])

    for isoform in gene.isoforms.values() :
        exons = isoform.sortedExons() # exons sorted according to strand
        for i in range(1,len(exons)) :
            intron = InferredIntron(chr, strand, exons[i-1].donor(), exons[i].acceptor())
            if len(intron) >= min_intron : result.add(intron)

    for mrna in gene.mrna.values() :
        cds = mrna.sortedExons() # exons inferred from CDS and sorted according to strand
        for i in range(1,len(cds)) :
            intron = InferredIntron(chr, strand, cds[i-1].donor(), cds[i].acceptor())
            if len(intron) >= min_intron : result.add(intron)

    return result

def getGraphIntrons(graph, min_intron=4) :
    """Uses a splice graph to infer all the introns in a gene."""
    chr    = graph.chromosome
    strand = graph.strand
    result = set([])
    for node in graph.nodeDict.values() :
        for c in node.children :
            intron = InferredIntron(chr, strand, SpliceGraph.donor(node), SpliceGraph.acceptor(c))
            if len(intron) >= min_intron : result.add(intron)
    return result

USAGE="""%prog [options]

Generates sequences for known splice junctions."""

parser = OptionParser(usage=USAGE)
parser.add_option('-C', dest='canonical',  default=False,         help='Canonical GT/AG and GC/AG sites only [default: %default]', action='store_true')
parser.add_option('-D', dest='add_dimers', default=False,         help='Include dimer in sequence output [default: %default]', action='store_true')
parser.add_option('-i', dest='minintron',  default=4,             help='Minimum allowed intron size [default: %default]', type='int')
parser.add_option('-f', dest='fasta',      default=SG_FASTA_REF,  help='FASTA reference file [default: %default]')
parser.add_option('-m', dest='model',      default=SG_GENE_MODEL, help='Gene model GFF file [default: %default]')
parser.add_option('-o', dest='outfile',    default=None,          help='Output file [default: %default]')
parser.add_option('-S', dest='sg_list',    default=None,          help='List of splice graph files to augment gene model [default: %default]')
parser.add_option('-P', dest='predfile',   default=None,          help='Store splice sites using prediction file format [default: %default]')
parser.add_option('-v', dest='verbose',    default=False,         help='Verbose mode [default: %default]', action='store_true')
parser.add_option('-w', dest='window',     default=WINDOW_SIZE,   help='Length for exon portions around splice site [default: %default]', type='int')

#-------------------------------------------------------
# Main program
#-------------------------------------------------------
opts, args = parser.parse_args(sys.argv[1:])

errStrings = []
if not (opts.model or opts.sg_list) :
    errStrings.append('** No GFF gene model specified.  Use the -m option or set SG_GENE_MODEL in your SpliceGrapher configuration.')
if not opts.fasta :
    errStrings.append('** No FASTA reference specified.  Use the -f option or set SG_FASTA_REF in your SpliceGrapher configuration.')
if errStrings :
    parser.print_help()
    sys.stderr.write('\n%s\n' % '\n'.join(errStrings))
    sys.exit(1)

if opts.model : validateFile(opts.model)
if opts.sg_list : validateFile(opts.sg_list)
validateFile(opts.fasta)
writeStartupMessage()

#-------------------------------------------------------
# Load reference data:
geneIds    = []
geneChrSet = set([])
if opts.model :
    geneModel  = loadGeneModels(opts.model, verbose=opts.verbose)
    geneIds    = geneModel.getAllGeneIds(geneFilter=gene_type_filter)
    geneChrSet = set([chrom.lower() for chrom in geneModel.getChromosomes()])

spliceGraphs = {}
if opts.sg_list :
    ctr=0
    for line in ezopen(opts.sg_list) :
        f = line.strip()
        graph = SpliceGraph.getFirstGraph(f)
        spliceGraphs[graph.getName()] = graph
        if not opts.model :
            geneIds.append(graph.getName())
            geneChrSet.add(graph.chromosome.lower())

geneIds.sort()

loader = FastaLoader(opts.fasta, verbose=opts.verbose)

# Validate chromosome names
fastaChrSet = set([k.lower() for k in loader.keys()])
commonChr   = geneChrSet.intersection(fastaChrSet)
if not commonChr :
    sys.stderr.write('No chromosomes from GFF3 gene model were found in FASTA reference:\n')
    if opts.model :
        sys.stderr.write('  Gene model chromosomes: %s\n' % ','.join(geneChrSet))
    else :
        sys.stderr.write('  Graph chromosomes: %s\n' % ','.join(geneChrSet))
    sys.stderr.write('  FASTA sequence names:   %s\n' % ','.join(loader.keys()))
    raise ValueError('GFF3 annotations incompatible with FASTA reference.\n')

#-------------------------------------------------------
# Generate splice junction sequences
outStream   = open(opts.outfile,'w') if opts.outfile else sys.stdout
if opts.predfile : spliceStream = open(opts.predfile,'w')
knownRecs   = 0
unknownRecs = 0
chrCounter  = {}
indicator   = ProgressIndicator(10000, description=' genes', verbose=opts.verbose)
for gid in geneIds :
    indicator.update()

    introns = set([])
    if opts.model :
        gene   = geneModel.getGeneByName(gid)
        chrom  = gene.chromosome.lower()
        strand = gene.strand
        name   = gene.id
        # Use introns as references:
        introns   = getGeneIntrons(gene, min_intron=opts.minintron)
    
    # Augment known sites with splice graphs (if any)
    try :
        sg       = spliceGraphs[gid]
        if not opts.model :
            chrom  = sg.chromosome.lower()
            strand = sg.strand
            name   = sg.getName()
        introns |= getGraphIntrons(sg, min_intron=opts.minintron)
    except KeyError :
        pass # Not all genes may have graphs
    
    chrCounter.setdefault(chrom,0)

    # Known junctions
    knownJct = {}
    for x in introns :
        try :
            (donExon, ignore, d) = getSpliceSiteParts(chrom, x.donor(), DONOR_SITE, strand, getSequence, exonWindow=opts.window, intronWindow=2)
            (accExon, ignore, a) = getSpliceSiteParts(chrom, x.acceptor(), ACCEPTOR_SITE, strand, getSequence, exonWindow=opts.window, intronWindow=2)

            if opts.canonical and not (d in KNOWN_DIMERS[DONOR_SITE] and a in KNOWN_DIMERS[ACCEPTOR_SITE]) : continue

            chrCounter[chrom] += 1
            header = makeHeader(chrom, name, strand, chrCounter[chrom], x.donor(), x.acceptor(), KNOWN_JCT)

            knownJct[x.donor()] = knownJct.setdefault(x.donor(), {})
            knownJct[x.donor()][x.acceptor()] = True

            if opts.add_dimers :
                seq = '%s|%s..%s|%s' % (donExon, d, a, accExon)
            else :
                seq = '%s%s' % (donExon, accExon)
            outStream.write(str(FastaRecord(header,seq)))
            if opts.predfile :
                spliceStream.write('%s\t%s\t%d\t%s\t1.0\td\n' % (chrom, strand, x.donor(), d.lower()))
                spliceStream.write('%s\t%s\t%d\t%s\t1.0\ta\n' % (chrom, strand, x.acceptor(), a.lower()))
            knownRecs += 1
        except KeyError :
            continue

    # Putative junctions from known sites
    donors = [x.donor() for x in introns]
    for donor in donors :
        # Only consider acceptors downstream from donor:
        if strand == '+' :
            acceptors = [x.acceptor() for x in introns if x.acceptor() > donor]
        else :
            acceptors = [x.acceptor() for x in introns if x.acceptor() < donor]
        if not acceptors : continue

        try :
            (donExon, ignore, d) = getSpliceSiteParts(chrom, donor, DONOR_SITE, strand, getSequence, exonWindow=opts.window, intronWindow=2)
            if opts.canonical and not (d in KNOWN_DIMERS[DONOR_SITE]) : continue
        except KeyError :
            continue

        for acceptor in acceptors :
            try :
                ignore = knownJct[donor][acceptor]
                continue
            except KeyError :
                pass

            try :
                (accExon, ignore, a) = getSpliceSiteParts(chrom, acceptor, ACCEPTOR_SITE, strand, getSequence, exonWindow=opts.window, intronWindow=2)
                if opts.canonical and not (a in KNOWN_DIMERS[ACCEPTOR_SITE]) : continue
                chrCounter[chrom] += 1
                header = makeHeader(chrom, name, strand, chrCounter[chrom], donor, acceptor, UNKNOWN_JCT)
                if opts.add_dimers :
                    seq = '%s|%s..%s|%s' % (donExon, d, a, accExon)
                else :
                    seq = '%s%s' % (donExon, accExon)
                outStream.write(str(FastaRecord(header,seq)))
                if opts.predfile :
                    spliceStream.write('%s\t%s\t%d\t%s\t1.0\td\n' % (chrom, strand, donor, d.lower()))
                    spliceStream.write('%s\t%s\t%d\t%s\t1.0\ta\n' % (chrom, strand, acceptor, a.lower()))
                unknownRecs += 1
            except KeyError :
                continue

indicator.finish()
sys.stderr.write("Wrote %d known and %d recombined junctions using known splice sites.\n" % (knownRecs, unknownRecs))
