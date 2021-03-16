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
Module that contains data structures and methods related to splice sites.
"""
ACCEPTOR_SITE     = 'acceptor'
DONOR_SITE        = 'donor'
SITE_TYPES        = [ACCEPTOR_SITE, DONOR_SITE]

ACCEPTOR_CODE     = 'a'
DONOR_CODE        = 'd'
SITE_CODES        = [ACCEPTOR_CODE, DONOR_CODE]
SITE_TYPE_NAME    = {ACCEPTOR_CODE:ACCEPTOR_SITE, DONOR_CODE:DONOR_SITE}

CANONICAL_ACCEPTOR = 'AG'
CANONICAL_DONOR    = 'GT'
NONCANONICAL_DONOR = 'GC'

# Strictly canonical dimers; this dict is designed to
# be interchangeable with KNOWN_DIMER (see below)
CANONICAL_DIMERS = {ACCEPTOR_SITE:[CANONICAL_ACCEPTOR], DONOR_SITE:[CANONICAL_DONOR]}

# Add others as they become known:
KNOWN_DIMERS     = {ACCEPTOR_SITE:[CANONICAL_ACCEPTOR], DONOR_SITE:[CANONICAL_DONOR, NONCANONICAL_DONOR]}

# Specific dimers to pick up for non-canonical sites
NC_DONORS        = ['AT']
NC_ACCEPTORS     = ['AC']
NC_DIMERS        = {ACCEPTOR_SITE:NC_ACCEPTORS, DONOR_SITE:NC_DONORS}

DEFAULT_EXON_WINDOW   = 100
DEFAULT_INTRON_WINDOW = 100

# Splice junction codes:
KNOWN_JCT     = 'K'
UNKNOWN_JCT   = 'U'
PREDICTED_JCT = 'P'

#
#  Splice site methods
#
def canonicalAcceptor(seq, dimers=CANONICAL_DIMERS) :
    """
    Returns true if the sequence has a canonical acceptor site (i.e., AG dimer); false otherwise.
    """
    return seq[-2:] in dimers[ACCEPTOR_SITE]

def canonicalDonor(seq, dimers=CANONICAL_DIMERS) :
    """
    Returns true if the sequence has a canonical donor site (i.e., GT dimer); false otherwise.
    """
    return seq[:2] in dimers[DONOR_SITE]

def canonicalIntron(seq, dimers=CANONICAL_DIMERS) :
    """
    Returns true if the sequence represents a canonical intron (i.e., with GT/AG dimers); false otherwise.
    """
    return canonicalDonor(seq, dimers) and canonicalAcceptor(seq, dimers)

# Example of getSequenceFunction that works with getSpliceSiteParts:
# NOTE: relies on global FastaLoader instance
# def getSequence(chrom, pos1, pos2, strand) :
#     return fastaLoader.subsequence(chr, pos1, pos2, reverse=(strand=='-'))
def getSpliceSiteParts(chrom, spliceSite, siteType, strand, getSequenceFunction, exonWindow=DEFAULT_EXON_WINDOW, intronWindow=DEFAULT_INTRON_WINDOW, verbose=False) :
    """
    Returns the exon sequence and intron sequence flanking a splice site
    along with the splice site dimer (usually 'GT' or 'AG').

    Parameters:
        chrom               - chromosome where sequence resides
        spliceSite          - position of the splice site within a genomic sequence
        siteType            - 'donor' or 'acceptor'
        strand              - '+' or '-'
        getSequenceFunction - function that returns a sequence given the following parameters:
                              (chromosome, fromPos, toPos, strand)
        exonWindow          - size of exon portion to return (default=100nt)
        intronWindow        - size of exon portion to return (default=100nt)
    """
    # Permit a range of possible donor site identifiers:
    ssType  = siteType.lower()
    isDonor = (ssType in ['gt','gc'] or DONOR_SITE.startswith(ssType))

    # Obtain the full-length sequence
    if (strand == '+' and isDonor) or (strand == '-' and not isDonor) :
        fullSeq = getSequenceFunction(chrom, spliceSite-exonWindow, spliceSite+intronWindow+1, strand)
    else : # + acceptor/- donor
        fullSeq = getSequenceFunction(chrom, spliceSite-intronWindow-1, spliceSite+exonWindow, strand)

    requiredLen = intronWindow+exonWindow+2

    # Locate the splice site within the sequence
    splicePos = exonWindow if isDonor else intronWindow

    dimer = fullSeq[splicePos:splicePos+2]

    if isDonor :
        exonPart   = fullSeq[:splicePos]
        intronPart = fullSeq[splicePos+2:]
    else : # acceptor
        intronPart = fullSeq[:splicePos]
        exonPart   = fullSeq[splicePos+2:]

    return (exonPart, intronPart, dimer)

def truncateSequences(fastaFile, exonSize, intronSize, outFile, acceptor=False, verbose=False) :
    """
    Method that truncates FASTA sequences based on given intron and exon sizes.
    """
    import sys
    from SpliceGrapher.formats.fasta import fasta_itr, FastaRecord
    fiter     = fasta_itr(fastaFile)
    outStream = file(outFile, 'w')

    pos = 0
    neg = 0
    tot = 0
    for rec in fiter :
        tot  += 1
        midpt = (len(rec.sequence)/2)
        seq   = rec.sequence[midpt-intronSize:midpt+exonSize] if acceptor else rec.sequence[midpt-exonSize:midpt+intronSize]
        if rec.header.find("label=0") > 0 : neg += 1
        if rec.header.find("label=1") > 0 : pos += 1
        rec   = FastaRecord(rec.header, seq)
        outStream.write(str(rec))
    if verbose : sys.stderr.write('truncateSequences found %d positive/%d negative in %d total sequences from %s\n' % (pos,neg,tot,fastaFile))

def validIntron(seq) :
    """
    Returns true if the sequence represents a valid intron (i.e., with known dimers); false otherwise.
    """
    return (seq[:2] in KNOWN_DIMERS[DONOR_SITE] and seq[-2:] in KNOWN_DIMERS[ACCEPTOR_SITE])
