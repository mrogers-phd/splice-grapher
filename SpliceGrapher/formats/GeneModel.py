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
Stores gene annotation information from a GFF3 annotation
file and provides methods for searching on the data.
"""
from SpliceGrapher.shared.utils import ProgressIndicator, ezopen, getAttribute, commaFormat

# TRYING A NEW WAY TO IDENTIFY GENES:
import difflib

import os, sys

# GFF record types
CDS_TYPE        = 'cds'
CHR_TYPE        = 'chromosome'
EXON_TYPE       = 'exon'
FP_UTR_TYPE     = 'five_prime_utr'
GENE_TYPE       = 'gene'
INTRON_TYPE     = 'intron'
MRNA_TYPE       = 'mrna'
MRNA_TE_TYPE    = 'mRNA_TE_gene'
NONUNIQUE_TYPE  = 'nonunique'
PROTEIN_RECORD  = 'protein'
PREDCDS_TYPE    = 'cds_predicted'
PREDGENE_TYPE   = 'predicted_gene'
TP_UTR_TYPE     = 'three_prime_utr'
TRANS_ELE_TYPE  = 'transposable_element_gene'

PSEUDOGENE_TYPE  = 'pseudogene'
PSEUDOTRANS_TYPE = 'pseudogenic_transcript'
PSEUDOEXON_TYPE  = 'pseudogenic_exon'

KNOWN_RECTYPES  = [CDS_TYPE, CHR_TYPE, EXON_TYPE, FP_UTR_TYPE, GENE_TYPE, INTRON_TYPE, MRNA_TYPE, MRNA_TE_TYPE,
                   NONUNIQUE_TYPE, PROTEIN_RECORD, PREDCDS_TYPE, PREDGENE_TYPE,
                   PSEUDOGENE_TYPE, PSEUDOEXON_TYPE, PSEUDOTRANS_TYPE,
                   TP_UTR_TYPE, TRANS_ELE_TYPE]
IGNORE_RECTYPES = [PROTEIN_RECORD, INTRON_TYPE, MRNA_TE_TYPE, TRANS_ELE_TYPE, NONUNIQUE_TYPE]
CDS_TYPES       = [FP_UTR_TYPE, TP_UTR_TYPE, CDS_TYPE]

# Special GTF types for conversions
GTF_GENE_ID    = 'gene_id'
GTF_GENE_NAME  = 'gene_name'
GTF_TRANSCRIPT = 'transcript_id'
GTF_TRANSNAME  = 'transcript_name'
GTF_SOURCE     = 'gene_biotype'
GTF_EXON_ID    = 'exon_number'
GTF_PROTEIN_ID = 'protein_id'

# Record type map allows mapping unusual names to known types:
RECTYPE_MAP                = dict([(s,s) for s in KNOWN_RECTYPES])
RECTYPE_MAP[PREDGENE_TYPE] = GENE_TYPE
RECTYPE_MAP[PREDCDS_TYPE]  = CDS_TYPE

# Virtual types (don't appear in GFF):
ISOFORM_TYPE = 'isoform'

# Annotation fields:
ID_FIELD     = 'ID'
NAME_FIELD   = 'Name'
NOTE_FIELD   = 'Note'
PARENT_FIELD = 'Parent'

# There seems to be no consensus about how these tags are used in UCSC, ENSEMBL
# and other forms of gene models, so we must try each kind:
POSSIBLE_GENE_FIELDS = [PARENT_FIELD, GTF_GENE_ID, GTF_GENE_NAME]
POSSIBLE_FORM_FIELDS = [PARENT_FIELD, GTF_TRANSCRIPT]

# ID for GFF output:
GFF_ID        = 'SpliceGrapher'
MAX_BAD_LINES = 3

# Some files contain '.' for an unassigned strand
VALID_STRANDS = ['-','+','.']

FORM_DELIMITERS = ['.','-','_',',']

# Genes may overlap or even contain one another, so binary search
# is performed only at a coarse resolution to get within this number
# of genes.  Then a linear search is performed with a slight margin.
# This makes it approximately O(lg n) + C
GENE_SEARCH_RANGE  = 16
GENE_SEARCH_MARGIN = 8

def gene_type_filter(g) :
    """Convenience filter for getting only 'gene' records."""
    return g.featureType == GENE_TYPE

def defaultGeneFilter(g) :
    """Default function for filtering genes from a list."""
    return True

def dictToGFF(d) :
    """Returns a string representation of a dictionary based on the GFF3 annotation format."""
    return ';'.join(['%s=%s' % (k,v) for k,v in sorted(d.items()) if k != 'parent'])

def dictToGTF(d) :
    """Returns a string representation of a dictionary based on the GTF annotation format."""
    return '; '.join(['%s "%s"' % (k,v) for k,v in sorted(d.items()) if k != 'parent'])

def cdsFactory(recType, startPos, endPos, chrName, strand, attr={}) :
    """Simple factory method for creating CDS-type records."""
    if recType == CDS_TYPE :
        return CDS(startPos, endPos, chrName, strand, attr)
    elif recType == FP_UTR_TYPE :
        return FP_UTR(startPos, endPos, chrName, strand, attr)
    elif recType == TP_UTR_TYPE :
        return TP_UTR(startPos, endPos, chrName, strand, attr)
    else :
        raise ValueError('Illegal CDS record type: %s' % recType)

class Chromosome(object) :
    """Class that encapsulates a chromosome GFF record."""
    def __init__(self, start, end, name) :
        self.minpos = start
        self.maxpos = end
        self.name   = name

    def __len__(self) :
        return self.maxpos-self.minpos+1

    def __str__(self) :
        return "%s: %d-%d" % (self.name, self.minpos, self.maxpos)

    def contains(self, pos) :
        return self.minpos <= pos <= self.maxpos

    def end(self) :
        return self.maxpos

    def gffString(self) :
        # Example: Chr1    TAIR9   chromosome  1   30427671    .   .   .   ID=Chr1;Name=Chr1
        nameStr = self.name.capitalize()
        return "%s\t%s\tchromosome\t%d\t%d\t.\t.\t.\tID=%s;Name=%s" % (nameStr, GFF_ID, self.start(), self.end(), nameStr, nameStr)

    def start(self) :
        return self.minpos

    def update(self, feature) :
        """
        Many species do not have 'chromosome' entries in their annotations, so we must infer the chromosome
        boundaries from the features found therein.
        """
        if feature.chromosome.lower() != self.name.lower() :
            raise ValueError('Cannot use feature from %s to update %s' % (feature.chromosome, self.name))
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)

def featureCmp(a, b) :
    """
    General comparison function for sorting any features that have 'minpos'
    and 'maxpos' attributes.
    """
    if a.minpos == b.minpos :
        return a.maxpos - b.maxpos
    else :
        return a.minpos - b.minpos
    
def featureOverlaps(a, b) :
    """
    General function for determining whether feature 'a' and feature 'b' overlap.
    """
    if not (a and b) : return False
    return a.minpos <= b.minpos <= a.maxpos or b.minpos <= a.minpos <= b.maxpos

def featureContains(a, b) :
    """
    General function for determining whether feature 'a' contains
    feature 'b'.  Note that both features must have 'minpos'
    and 'maxpos' attributes.
    """
    if not (a and b) : return False
    return a.minpos <= b.minpos and a.maxpos >= b.maxpos

def featureSearch(features, query, lo=0, hi=None) :
    """
    Binary search through a list of sorted features, each of which has
    'minpos' and 'maxpos' attributes.  'query' is an example of a feature
    (such as an exon) to be found in the list.

    Returns either the feature in the list that contains the query feature,
    or the one that would immediately precede it in the list.
    """
    if hi is None :
        hi = len(features)-1

    if lo+1 >= hi :
        if contains(features[lo], query) :
            return features[lo]
        elif contains(features[hi], query) :
            return features[hi]
        elif featureCmp(features[lo], features[hi]) >= 0 :
            return features[lo]
        else :
            return features[hi]
    else :
        midpt = (lo+hi)/2
        if features[midpt].minpos <= query.minpos :
            return featureSearch(features, query, midpt, hi)
        elif features[midpt].maxpos >= query.maxpos :
            return featureSearch(features, query, lo, midpt)
        else : # query contains midpoint
            return featureSearch(features, query, midpt-1, midpt)

class BaseFeature(object) :
    def __init__(self, featureType, start, end, chromosome, strand, attr={}) :
        self.chromosome  = chromosome
        self.strand      = strand
        self.parent      = None
        self.minpos      = min(start,end)
        self.maxpos      = max(start,end)
        self.featureType = featureType
        self.attributes  = attr

    def acceptor(self) :
        """
        Returns the location where the acceptor dimer begins: 2nt
        upstream of the exon start position.  This is 2 before the
        start on the + strand and the exact start on the - strand.
          + Example: AG|CGTATTC
          - Example: GAATACG|CT (reverse-complement)
        """
        return self.start()-2 if self.strand == '+' else self.start()

    def __cmp__(self, other) :
        if self.chromosome != other.chromosome :
            return 2*(self.chromosome > other.chromosome) - 1
        elif self.minpos == other.minpos :
            return self.maxpos - other.maxpos
        else :
            return self.minpos - other.minpos

    def contains(self, pos, strand) :
        return (strand == self.strand and self.minpos <= pos <= self.maxpos)

    def detailString(self) :
        return "Feature: %s\nChromosome: %s\nStart: %d; End: %d; Strand: '%s'" % \
                (self.featureType, self.chromosome, self.start(), self.end(), self.strand)

    def donor(self) :
        """
        Returns the location where the donor dimer begins: at the
        end of the exon.  This is 2 nt before the start on the - strand.
          + strand example: CGTATTC|GT
          - strand example: AC|GAATACG
        """
        return self.end() if self.strand == '+' else self.end()-2

    def end(self) :
        return self.maxpos if self.strand == '+' else self.minpos

    def __eq__(self, other) :
        return self.minpos == other.minpos and self.maxpos == other.maxpos and self.strand == other.strand and self.chromosome == other.chromosome

    def gffString(self, altAttributes={}) :
        """Returns a GFF-formatted representation of the feature."""
        attrs = dict(self.attributes)
        attrs.update(altAttributes)

        return "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s" % \
                (self.chromosome.capitalize(), GFF_ID, self.featureType, self.minpos, self.maxpos, self.strand, dictToGFF(attrs))

    def gtfString(self, transcript, genePtr, exonId) :
        """Returns a GTF-formatted representation of the feature."""
        attrs = dict([(k.lower(),self.attributes[k]) for k in self.attributes])
        attrs[GTF_GENE_ID]    = genePtr.id
        attrs[GTF_GENE_NAME]  = genePtr.name
        attrs[GTF_TRANSCRIPT] = transcript
        attrs[GTF_EXON_ID]    = exonId
        try :
            source = attrs[GTF_SOURCE]
        except KeyError :
            source = GFF_ID
        return '%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s' % \
                (self.chromosome.capitalize(), source, self.featureType, self.minpos, self.maxpos, self.strand, dictToGTF(attrs))

    def __hash__(self) :
        return self.__str__().__hash__()

    def __len__(self) :
        return self.maxpos-self.minpos+1

    def setParent(self, id) :
        self.parent = id

    def start(self) :
        return self.minpos if self.strand == '+' else self.maxpos

    def __str__(self) :
        return "%s %s: %d-%d (%s)" % (self.chromosome, self.featureType, self.start(), self.end(), self.strand)

class Exon(BaseFeature) :
    def __init__(self, start, end, chromosome, strand, attr={}) :
        BaseFeature.__init__(self, EXON_TYPE, start, end, chromosome, strand, attr)
        self.parents = []

    def __str__(self) :
        if self.parents :
            return "%s %d-%d(%s) (isoforms: %s)" % (self.featureType, self.start(), self.end(), self.strand, ','.join([x.id for x in self.parents]))
        else :
            return "%s %d-%d(%s)" % (self.featureType, self.start(), self.end(), self.strand)

    def addParent(self, isoform) :
        if isoform not in self.parents :
            self.parents.append(isoform)

class Isoform(BaseFeature) :
    def __init__(self, id, start, end, chromosome, strand, attr={}) :
        BaseFeature.__init__(self, ISOFORM_TYPE, start, end, chromosome, strand, attr)
        self.id       = id
        self.features = []
        self.exons    = []
        self.exonMap  = {}

    def __eq__(self, other) :
        return self.featureType == other.featureType and self.id == other.id

    def acceptorList(self) :
        """
        Returns a list of acceptor positions for the isoform.
        """
        self.exons.sort(reverse=(self.strand=='-'))
        return [e.acceptor() for e in self.exons[1:]]

    def addExon(self, exon) :
        """
        Adds an exon to the isoform if it's unique.
        Returns True if the exon was added; false otherwise.
        """
        if exon.strand != self.strand :
            raise Exception("Exon strand '%s' does not match isoform strand '%s' for %s" % (exon.strand, self.strand, self.id))
        if exon.chromosome != self.chromosome :
            raise Exception("Exon chromosome '%s' does not match isoform chromosome '%s' for %s" % (exon.chromosome, self.chromosome, self.id))

        exonTuple = (exon.minpos,exon.maxpos)
        try :
            ignore = self.exonMap[exonTuple]
            return False
        except KeyError :
            self.exons.append(exon)
            self.exonMap[exonTuple] = exon
            if self not in exon.parents :
                exon.parents.append(self)
            self.minpos = min(self.minpos, exon.minpos)
            self.maxpos = max(self.maxpos, exon.maxpos)
        return True

    def addFeature(self, feature) :
        if feature.strand != self.strand :
            raise Exception("ERROR: feature strand '%s' does not match form strand '%s'" % (feature.strand, self.strand))

        if feature.chromosome != self.chromosome :
            raise Exception("ERROR: feature chromosome '%s' does not match form chromosome '%s'" % (feature.chromosome, self.chromosome))

        feature.setParent(self.id)
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)
        self.features.append(feature)

    def detailString(self) :
        return "Isoform %s\nStart: %d; End: %d; Strand: %s\nExons: [%s]\n" % \
                (self.id, self.start(), self.end(), self.strand, self.exonString())

    def donorList(self) :
        """
        Returns a list of donor positions for the isoform.
        """
        self.exons.sort(reverse=(self.strand=='-'))
        return [e.donor() for e in self.exons[:-1]]

    def exonString(self) :
        return ','.join([str(e) for e in self.exons])

    def getFeatureList(self, featureType) :
        return [f for f in self.features if f.featureType == featureType]

    def gffStrings(self) :
        result  = []
        # Attributes depend on where data originated
        gffAttr = {PARENT_FIELD:self.id}
        gtfAttr = {PARENT_FIELD:self.id, GTF_TRANSCRIPT:self.id, GTF_TRANSNAME:self.id, GTF_PROTEIN_ID:self.id}
        for exon in self.exons :
            if GTF_TRANSCRIPT in exon.attributes :
                result.append(exon.gffString(altAttributes=gtfAttr))
            else :
                result.append(exon.gffString(altAttributes=gffAttr))
        return result

    def gtfStrings(self) :
        result   = []
        # Always sort in ascending order by position
        exonList = sorted(self.exons)
        for i in range(len(exonList)) :
            exon = exonList[i]
            result.append(exon.gtfString(self.id, self.parent, i+1))
        return result

    def sortedExons(self) :
        """
        Sorts the exons in an isoform based on its strand and
        returns the sorted list of exon objects.
        """
        self.exons.sort(reverse=(self.strand=='-'))
        return self.exons

    def sortedIntrons(self) :
        """
        Returns a list of intron (donor,acceptor) tuples sorted based on strand.
        """
        result = []
        exons  = self.sortedExons()
        prev   = exons[0]
        for e in exons[1:] :
            result.append((prev.donor(),e.acceptor()))
            prev = e
        return result

    def __str__(self) :
        return "%s (%s): %d-%d (len=%d, strand=%s), %d exons/cds, range %d to %d" % \
                (self.id, self.chromosome, self.start(), self.end(), len(self), self.strand, len(self.exons), self.minpos, self.maxpos)

# Just as exons are part of an isoform, so CDS elements are part of an mRNA sequence:
class CDS(Exon) :
    def __init__(self, start, end, chromosome, strand, attr={}) :
        BaseFeature.__init__(self, CDS_TYPE, start, end, chromosome, strand, attr)
        self.parents = []

    def __cmp__(self, o) :
        """Special for CDS records, as two records may be the same in other regards but different types."""
        result = BaseFeature.__cmp__(self,o)
        return cmp(self.featureType, o.featureType) if result == 0 else result

    def __eq__(self, o) :
        """Special for CDS records, as two records may have the same locations but different types."""
        return self.featureType == o.featureType and self.minpos == o.minpos and self.maxpos == o.maxpos \
                and self.strand == o.strand and self.chromosome == o.chromosome

    def __str__(self) :
        return "%s %d-%d (isoforms: %s)" % (self.featureType, self.start(), self.end(), ','.join([x.id for x in self.parents]))

# We treat UTR records the same way as CDS records
class FP_UTR(CDS) :
    def __init__(self, start, end, chromosome, strand, attr={}) :
        BaseFeature.__init__(self, FP_UTR_TYPE, start, end, chromosome, strand, attr)
        self.parents = []

class TP_UTR(CDS) :
    def __init__(self, start, end, chromosome, strand, attr={}) :
        BaseFeature.__init__(self, TP_UTR_TYPE, start, end, chromosome, strand, attr)
        self.parents = []

class mRNA(Isoform) :
    """
    An mRNA acts like an isoform in that it is associated with a parent gene
    and contains a number of coding sequences (CDS).
    """
    def __init__(self, id, start, end, chromosome, strand, attr={}) :
        BaseFeature.__init__(self, MRNA_TYPE, start, end, chromosome, strand, attr)
        self.id          = id
        self.exons       = []
        self.features    = []
        self.cds         = []
        self.cdsMap      = {}
        self.start_codon = None
        self.end_codon   = None

    def acceptorList(self) :
        """Returns a list of acceptor positions for the mRNA."""
        self.cds.sort(reverse=(self.strand=='-'))
        return [c.acceptor() for c in self.cds[1:]]

    def addCDS(self, cds) :
        """
        Adds a CDS to the mRNA if it's unique.
        Returns True if the CDS was added; false otherwise.
        """
        if cds.strand != self.strand :
            raise Exception("ERROR: CDS strand '%s' does not match gene strand '%s'" % (cds.strand, self.strand))
        if cds.chromosome != self.chromosome :
            raise Exception("ERROR: CDS chromosome '%s' does not match gene chromosome '%s'" % (cds.chromosome, self.chromosome))

        cdsTuple = (cds.minpos,cds.maxpos)
        try :
            ignore = self.cdsMap[cdsTuple]
            return False
        except KeyError :
            self.cdsMap[cdsTuple] = cds

        if self not in cds.parents :
            cds.parents.append(self)
        self.cds.append(cds)
        return True

    def donorList(self) :
        """Returns a list of donor positions for the mRNA."""
        self.cds.sort(reverse=(self.strand=='-'))
        return [c.donor() for c in self.cds[:-1]]

    def endCodon(self) :
        """Returns the end codon for this splice form as a duple of (start,end)
        positions, or None if there are none found.  Positions are relative to
        the start of the chromosome (1-based)."""
        if not self.end_codon : self.findCodons()
        return self.end_codon

    def findCodons(self) :
        """Infers a transcript's start and end codon positions based on
        the relative positions of UTR and CDS records."""
        if self.end_codon and self.start_codon : return
        if not self.cds : return
        self.cds.sort(reverse=(self.strand=='-'))
        prev = self.cds[0]
        for c in self.cds[1:] :
            if not self.start_codon and prev.featureType == FP_UTR_TYPE and c.featureType == CDS_TYPE :
                self.start_codon = (c.minpos,c.minpos+2) if self.strand == '+' else (c.maxpos-2,c.maxpos)
            elif not self.end_codon and prev.featureType == CDS_TYPE and c.featureType == TP_UTR_TYPE :
                self.end_codon = (prev.maxpos-2,prev.maxpos) if self.strand == '+' else (prev.minpos,prev.minpos+2)
            prev = c

    def inferCodons(self) :
        """This method will infer start and end codons even when there are no
        UTR records for a transcript.  It is best to use this only after all
        data have been loaded for a gene."""
        # No CDS records --> no way to infer codons
        if not self.cds : return

        # First try UTR inference:
        self.findCodons()

        # If either codon is missing, assume the CDS endpoints
        # represent start/stop codons.
        if not self.start_codon :
            c = self.cds[0]
            self.start_codon = (c.minpos,c.minpos+2) if self.strand == '+' else (c.maxpos-2,c.maxpos)

        if not self.end_codon :
            c = self.cds[-1]
            self.end_codon = (c.maxpos-2,c.maxpos) if self.strand == '+' else (c.minpos,c.minpos+2)

    def getUTRs(self) :
        """Returns a list of all UTR records in the mRNA object."""
        return [c for c in self.cds if c.featureType in [FP_UTR_TYPE, TP_UTR_TYPE]]

    def gffStrings(self) :
        result  = [self.gffString()]
        gffAttr = {PARENT_FIELD:self.id}
        gtfAttr = {PARENT_FIELD:self.id, GTF_TRANSCRIPT:self.id, GTF_TRANSNAME:self.id, GTF_PROTEIN_ID:self.id}
        for c in self.cds :
            if GTF_TRANSCRIPT in c.attributes :
                result.append(c.gffString(altAttributes=gtfAttr))
            else :
                result.append(c.gffString(altAttributes=gffAttr))
        return result

    def gtfStartCodon(self) :
        if self.start_codon :
            return '%s\t%s\tstart_codon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"' % \
                (self.chromosome, GFF_ID, self.start_codon[0], self.start_codon[1], self.strand, self.parent.id, self.id)

    def gtfStopCodon(self) :
        if self.end_codon :
            return '%s\t%s\tstop_codon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"' % \
                (self.chromosome, GFF_ID, self.end_codon[0], self.end_codon[1], self.strand, self.parent.id, self.id)

    def gtfStrings(self) :
        """
        Returns GTF strings for all elements of the mRNA transcript, including
        start/stop codon locations.
        """
        result  = []
        codonString = self.gtfStartCodon() if self.strand == '+' else self.gtfStopCodon()
        if codonString : result.append(codonString)

        # Always sort in ascending order by position
        cdsList = sorted(self.cds)
        for i in range(len(cdsList)) :
            c = cdsList[i]
            result.append(c.gtfString(self.id, self.parent, i+1))

        codonString = self.gtfStopCodon() if self.strand == '+' else self.gtfStartCodon()
        if codonString : result.append(codonString)

        return result

    def sortedCDS(self) :
        """
        Sorts the CDS in an mRNA based on its strand and returns the sorted list of CDS objects.
        """
        self.cds.sort(reverse=(self.strand=='-'))
        return self.cds

    def sortedExons(self, **args) :
        """
        Infers an exon list from a CDS list and sorts the list based on strand.  Usually
        this means a 5' UTR record abuts a CDS record or a CDS record abuts a 3' UTR record,
        in which case an expanded exon represents both.
        """
        minintron = getAttribute('minintron', 2, **args)
        cdsList   = self.sortedCDS()
        result    = []
        for i in range(len(cdsList)) :
            cds  = cdsList[i]
            if len(result) > 0 and abs(cds.start() - result[-1].end()) < minintron :
                prev       = result[-1]
                result[-1] = Exon(prev.start(), cds.end(), self.chromosome, self.strand)
            else :
                result.append(Exon(cds.start(), cds.end(), self.chromosome, self.strand))
        # Revise feature types to indicate CDS instead of exon
        for e in result :
            e.featureType = CDS_TYPE
        return result

    def startCodon(self) :
        """Returns the start codon for this splice form as a duple of (start,end)
        positions, or None if there are none found.  Positions are relative to
        the start of the chromosome (1-based)."""
        if not self.start_codon : self.findCodons()
        return self.start_codon

    def __str__(self) :
        return "%s (%s): %d-%d (len=%d, strand=%s), %d exons/cds, range %d to %d" % \
                (self.id, self.chromosome, self.start(), self.end(), len(self), self.strand, len(self.cds), self.minpos, self.maxpos)

class Gene(BaseFeature) :
    def __init__(self, id, note, start, end, chromosome, strand, name=None, attr={}) :
        BaseFeature.__init__(self, GENE_TYPE, start, end, chromosome, strand, attr)
        self.id       = id
        self.name     = name if name is not None else id
        self.note     = note
        self.isoforms = {}
        self.mrna     = {}
        self.exons    = []
        self.cds      = []
        self.exonMap  = {}
        self.cdsMap   = {}
        self.string   = ''

        # Codons for all transcripts/mRNA, one entry per transcript
        self.start_codons = {}
        self.end_codons   = {}

        # All other features associated with genes in an annotation file, such as:
        #    3'/5' UTRs, mRNA, miRNA, siRNA, tRNA, rRNA, ncRNA, snRNA, snoRNA
        self.features = []

    def acceptorList(self) :
        """
        Returns a list of acceptors for this gene.
        """
        acceptorSet = set()
        for k in list(self.isoforms.keys()) :
            acceptorSet.update(self.isoforms[k].acceptorList())
        for k in list(self.mrna.keys()) :
            acceptorSet.update(self.mrna[k].acceptorList())
        return sorted(list(acceptorSet), reverse=(self.strand=='-'))

    def addCDS(self, newmRNA, newCDS) :
        """
        Adds a CDS to the gene if it's unique.
        Returns True if the CDS was added; false otherwise.
        """
        result = False
        cdsTuple = (newCDS.featureType,newCDS.minpos,newCDS.maxpos)
        try :
            cds = self.cdsMap[cdsTuple]
        except KeyError :
            cds = newCDS
            self.cds.append(cds)
            self.cdsMap[cdsTuple] = cds
            self.minpos = min(self.minpos, cds.minpos)
            self.maxpos = max(self.maxpos, cds.maxpos)
            result = True

        mrna = self.addmRNA(newmRNA)
        mrna.addCDS(cds)
        self.start_codons[mrna.id] = mrna.startCodon()
        self.end_codons[mrna.id]   = mrna.endCodon()
        return result

    def addExon(self, newIsoform, newExon) :
        """
        Adds an exon to the gene if it's unique.
        Returns True if the exon was added; false otherwise.
        """
        if newIsoform is None : raise ValueError('Illegal null isoform in exon %s' % newExon)
        result = False
        posTuple = (newExon.minpos, newExon.maxpos)
        try :
            exon = self.exonMap[posTuple]
        except KeyError :
            exon   = newExon
            self.exons.append(exon)
            self.exonMap[posTuple] = exon
            self.minpos = min(self.minpos, exon.minpos)
            self.maxpos = max(self.maxpos, exon.maxpos)
            result = True

        isoform = self.addIsoform(newIsoform)
        isoform.addExon(exon)
        return result

    def addFeature(self, feature) :
        if feature.strand != self.strand :
            raise Exception("ERROR: feature strand '%s' does not match gene strand '%s'" % (feature.strand, self.strand))

        if feature.chromosome != self.chromosome :
            raise Exception("ERROR: feature chromosome '%s' does not match gene chromosome '%s'" % (feature.chromosome, self.chromosome))

        feature.setParent(self.id)
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)
        self.features.append(feature)

    def addIsoform(self, isoform) :
        isoform.setParent(self)
        return self.isoforms.setdefault(isoform.id, isoform)

    def addmRNA(self, mrna) :
        mrna.parent = self
        return self.mrna.setdefault(mrna.id, mrna)

    def detailString(self) :
        result = "%s (%s): %d-%d (len=%d, strand=%s), %d exons/cds, range %d to %d\n" % \
                (self.id, self.chromosome, self.start(), self.end(), len(self), self.strand, (len(self.exons)+len(self.cds)), self.minpos, self.maxpos)
        result += '\n  '.join([x.detailString() for x in list(self.isoforms.values())])
        return result

    def donorList(self) :
        """
        Returns a list of donors for this gene.
        """
        donorSet = set()
        for k in list(self.isoforms.keys()) :
            donorSet.update(self.isoforms[k].donorList())
        for k in list(self.mrna.keys()) :
            donorSet.update(self.mrna[k].donorList())
        return sorted(list(donorSet), reverse=(self.strand=='-'))

    def endCodons(self) :
        """Returns a dictionary of splice forms and their associated end codon locations.
        A codon location is given as a duple of (start,end) positions."""
        return self.end_codons

    def getFeatureList(self, featureType) :
        return [f for f in self.features if f.featureType == featureType]

    def getIntrons(self) :
        """
        Returns a list of duples containing start/end positions of introns in this gene.
        """
        result = {}
        for iid in list(self.isoforms.keys()) :
            exons = self.isoforms[iid].sortedExons()
            for i in range(1,len(exons)) :
                key = (exons[i-1].end(), exons[i].start())
                result[key] = 1

        for mid in list(self.mrna.keys()) :
            cds = self.mrna[mid].sortedExons()
            for i in range(1,len(cds)) :
                key = (cds[i-1].end(), cds[i].start())
                result[key] = 1

        return list(result.keys())

    def getIsoform(self, id) :
        return self.isoforms[id]

    def getJunctions(self) :
        """
        Returns a list of all known splice junctions for this gene
        based on its isoform list.  List contains only unique
        donor-acceptor duples.
        """
        result = {}
        for iid in list(self.isoforms.keys()) :
            iso   = self.isoforms[iid]
            exons = sorted(iso.exons, reverse=(self.strand=='-'))
            for i in range(1,len(exons)) :
                duple = (exons[i-1].donor(), exons[i].acceptor())
                result[duple] = 1
        return list(result.keys())

    def gffStrings(self) :
        """Returns a GFF string representation of the gene record plus
        all elements within the gene."""
        stringList = [self.gffString(altAttributes={ID_FIELD:self.id})]
        isoSet     = set(self.isoforms.keys())
        mrnaSet    = set(self.mrna.keys())
        commonKeys = isoSet & mrnaSet
        isoKeys    = isoSet-mrnaSet
        mrnaKeys   = mrnaSet-isoSet
        allKeys    = sorted(isoSet | mrnaSet)
        for k in allKeys :
            if k in commonKeys :
                allExons = self.isoforms[k].exons + self.mrna[k].cds
                allExons.sort(reverse=(self.strand=='-'))
                stringList.append(self.mrna[k].gffString())
                gffAttr = {PARENT_FIELD:k}
                gtfAttr = {PARENT_FIELD:k, GTF_TRANSCRIPT:k, GTF_TRANSNAME:k, GTF_PROTEIN_ID:k}
                for e in allExons :
                    if GTF_TRANSCRIPT in e.attributes :
                        stringList.append(e.gffString(altAttributes=gtfAttr))
                    else :
                        stringList.append(e.gffString(altAttributes=gffAttr))

            elif k in isoKeys :
                stringList += self.isoforms[k].gffStrings()
            elif k in mrnaKeys :
                stringList += self.mrna[k].gffStrings()

        return '\n'.join(stringList)

    def gtfStrings(self) :
        """Returns a GTF string representation of the gene record plus
        all elements within the gene."""
        stringList = []
        isoSet     = set(self.isoforms.keys())
        mrnaSet    = set(self.mrna.keys())
        commonKeys = isoSet & mrnaSet
        isoKeys    = isoSet-mrnaSet
        mrnaKeys   = mrnaSet-isoSet
        allKeys    = sorted(isoSet | mrnaSet)
        for k in allKeys :
            if k in commonKeys :
                allExons = self.isoforms[k].exons + self.mrna[k].cds
                allExons.sort()

                # Ensure that there are codons to write
                self.mrna[k].inferCodons()

                codonString = self.mrna[k].gtfStartCodon() if self.strand == '+' else self.mrna[k].gtfStopCodon()
                if codonString : stringList.append(codonString)

                eCtr = 0
                cCtr = 0
                for i in range(len(allExons)) :
                    item = allExons[i]
                    if item.featureType == EXON_TYPE :
                        eCtr += 1
                        stringList.append(item.gtfString(k,self,eCtr))
                    elif item.featureType == CDS_TYPE :
                        cCtr += 1
                        stringList.append(item.gtfString(k,self,cCtr))

                codonString = self.mrna[k].gtfStopCodon() if self.strand == '+' else self.mrna[k].gtfStartCodon()
                if codonString : stringList.append(codonString)

            elif k in isoKeys :
                stringList += self.isoforms[k].gtfStrings()
            elif k in mrnaKeys :
                stringList += self.mrna[k].gtfStrings()
        return '\n'.join(stringList)

    def isSingleExon(self) :
        return (len(self.exons) == 1)

    def sortedExons(self) :
        """Returns a list of all exons inferred by the gene model, sorted 5' to 3'."""
        tmpset = set()
        for iid in list(self.isoforms.keys()) :
            tmpset.update(self.isoforms[iid].sortedExons())

        # Avoid returning duplicates:
        stored = set([(e.minpos,e.maxpos) for e in tmpset])
        for mid in list(self.mrna.keys()) :
            tmpset.update([e for e in self.mrna[mid].sortedExons() if (e.minpos,e.maxpos) not in stored])

        result = list(tmpset)
        result.sort(reverse=(self.strand=='-'))
        return result

    def startCodons(self) :
        """Returns a dictionary of splice forms and their associated start codon locations.
        A codon location is given as a duple of (start,end) positions."""
        return self.start_codons

    def __str__(self) :
        if not self.string :
            self.string = "%s (%s): %d-%d (len=%d, strand=%s), %d exons/cds, range %d to %d" % \
                (self.id, self.chromosome, self.start(), self.end(), len(self), self.strand, (len(self.cds)+len(self.exons)), self.minpos, self.maxpos)
        return self.string

class PseudoGene(Gene) :
    def __init__(self, id, note, start, end, chromosome, strand, name=None, attr={}) :
        BaseFeature.__init__(self, PSEUDOGENE_TYPE, start, end, chromosome, strand, attr)
        self.id       = id
        self.name     = name
        self.note     = note
        self.features = []
        self.exons    = []
        self.cds      = []
        self.mrna     = {}
        self.isoforms = {}
        self.exonMap  = {}
        self.cdsMap   = {}
        self.string   = ''

        self.start_codons = {}
        self.end_codons   = {}

    def detailString(self) :
        return self.__str__()

class GeneModel(object) :
    def __init__(self, gffPath, **args) :
        """Instantiates a GeneModel object.  If a path is provided, this will
        load gene models from the given file."""
        # 2-dimensional model indexed by chromosome then gene
        self.allGenes   = {}
        self.allChr     = {}
        self.foundTypes = {}
        self.model      = {}
        self.mRNAforms  = {}
        self.sorted     = {}

        if gffPath : # Load gene models from a file
            if type(gffPath) == str and not os.path.exists(gffPath) :
                raise ValueError('Gene model file not found: %s' % gffPath)

            self.loadGeneModel(gffPath, **args)

            if not self.model :
                raise ValueError('No gene models found in %s' % gffPath)
            self.makeSortedModel()

    def __contains__(self, gene) :
        """Returns true if a gene is in the model; false otherwise."""
        return str(gene) in self.allGenes

    def addChromosome(self, start, end, name) :
        """Adds a chromosome to a gene model or updates the end points
        if the record already exists."""
        key = name.lower()
        try :
            rec        = self.allChr[key]
            rec.minpos = min(rec.minpos, start)
            rec.maxpos = max(rec.maxpos, end)
        except KeyError :
            self.allChr[key] = Chromosome(start, end, name)
            self.model.setdefault(key,{})

    def addGene(self, gene) :
        """Adds a gene to a gene model.  Raises a ValueError if the
        gene has already been added."""
        try :
            ##ignore = (self.model[gene.chromosome][gene.id], self.allGenes[gene.id])
            ignore = (self.model[gene.chromosome][gene.id], self.allGenes[str(gene)])
            raise ValueError('Gene %s already stored in gene model' % gene.id)
        except KeyError :
            self.model[gene.chromosome][gene.id] = gene
            ##self.allGenes[gene.id]       = gene
            self.allGenes[str(gene)]       = gene

    def binarySearchGenes(self, geneList, loc, lo, hi, pfx='') :
        """Quasi-binary search through a gene list.  Performs binary search to a point,
        then performs linear search within a refined gene range.  Real binary search may
        not be possible since gene models may overlap or contain one within another."""
        if (hi-lo) <= GENE_SEARCH_RANGE :
            loBound = max(0, lo-GENE_SEARCH_MARGIN)
            hiBound = min(hi+GENE_SEARCH_MARGIN, len(geneList))
            for gene in geneList[loBound:hiBound] :
                if gene.contains(loc, gene.strand) :
                    return (gene,gene)
            return (None,None)
        else :
            midpt   = (hi+lo)/2
            midGene = geneList[midpt]
            if midGene.contains(loc, midGene.strand) :
                return (midGene,midGene)
            elif loc > midGene.maxpos :
                return self.binarySearchGenes(geneList, loc, midpt+1, hi, pfx+' ')
            elif loc < midGene.minpos :
                return self.binarySearchGenes(geneList, loc, lo, midpt-1, pfx+' ')

    def cleanName(self, s) :
        """
        Some feature names include URL characters that we may wish to fix.
        """
        import urllib.request, urllib.parse, urllib.error
        revised = urllib.parse.unquote(s)
        revised = revised.replace(',','')
        return revised.replace(' ','-')

    def getAllAcceptors(self, geneFilter=defaultGeneFilter) :
        """
        Returns a dictionary of all known acceptor sites in the gene model,
        indexed by chromosome and strand.
        """
        result = {}
        for chrom in list(self.model.keys()) :
            result[chrom] = self.getKnownAcceptors(chrom, geneFilter)
        return result

    def getAllDonors(self, geneFilter=defaultGeneFilter) :
        """
        Returns a dictionary of all known donor sites in the gene model,
        indexed by chromosome and strand.
        """
        result = {}
        for chrom in list(self.model.keys()) :
            result[chrom] = self.getKnownDonors(chrom, geneFilter)
        return result

    def getAllGeneIds(self, geneFilter=defaultGeneFilter) :
        """Returns a list of ids for all genes stored."""
        return [g.id for g in list(self.allGenes.values()) if geneFilter(g)]

    def getAllGenes(self, geneFilter=defaultGeneFilter, **args) :
        """Returns a list of all genes stored."""
        verbose   = getAttribute('verbose', False, **args)
        indicator = ProgressIndicator(10000, verbose=verbose)
        result    = []
        for g in list(self.allGenes.values()) :
            indicator.update()
            if geneFilter(g) : result.append(g)
        indicator.finish()
        return result

    def getAnnotation(self, key, annotDict, default=None) :
        """
        Convenience method for retrieving a value from an annotation dictionary
        """
        try :
            return annotDict[key]
        except KeyError :
            return default
        
    def getAnnotationDict(self, s) :
        """
        Parses a ';'-separated annotation string containing key-value pairs
        and returns them as a dictionary.
        """
        valStr = s.replace(' ','')
        parts  = valStr.split(';')
        result = {}
        for p in parts :
            if p.find('=') >= 0 :
                keyval = p.split('=')
                result[keyval[0]] = keyval[1]
        return result
        
    def getChromosome(self, chrName) :
        """Returns a simple record with basic chromosome information."""
        try :
            return self.allChr[chrName]
        except KeyError :
            return None

    def getChromosomes(self) :
        """Returns a list of all chromosomes represented in the model."""
        return list(self.model.keys())

    def getFeatureList(self, featureType) :
        """
        Returns a list of all features of the given type found in all genes.
        """
        result = []
        for gene in list(self.allGenes.values()) :
            fList = gene.getFeatureList(featureType)
            if fList :
                result += fList
        return result

    def getGene(self, chrom, geneId) :
        """
        Returns a gene from within a chromosome.  The gene will contain
        information on all exons within it.
        """
        try :
            return self.model[chrom.lower()][geneId]
        except KeyError :
            return None

    def getGeneByName(self, id) :
        """
        Returns a gene with the given id if it exists.
        """
        for k in list(self.model.keys()) :
            try :
                return self.model[k][id.upper()]
            except KeyError :
                pass
        return None

    def getGeneFromLocations(self, chrom, startPos, endPos, strand) :
        """
        Finds the gene within the given chromosome that contains the given start
        and end positions.  Uses quasi-binary search through the sorted list of genes.
        """
        # Only get genes on the correct strand:
        try :
            geneList  = self.sorted[chrom.lower()][strand]
        except KeyError :
            raise KeyError('Key %s not found in %s' % (chrom.lower(), ','.join(list(self.sorted.keys()))))

        (logene,higene) = self.binarySearchGenes(geneList, startPos, 0, len(geneList)-1)

        if not (logene or higene) :
            return None
        elif logene.contains(startPos, strand) or logene.contains(endPos, strand) :
            return logene
        elif higene.contains(startPos, strand) or higene.contains(endPos, strand) :
            return higene

        return None

    def getGeneRecords(self, chrom, geneFilter=defaultGeneFilter, **args) :
        """
        Returns a list of all gene instances represented within a given chromosome.
        The gene list may be filtered by changing the geneFilter function.
        """
        verbose = getAttribute('verbose', False, **args)
        try :
            #return [g for g in self.model[chrom.lower()].values() if geneFilter(g)]
            indicator = ProgressIndicator(10000, verbose=verbose)
            result    = []
            for g in list(self.model[chrom.lower()].values()) :
                indicator.update()
                if geneFilter(g) : result.append(g)
            indicator.finish()
            return result
        except KeyError :
            return []

    def getGenes(self, chrom) :
        """
        Returns a list of all genes represented within a given chromosome.
        """
        try :
            return list(self.model[chrom.lower()].keys())
        except KeyError :
            return []

    def getGenesInRange(self, chrom, minpos, maxpos, strand=None) :
        """
        Returns a list of all gene instances represented within a given chromosome
        that overlap the given range.  If no strand is specified, this will
        return all genes on both strands that overlap the range.
        """
        result = []
        for g in self.getGeneRecords(chrom) :
            if g.maxpos < minpos or g.minpos > maxpos : continue
            if strand and g.strand != strand : continue
            result.append(g)
        return result

    def getKnownAcceptors(self, chrom, geneFilter=defaultGeneFilter) :
        """
        Returns a dictionary of all known acceptor sites for the chromosome,
        indexed by strand.
        """
        result = {'-':set([]), '+':set([])}
        for g in self.getGeneRecords(chrom, geneFilter) :
            result[g.strand].update(g.acceptorList())
        return result

    def getKnownDonors(self, chrom, geneFilter=defaultGeneFilter) :
        """
        Returns a dictionary of all known donor sites for the chromosome,
        indexed by strand.
        """
        result = {'-':set([]), '+':set([])}
        for g in self.getGeneRecords(chrom, geneFilter) :
            result[g.strand].update(g.donorList())
        return result

    def getParent(self, s, chrom, searchGenes=True, searchmRNA=True) :
        """
        Parent identifiers are not stored in a consistent manner.  We may have
        'AT1G01160', 'AT1G01160.1' or '12345.AT1G01160' or possibly something else.
        This method looks for the most specific candidate name to identify a record's parent.
        """
        def getSubnames(fullString,d) :
            parts  = fullString.split(d)
            result = []
            # most-to-least specific: 'a-b-c-d' --> ['a-b-c', 'a-b', 'a'] or 'a,b,c' --> ['a,b', 'a']
            for i in range(len(parts)-1,0,-1) :
                result.append(d.join(parts[:i]))
            return result

        parentString = s.upper()
        # First try the easiest method:
        if searchmRNA :
            try :
                return self.mRNAforms[chrom][parentString]
            except KeyError :
                pass

        if searchGenes :
            try :
                return self.model[chrom][parentString]
            except KeyError :
                pass

        # Iterate over a set of known delimiters 
        delim      = [c for c in FORM_DELIMITERS if c in parentString]
        candidates = [parentString]
        for c in delim :
            candidates += getSubnames(parentString,c)

        # First try known mRNA records
        if searchmRNA and chrom in self.mRNAforms :
            commasep   = parentString.split(',')
            additional = list(commasep)
            for c in commasep :
                additional += c.split('.')
            candidates += additional
            for c in candidates :
                try :
                    return self.mRNAforms[chrom][c]
                except KeyError :
                    pass

        # Next try known genes; 
        if searchGenes and chrom in self.model :
            for c in candidates :
                if c in self.model[chrom] :
                    return self.model[chrom][c]

    def getRecordTypes(self) :
        """Returns a list of all record types found in the input file."""
        return [k for k in self.foundTypes if self.foundTypes[k]]

    def isoformDict(self, **args) :
        """Returns a dictionary that maps gene names to their corresponding
        isoform identifiers.  Each gene is associated with a set of isoform ids."""
        result = {}
        for g in self.getAllGenes(**args) :
            result[g.id] = set(list(g.mrna.keys()) + list(g.isoforms.keys()))
        return result

    def loadGeneModel(self, gffRecords, **args) :
        """
        Reads a tab-delimited gene annotation GFF file and stores information
        on chromosomes, the genes within each chromosome and exons within each gene.

        Parameters:
          'gffRecords'   - source of GFF records; may be a file path, a file stream
                           or a list/set of strings
          'requireNotes' - require annotations for all gene records (default=False)
          'chromosomes'  - chromosome name or list of chromosomes to store (default=all)
          'verbose'      - provide verbose feedback (default=False)
          'ignoreErrors' - ignore error conditions (default=False)
        """
        verbose      = getAttribute('verbose', False, **args)
        requireNotes = getAttribute('requireNodes', False, **args)
        ignoreErrors = getAttribute('ignoreErrors', False, **args)
        chromosomes  = None

        # Convenience method for handling exceptions that could be ignored:
        def conditionalException(s) :
            if not ignoreErrors : raise Exception(s)

        if type(gffRecords) == str :
            if verbose : sys.stderr.write('Loading and validating gene models in %s\n' % gffRecords)
            instream = ezopen(gffRecords)
        elif type(gffRecords) in [list, set, file] :
            if verbose : sys.stderr.write('Loading and validating gene models in %d records\n' % len(gffRecords))
            instream = gffRecords
        else :
            raise ValueError('Unrecognized GFF record source (%s); must be file path or a list/set of strings.')

        if 'chromosomes' in args :
            chrParam = args['chromosomes']
            if (type(chrParam) == type([])) :
                chromosomes = [s.lower() for s in chrParam]
            elif chrParam :
                chromosomes = [chrParam]
            if verbose : sys.stderr.write("GeneModel loading chromosomes %s\n" % ','.join(chromosomes))

        self.model      = {}
        self.mRNAforms  = {}
        self.allGenes   = {}
        self.foundTypes = {}
        self.allChr     = {}
        lineCtr         = 0

        geneAlias = {}
        geneCount = 0
        exonCount = 0
        isoCount  = 0
        mrnaCount = 0
        cdsCount  = 0
        badLines  = 0
        indicator = ProgressIndicator(1000000, verbose=verbose)

        for line in instream :
            lineCtr  += 1
            indicator.update()
            s     = line.rstrip()
            if not s or s[0] == '#' : continue # GFF comments:

            parts = s.split('\t')

            if len(parts) < 7 :
                badLines += 1
                if verbose : sys.stderr.write("line %d: invalid GFF format (not enough columns); file may be corrupt\n" % lineCtr)
                if badLines >= MAX_BAD_LINES :
                    sys.stderr.write('\nInput GFF file appears to be invalid; aborting.\n')
                    raise ValueError('Invalid GFF input file')
                continue

            annots    = self.getAnnotationDict(parts[-1])
            # Columns in GFF file
            # Since chromosome name and record type may be used in
            # comparisions, make them all lowercase.
            chrName   = parts[0].lower()
            if chromosomes and chrName not in chromosomes : continue

            # Try to map record types to known types, but default to string if necessary
            recType = parts[2].lower()
            try :
                recType = RECTYPE_MAP[recType]
            except KeyError :
                pass

            startPos = int(parts[3])
            endPos   = int(parts[4])
            strand   = parts[6]

            # Track all known record types and which were stored
            self.foundTypes[recType] = (recType in KNOWN_RECTYPES and recType not in IGNORE_RECTYPES)
            if not self.foundTypes[recType] : continue

            # Many gene types possible, usually ending in '_gene'
            if recType in [GENE_TYPE, PSEUDOGENE_TYPE] :
                # Id is used for comparision, so use all uppercase
                try :
                    gid  = annots[ID_FIELD].upper()
                except KeyError :
                    raise ValueError('line %d: %s record has no ID field:\n%s\n' % (lineCtr,recType,line))

                name = self.getAnnotation(NAME_FIELD, annots, None)
                if name : name = self.cleanName(name)

                note = self.getAnnotation(NOTE_FIELD, annots)
                if not note and requireNotes : continue

                if strand not in VALID_STRANDS :
                    conditionalException("line %d: %s record with unknown strand" % (lineCtr,recType))

                if not chrName in self.model :
                    self.model[chrName] = {}
                    self.addChromosome(1,endPos,chrName)

                if recType == PSEUDOGENE_TYPE :
                    gene = PseudoGene(gid, note, startPos, endPos, chrName, strand, name, annots)
                else : # genes and predicted_genes
                    gene = Gene(gid, note, startPos, endPos, chrName, strand, name, annots)

                try :
                    other = self.allGenes[str(gene)]
                    conditionalException("line %d: gene %s associated with multiple loci: %d-%d and %d-%d" % \
                            (lineCtr, gene.id, other.minpos, other.maxpos, startPos, endPos))
                except KeyError :
                    pass

                self.allChr[chrName].update(gene)
                self.addGene(gene)
                # Map both relationships
                geneAlias[gene.name.upper()] = gene.id.upper()
                geneAlias[gene.id.upper()] = gene.name.upper()
                geneCount += 1

            elif recType in [EXON_TYPE, PSEUDOEXON_TYPE] :
                if chrName not in self.model : continue

                # Parent identifiers for exons
                parent  = None
                tried   = set()
                for k in POSSIBLE_GENE_FIELDS :
                    if (k not in annots) or (annots[k] in tried) : continue
                    isoName = annots[k]
                    parent  = self.getParent(isoName, chrName)
                    if parent : break
                    tried.add(isoName)
                if not parent : continue

                gene = parent.parent if parent.featureType == MRNA_TYPE else parent

                # Next get the isoform if it already exists:
                isoform = None
                isoName = ''
                tried   = set()
                # Look through known transcript ID keys for a suitable attribute:
                for k in POSSIBLE_FORM_FIELDS :
                    if k not in annots :
                        continue
                    if annots[k] in tried : continue
                    isoName = annots[k]
                    if isoName in gene.isoforms :
                        isoform = gene.isoforms[isoName]
                        break
                    else :
                        tried.add(isoName)

                # If there is no form or name to trace we're out of luck:
                if not (isoform or isoName) : continue

                if strand in VALID_STRANDS and strand != gene.strand :
                    conditionalException("line %d: exon strand (%s) != gene strand (%s) for %s" % (lineCtr,strand,gene.strand,gene.id))
                else :
                    strand = gene.strand

                # Must add isoform to gene before adding exon
                if not isoform :
                    isoAttr   = {PARENT_FIELD:gene.id, NAME_FIELD:isoName, ID_FIELD:isoName}
                    isoform   = Isoform(isoName, startPos, endPos, chrName, strand, attr=isoAttr)
                    isoCount += 1

                exon      = Exon(startPos, endPos, chrName, strand, annots)
                exonCount = exonCount+1 if gene.addExon(isoform, exon) else exonCount

            elif recType in [MRNA_TYPE, PSEUDOTRANS_TYPE] :
                if chrName not in self.model :
                    conditionalException("line %d: mRNA with missing chromosome dictionary %s (known: %s)" \
                            % (lineCtr, chrName, ','.join(list(self.model.keys()))))

                if not ID_FIELD in annots :
                    conditionalException("line %d: mRNA with missing ID" % lineCtr)

                id     = annots[ID_FIELD].upper()
                parent = self.getAnnotation(PARENT_FIELD, annots)
                if not parent : continue

                parent = parent.upper()
                gene   = self.getParent(parent, chrName)
                if not gene :
                    try :
                        alias = geneAlias[parent]
                        gene  = self.getParent(alias, chrName)
                    except KeyError :
                        if verbose :
                            if alias == parent :
                                sys.stderr.write("line %d: no gene %s found for %s\n" % (lineCtr, parent, recType))
                            else :
                                sys.stderr.write("line %d: no gene '%s' or '%s' found for %s\n" % (lineCtr, parent, alias, recType))
                        continue

                if strand in VALID_STRANDS and strand != gene.strand :
                    conditionalException("line %d: mRNA strand (%s) does not match gene strand (%s)" % (lineCtr,strand,gene.strand))
                else :
                    strand = gene.strand

                mrnaAttr = {PARENT_FIELD:gene.id, NAME_FIELD:id, ID_FIELD:id}
                mrna     = mRNA(id, startPos, endPos, chrName, strand, attr=mrnaAttr)
                gene.addmRNA(mrna)
                self.mRNAforms[chrName] = self.mRNAforms.setdefault(chrName, {})
                self.mRNAforms[chrName][id] = mrna
                mrnaCount += 1

            elif recType in CDS_TYPES :
                if chrName not in self.model :
                    conditionalException("line %d: %s has unrecognized chromosome: %s (known: %s)" % \
                            (lineCtr, recType, chrName, ','.join(list(self.model.keys()))))
                if chrName not in self.mRNAforms :
                    conditionalException("line %d: %s has unrecognized chromosome: %s (known: %s)" % \
                            (lineCtr, recType, chrName, ','.join(list(self.mRNAforms.keys()))))

                mrna = self.getParent(annots[PARENT_FIELD], chrName, searchGenes=False)
                if not mrna :
                    if verbose : sys.stderr.write("line %d: no mRNA %s found for %s\n" % (lineCtr, annots[PARENT_FIELD], recType))
                    continue

                if strand in VALID_STRANDS and strand != mrna.strand :
                    conditionalException("line %d: CDS strand (%s) does not match mRNA strand (%s)" % (lineCtr,strand,mrna.strand))
                else :
                    strand = mrna.strand

                geneId = mrna.parent.id
                cds    = cdsFactory(recType, startPos, endPos, chrName, strand, annots)

                if not gene : conditionalException("line %d: No gene loaded for CDS record" % lineCtr)
                if gene.addCDS(mrna, cds) : cdsCount += 1

            elif recType == CHR_TYPE :
                self.addChromosome(startPos, endPos, chrName)

            elif recType in [FP_UTR_TYPE, TP_UTR_TYPE] :
                if chrName not in self.model : continue
                if PARENT_FIELD not in annots : continue

                parentId = annots[PARENT_FIELD]
                parent   = self.getParent(parentId, chrName)
                if parent :
                    if strand in VALID_STRANDS and strand != parent.strand :
                        conditionalException("line %d: %s strand (%s) does not match parent strand (%s)" % (lineCtr,recType,strand,parent.strand))
                    else :
                        strand = parent.strand
                    parent.addFeature(BaseFeature(recType, startPos, endPos, chrName, strand, annots))
                else :
                    if verbose : sys.stderr.write("line %d: no parent %s found for %s at line %d: %s\n" % (lineCtr, parentId, recType, e))
                    continue

            elif recType not in IGNORE_RECTYPES :
                if chrName not in self.model : continue
                if PARENT_FIELD not in annots : continue

                parent = annots[PARENT_FIELD].split('.')[0]
                if parent not in self.model[chrName] :
                    continue

                gene = self.model[chrName][parent]
                try :
                    gene.addFeature(BaseFeature(recType, startPos, endPos, chrName, strand, annots))
                except Exception as e :
                    conditionalException("line %d: %s" % (lineCtr, e))

        indicator.finish()
        if verbose :
            if geneCount > 0 :
                sys.stderr.write("Loaded %s genes with %s isoforms, %s exons (avg. %.1f/gene), %s mRNA, %s CDS (avg. %.1f/gene)\n" % \
                    (commaFormat(geneCount), commaFormat(isoCount), commaFormat(exonCount), float(exonCount)/geneCount, commaFormat(mrnaCount), commaFormat(cdsCount), float(cdsCount/geneCount)))
            else :
                sys.stderr.write("** Warning: no genes loaded!\n")

    def makeSortedModel(self) :
        self.sorted = {}
        for chrom in list(self.model.keys()) :
            self.sorted[chrom] = {}
            # Get a list of gene objects sorted by position
            geneList = list(self.model[chrom].values())
            geneList.sort()
            # Note: we accept unknown ('.') strands, but it's
            # up to calling routines to ask for them specifically
            self.sorted[chrom] = {'-':[], '+':[], '.':[]}
            for g in geneList :
                self.sorted[chrom][g.strand].append(g)

    def writeGFF(self, gffPath, **args) :
        """Writes a complete gene model out to a GFF file."""
        geneFilter = getAttribute('geneFilter', defaultGeneFilter, **args)
        geneSubset = getAttribute('geneSet', None, **args)
        verbose    = getAttribute('verbose', False, **args)

        outStream = open(gffPath, 'w') if type(gffPath) == type('') else gffPath
        chromList = sorted(self.allChr.keys())
        indicator = ProgressIndicator(10000, verbose=verbose)
        for c in chromList :
            chrom = self.getChromosome(c)
            outStream.write('%s\n' % chrom.gffString())
            genes = self.getGeneRecords(c, geneFilter)
            if geneSubset :
                genes = [g for g in genes if g.id in geneSubset or g.name in geneSubset]
            genes.sort()
            for g in genes :
                indicator.update()
                strings = g.gffStrings()
                if strings : outStream.write('%s\n' % strings)
        indicator.finish()

    def writeGTF(self, gtfPath, geneFilter=defaultGeneFilter, **args) :
        """Writes a complete gene model out to a GTF file."""
        verbose   = getAttribute('verbose', False, **args)
        outStream = open(gtfPath, 'w') if type(gtfPath) == type('') else gtfPath
        chromList = sorted(self.allChr.keys())
        indicator = ProgressIndicator(10000, verbose=verbose)
        for c in chromList :
            chrom = self.getChromosome(c)
            genes = self.getGeneRecords(c, geneFilter)
            genes.sort()
            for g in genes :
                indicator.update()
                strings = g.gtfStrings()
                if strings : outStream.write('%s\n' % strings)
        indicator.finish()
