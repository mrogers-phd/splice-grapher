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
from SpliceGrapher.shared.utils      import *
from SpliceGrapher.formats.gtf       import *
from SpliceGrapher.formats.GeneModel import *

import sys

# Types found in feature-type column 3 (ENSEMBL/UCSC)
KNOWN_TYPES = [GTF_CDS, GTF_EXON, GTF_START_CODON, GTF_STOP_CODON] = ['cds', 'exon', 'start_codon', 'stop_codon']

# Attributes that may appear in GTF records but are not appropriate for a gene:
NON_GENE_ATTRIBUTES = ['exon_number', 'transcript_id', 'transcript_name', 'parent']

class GTFException(Exception) :
    pass

class GTFLoader(object) :
    """Class for loading gene models from a GTF file into memory."""

    def __init__(self, path, **args) :
        self.model   = self.loadGeneModels(path, **args)
        self.unknown = set()

    def addGeneToModel(self, model, rec, start, end) :
        """Adds a new gene to the model if it is not already there.
        Returns the new gene or the original gene."""
        minpos = min(start,end)
        maxpos = max(start,end)
        chrom  = rec.parts[SEQNAME_INDEX]
        gid    = rec.attrs[GENE_ID_ATTR].upper()
        try :
            gene        = model.model[chrom][gid]
            gene.minpos = min(minpos, gene.minpos)
            gene.maxpos = max(maxpos, gene.maxpos)
        except KeyError :
            attrs = dict([(k,rec.attrs[k]) for k in rec.attrs if k not in NON_GENE_ATTRIBUTES])
            gene = Gene(gid, None, start, end, chrom, rec.parts[STRAND_INDEX], rec.attrs[GENE_NAME_ATTR], attrs)
            model.model[chrom][gid] = gene
            model.allGenes[gid] = model.model[chrom][gid]
        return gene

    def addExonToGene(self, gene, rec, verbose=False) :
        """Adds an exon to the gene along with its isoform."""
        isoName = rec.attrs[TRANS_ID_ATTR]
        chrom   = rec.parts[SEQNAME_INDEX]
        strand  = rec.parts[STRAND_INDEX]

        try :
            newIsoform = gene.isoforms[isoName]
        except KeyError :
            isoAttr    = {PARENT_FIELD:gene.id, NAME_FIELD:isoName, ID_FIELD:isoName}
            newIsoform = Isoform(isoName, rec.start(), rec.end(), chrom, strand, attr=isoAttr)

        newExon  = Exon(rec.start(), rec.end(), chrom, strand, rec.attrs)
        posTuple = (newExon.minpos, newExon.maxpos)
        try :
            exon = gene.exonMap[posTuple]
        except KeyError :
            exon = newExon
            gene.exons.append(exon)
            gene.exonMap[posTuple] = exon

        isoform = gene.addIsoform(newIsoform)
        isoform.addExon(exon)
        return exon

    def addmRNA(self, gene, rec) :
        """Adds an mRNA record to a gene, if necessary."""
        chrom    = rec.parts[SEQNAME_INDEX]
        strand   = rec.parts[STRAND_INDEX]
        mrnaName = rec.attrs[TRANS_ID_ATTR]
        try :
            return gene.mrna[mrnaName]
        except KeyError :
            isoAttr = {PARENT_FIELD:gene.id, NAME_FIELD:mrnaName, ID_FIELD:mrnaName}
            gene.addmRNA(mRNA(mrnaName, rec.start(), rec.end(), chrom, strand, attr=isoAttr))
            return gene.mrna[mrnaName]

    def addCDSToGene(self, gene, rec) :
        """Adds a coding sequence the gene along with its mRNA."""
        mrnaName = rec.attrs[TRANS_ID_ATTR]
        chrom    = rec.parts[SEQNAME_INDEX]
        strand   = rec.parts[STRAND_INDEX]
        newCDS   = CDS(rec.start(), rec.end(), chrom, strand, rec.attrs)
        posTuple = (newCDS.minpos, newCDS.maxpos)
        try :
            cds = gene.cdsMap[posTuple]
        except KeyError :
            cds = newCDS
            gene.cds.append(cds)
            gene.cdsMap[posTuple] = cds

        mrna = self.addmRNA(gene, rec)
        mrna.addCDS(cds)
        if mrnaName in gene.start_codons : mrna.start_codon = gene.start_codons[mrnaName]
        if mrnaName in gene.end_codons   : mrna.end_codon   = gene.end_codons[mrnaName]

    def addUTRToGene(self, gene, rec, cds) :
        """This will see if the CDS start or end position matches a codon
        record start/end.  If such a codon exists, then a new 5' or 3' UTR
        record may be added to the gene."""
        transId = rec.attrs[TRANS_ID_ATTR]
        start   = gene.start_codons[transId] if transId in gene.start_codons else None
        stop    = gene.end_codons[transId] if transId in gene.end_codons else None
        if not (start or stop) : return

        def cdsContains(e, duple) :
            return duple and e.minpos < min(duple) and max(duple) < e.maxpos

        def cdsDownstream(e, duple) :
            if not duple : return False
            return bool(e.minpos > max(duple)) if e.strand == '+' else bool(e.maxpos < min(duple))

        def cdsUpstream(e, duple) :
            if not duple : return False
            return bool(e.maxpos < min(duple)) if e.strand == '+' else bool(e.minpos > max(duple))

        fpUTR = None
        if cdsContains(cds, start) :
            (minpos,maxpos) = (cds.minpos, min(start)-1) if gene.strand == '+' else (max(start)+1, cds.maxpos)
            fpUTR = FP_UTR(minpos, maxpos, gene.chromosome, gene.strand, rec.attrs)
        elif cdsUpstream(cds, start) :
            fpUTR = FP_UTR(cds.minpos, cds.maxpos, gene.chromosome, gene.strand, rec.attrs)

        tpUTR = None
        if cdsContains(cds, stop) :
            (minpos,maxpos) = (max(stop)+1, cds.maxpos) if gene.strand == '+' else (cds.minpos, min(stop)-1)
            tpUTR = TP_UTR(minpos, maxpos, gene.chromosome, gene.strand, rec.attrs)
        elif cdsDownstream(cds, stop) :
            tpUTR = TP_UTR(cds.minpos, cds.maxpos, gene.chromosome, gene.strand, rec.attrs)

        for utr in [fpUTR, tpUTR] :
            if not utr : continue
            ignore   = self.addmRNA(gene, rec)
            posTuple = (utr.minpos, utr.maxpos)
            try :
                c = gene.cdsMap[posTuple]
                c.featureType = utr.featureType
            except KeyError :
                gene.cds.append(utr)
                gene.cdsMap[posTuple] = utr
                gene.mrna[transId].addCDS(utr)

    def loadGeneModels(self, path, **args) :
        """Loads gene models from a GTF file and returns a GeneModel object."""
        verbose   = getAttribute('verbose', False, **args)
        result    = GeneModel(None)
        unknown   = set()
        badGenes  = set()
        oldGene   = None
        gene      = None
        skipped   = 0
        if verbose : sys.stderr.write('Loading and validating gene models from %s\n' % path)
        indicator = ProgressIndicator(1000000, verbose=verbose)
        for line in ezopen(path) :
            indicator.update()
            try :
                rec   = GTF_Line(line.strip())
            except ValueError, ve :
                skipped += 1
                continue

            ftype = rec.parts[FEATURE_INDEX].lower()
            if ftype not in KNOWN_TYPES :
                unknown.add(ftype)
                continue
    
            # Ensure that chromosomes are always in lower case
            rec.parts[SEQNAME_INDEX] = rec.parts[SEQNAME_INDEX].lower()
            start  = rec.start()
            end    = rec.end()
            strand = rec.strand()
            result.addChromosome(1, end, rec.parts[SEQNAME_INDEX])
            gene  = self.addGeneToModel(result, rec, start, end)

            # When we switch from one gene to the next, 
            # make sure UTRs are created for the old gene:
            if oldGene and gene != oldGene :
                for cds in oldGene.cds :
                    self.addUTRToGene(oldGene, rec, cds)

            # May not want this to be a catastrophic error
            # and accept some coincidental correctness
            if strand != gene.strand :
                badGenes.add(gene.id)
                continue

            tid = rec.attrs[TRANS_ID_ATTR]
            if ftype == GTF_EXON :
                exon = self.addExonToGene(gene, rec, verbose=verbose)
            elif ftype == GTF_CDS :
                self.addCDSToGene(gene, rec)
            elif ftype == GTF_START_CODON :
                gene.start_codons[tid] = (start,end)
                if tid in gene.mrna : gene.mrna[tid].start_codon = (start,end)
            elif ftype == GTF_STOP_CODON :
                gene.end_codons[tid] = (start,end)
                if tid in gene.mrna : gene.mrna[tid].end_codon = (start,end)

            oldGene = gene

        # Create any missing UTRs in the last gene:
        if gene :
            for exon in gene.exons :
                self.addUTRToGene(oldGene, rec, exon)
        indicator.finish()

        # Remove all invalid genes from the model
        if verbose and skipped > 0 :
            sys.stderr.write('Skipped %s records missing transcript_id attributes\n' % commaFormat(skipped))
        if verbose and badGenes :
            sys.stderr.write('Removing %s genes with missing or inconsistent strands\n' % commaFormat(len(badGenes)))
        for gid in badGenes :
            del result.allGenes[gid]

        # Final step required for any location-based searches
        result.makeSortedModel()
        return result
