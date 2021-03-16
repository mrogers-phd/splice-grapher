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
from SpliceGrapher.SpliceGraph  import SpliceGraph

from sys import maxint as MAXINT
import sys

# GTF columns (summarized from GTF 2.2 spec):
#  0 - <seqname> chromosome ID or contig ID
#  1 - <source> annotation source --- typically the name of either a prediction program or a public database.
#  2 - <feature> Valid feature types (must have correct capitalization):
#      CDS           (required)
#      start_codon   (required)
#      stop_codon    (required)
#      5UTR          (optional)
#      3UTR          (optional)
#      inter         (optional) intergenic
#      inter_CNS     (optional) intergenic conserved noncoding sequence
#      intron_CNS    (optional) intronic conserved noncoding sequence
#      exon          (optional)
#    Additional feature types found in gencode files:
#      gene          (optional)
#      transcript    (optional)
#      UTR           (optional)
#      Selenocysteine
#  3 - <start> (1-based start, end positions; start <= end)
#  4 - <end> 
#  5 - <score>  confidence score
#  6 - <strand> 
#  7 - <frame>  0=5' most base; 1=next one; 2 = next one after that
#  8 - [attributes] All features have the same two mandatory attributes:
#      gene_id <value>;
#      transcript_id <value>; 
#  [comments] Comments begin with a hash ('#') and continue to the end of the line.

# All known ENSEMBL GTF gene types:
ALL_ENSEMBL_SOURCES = [ "3prime_overlapping_ncrna", "ambiguous_orf", "antisense", "disrupted_domain",
                        "IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_J_gene", "IG_J_pseudogene",
                        "IG_M_gene", "IG_V_gene", "IG_V_pseudogene", "IG_Z_gene", "lincRNA",
                        "miRNA", "miRNA_pseudogene", "misc_RNA", "misc_RNA_pseudogene", "Mt_rRNA",
                        "Mt_tRNA", "Mt_tRNA_pseudogene", "ncRNA", "ncrna_host", "non_coding",
                        "nonsense_mediated_decay", "non_stop_decay", "polymorphic_pseudogene", "processed_pseudogene", "processed_transcript",
                        "protein_coding", "pseudogene", "retained_intron", "retrotransposed", "rRNA",
                        "rRNA_pseudogene", "scRNA_pseudogene", "sense_intronic", "sense_overlapping", "snlRNA",
                        "snoRNA", "snoRNA_pseudogene", "snRNA", "snRNA_pseudogene", "TEC",
                        "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "TR_C_gene", "TR_D_gene", "TR_J_gene",
                        "TR_J_pseudogene", "tRNA", "tRNA_pseudogene", "TR_V_gene", "TR_V_pseudogene",
                        "unitary_pseudogene", "unprocessed_pseudogene",]

(SEQNAME_INDEX, SOURCE_INDEX, FEATURE_INDEX, START_INDEX, END_INDEX, SCORE_INDEX, STRAND_INDEX, FRAME_INDEX, ATTR_INDEX, COMMENT_INDEX) = range(10)

# EMBL appears to store chromosome ids as integers or X, Y or M:
KNOWN_CHROMOSOMES = ['X', 'Y', 'M']

# EMBL sources we care about:
PROTEIN_CODING  = 'protein_coding'

# Standard GTF feature types:
CDS_TYPE        = 'CDS'
START_TYPE      = 'start_codon'
STOP_TYPE       = 'stop_codon'
FPU_TYPE        = '5UTR'
TPU_TYPE        = '3UTR'
INTER_TYPE      = 'inter'
INTER_CNS_TYPE  = 'inter_CNS'
INTRON_CNS_TYPE = 'intron_CNS'
EXON_TYPE       = 'exon'
# Additional GTF feature types:
GENE_TYPE       = 'gene'
TRANS_TYPE      = 'transcript'
UTR_TYPE        = 'UTR'
SEC_TYPE        = 'Selenocysteine'
GENE_BIOTYPE    = 'gene_biotype'
# Types given in the GTF 2.2 spec:
VALID_TYPES     = [CDS_TYPE, START_TYPE, STOP_TYPE, FPU_TYPE, TPU_TYPE, INTER_TYPE, INTER_CNS_TYPE, INTRON_CNS_TYPE, EXON_TYPE]
EXTRA_TYPES     = [GENE_TYPE, TRANS_TYPE, UTR_TYPE, SEC_TYPE]
# Types found in EMBL files that we actually care about:
KNOWN_TYPES     = [CDS_TYPE, EXON_TYPE]

# Cufflinks GTF attribute:
TRANSCRIPT      = 'transcript'

# Required GTF 2.2 attributes
GENE_ID_ATTR    = 'gene_id'
TRANS_ID_ATTR   = 'transcript_id'

# EMBL/optional GTF 2.2 attributes
EXON_NO_ATTR    = 'exon_number'
GENE_NAME_ATTR  = 'gene_name'
PROTEIN_ID_ATTR = 'protein_id'
TRANS_NAME_ATTR = 'transcript_name'
CASE_INSENSITIVE = [GENE_ID_ATTR, GENE_NAME_ATTR, PROTEIN_ID_ATTR, TRANS_ID_ATTR, TRANS_NAME_ATTR]

# Attributes we care about:
REQD_ATTRS      = [GENE_ID_ATTR, TRANS_ID_ATTR]
OPT_ATTRS       = [GENE_NAME_ATTR, TRANS_NAME_ATTR]
KNOWN_ATTRS     = REQD_ATTRS + OPT_ATTRS

# Synonyms for optional attributes:
SYNONYM         = {GENE_NAME_ATTR:GENE_ID_ATTR, TRANS_NAME_ATTR:TRANS_ID_ATTR}

COMMENT_START   = '#'

#=====================================================================
# Generic GTF Classes
#=====================================================================
class GTF_Chromosome(object) :
    """Encapsulates a chromosome inferred from GTF record."""
    def __init__(self, name) :
        self.name   = name
        self.genes  = {}
        self.minpos = 1
        self.maxpos = 0

    def append(self, rec) :
        self.maxpos = max(self.maxpos, rec.end())
        self.genes.setdefault(rec.gene_name(), GTF_Gene(rec.gene_name(), rec.gene_id(), self.name, rec.strand()))
        self.genes[rec.gene_name()].append(rec)

    def keys(self) :
        return self.genes.keys()

    def __getitem__(self,k) :
        return self.genes[k]

    def __len__(self) :
        return len(self.genes)

    def gff3String(self) :
        if self.name is None : raise ValueError('Illegal null name for chromosome record')
        # Chr1    TAIR9   chromosome  1   30427671    .   .   .   ID=Chr1;Name=Chr1
        return "%s\tGTF2GFF\tchromosome\t1\t%d\t.\t.\t.\tID=%s;Name=%s" % (self.name,self.maxpos,self.name,self.name)

    def writeGFF3(self, stream) :
        stream.write('%s\n' % self.gff3String())
        genes = sorted(self.genes.values())
        for g in genes :
            if g.valid : g.writeGFF3(stream)

class GTF_Gene(object) :
    def __init__(self, name, id, chrom, strand) :
        self.name     = name
        self.id       = id
        self.chrom    = chrom
        self.strand   = strand
        self.source   = None
        self.geneType = None 
        self.minpos   = MAXINT
        self.maxpos   = 0
        self.valid    = True
        self.xcripts  = {}

    def append(self, rec) :
        ## # For now, ignore non-GTF2.2 types:
        ## if rec.feature() not in VALID_TYPES : return
        self.minpos = min(self.minpos, rec.start())
        self.maxpos = max(self.maxpos, rec.end())
        if not self.source :
            self.source = rec.source()
        elif self.source.lower() != rec.source().lower() :
            self.valid = False
        if rec.strand() != self.strand :
            self.valid = False
        if not self.geneType :
            self.geneType = rec.gene_type()
        tname = rec.transcript_id()
        self.xcripts.setdefault(tname,GTF_Transcript(tname))
        self.xcripts[tname].append(rec)

    def __cmp__(self,other) :
        if self.minpos == other.minpos :
            return self.maxpos - other.maxpos
        else :
            return self.minpos - other.minpos

    def __len__(self) :
        return len(self.records)

    def gff3String(self) :
        if self.name is None : raise ValueError('Illegal null name for gene record')
        # Chr1  TAIR9   gene    3631    5899    .   +   .   ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
        # Note: ENSEMBL uses the GTF source column to store the type of gene (e.g., protein_coding, ncRNA, etc.)
        #       These may not always make sense in the GFF3 "type" column, so we use 'gene' for all types
        #       and use 'Note' to record the ENSEMBL type.
        srcStr = 'protein_coding_gene' if not self.source else self.source
        # Want to use this later to discriminate between 'gene', 'pseudogene', 'snoRNA', etc.
        ##geneStr = self.geneType if self.geneType else 'gene'
        geneStr = 'gene'
        return "%s\tGTF2GFF\t%s\t%d\t%d\t.\t%s\t.\tID=%s;Note=%s;Name=%s" % (self.chrom, geneStr, self.minpos, self.maxpos, self.strand, self.id, srcStr, self.name)

    def writeGFF3(self, stream) :
        stream.write('%s\n' % self.gff3String())
        keys = sorted(self.xcripts.keys())
        for k in keys :
            self.xcripts[k].writeGFF3(stream)

class GTF_Line(object) :
    """Encapsulates a single line/record in a GTF file."""
    def __init__(self, rawstring) :
        ignorepos  = rawstring.find(COMMENT_START)
        s          = rawstring[:ignorepos] if ignorepos >= 0 else rawstring.strip()
        self.parts = s.split('\t')
        if len(self.parts) != ATTR_INDEX+1 :
            raise ValueError('Bad GTF record: had %d columns where %d are required' % (len(self.parts), ATTR_INDEX+1))
        self.attrs = self.setAttributes()

        # Account for cases where either the gene name or gene ID are missing.
        if not (self.attrs[GENE_NAME_ATTR] or self.attrs[GENE_ID_ATTR]) :
            raise ValueError('Missing gene name and gene id')
        elif not self.attrs[GENE_NAME_ATTR] :
            self.attrs[GENE_NAME_ATTR] = self.attrs[GENE_ID_ATTR]
        elif not self.attrs[GENE_ID_ATTR] :
            self.attrs[GENE_ID_ATTR] = self.attrs[GENE_NAME_ATTR]
        # Ensure chromosome is case-insensitive:
        self.parts[SEQNAME_INDEX] = self.parts[SEQNAME_INDEX].lower()

    # Special cases, as these values are guaranteed to exist (see __init__)
    def gene_id(self) :   return self.attrs[GENE_ID_ATTR]
    def gene_name(self) : return self.attrs[GENE_NAME_ATTR]

    # Main column values
    def attributes(self) :      return self.attrs
    def end(self) :             return int(self.parts[END_INDEX])
    def feature(self) :         return self.parts[FEATURE_INDEX]
    def frame_index(self) :     return self.parts[FRAME_INDEX]
    def raw_attributes(self) :  return self.parts[ATTR_INDEX]
    def score(self) :           return self.parts[SCORE_INDEX]
    def seqname(self) :         return self.parts[SEQNAME_INDEX]
    def source(self) :          return self.parts[SOURCE_INDEX]
    def start(self) :           return int(self.parts[START_INDEX])
    def strand(self) :          return self.parts[STRAND_INDEX]

    # Attribute values
    #def gene_type(self) :       return self.getAttribute(GENE_BIOTYPE)
    def gene_type(self) :       pass

    def transcript_id(self) :   return self.getAttribute(TRANS_ID_ATTR)
    def transcript_name(self) : return self.getAttribute(TRANS_NAME_ATTR)

    def attributeString(self) :
        return '; '.join(['%s=%s'%(k,self.attrs[k]) for k in self.attrs])

    def getAttribute(self, key) :
        try :
            return self.attrs[key]
        except KeyError :
            pass

    def setAttributes(self) :
        # From the GTF specification:
        #   "Attributes must end in a semicolon which must then be separated
        #    from the start of any subsequent attribute by exactly one space
        #    character (NOT a tab character)."
        keyvalPairs = self.parts[ATTR_INDEX].split(';')
        result      = {}
        for pair in keyvalPairs :
            if not pair : continue
            parts = pair.split()
            key   = parts[0].strip()
            val   = parts[-1].strip().replace('"','')
            # HACK: semicolons may still appear within
            # values but can wreak havoc on Unix/GFF files
            val   = val.replace(';','-')
            result[key] = val.upper() if key in CASE_INSENSITIVE else val

        # validate attributes:
        missing = set()
        for attr in REQD_ATTRS :
            if attr not in result :
                missing.add(attr)
        if missing :
            raise ValueError("GTF record is missing required attributes '%s'" % "', '".join(missing))

        # set defaults for optional attributes (requires REQD_ATTRS):
        for attr in OPT_ATTRS :
            if (attr not in result) or (result[attr] is None) :
                result[attr] = result[SYNONYM[attr]]

        return result

    def __str__(self) :
        attrString = '; '.join(['%s "%s"'%(k,self.attrs[k]) for k in sorted(self.attrs.keys())])
        return '%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s' % (self.seqname(), self.source(), self.feature(), self.start(), self.end(), self.strand(), attrString)

    def writeGFF3(self, stream) :
        if self.transcript_id() is None : raise ValueError('Illegal null parent name for %s record'%self.feature())
        # Chr1  TAIR9   exon    3631    3913    .   +   .   Parent=AT1G01010.1
        attrString = self.attributeString()
        if attrString : attrString = '; '+attrString
        stream.write("%s\tGTF2GFF\t%s\t%d\t%d\t.\t%s\t.\tParent=%s%s\n" % (self.seqname(), self.feature(), self.start(), self.end(), self.strand(), self.transcript_id(),attrString))

class GTF_Transcript(object) :
    """Encapsulates a single mRNA transcript gleaned from GTF data."""
    def __init__(self, name) :
        self.name    = name
        self.minpos  = MAXINT
        self.maxpos  = 0
        self.records = []
        self.hasCDS  = False
        self.chrom   = None
        self.strand  = None
        self.gene    = None
        self.id      = None

    def append(self, rec) :
        """Appends a feature (exon/codon/etc.)  to the current transcript."""
        if not self.chrom :
            self.chrom  = rec.seqname()
            self.strand = rec.strand()
            self.gene   = rec.gene_name()
            self.id     = rec.gene_id()
        self.minpos = min(self.minpos, rec.start())
        self.maxpos = max(self.maxpos, rec.end())
        if rec.feature() == CDS_TYPE : self.hasCDS = True
        self.records.append(rec)

    def gff3String(self) :
        """Returns a string representation of the transcript CDS record."""
        if self.name is None :
            raise ValueError('Illegal null name for mRNA record %s (%d-%d) %s in gene %s' % (self.chrom, self.minpos, self.maxpos, self.strand, self.gene))

        # Chr1  TAIR9   gene    3631    5899    .   +   .   ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
        return "%s\tGTF2GFF\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s;Name=%s" % \
                (self.chrom, self.minpos, self.maxpos, self.strand, self.name, self.gene, self.name)

    def writeGFF3(self, stream) :
        """Writes the transcript as a sequence of GFF3-formatted records."""
        if self.hasCDS : stream.write('%s\n' % self.gff3String())
        for r in self.records :
            r.writeGFF3(stream)

#=====================================================================
# Cufflinks GTF Classes
#=====================================================================
class GTFRecord(object) :
    """Encapsulates a single record in a Cufflinks GTF file."""
    def __init__(self, s, **args) :
        self.verbose = getAttribute('verbose', False, **args)
        self.cols    = s.strip().split('\t')
        self.attrs   = self.getAttributes(self.cols[ATTR_INDEX])
        for id in REQD_ATTRS : # simply raises a KeyError if a required attribute is missing:
            ignore = self.attrs[id]

    def chromosome(self) :
        return self.cols[SEQNAME_INDEX]

    def __eq__(self, other) :
        return self.chromosome() == other.chromosome() \
                and self.strand() == other.strand() \
                and self.startpos() == other.startpos() \
                and self.endpos() == other.endpos()

    def endpos(self) :
        return int(self.cols[END_INDEX])

    def getAttributes(self, field) :
        parts  = [s.strip() for s in field.split(';')]
        result = dict([p.split() for p in parts if len(p) > 0])
        if self.verbose : fpkm = False
        for k in result.keys() :
            if self.verbose : fpkm |= (k == 'FPKM')
            result[k] = result[k].replace('"','')
        if self.verbose and not fpkm : sys.stderr.write('%s HAS NO FPKM ATTRIBUTE\n' % result[GENE_ID_ATTR])
        return result

    def startpos(self) :
        return int(self.cols[START_INDEX])

    def __str__(self) :
        return '\t'.join(self.cols)

    def strand(self) :
        return self.cols[STRAND_INDEX]

class GTFParser(object) :
    """
    Parses a GTF file filled with isoform transcripts to produce a
    splice graph iterator through a Cufflinks GTF file.
    """
    def __init__(self, path, **args) :
        if type(path) == type('') :
            self.instream = open(path, 'r')
        else :
            self.instream = path

        self.verbose   = getAttribute('verbose', False, **args)
        self.graphDict = {}
        self.graphIds  = []
        self.loadFromFile()
        self.finishGraphs(**args)

    def __iter__(self) :
        return self

    def next(self) :
        """Iterator implementation."""
        try :
            return self.graphDict[self.graphIter.next()]
        except Exception :
            raise StopIteration

    def finishGraphs(self, **args) :
        """Completes all graphs and annotates them."""
        graphs  = self.graphDict.values()
        for g in graphs :
            g.annotate()

    def loadFromFile(self) :
        """
        Loads a splice graph from a GTF file.
        """
        lineNo        = 0
        graph         = None
        currentGene   = None
        currentForm   = None
        currentNodeId = None
        exonDict      = {}
        nodeCtr       = 0
        edges         = set([])
        indicator = ProgressIndicator(100000, verbose=self.verbose)
        for line in self.instream :
            lineNo += 1
            indicator.update()
            if line.startswith('#') : continue
            s       = line.rstrip()
            rec     = GTFRecord(s, verbose=self.verbose)
            feature = rec.cols[FEATURE_INDEX].lower()
            geneId  = rec.attrs[GENE_ID_ATTR]
            formId  = rec.attrs[TRANS_ID_ATTR]

            # Infer new gene name when it changes (no gene records in GTF)
            if geneId.lower() != currentGene :
                if graph :
                    for edge in edges :
                        graph.addEdge(edge[0], edge[1])
                graph        = SpliceGraph(name=geneId, chromosome=rec.chromosome(), strand=rec.strand())
                graph.minpos = min(rec.startpos(), rec.endpos())
                graph.maxpos = max(rec.startpos(), rec.endpos())
                self.graphIds.append(geneId)
                self.graphDict[geneId] = graph
                currentNodeId = None
                currentGene  = geneId.lower()
                nodeCtr      = 0
                edges        = set([])
            elif graph is None :
                raise ValueError("Graph feature found before graph header at line %d" % lineNo)

            # Start a new transcript when necessary
            if feature == TRANSCRIPT :
                currentNodeId = None
                currentForm   = formId
                continue

            nodeCtr += 1
            nodeId   = "%s-%d" % (currentForm,nodeCtr)
            node     = graph.addNode(nodeId, rec.startpos(), rec.endpos())
            node.addIsoform(currentForm)

            for k in rec.attrs :
                node.addAttribute(k, rec.attrs[k])

            # All records represent a node, but some will be duplicates
            if currentNodeId :
                if rec.strand() == '+' :
                    edges.add((currentNodeId, node.id))
                else :
                    edges.add((node.id, currentNodeId))
            currentNodeId = node.id
            
        indicator.finish()
        self.graphIter = iter(self.graphIds)

#=====================================================================
# Functions
#=====================================================================
def keyString(rec) :
    """Create a key that is unique for chromosome, strand, start and end position. """
    return "%s;%s;%d;%d" % (rec.chromosome(), rec.strand(), rec.startpos(), rec.endpos())

def getFirstCufflinksGraph(f) :
    """Returns the first graph in a Cufflinks GTF file."""
    return GTFParser(f).next()
