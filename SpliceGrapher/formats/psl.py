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
Module that loads and converts PSL-formatted EST alignments.
"""
from SpliceGrapher.SpliceGraph        import SpliceGraph, updateRoot, updateLeaf, donor, acceptor
from SpliceGrapher.formats.GeneModel  import Exon
from SpliceGrapher.shared.utils       import *
from SpliceGrapher.predict.SpliceSiteValidator import *
from sys import maxsize

def estsToSpliceGraph(geneName, recs, chromosome, **args) :
    """Method that creates a splice graph from a list of PSLRecord instances.  Single-exon
    alignments may have a bearing on intron retention events, but are excluded by default."""
    geneBounds     = getAttribute('geneBounds', None, **args)
    includeSingles = getAttribute('includeSingles', False, **args)
    validator      = getAttribute('validator', None, **args)
    min_intron     = getAttribute('min_intron', 4, **args)
    strand         = getAttribute('strand', None, **args)
    verbose        = getAttribute('verbose', False, **args)

    if verbose : import sys
    
    exonIds  = idFactory('%s_' % geneName)
    if not strand : strand = inferStrand(recs)
    graph    = SpliceGraph(geneName, chromosome, strand)
    geneMin  = min(geneBounds) if geneBounds else 0
    geneMax  = max(geneBounds) if geneBounds else maxint
    paths    = []
    exonDict = {}
    ctr      = 0
    for r in recs :
        exonList = psl_to_exons(r, min_intron=min_intron)

        # Only accept exons that overlap gene boundaries:
        exons = [e for e in exonList if e.maxpos > geneMin and e.minpos < geneMax]

        if len(exons) == 1 and not includeSingles : continue

        if not exons :
            if verbose :
                sys.stderr.write('None of these exons from %s overlap %s (%d-%d):\n' % (r.Qname, geneName, geneMin, geneMax))
                for e in exonList :
                    sys.stderr.write('   [%d-%d]\n' % (e.minpos, e.maxpos))
            continue

        # reverse exon order to match graph strand:
        if r.strand != strand :
            exons = exons[::-1]
            for e in exons : e.strand = strand

        subgraph = SpliceGraph('tmp', chromosome, strand)
        prev     = None
        for i in range(len(exons)) :
            exon = exons[i]
            eid  = next(exonIds)
            node = subgraph.addNode(eid, exon.minpos, exon.maxpos)
            #assert(node.id == eid)
            node.addIsoform(r.Qname)
            if prev : subgraph.addEdge(prev, eid)
            prev = eid

        # If a FASTA reference is provided, verify that splice sites
        # are valid according to the given dimer dictionary
        goodSpliceSites = True
        if validator :
            for e in list(subgraph.nodeDict.values()) :
                if not e.isLeaf() :
                    goodSpliceSites &= validator.validDonor(chromosome, donor(e), strand)
                if not e.isRoot() :
                    goodSpliceSites &= validator.validAcceptor(chromosome, acceptor(e), strand)
                if not goodSpliceSites : break

        if not goodSpliceSites :
            if verbose : sys.stderr.write('Omitting subgraph for %s: contains invalid splice site\n' % r.Qname)
            continue

        if verbose :
            ctr += 1
            sys.stderr.write('Graph %d (%s strand) from %s:\n' % (ctr,r.strand,r.Qname))
            sys.stderr.write('    %s\n' % '\n    '.join(['%s'%x for x in list(subgraph.nodeDict.values())]))

        if geneBounds :
            if subgraph.maxpos < geneMin or subgraph.minpos > geneMax :
                if verbose :
                    sys.stderr.write('Omitting subgraph for %s: outside gene boundaries (%d < %d or %d > %d)\n' % (r.Qname,subgraph.maxpos, geneMin, subgraph.minpos, geneMax))
                continue
            
        # Must merge nodes in both directions since
        # PSL records are in no particular order:
        updateRoot(subgraph, graph)
        updateLeaf(subgraph, graph)
        updateRoot(graph, subgraph)
        updateLeaf(graph, subgraph)

        graph = graph.union(subgraph, keepName=True, verbose=True)
        graph.validate(halt=True)

    return graph

def inferStrand(recs) :
    """Method to infer the correct strand from a list of PSL records using a simple majority."""
    plusCount  = sum([1 for r in recs if r.strand == '+'])
    minusCount = sum([1 for r in recs if r.strand == '-'])
    return ['+','-'][minusCount>plusCount]

def loadPSLRecords(f) :
    """Loads PSL records from a file and returns them in a list."""
    hstring   = None
    result    = []
    pslStream = ezopen(f)
    for s in pslStream :
        if not s.startswith('@') :
            result.append(PSLRecord(s))
    return result

def loadPSLDict(f) :
    """Loads PSL records from a file and returns them in a dictionary keyed by target name."""
    records = loadPSLRecords(f)
    result  = {}
    for rec in records :
        result[rec.Tname] = result.setdefault(rec.Tname, [])
        result[rec.Tname].append(rec)
    return result

#
# Original class and PSL parsing written by Adam Labadorf
# Modified for SpliceGrapher by Mark Rogers
#
class PSLRecord(object) :
    """Class representing the data object of a line in a PSL-formatted file.  Data fields follow field
    names specified in `<http://genome.ucsc.edu/FAQ/FAQformat#format2>`_ .  If no argument is supplied
    to the constructor all the fields are initialized to empty strings.  This is useful for constructing
    a PSL file from some other data format."""

    FIELDS = ['match','mismatch','repmatch','Ns', \
              'Qgapcnt','Qgapbas','Tgapcnt','Tgapbas','strand','Qname', \
              'Qsize','Qstart','Qend','Tname','Tsize','Tstart','Tend', \
              'blkcnt','blksize','qstarts','tstarts']#,'splicescores']

    def __init__(self,line=None) :
        if line is None :
            line = [''] * len(PSLRecord.FIELDS)

        if type(line) == str :
            line = line.split()

        # indices for easy indexing
        match,mismatch,repmatch,Ns, \
        Qgapcnt,Qgapbas,Tgapcnt,Tgapbas,strand,Qname, \
        Qsize,Qstart,Qend,Tname,Tsize,Tstart,Tend, \
        blkcnt,blksize,qstarts,tstarts = list(range(len(PSLRecord.FIELDS)))

        try :
            for i,f in enumerate(PSLRecord.FIELDS) :
                if i in [strand,Qname,Tname] :
                    self.__dict__[f] = line[i]
                elif i in [blksize,qstarts,tstarts] :
                    if type(line[i]) in [list,tuple] :
                        self.__dict__[f] = line[i]
                    elif type(line[i]) is str :
                        self.__dict__[f] = [to_numeric(x) for x in line[i].split(',')[0:-1]]
                    else :
                        raise Exception("PSLRecord constructor received invalid input for %s field: %s"%(f,str(line[i])))
                else :
                    self.__dict__[f] = to_numeric(line[i])
        except :
            raise Exception("Error initializing PSLRecord object from input: %s" % repr(line))

    def tolist(self) :
        """Returns datafields as a  list"""
        return [self.__dict__[f] for f in PSLRecord.FIELDS]

    def __repr__(self) :
        """Returns string representation of the GFFLine, usable with `eval()`"""
        return "PSLRecord(%s)" % self.tolist()

    def __eq__(self,o) :
        try :
            return all([self.__dict__[f]==o.__dict__[f] for f in PSLRecord.FIELDS])
        except AttributeError as e :
            raise Exception('Argument to PSLRecord.__cmp__(o) must be a PSLRecord instance. Exception: %s'%e)

    def __cmp__(self,o) :
        try :
            c = cmp(self.Tname,o.Tname)
            if c == 0 :
                return cmp(self.Tstart,o.Tstart)
            return c
        except AttributeError as e :
            raise Exception('Argument to PSLRecord.__cmp__(o) must be a PSLRecord instance. Exception: %s'%e)

    def __hash__(self) :
        try :
            res = self._hash
        except AttributeError :
            res = self._hash = hash(self.output_format())
        return res

    def list_format(self) :
        """Returns string representation in compact list format"""
        ret = ""
        for f in PSLRecord.FIELDS :
            ret += "%s:%s, "%(f,self.__dict__[f])
        return ret

    def print_format(self) :
        """Returns string representation in human-readable(ish) format"""
        ret = ""
        for f in PSLRecord.FIELDS :
            ret += "%s:%s\n"%(f,self.__dict__[f])
        return ret

    def output_format(self) :
        """Returns string representation as it would exist in a PSL file"""
        line = []
        for i,f in enumerate(PSLRecord.FIELDS) :
            if f in ['blksize','qstarts','tstarts'] :
                line.append("".join(["%s,"%x for x in self.__dict__[f]]))
            else :
                line.append(str(self.__dict__[f]))
        return "\t".join(line)

#
# Original method written by Adam Labadorf to convert PSL records
# into GFF output.  Modified for SpliceGrapher by Mark Rogers
#
def psl_to_exons(psl, min_intron=0) :
    """Convert a PSLRecord to a list of exons.  We defer to the reference
    genome, so we ignore all insertions in the query sequence.  Exons that
    that flank introns smaller than the min_intron value will be merged."""

    qry_coords = list(zip(psl.qstarts, [x+b for x,b in zip(psl.qstarts,psl.blksize)] ))
    qry_coords.extend(list(zip([x+b for x,b in zip(psl.qstarts[:-1],psl.blksize)], psl.qstarts[1:])))
    qry_coords.sort()

    genm_coords = list(zip(psl.tstarts, [x+b for x,b in zip(psl.tstarts,psl.blksize)]))
    genm_coords.extend(list(zip([x+b for x,b in zip(psl.tstarts[:-1],psl.blksize)], psl.tstarts[1:])))
    genm_coords.sort()

    if len(qry_coords) != len(genm_coords) :
        raise Exception('Illegal PSL format: length of qstarts and tstarts must agree in PSL record %s'%repr(psl))

    exons = []
    for i in range(len(qry_coords)) :
        qry_span  = qry_coords[i][1]-qry_coords[i][0]
        genm_span = genm_coords[i][1]-genm_coords[i][0]
        delta     = genm_span-qry_span

        if qry_span < 0  : raise Exception('Illegal query alignment (negative span) from %d to %d in PSL file' % (qry_coords[i][0],qry_coords[i][1]))
        if qry_span == 0 or delta < 0 : continue

        startpos = genm_coords[i][0]
        exons.append(Exon(startpos+1, startpos+qry_span, psl.Tname, psl.strand))

    # put exons in 5'-to-3' order prior to collapsing
    exons.sort(reverse=(psl.strand=='-'))

    # collapse consecutive match elements
    if min_intron > 0 :
        for i in range(1,len(exons)) :
            if abs(exons[i].start() - exons[i-1].end()) < min_intron :
                exons[i]   = Exon(exons[i-1].start(), exons[i].end(), psl.Tname, psl.strand)
                exons[i-1] = None # don't mutate list, just mark for removal

    return [e for e in exons if e is not None]
