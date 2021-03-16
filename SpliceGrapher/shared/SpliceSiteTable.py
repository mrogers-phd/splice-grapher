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
from SpliceGrapher.formats.GeneModel import *

class SpliceSiteTable(object):
    """A table that provides quick access to all splice sites in a gene model."""
    def __init__(self, geneModel):
        self.donors    = {}
        self.acceptors = {}
        self.junctions = {}
        allGenes       = geneModel.getAllGenes(geneFilter=gene_type_filter)
        for g in allGenes:
            c = g.chromosome.lower()
            s = g.strand.lower()
            self.donors.setdefault(c,{})
            self.acceptors.setdefault(c,{})
            self.junctions.setdefault(c, {})

            self.donors[c].setdefault(s,{})
            self.acceptors[c].setdefault(s,{})
            self.junctions[c].setdefault(s,{})

            for iid in list(g.isoforms.keys()):
                exons = g.isoforms[iid].sortedExons()
                for i in range(1,len(exons)):
                    d = exons[i-1].donor()
                    a = exons[i].acceptor()
                    self.donors[c][s].setdefault(d,set([]))
                    self.acceptors[c][s].setdefault(a,set([]))
                    self.junctions[c][s].setdefault(d,{})

                    self.donors[c][s][d].add(g)
                    self.acceptors[c][s][a].add(g)
                    self.junctions[c][s][d][a] = True

            for mid in list(g.mrna.keys()):
                cds = g.mrna[mid].sortedExons()
                for i in range(1,len(cds)):
                    d = cds[i-1].donor()
                    a = cds[i].acceptor()
                    self.donors[c][s].setdefault(d,set([]))
                    self.acceptors[c][s].setdefault(a,set([]))
                    self.junctions[c][s].setdefault(d,{})

                    self.donors[c][s][d].add(g)
                    self.acceptors[c][s][a].add(g)
                    self.junctions[c][s][d][a] = True

    def acceptorGenes(self, chrom, strand, a):
        """Returns all genes that contain the given acceptor site."""
        try:
            return self.acceptors[chrom.lower()][strand.lower()][a]
        except KeyError:
            return set([])

    def acceptorSites(self, chrom, strand):
        """Returns an unsorted list of all acceptors for the given chromosome and strand."""
        return list(self.acceptors[chrom.lower()][strand.lower()].keys())

    def donorSites(self, chrom, strand):
        """Returns an unsorted list of all donors for the given chromosome and strand."""
        return list(self.donors[chrom.lower()][strand.lower()].keys())

    def donorGenes(self, chrom, strand, d):
        """Returns all genes that contain the given donor site."""
        try:
            return self.donors[chrom.lower()][strand.lower()][d]
        except KeyError:
            return set([])

    def isAcceptor(self, chrom, strand, pos):
        """Returns True if the given position, chromosome and strand
        represents an acceptor site; False otherwise."""
        try:
            return len(self.acceptorGenes(chrom,strand,pos) > 0)
        except KeyError:
            return False

    def isDonor(self, chrom, strand, pos):
        """Returns True if the given position, chromosome and strand
        represents a donor site; False otherwise."""
        try:
            return len(self.donorGenes(chrom,strand,pos) > 0)
        except KeyError:
            return False

    def knownJunction(self, chrom, strand, d, a):
        """Returns True if the given donor and acceptor represent
        a known junction for the given chromosome and strand; False otherwise."""
        try:
            return self.junctions[chrom.lower()][strand.lower()][d][a]
        except KeyError:
            return False

    def recombinedJunction(self, chrom, strand, d, a):
        """Returns True if the given donor and acceptor are found
        in the same gene; False otherwise."""
        try:
            D = self.donorGenes(chrom,strand,d)
            A = self.acceptorGenes(chrom,strand,a)
            correctOrder = (d < a) if strand == '+' else (d > a)
            return (len(D.intersection(A)) > 0) and correctOrder
        except KeyError:
            return False
