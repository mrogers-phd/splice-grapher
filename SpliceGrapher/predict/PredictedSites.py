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
from SpliceGrapher.shared.utils       import ProgressIndicator, getAttribute, ezopen
from SpliceGrapher.predict.SpliceSite import *
from SpliceGrapher.formats.GeneModel  import *

import os, sys

ACCEPTOR = 'a'
DONOR    = 'd'

DEFAULT_POSITIVE = 99.0
DEFAULT_NEGATIVE = -99.0

class PredictedSites(object) :
    """
    Class that encapsulates a set of predicted splice sites loaded from a
    file of prediction values.  The file must have the following format as
    output by SpliceGrapher prediction software:

          chromosome	strand	position	dimer	score	[a/d]
    """
    def __init__(self, dataFile, **args) :
        self.verbose = getAttribute('verbose', False, **args)
        if self.verbose : sys.stderr.write('Loading predicted splice site database from %s\n' % dataFile)
        indicator = ProgressIndicator(1000000, verbose=self.verbose)

        self.sites     = {}
        lineCtr    = 0
        for line in ezopen(dataFile) :
            indicator.update()
            lineCtr += 1
            # 0             1       2           3       4       5
            # chromosome	strand	position	dimer	score	[a/d]
            parts    = line.strip().split('\t')
            chrom    = parts[0].lower()
            strand   = parts[1]
            pos      = int(parts[2])
            score    = float(parts[4])

            siteType = parts[-1].lower()
            if not siteType in [ACCEPTOR,DONOR] :
                raise Exception('Site type %s is not a or d at line %d' % (siteType,lineCtr))

            self.sites.setdefault(chrom,{})
            self.sites[chrom].setdefault(strand,{})
            self.sites[chrom][strand].setdefault(siteType,{})
            self.sites[chrom][strand][siteType][pos] = score
        indicator.finish()

        # Sort indexes to permit searching
        self.donors    = {}
        self.acceptors = {}
        for c in self.sites :
            self.donors[c]    = {}
            self.acceptors[c] = {}
            for s in ['-','+'] :
                try :
                    self.donors[c][s] = sorted(self.sites[c][s][DONOR].keys())
                except KeyError :
                    self.acceptors[c][s] = []
                try :
                    self.acceptors[c][s] = sorted(self.sites[c][s][ACCEPTOR].keys())
                except KeyError :
                    self.donors[c][s] = []

    def acceptorScore(self, chrom, strand, pos) :
        """
        Returns the classifier score for the given acceptor site, if it was
        predicted, or -99.0 if it was not.
        """
        return self.siteScore(chrom, strand, pos, ACCEPTOR)

    def addGeneModel(self, geneModel, **args) :
        """
        Allows a client to add known sites from a gene model.
        """
        geneFilter = getAttribute('geneFilter', gene_type_filter, **args)
        siteDict   = { DONOR:geneModel.getAllDonors(geneFilter=geneFilter),
                       ACCEPTOR:geneModel.getAllAcceptors(geneFilter=geneFilter) }
        for siteType in [DONOR,ACCEPTOR] :
            for chrom in siteDict[siteType] :
                self.sites.setdefault(chrom,{})
                for strand in siteDict[siteType][chrom] :
                    self.sites[chrom].setdefault(strand,{})
                    self.sites[chrom][strand].setdefault(siteType,{})
                    for pos in siteDict[siteType][chrom][strand] :
                        if type(pos) != int : raise ValueError('Invalid splice site %s (type %s) in gene model; should be type int' % (str(pos), type(pos)))
                        self.sites[chrom][strand][siteType][pos] = DEFAULT_POSITIVE

    def donorScore(self, chrom, strand, pos) :
        """
        Returns the classifier score for the given donor site, if it was
        predicted, or -99.0 if it was not.
        """
        return self.siteScore(chrom, strand, pos, DONOR)

    def isPredictedAcceptor(self, chrom, strand, pos, threshold=0.0) :
        """
        Returns true if the given acceptor site was predicted; false otherwise.
        """
        return self.isPredicted(chrom, strand, pos, ACCEPTOR, threshold)

    def isPredictedDonor(self, chrom, strand, pos, threshold=0.0) :
        """
        Returns true if the given donor site was predicted; false otherwise.
        """
        return self.isPredicted(chrom, strand, pos, DONOR, threshold)

    def isPredicted(self, chrom, strand, pos, siteType, threshold=0.0) :
        """
        Returns true if the given site was predicted; false otherwise.
        """
        try :
            score = self.sites[chrom.lower()][strand][siteType][pos]
            return (score > threshold)
        except KeyError :
            pass
        return False

    def siteScore(self, chrom, strand, pos, siteType) :
        """
        Returns the classifier score for the given site, if it was
        predicted, or -99.0 if it was not.
        """
        try :
            return self.sites[chrom.lower()][strand][siteType][pos]
        except KeyError :
            pass
        return DEFAULT_NEGATIVE

class ChromosomeSiteMap(object) :

    def __init__(self, path, **args) :
        """
        Stores a 3-D dictionary of splice sites for a single chromsome
        indexed by site type, strand and position.

        :parameters:
            'offset'  - offset for positions (default=0)
        """
        self.chromosome  = None
        self.path        = path
        self.siteTypes   = {}
        self.siteDict    = {}
        self.entryCount  = 0
        self.verbose     = getAttribute('verbose', False, **args)

        offset = getAttribute('offset', 0, **args)
        self.loadSites(path, offset=offset)

    def loadSites(self, path, offset=0) :
        """Loads positions from a file in which each record has the format:
            chromosome	strand	position	dimer	score	site-code
        """
        ssCount = 0
        indicator = ProgressIndicator(1000000, verbose=self.verbose)
        for line in ezopen(path) :
            vals     = line.strip().split('\t')
            ssCount += 1
            if len(vals) != 6 :
                raise ValueError('Expected 6 columns in file but found %d at line %d' % (len(vals),ssCount))

            indicator.update()

            chrom = vals[0].lower()
            if self.chromosome :
                if chrom != self.chromosome :
                    raise ValueError('All records in the file must have the same chromosome: %s != %s at line %d' % (chrom, self.chromosome, ssCount))
            else :
                self.chromosome = chrom

            strand = vals[1]
            try :
                siteType = SITE_TYPE_NAME[vals[5]]
            except KeyError :
                raise ValueError('Invalid site-type code at line %d (was %s, should be %s)' % (ssCount, vals[5], '/'.join(SITE_TYPE_NAME.keys())))

            pos    = int(vals[2]) + offset
            score  = float(vals[4])

            self.siteTypes[siteType] = 1
            self.siteDict[siteType]  = self.siteDict.setdefault(siteType, {'-':{}, '+':{}})
            self.siteDict[siteType][strand][pos] = score
            self.entryCount += 1

        indicator.finish()
        if self.verbose : sys.stderr.write('Loaded %d predicted splice sites from %s\n' % (ssCount, path))

    def getAcceptors(self, strand, minpos=0, maxpos=10**10, minScore=0.0) :
        """Convenience method that returns just acceptor sites."""
        return self.getSites(ACCEPTOR, strand, minpos, maxpos, minScore)

    def getAllsites(self, minpos=0, maxpos=10**10, minScore=0.0) :
        """Convenience method that returns all sites on both strands."""
        result = []
        for sType in self.siteTypes.keys() :
            for strand in self.siteTypes[sType].keys() :
                result += self.getSites(sType, strand, minpos, maxpos, minScore)
        return result

    def getDonors(self, strand, minpos=0, maxpos=10**10, minScore=0.0) :
        """Convenience method that returns just donor sites."""
        return self.getSites(DONOR, strand, minpos, maxpos, minScore)

    def getScores(self) :
        """Returns a list of all scores in the dictionary."""
        result = []
        for a in self.siteDict :
            for b in self.siteDict[a] :
                for c in self.siteDict[a][b] :
                    result.append(self.siteDict[a][b][c])
        return result

    def getSites(self, siteType, strand, minpos=0, maxpos=10**10, minScore=0.0) :
        """Return a list of predicted sites for the given site type and strand."""
        sType = siteType.lower()
        if sType in SITE_CODES :
            sType = SITE_TYPE_NAME[sType]
        elif sType not in SITE_TYPES :
            raise ValueError("Unrecognized site type '%s' -- must be a recognized site code (%s) or site type (%s)" % (siteType, '/'.join(SITE_CODES), '/'.join(SITE_TYPES)))

        result = []
        if sType in self.siteDict :
            result = [x for x in self.siteDict[sType][strand].keys() if minpos <= x <= maxpos and self.siteDict[sType][strand][x] >= minScore]
        result.sort()
        return result

    def getSiteTypes(self) :
        """Returns the site types represented in the container."""
        return sorted(self.siteTypes.keys())

    def getSiteDict(self) :
        """
        Returns a 3-D dictionary of splice sites indexed by
        site type, strand and position.
        """
        return self.siteDict

    def size(self) :
        """Returns the total number of entries in the dictionary."""
        return self.entryCount
