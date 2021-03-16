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
from SpliceGrapher.formats.FastaLoader import FastaLoader
from SpliceGrapher.predict.SpliceSite  import *

class SpliceSiteValidator(object) :
    """
    Class that validates splice sites based on the given site type.
    """
    def __init__(self, dbPath, **args) :
        self.verbose = ('verbose' in args and args['verbose'])
        self.loader  = FastaLoader(dbPath, verbose=self.verbose)

    def getAcceptorDimer(self, chrom, pos, strand, window=0) :
        """Convenience method that returns the acceptor dimer at the given location."""
        return self.getDimer(chrom, pos, strand, ACCEPTOR_SITE, window=window)

    def getChromosomes(self) :
        """Returns the chromosomes found in the FASTA file."""
        return self.loader.keys()

    def getDimer(self, chrom, pos, strand, siteType, window=0) :
        """Returns the given splice site dimer."""
        if not siteType in SITE_TYPES :
            raise ValueError('Given site type was %s; must be one of %s' % (siteType,'/'.join(SITE_TYPES)))

        (estr, istr, dimer) = getSpliceSiteParts(chrom, pos, siteType, strand, self.getSequence, exonWindow=window, intronWindow=window)
        if window == 0 :
            return dimer
        else :
            values = (estr,dimer,istr) if siteType == DONOR_SITE else (istr,dimer,estr)
            return '%s|%s|%s' % values

    def getDonorDimer(self, chrom, pos, strand, window=0) :
        """Convenience method that returns the donor dimer at the given location."""
        return self.getDimer(chrom, pos, strand, DONOR_SITE, window=window)

    def getSequence(self, chrom, pos1, pos2, strand) :
        """Method used by getSpliceSiteParts (see SpliceSite module)"""
        return self.loader.subsequence(chrom, pos1, pos2, reverse=(strand=='-'))

    def validAcceptor(self, chrom, pos, strand, dimer=None) :
        """Returns true if the given splice site is a valid acceptor; false otherwise."""
        return self.validSite(chrom, pos, strand, ACCEPTOR_SITE, dimer=dimer)

    def validDonor(self, chrom, pos, strand, dimer=None) :
        """Returns true if the given splice site is a valid donor with the given dimer; false otherwise."""
        return self.validSite(chrom, pos, strand, DONOR_SITE, dimer=dimer)

    def validSite(self, chrom, pos, strand, siteType, dimer=None) :
        """Returns true if the given splice site has the given dimer; false otherwise."""
        if not siteType in SITE_TYPES :
            raise ValueError('Given site type was %s; must be one of %s' % (siteType,'/'.join(SITE_TYPES)))

        (estr, istr, dstr) = getSpliceSiteParts(chrom, pos, siteType, strand, self.getSequence, exonWindow=1, intronWindow=1)

        if dimer :
            return dstr.upper() == dimer.upper()
        elif siteType in KNOWN_DIMERS :
            return dstr.upper() in KNOWN_DIMERS[siteType]
        else :
            return False
