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
Module that encapsulates a configuration for a splice-site dimer classifier.
"""
from SpliceGrapher.shared.utils import getAttribute
from configparser import ConfigParser
import sys

# Configuration section name
CONFIG_SECTION    = 'SpliceGrapherConfig'

# Individual configuration parameters
ACCEPTOR_TAG      = 'acceptor'
C_TAG             = 'c'
FASTA_PATH        = 'fasta_path'
DIMER_TAG         = 'dimer'
EXON_TAG          = 'exon_size'
INTRON_TAG        = 'intron_size'
MINK_TAG          = 'mink'
MAXK_TAG          = 'maxk'
MISMATCH_TAG      = 'mismatches'
MPROFILE_TAG      = 'mismatch_profile'
NOSHIFT_START_TAG = 'noshift_start'
NOSHIFT_END_TAG   = 'noshift_end'
NORMALIZE_TAG     = 'normalize'
SHIFT_TAG         = 'shift_size'
SVM_PATH          = 'svm_path'
THRESHOLD_TAG     = 'threshold'

# Configuration section name
PERFORMANCE_SECTION = 'ClassifierPerformance'
ROC_TAG             = 'roc_score'
PERFORMANCE_TAGS    = [ROC_TAG]

# Parameter lists
STRING_TAGS  = [FASTA_PATH, DIMER_TAG, SVM_PATH]
FLOAT_TAGS   = [C_TAG, THRESHOLD_TAG]
INTEGER_TAGS = [EXON_TAG, INTRON_TAG, MINK_TAG, MAXK_TAG, SHIFT_TAG, NOSHIFT_START_TAG, NOSHIFT_END_TAG, MISMATCH_TAG]
BOOLEAN_TAGS = [ACCEPTOR_TAG, MPROFILE_TAG, NORMALIZE_TAG]
CONFIG_TAGS  = STRING_TAGS + FLOAT_TAGS + INTEGER_TAGS + BOOLEAN_TAGS
ALL_TAGS     = {CONFIG_SECTION:CONFIG_TAGS, PERFORMANCE_SECTION:PERFORMANCE_TAGS}

class ClassifierConfig :
    """
    Class for reading/writing classifier configurations that contain all the
    parameters for a splice site dimer classifier.
    """
    def __init__(self, **args) :
        self.comments = []
        self.verbose  = getAttribute('verbose', False, **args)
        self.config   = ConfigParser()

        # Internal work-around since ConfigParser fails ungracefully
        # if users try to retrieve values they've just set.
        self.map      = {}

        # Initialize configuration before anything else,
        # to ensure all values are represented.
        self.config.add_section(CONFIG_SECTION)
        self.map[CONFIG_SECTION] = {}
        for tag in STRING_TAGS  : self.setValue(tag,'')
        for tag in FLOAT_TAGS   : self.setValue(tag,0.0)
        for tag in INTEGER_TAGS : self.setValue(tag,0)
        for tag in BOOLEAN_TAGS : self.setValue(tag,False)

        self.config.add_section(PERFORMANCE_SECTION)
        self.map[PERFORMANCE_SECTION] = {}
        self.setROC(0.0)

        # Note that this will not update the local map, but in this
        # case the workaround is not needed.
        if 'file' in args :
            self.config.read(args['file'])

    # Convenience getter methods:
    def acceptor(self)        : return self.getValue(ACCEPTOR_TAG, default=False)
    def fastaPath(self)       : return self.getValue(FASTA_PATH)
    def dimer(self)           : return self.getValue(DIMER_TAG)
    def C(self)               : return self.getValue(C_TAG)
    def exonSize(self)        : return self.getValue(EXON_TAG)
    def intronSize(self)      : return self.getValue(INTRON_TAG)
    def mink(self)            : return self.getValue(MINK_TAG)
    def maxk(self)            : return self.getValue(MAXK_TAG)
    def maxShift(self)        : return self.getValue(SHIFT_TAG)
    def mismatches(self)      : return self.getValue(MISMATCH_TAG)
    def mismatchProfile(self) : return self.getValue(MPROFILE_TAG)
    def normalize(self)       : return self.getValue(NORMALIZE_TAG, default=False)
    def noShiftStart(self)    : return self.getValue(NOSHIFT_START_TAG)
    def noShiftEnd(self)      : return self.getValue(NOSHIFT_END_TAG)
    def roc(self)             : return self.getValue(ROC_TAG, default=0.0, section=PERFORMANCE_SECTION)
    def SVMPath(self)         : return self.getValue(SVM_PATH)
    def threshold(self)       : return self.getValue(THRESHOLD_TAG, default=0.0)

    # Convenience setter methods:
    def setAcceptor(self, b)        : self.setValue(ACCEPTOR_TAG, b)
    def setC(self, x)               : self.setValue(C_TAG, x)
    def setFastaPath(self, s)       : self.setValue(FASTA_PATH, s)
    def setExonSize(self, k)        : self.setValue(EXON_TAG, k)
    def setIntronSize(self, k)      : self.setValue(INTRON_TAG, k)
    def setMink(self, k)            : self.setValue(MINK_TAG, k)
    def setMaxk(self, k)            : self.setValue(MAXK_TAG, k)
    def setMaxShift(self, k)        : self.setValue(SHIFT_TAG, k)
    def setMismatches(self, k)      : self.setValue(MISMATCH_TAG, k)
    def setNoShiftStart(self, k)    : self.setValue(NOSHIFT_START_TAG, k)
    def setNoShiftEnd(self, k)      : self.setValue(NOSHIFT_END_TAG, k)
    def setMismatchProfile(self, b) : self.setValue(MPROFILE_TAG, b)
    def setNormalize(self, b)       : self.setValue(NORMALIZE_TAG, b)
    def setROC(self,x)              : self.setValue(ROC_TAG, x, section=PERFORMANCE_SECTION)
    def setSVMPath(self, s)         : self.setValue(SVM_PATH, s)
    def setThreshold(self, x)       : self.setValue(THRESHOLD_TAG, x)

    def addComment(self, s) :
        """Adds a comment to a new configuration.  Note that these will be
        initialized to an empty list when an existing config file is read."""
        self.comments.append(s)

    def getComments(self) :
        """Returns any comments added during a session.  Note that the ConfigParser
        interface does not provide access to comments in a file, so these will only
        be the comments added during creation of a config file."""
        return self.comments

    def setDimer(self, s) :
        """Set the splice-site dimer for the configuration.  Raises a ValueError
        exception if the dimer is anything other than 2 characters."""
        if len(s) != 2 : raise ValueError('Attempted to set a dimer that was not 2nt: %s (%dnt)' % (s,len(s)))
        self.setValue(DIMER_TAG, s)

    # Workhorse methods:
    def getMapValue(self, name, default=None, section=CONFIG_SECTION) :
        """Returns a value stored in the local map.  This is a workaround
        for defects in the ConfigParser module."""
        try :
            return self.map[section][name]
        except KeyError :
            return default

    def getValue(self, name, default=None, section=CONFIG_SECTION) :
        """Returns a value from the configuration."""
        tag = name.lower()
        if name not in ALL_TAGS[section] :
            raise ValueError("Unrecognized %s tag %s" % (section,name))

        try :
            # When ConfigParser fails, it usually throws a TypeError
            # from 'way down in its bowels, but we catch them all.
            value = self.config.get(section,name)
        except Exception :
            value = self.getMapValue(name)

        try :
            if tag in FLOAT_TAGS :
                return float(value)
            elif tag in INTEGER_TAGS :
                return int(value)
            elif tag in BOOLEAN_TAGS :
                return str(value).lower() in ['t','true','0']
            else : # STRING_TAGS
                return value
        except ValueError :
            return default

    def save(self, fileName) :
        """Writes a classifier configuration to the specified file."""
        if self.verbose : sys.stderr.write('Writing configuration to %s\n' % fileName)
        outstream = file(fileName, 'w')
        self.write(outstream)

    def setValue(self, name, value, section=CONFIG_SECTION) :
        """Sets a classifier configuration value."""
        tag = name.lower()
        if tag in ALL_TAGS[section] :
            self.config.set(section,name,value)
            self.map[section][name] = value
        else :
            raise ValueError("Unrecognized configuration tag %s" % name)

    def write(self, outstream=sys.stdout) :
        """Writes a classifier configuration to the specified output stream."""
        for c in self.comments :
            outstream.write('# %s\n' % c)
        self.config.write(outstream)

## if __name__ == '__main__' :
##     print "Testing ClassifierConfig:"
##     # Create a configuration from scratch:
##     config = ClassifierConfig()
##     config.setDimer('gt')
##     config.setC(0.01)
##     config.setExonSize(20)
##     config.setIntronSize(20)
##     config.setMink(1)
##     config.setMaxk(12)
##     config.setMaxShift(4)
##     config.setNoShiftStart(8)
##     config.setNoShiftEnd(12)
##     config.setMismatchProfile(False)
##     config.setNormalize(True)
##     config.setROC(0.8765)
##     config.setFastaPath('/a/b/c/gt.fa')
##     config.setSVMPath('/a/b/c/gt.svm')
##     config.addComment('This is a comment added to the file.')
##     config.addComment('This is a second comment added to the file.')
##     config.save('ClassifierConfigTest.cfg')
## 
##     # Load a configuration from a file
##     config = ClassifierConfig(file='ClassifierConfigTest.cfg')
##     print "Values:"
##     for tag in CONFIG_TAGS :
##         value = config.getValue(tag)
##         if value : print "  %s = %s" % (tag, str(value))
##     for tag in PERFORMANCE_TAGS :
##         value = config.getValue(tag, section=PERFORMANCE_SECTION)
##         if value : print "  %s = %s" % (tag, str(value))
