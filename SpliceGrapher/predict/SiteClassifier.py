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
Python script for classifying putative splice sites within a FASTA file.
"""
from SpliceGrapher.shared.config            import *
from SpliceGrapher.shared.utils             import *
from SpliceGrapher.predict.ClassifierConfig import *

from PyML.containers import SequenceData
from PyML.classifiers.svm import loadSVM

from zipfile import *

def positionalKmerData(config, dataSet, **args) :
    """
    Code copied and modified from PyML ker.py
    """
    verbose    = getAttribute('verbose', False, **args)
    exonSize   = config.exonSize()
    intronSize = config.intronSize()

    # if no mismatches are allowed the mismatch profile needs to be all 0
    if config.mismatches() == 0 :
        mismatchProfile = [0 for i in range(config.mink(), config.maxk()+1)]
    elif config.acceptor() :
        zeroStart        = intronSize-10
        mismatchProfile  = [1 for i in range(zeroStart)] + [0 for i in range(zeroStart,intronSize)]
        zeroEnd          = 3
        mismatchProfile += [0 for i in range(zeroEnd)] + [1 for i in range(zeroEnd,exonSize)]
    else : # donor
        zeroStart        = exonSize-3
        mismatchProfile  = [1 for i in range(zeroStart)] + [0 for i in range(zeroStart,exonSize)]
        zeroEnd          = 10
        mismatchProfile += [0 for i in range(zeroEnd)] + [1 for i in range(zeroEnd,intronSize)]

    if verbose : logMsg('Mismatch profile for kernel: %s' % mismatchProfile)

    return SequenceData(dataSet,
                        mink=config.mink(),
                        maxk=config.maxk(),
                        mismatches=config.mismatches(),
                        mismatchProfile=mismatchProfile,
                        maxShift=config.maxShift(),
                        noShiftStart=config.noShiftStart(),
                        noShiftEnd=config.noShiftEnd())

def unzipClassifiers(path, **args) :
    """
    Unpacks all classifier files from the given zip file.
    Returns a list of configuration files available in the current
    working directory.  Note: this will over-write any existing
    SVM-related files (xx_xxx.fa, xx_xxx.svm, xx_xxx.cfg) in the directory.
    """
    verbose = getAttribute('verbose', False, **args)
    if not is_zipfile(path) :
        raise ValueError('%s is not a zip file\n' % path)

    # File extensions that must be present for an SVM
    requiredExts = ['.cfg', '.fa', '.svm']
    zipFile      = ZipFile(path)
    names        = zipFile.namelist()
    result       = []
    cfgFiles     = [s for s in names if s.endswith('.cfg')]

    for cfg in cfgFiles :
        # Get the classifier name, e.g., 'ag_acc.cfg' --> 'ag_acc'
        pfx,ignore = os.path.splitext(cfg)
        # Create list of expected files
        required   = set(['%s%s'%(pfx,ext) for ext in requiredExts])
        # Find all files with given prefix
        others     = set([s for s in names if s.startswith(pfx+'.')])
        # See if any are missing
        missing    = required-others
        if missing :
            if verbose : sys.stderr.write("** Warning: missing %s for classifier '%s'\n" % (','.join(missing), pfx))
            continue

        # Extract the files and record the classifier
        for r in required :
            zipFile.extract(r)
        if verbose : sys.stderr.write("Added classifier '%s'\n" % pfx)
        result.append(cfg)

    return result

class SiteClassifier :
    """Wrapper class for a splice site classifier."""
    def __init__(self, filePath, **args) :
        self.verbose   = getAttribute('verbose', False, **args)
        self.config    = ClassifierConfig(file=filePath)
        self.svmPath   = self.config.SVMPath()
        self.fastaPath = self.config.fastaPath()

        if self.config.mismatchProfile() :
            self.data = positionalKmerData(self.config, self.fastaPath)
        else :
            self.data = SequenceData(self.fastaPath,
                                     mink=self.config.mink(), maxk=self.config.maxk(),
                                     maxShift=self.config.maxShift(),
                                     headerHandler=process_labeled_fasta_header)

        if self.config.normalize() :
            self.data.attachKernel('cosine')
        self.svm = loadSVM(self.svmPath, self.data)

    def classify(self, itemList, k) :
        """
        Classifies the given list item and applies a threshold given in the configuration.
        """
        (ssClass,value) = self.svm.classify(itemList,k)
        if value >= self.config.threshold() :
            return (1,value)
        else :
            return (0,value)

    def runCV(self) :
        """
        Classifies the given list item and applies a threshold given in the configuration.
        """
        return self.svm.stratifiedCV(self.data)
