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
from SpliceGrapher.shared.utils import ezopen, bsearch, getAttribute
from SpliceGrapher.predict.PredictedSites import ACCEPTOR, DONOR
from PyML import *
from PyML.classifiers.svm import loadSVM
from sys import maxint as MAXINT
import sys,os

#==================
# Global constants
#==================
# Gap classifier features:
CLASS         = 'class'
GENE_FEATURES = [CLUST_MEAN, CLUST_VAR, GENE_MEAN, INTRON_COUNT, INTRON_MEAN] = \
        ['clustmean', 'clustvar', 'genemean', 'introncount', 'intronmean']
# US = upstream / DS = downstream
US_FEATURES   = [US_GAP, US_MEAN, US_VAR, US_ACC, US_DON ]         = ['usize', 'umean', 'uvar', 'uacc', 'udon']
GAP_FEATURES  = [GAP_SIZE, GAP_ACC, GAP_DON, INTRON_SPAN, NONUNIQ] = ['gap', 'gacc', 'gdon', 'intronspan', 'nonuniq']
DS_FEATURES   = [DS_GAP, DS_MEAN, DS_VAR, DS_ACC, DS_DON ]         = ['dsize', 'dmean', 'dvar', 'dacc', 'ddon']
FEATURES      = GENE_FEATURES + US_FEATURES + GAP_FEATURES + DS_FEATURES
FIELDS        = [CLASS] + FEATURES
FIELD_STRING  = ','.join(FIELDS)
FIELD_FMT     = ','.join(['%f' for x in FEATURES])
# Indexes to all features (add 1 to each since label is at 0)
ALL_FEATURES  = range(1,len(FEATURES)+1)

# Feature type: discrete or continuous (directs representation)
DISCRETE_TAG   = 'discrete'
CONTINUOUS_TAG = 'continuous'
FEATURE_TYPE   = {CLUST_MEAN   : CONTINUOUS_TAG,
                  CLUST_VAR    : CONTINUOUS_TAG,
                  GENE_MEAN    : CONTINUOUS_TAG,
                  INTRON_COUNT : DISCRETE_TAG,
                  INTRON_MEAN  : CONTINUOUS_TAG,
                  US_GAP       : DISCRETE_TAG,
                  US_MEAN      : CONTINUOUS_TAG,
                  US_VAR       : CONTINUOUS_TAG,
                  US_ACC       : DISCRETE_TAG,
                  US_DON       : DISCRETE_TAG,
                  GAP_SIZE     : DISCRETE_TAG,
                  GAP_ACC      : DISCRETE_TAG,
                  GAP_DON      : DISCRETE_TAG,
                  INTRON_SPAN  : CONTINUOUS_TAG,
                  NONUNIQ      : DISCRETE_TAG,
                  DS_GAP       : DISCRETE_TAG,
                  DS_MEAN      : CONTINUOUS_TAG,
                  DS_VAR       : CONTINUOUS_TAG,
                  DS_ACC       : DISCRETE_TAG,
                  DS_DON       : DISCRETE_TAG,
                  }

DEFAULT_DIST = 10000

# Gap classifier configuration:
TAG_DELIM    = '='
KNOWN_TAGS   = [SVM_TAG, C_TAG, CSV_TAG, NAME_TAG, GAMMA_TAG, DEGREE_TAG, COSINE_TAG, SPARSE_TAG, ROC_TAG, ACC_TAG, FEAT_TAG] \
             = ['svm', 'c', 'csv', 'name', 'gamma', 'degree', 'cosine', 'sparse', 'roc', 'acc', 'features']
TAG_DEFAULTS = {SVM_TAG    : None,
                C_TAG      : 0.0,
                CSV_TAG    : None,
                NAME_TAG   : '',
                GAMMA_TAG  : -1,
                DEGREE_TAG : 1,
                SPARSE_TAG : False,
                COSINE_TAG : False,
                ROC_TAG    : 0.0,
                ACC_TAG    : 0.0,
                FEAT_TAG   : None
                }

#==================
# Classes
#==================
class Gap(object) :
    """Simple class that encapsulates a gap between two Clusters."""
    def __init__(self, cluster1, cluster2) :
        self.minpos = cluster1.maxpos+1 if cluster1 else 0
        self.maxpos = cluster2.minpos-1 if cluster2 else MAXINT
        self.prior  = cluster1
        self.after  = cluster2
        if self.minpos > self.maxpos : raise ValueError('Received minpos %d > maxpos %d' % (minpos, maxpos))

    def contains(self, o) : return self.minpos < o.minpos and o.maxpos < self.maxpos 
    def overlaps(self, o) : return self.minpos < o.maxpos and o.minpos < self.maxpos 
    def __len__(self)     : return self.maxpos-self.minpos+1 
    def __str__(self)     : return 'gap %d-%d (%d)' % (self.minpos, self.maxpos, len(self))

class GapFeatures(object) :
    """Encapsulates all features for one example of a gap."""
    def __init__(self) :
        self.values = {}.fromkeys(FEATURES,0.0)
        self.label  = None

    def featureString(self, features=FEATURES) :
        fmtString   = ','.join(['%f' for x in features])
        fieldString = fmtString % tuple([self.value(k) for k in features])
        return '%d,%s' % (self.label, fieldString)

    def featureList(self, features=FEATURES) :
        return [self.values[x] for x in features]

    def fieldString(self) :
        return FIELD_FMT % tuple([self.value(k) for k in FEATURES])

    def set(self, k, v) : self.values[k] = float(v)
    def __str__(self)   : return self.featureString()
    def value(self, k)  : return self.values[k]

class GapSVMConfig(object) :
    """
    Encapsulates all the information necessary to completely describe
    a gap-merge classifier.
    """
    def __init__(self, path=None) :
        self.__dict__.update(TAG_DEFAULTS)
        if path : self.load(path)

    def load(self, path) :
        """Loads a gap SVM configuration from a file."""
        if type(path) != type('') :
            raise ValueError('Path must be a string')
        ctr = 0
        for line in ezopen(path) :
            ctr += 1
            if line.startswith('#') : continue # Comment line
            s     = line.strip()
            parts = s.replace(' ','').split(TAG_DELIM)
            key   = parts[0].lower()
            if key not in KNOWN_TAGS :
                sys.stderr.write("Warning: unrecognized key '%s' in configuration file, line %d.\n" % (key,ctr))
            if key != NAME_TAG and len(parts) != 2 :
                raise ValueError('Too many keys/values on same line: %s' % s)

            value = parts[-1]
            if key == GAMMA_TAG :
                value = float(parts[-1])
            elif key == DEGREE_TAG :
                value = int(parts[-1])
            elif key == SPARSE_TAG :
                value = (parts[-1][0].lower() == 't')
            elif key == FEAT_TAG :
                value = [int(x) for x in parts[-1].split(',')]
            elif key == NAME_TAG :
                value = TAG_DELIM.join(s.split(TAG_DELIM)[1:]).strip()
            self.__dict__[key] = value
        self.validate()

    def isGaussian(self)   : return self.gamma > 0
    def isLinear(self)     : return self.gamma < 0 and self.degree < 2
    def isPolynomial(self) : return self.degree > 1

    def save(self, path) :
        outStream = open(path,'w')
        for k in KNOWN_TAGS :
            outStream.write('%-8s = %s\n' % (k, str(self.__dict__[k])))
        outStream.close()

    def __str__(self) :
        return str(self.__dict__)

    def validate(self) :
        # Check for invalid parameter combinations:
        if self.gamma > 0 and self.degree > 1 :
            raise ValueError('Both gamma (%f) and polynomial degree (%d) specified; only one may be used.' % (self.gamma, self.degree))

#==================
# Functions
#==================
def clusterDepth(c) :
    return c.avgDepth() if c else 0.0

def clusterLength(c) :
    return len(c) if c else 0.0

def clusterVariance(c) :
    if not c : return MAXINT
    mean    = c.avgDepth()
    N       = 0
    depthSS = 0.0
    for k in c.depths.keys() :
        delta    = c.depths[k] - mean
        depthSS += delta*delta
        N       += 1
    return depthSS/N if N > 0 else MAXINT

def clusterListDepth(clusters) :
    if not clusters : return 0.0
    allDepths = {}
    for c in clusters :
        allDepths.update(c.depths)
    return float(sum(allDepths.values()))/len(allDepths.keys())

def clusterListVariance(clusters) :
    if not clusters : return 0.0
    mean      = clusterListDepth(clusters)
    allDepths = {}
    for c in clusters :
        allDepths.update(c.depths)
    SS = 0.0
    for k in allDepths.keys() :
        delta = allDepths[k]-mean
        SS   += delta*delta
    return float(SS)/len(allDepths)

def gapStatistics(gap, sitePredictor, nonunique, chrom, strand) :
    """Returns a list of values associate with a gap.  This function
    provides a unique and common method for obtaining gap information."""
    rec   = gapStatisticsToFeatures(gap, sitePredictor, nonunique, chrom, strand)
    order = [GAP_SIZE, US_GAP, US_MEAN, US_VAR, DS_GAP, DS_MEAN, DS_VAR, US_ACC, US_DON, DS_ACC, DS_DON, GAP_ACC, GAP_DON, INTRON_SPAN, NONUNIQ]
    return rec.featureList(features=order)

def gapStatisticsString(gap, sitePredictor, nonunique, chrom, strand) :
    """Returns a comma-delmited string of values associate with a gap."""
    values       = gapStatistics(gap, sitePredictor, nonunique, chrom, strand)
    formats      = ['%f' for i in range(len(values))]
    formatString = ','.join(formats)
    return formatString % tuple(values)

def gapStatisticsToFeatures(gap, sitePredictor, nonunique, chrom, strand, **args) :
    """Computes statistics relevant to the given gap and stores the results
    in a GapFeatures object."""

    #DEBUG = getAttribute('debug', False, **args)
    DEBUG = False

    # Splice-site and non-uniquely-alignable statistics:
    gdCount = gaCount = nuCount = 0
    udDist = uaDist = ddDist = daDist = MAXINT

    # Upstream gap donor/acceptor counts
    if sitePredictor :
        # Find minimum distance between gap's 5' end and nearest sites
        # and between 3' end and nearest sites
        donors    = sitePredictor.donors[chrom][strand]
        acceptors = sitePredictor.acceptors[chrom][strand]
        dmin      = bsearch(donors, gap.minpos-1)
        amin      = bsearch(acceptors, gap.minpos-1)
        dmax      = bsearch(donors, gap.maxpos+1, lo=dmin)
        amax      = bsearch(acceptors, gap.maxpos+1, lo=amin)

        # Adjust min/max position to be right at or below gap
        while dmin > 0 and donors[dmin] > gap.minpos : dmin -= 1
        while dmax < len(donors) and donors[dmax] < gap.maxpos : dmax += 1
        while amin > 0 and acceptors[amin] > gap.minpos : amin -= 1
        while amax < len(acceptors) and acceptors[amax] < gap.maxpos : amax += 1
        dmin = max(0,dmin)
        amin = max(0,amin)
        amax = min(len(acceptors)-1,amax)
        dmax = min(len(acceptors)-1,dmax)

        if strand == '+' :
            udDist = gap.minpos - donors[dmin]
            uaDist = gap.minpos - acceptors[amin]
            ddDist = donors[dmax] - gap.maxpos
            daDist = acceptors[amax] - gap.maxpos
            if DEBUG :
                sys.stderr.write('gapStatisticsToFeatures: DEBUG %s (%s strand):\n' % (gap, strand))
                sys.stderr.write('    donor    %d upstream of %d\n' % (donors[dmin], gap.minpos))
                sys.stderr.write('    acceptor %d upstream of %d\n' % (acceptors[amin], gap.minpos))
                sys.stderr.write('    donor    %d downstream of %d\n' % (donors[dmax], gap.maxpos))
                sys.stderr.write('    acceptor %d downstream of %d\n' % (acceptors[amax], gap.maxpos))
                sys.stderr.write('    upstream donor: %d / downstream acceptor: %d\n' % (donors[dmin], acceptors[amax]))
        else :
            ddDist = gap.minpos - donors[dmin]
            daDist = gap.minpos - acceptors[amin]
            udDist = donors[dmax] - gap.maxpos
            uaDist = acceptors[amax] - gap.maxpos
            if DEBUG :
                sys.stderr.write('gapStatisticsToFeatures: DEBUG %s (%s strand):\n' % (gap, strand))
                sys.stderr.write('    donor    %d upstream of %d\n' % (donors[dmax], gap.minpos))
                sys.stderr.write('    acceptor %d upstream of %d\n' % (acceptors[amax], gap.minpos))
                sys.stderr.write('    donor    %d downstream of %d\n' % (donors[dmin], gap.maxpos))
                sys.stderr.write('    acceptor %d downstream of %d\n' % (acceptors[amin], gap.maxpos))
                sys.stderr.write('    upstream donor: %d / downstream acceptor: %d\n' % (donors[dmax], acceptors[amin]))

        # Within-gap donor/acceptor counts
        while dmin < len(donors) and donors[dmin] < gap.maxpos :
            if dmin > gap.minpos : gdCount += 1
            dmin += 1

        while amin < len(acceptors) and acceptors[amin] < gap.maxpos :
            if amin > gap.minpos : gaCount += 1
            amin += 1

    if nonunique :
        for pos in range(gap.minpos, gap.maxpos) :
            if nonunique[pos] : nuCount += 1

    # Determine the minimum possible intron size that might contain a gap:
    # the gap size plus the distances to the nearest upstream donor and
    # downstream acceptor
    minIntronSize = len(gap) + udDist + daDist

    # Intron span proportion is the gap size divided by the smallest intron container
    intronSpan    = float(len(gap))/minIntronSize

    result = GapFeatures()
    result.set(GAP_SIZE, len(gap))
    result.set(US_GAP,  clusterLength(gap.prior))
    result.set(US_MEAN, clusterDepth(gap.prior))
    result.set(US_VAR,  clusterVariance(gap.prior))
    result.set(DS_GAP,  clusterLength(gap.after))
    result.set(DS_MEAN, clusterDepth(gap.after))
    result.set(DS_VAR,  clusterVariance(gap.after))
    result.set(US_ACC, uaDist)
    result.set(US_DON, udDist)
    result.set(DS_ACC, daDist)
    result.set(DS_DON, ddDist)
    result.set(GAP_ACC, gaCount)
    result.set(GAP_DON, gdCount)
    result.set(INTRON_SPAN, intronSpan)
    result.set(NONUNIQ, nuCount)
    return result

def loadGapData(config, std=None, dataFile=None) :
    """
    Returns a data set using the configuration file to initialize it.
    """
    if not dataFile :
        dataFile = config.csv

    data = VectorDataSet(dataFile, labelsColumn=0)

    if std :
        std.test(data)
    else :
        std  = preproc.Standardizer()
        std.skip_sparse = (not config.sparse)
        std.train(data)

    if config.cosine :
        data.attachKernel('cosine')

    if config.isPolynomial() :
        data.attachKernel('poly', degree=int(config.degree))
    elif config.isGaussian() :
        data.attachKernel('gaussian', gamma=float(config.gamma))
    return data, std

def loadGapSVM(config) :
    """
    Returns a classifier based on the given configuration.
    """
    data, std = loadGapData(config)
    svm       = loadSVM(config.svm, data)
    return svm, data, std

def makeGapVector(values, config, std, filterFeatures=False) :
    """
    Converts the given list of values into a data set for classifying.
    """
    data   = [values[k] for k in config.features] if filterFeatures else list(values)
    result = VectorDataSet([data])
    std.test(result)

    if config.cosine :
        result.attachKernel('cosine')

    if config.isPolynomial() :
        result.attachKernel('poly', degree=int(config.degree))
    elif config.isGaussian() :
        result.attachKernel('gaussian', gamma=float(config.gamma))
    return result

def makeFeatureVector(gapFeatures, config, std) :
    """
    Converts the given GapFeatures values into a data set for classifying.
    """
    values = [gapFeatures.value(k) for k in config.features]
    result = VectorDataSet([values])
    std.test(result)

    if config.cosine :
        result.attachKernel('cosine')

    if config.isPolynomial() :
        result.attachKernel('poly', degree=int(config.degree))
    elif config.isGaussian() :
        result.attachKernel('gaussian', gamma=float(config.gamma))
    return result
