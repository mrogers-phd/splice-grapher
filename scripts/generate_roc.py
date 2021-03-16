#! /usr/bin/env python
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
Script that generates an ROC curve for a splice site dimer classifier.
"""
from SpliceGrapher.shared.config            import *
from SpliceGrapher.shared.utils             import *
from SpliceGrapher.shared.streams           import *
from SpliceGrapher.predict.ClassifierConfig import *
from SpliceGrapher.predict.SiteClassifier   import *

from optparse import OptionParser
import sys, os.path, random, numpy

try :
    from PyML                 import *
    from PyML.containers      import SequenceData
    from PyML.containers      import Labels
    from PyML.classifiers     import SVM
    from PyML.classifiers.svm import loadSVM
except Exception :
    sys.stderr.write('\n** Unable to import PyML modules required for this script.\n')
    sys.exit(1)

#==================================================================================================
# Incompatibilities in some matplotlib versions may yield the following runtime warning:
# ".../matplotlib/backends/backend_gtk.py:639: DeprecationWarning: Use the new widget gtk.Tooltip"
import warnings
warnings.filterwarnings('ignore')
#==================================================================================================

DEFAULT_DIMERS       = ['gt', 'gc', 'ag']
DEFAULT_DIMER_STRING = ','.join(DEFAULT_DIMERS)

def roc(Y, givenY, decisionFunc, n = None, targetClass = 1, normalize = True) :
    """Compute the ROC curve and area under the curve for a two class problem

	:Parameters:
      - `Y` - the predicted labels (can put None instead)
      - `givenY` - the true labels
	  - `decisionFunc` - the values of the decision function
	  - `n` - the number of false positives to take into account (roc_n)
	  - `targetClass` - the "positive" class
      - `normalize` whether to normalize the roc curve (default: True)
        when this is set to False, TP/FP counts are output rather than TP/FP rates
            
    """
    # the number of false positives to take into account
    # note that n can be either an integer or a fraction
    if n is not None and n < 1 :
        n = int(n * numpy.sum(numpy.not_equal(givenY, targetClass)))

    I = range(len(decisionFunc))
    random.shuffle(I)
    decisionFunc = [decisionFunc[i] for i in I]
    givenY = [givenY[i] for i in I]
    f = numpy.array(decisionFunc)
    tp = [0.0]
    fp = [0.0]
    I = numpy.argsort(-f)
    
    rdf = [f[I[0]]]
    for patternIdx in I :
        if givenY[patternIdx] == targetClass :
            tp[-1] += 1
            rdf[-1] = f[patternIdx]
        else :
            tp.append(tp[-1])
            fp.append(fp[-1] + 1.0)
            rdf.append(f[patternIdx])
        if n is not None and fp[-1] >= n :
            break

    numTP = numpy.sum(numpy.equal(givenY, targetClass))
    
    if normalize : 
        for i in xrange(len(tp)):
            if tp[-1] > 0 : tp[i] /= float(numTP)
        for i in xrange(len(fp)) :
            if fp[-1] > 0 : fp[i] /= float(fp[-1])

        area = numpy.sum(tp) / len(tp)

    else :
        area = numpy.sum(tp) / (len(tp) * numTP)

    return tp,fp, area, rdf

def plotROC(res, fileName = None, **args) :
    """plot the ROC curve from a given Results (or Results-like) object

    :Parameters:
      - `res` - Results (or Container object that was made by saving a a
        Results object (note that if you have a Results object you can
        use this function as a method so there is no need to supply this
        argument).
      - `fileName` - optional argument - if given, the roc curve is saved
        in the given file name.  The format is determined by the extension.
        Supported extensions: .eps, .png, .svg
    
    :Keywords:
      - `rocN` - what type of ROC curve to plot (roc50, roc10 etc.) default is
                 full ROC curve
      - `normalize` - whether to normalize the ROC curve (default: True)
      - `plotStr` - which string to pass to matplotlib's plot function
                 default: 'ob'
      - `axis` - redefine the figure axes; takes a list of the form
                 [xmin,xmax,ymin,ymax]
      - `show` - whether to show the ROC curve (default: True)
                 useful when you just want to save the curve to a file.
                 The use of some file formats automatically sets this to False
                 (e.g. svg files).  This relates to quirks of matplotlib.
    """
    rocN         = getAttribute('rocN', None, **args)
    show         = getAttribute('show', True, **args)
    rocNormalize = getAttribute('normalize', True, **args)
    plotStr      = getAttribute('plotStr', 'ob', **args)
    numPoints    = getAttribute('numPoints', 200, **args)
    lineWidth    = getAttribute('lineWidth', 2, **args)
    keyPoints    = getAttribute('keyPoints', [], **args)
    zeroDF       = getAttribute('zeroDF', False, **args)
    targetClass  = 1

    if type(res) == type([]) :
        feature = res[0]
        givenY  = res[1]
        rocTP, rocFP, rocArea, rocDF = roc(None, givenY, feature, rocN, targetClass, rocNormalize)
    else :
        rocTP, rocFP, rocArea, rocDF = roc(res.Y, res.givenY, res.decisionFunc, rocN, targetClass, rocNormalize)    

    stride = int(max(1, float(len(rocTP)) / float(numPoints)))
    if stride > 1 :
        rocTP = [rocTP[i] for i in xrange(0,len(rocTP), stride)]
        rocFP = [rocFP[i] for i in xrange(0,len(rocFP), stride)]        
        rocDF = [rocDF[i] for i in xrange(0,len(rocDF), stride)]        

    import matplotlib
    if fileName is not None :
        if fileName.find('.svg') > 0 : matplotlib.use('SVG')
        if fileName.find('.eps') > 0 : matplotlib.use('PS')

    # Plot ROC curve
    from matplotlib import pylab
    lines = pylab.plot(rocFP, rocTP, plotStr, markersize = 8, linewidth = lineWidth)

    # Add diagonal to plot
    pylab.plot([0,max(rocFP)], [0,max(rocTP)], ':', linewidth = 1)

    if rocNormalize :
        pylab.xlabel('False positive rate', fontsize = 12)
        pylab.ylabel('True positive rate', fontsize = 12)
    else :
        pylab.xlabel('False positives', fontsize = 12)
        pylab.ylabel('True positives', fontsize = 12)

    if rocNormalize :
        pylab.axis([0, 1, 0, 1])

    if 'axis' in args :
        pylab.axis(args['axis'])

    axis = pylab.gca()
    if 'logscale' in args :
        axis.set_xscale('log',nonposx='clip')
        axis.set_xlim(min(rocFP),1.0)

    if 'title' in args :
        axis.set_title(args['title'])

    # Add 'L' markers to show key points:
    # Use first TP value that exceeds each key point
    if keyPoints :
        for val in keyPoints :
            higher = [i for i in xrange(len(rocTP)) if rocTP[i] >= val]
            if higher :
                i     = higher[0]
                lines = pylab.plot([0,rocFP[i],rocFP[i]], [rocTP[i],rocTP[i],0] , '--', linewidth = 1)
                color = lines[0].get_color()
                bbox  = {'facecolor':'white', 'edgecolor':color}
                axis.text(-0.12, rocTP[i], '%3.3f'%rocTP[i], size='x-small', weight='bold', color=color)
                axis.text(rocFP[i]+0.01, rocTP[i]-0.01, 'df=%.2f'%rocDF[i], va='top', size='x-small', weight='bold', bbox=bbox, color=color)

    if zeroDF :
        # Add a marker where DF value is 0
        smallDFidx = [i for i in xrange(len(rocTP)) if rocDF[i] <= 0]
        if smallDFidx :
            i     = smallDFidx[0]
            lines = pylab.plot([0,rocFP[i],rocFP[i]], [rocTP[i],rocTP[i],0] , '-k', linewidth = 1)
            axis.text(rocFP[i]+0.01, 0.5, '~zero: df=%.2f'%rocDF[i], va='center', size='x-small', weight='bold', color='black', rotation='vertical')


    if fileName is not None :
        pylab.savefig(fileName)
    else :
        pylab.show()


def validateFile(path) :
    if not os.path.exists(path) :
        raise Exception("File '%s' not found; exiting." % path)
        sys.exit(1)

USAGE="""%prog config-file [options]

Generates an ROC curve for the SVM for the given configuration file."""

parser    = OptionParser(usage=USAGE)
parser.add_option('-o', dest='output',    default=None,  help='Output file [default= <dimer>_roc.pdf]')
parser.add_option('-s', dest='species',   default=None,  help='Species name for title [default= %default]')
parser.add_option('-A', dest='auc',       default=False, help='Display AUC and accuracy values for curve [default: %default]', action='store_true')
parser.add_option('-D', dest='dfunc',     default=False, help='Display decision function thresholds at ROC=0.9,0.95,0.975 [default: %default]', action='store_true')
parser.add_option('-E', dest='thickness', default=2,     help='Adjust graph line thickness [default: %default]', type='int')
parser.add_option('-Z', dest='zero',      default=False, help='Display ROC score for decision function 0 value [default: %default]', action='store_true')
parser.add_option('-v', dest='verbose',   default=False, help='Verbose mode [default: %default]', action='store_true')

#=======================================================
# Main program
opts, args = parser.parse_args(sys.argv[1:])

if len(args) != 1 :
    parser.print_help()
    sys.exit(1)

validateFile(args[0])
writeStartupMessage()

#-------------------------------------------------------
# Required parameters:
if opts.verbose : sys.stderr.write('ROC script initiating SiteClassifier with %s\n' % args[0])
classifier = SiteClassifier(args[0], verbose=opts.verbose)
dimer      = classifier.config.dimer()
hideStdout()
results    = classifier.runCV()
showStdout()
if opts.auc :
    titleStr   = '%s Site Classifier ROC (AUC=%.2f, ACC=%.2f)' % (dimer.upper(), results[0].roc, results[0].successRate)
else :
    titleStr   = '%s Site Classifier ROC' % dimer.upper()

if opts.species :
    titleStr   = '%s %s' % (opts.species, titleStr)

if not opts.output : opts.output = '%s_roc.pdf' % dimer
keyPoints = [ 0.90, 0.95, 0.975] if opts.dfunc else []

plotROC(results[0], plotStr='-', keyPoints=keyPoints, zeroDF=opts.zero, fileName=opts.output, title=titleStr, lineWidth=opts.thickness)

if opts.verbose : sys.stderr.write('Output written to %s\n' % opts.output)
