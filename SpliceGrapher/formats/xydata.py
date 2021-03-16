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
Simple interface to X,Y data files.
"""
from SpliceGrapher.shared.utils import *
from sys import maxint as MAXINT

def getXYData(f, **args) :
    """Loads X,Y value pairs from the given file and returns them as two lists."""
    verbose = getAttribute('verbose', False, **args)
    minpos  = getAttribute('minpos', 0, **args)
    maxpos  = getAttribute('maxpos', MAXINT, **args)
    if verbose : sys.stderr.write('Loading X,Y pairs from %s\n' % f)
    indicator = ProgressIndicator(1000000, verbose=verbose)
    X = []
    Y = []
    linectr = 0
    for line in ezopen(f) :
        linectr += 1
        indicator.update()
        parts = line.strip().split(',')
        if len(parts) != 2 : raise Exception('Invalid X,Y pair in %s (line %d): %s' % (f, linectr, line.strip())) 
        x = float(parts[0])
        if x < minpos : continue
        if x > maxpos : break
        X.append(x)
        Y.append(float(parts[1]))
    indicator.finish()
    if verbose : sys.stderr.write('Loaded %d X,Y pairs from %s\n' % (len(X),f))
    return (X,Y)

def writeXYData(f, X, Y, **args) :
    """Writes X,Y value pairs to the given file."""
    if len(X) != len(Y) : raise ValueError('X (%d values) and Y (%d values) are not the same length' % (len(X), len(Y)))
    verbose = getAttribute('verbose', False, **args)
    if verbose : sys.stderr.write('Writing X,Y pairs to %s\n' % f)
    indicator = ProgressIndicator(1000000, verbose=verbose)
    outStream = open(f,'w')
    for i in range(len(X)) :
        indicator.update()
        outStream.write('%f,%f\n' % (X[i], Y[i]))
    indicator.finish()
    outStream.close()
