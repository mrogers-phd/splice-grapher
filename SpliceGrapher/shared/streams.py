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
Module that provides control over the sys.stderr and sys.out streams.
This is particularly useful for Python-wrapped C code that may not
provide control over these streams.
"""
import sys,os

SAVE_STDOUT = os.dup(sys.stdout.fileno())
SAVE_STDERR = os.dup(sys.stderr.fileno())
NULL_STREAM = open(os.devnull,'a')

def hideAll(verbose=False):
    """Hides output to both stdout and stderr."""
    hideStderr(verbose=verbose)
    hideStdout(verbose=verbose)

def hideStderr(verbose=False):
    """Hides output to sys.stderr until the next call to showStderr()."""
    if verbose:
        sys.stderr.write('Hiding stderr\n')
    os.dup2(NULL_STREAM.fileno(), 2)

def hideStdout(verbose=False):
    """Hides output to sys.stdout until the next call to showStdout()."""
    if verbose:
        sys.stdout.write('Hiding stdout\n')
    os.dup2(NULL_STREAM.fileno(), 1)

def showAll(verbose=False):
    """Reinstates output to both stdout and stderr."""
    showStderr(verbose=verbose)
    showStdout(verbose=verbose)

def showStderr(verbose=False):
    """Reinstates output to sys.stderr."""
    sys.stderr.flush()
    os.dup2(SAVE_STDERR, 2)
    if verbose:
        sys.stderr.write('Reinstated stderr\n')

def showStdout(verbose=False):
    """Reinstates output to sys.stdout."""
    sys.stdout.flush()
    os.dup2(SAVE_STDOUT, 1)
    if verbose:
        sys.stdout.write('Reinstated stdout\n')
