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
Module containing general utility methods.
"""
import configparser, gzip, locale, os, random, subprocess, sys, time
from glob import glob
from io import TextIOWrapper

# Python has no maxint equivalent, so we define our own useful one.
MAXINT = 10**10

def asList(value, delim=','):
    """Generic method for creating a list from a given value.
    The value may be a string of delimiter-separated values, 
    a list or a set."""
    if type(value) == str:
        return value.split(delim)
    elif type(value) == list:
        return value
    elif type(value) == set:
        return list(value)
    else:
        raise ValueError('Expected a string, a list or a set; received %s' % type(value))

def asSet(value, delim=','):
    """Generic method for creating a set from a given value.
    The value may be a string of delimiter-separated values, 
    a list or a set."""
    if type(value) == str:
        return set(value.split(delim))
    elif type(value) == list:
        return set(value)
    elif type(value) == set:
        return value
    else:
        raise ValueError('Expected a string, a list or a set; received %s' % type(value))

def assureDirectory(d, verbose=False):
    """Guarantee that the given directory exists, or raise an
    exception if there is a problem with it."""
    if not os.path.isdir(d):
        if verbose:
            sys.stderr.write('Creating directory %s\n' % d)
        os.makedirs(d)
    return d

def bsearch(X, target, getValue=lambda a:a, **args):
    """
    Generic binary search returns the index of the value in X that
    is closest to the target.  'getValue' provides access to the correct
    attribute of each element in X.  Assumes X is sorted such that
    getValue(a) < getValue(b) iff a < b.
    """
    if not X:
        raise ValueError('Cannot perform binary search on an empty list')
    lo  = getAttribute('low', 0, **args)
    hi  = getAttribute('high', len(X)-1, **args)
    mid = lo
    while lo < hi:
        mid    = (lo+hi)/2
        midval = getValue(X[mid])
        if midval < target:
            lo = mid+1
        elif midval > target:
            hi = mid-1
        else:
            return mid
    return lo

def commaFormat(d):
    """Formats integer values using commas.  For example, 123456789 becomes '123,456,789'"""
    # Establish user's default their OS:
    locale.setlocale(locale.LC_ALL, '')
    return locale.format("%d", d, grouping=True)

def configMap(cfgFile):
    """Reads a configuration file and returns a dictionary of all sections, options and values."""
    result = {}
    config = configparser.ConfigParser()
    # Activates case-sensitivity:
    config.optionxform = str

    try:
        config.read(cfgFile)
    except configparser.ParsingError:
        raise ValueError('Invalid configuration file %s' % cfgFile)

    for sect in config.sections():
        result[sect] = {}
        options = config.options(sect)
        for opt in options:
            try:
                result[sect][opt] = config.get(sect,opt)
            except:
                result[option] = None
    return result

def dictString(valDict, delim=','):
    """Returns a simple string representation of a list of values."""
    return delim.join(['%s -> %s' % (x,valDict[x]) for x in valDict])

def ezopen(fileName):
    """Allows clients to open files without regard for whether they're gzipped."""
    if not (os.path.exists(fileName) and os.path.isfile(fileName)):
        raise ValueError('file does not exist at %s' % fileName)
    
    fileHandle  = gzip.GzipFile(fileName)
    gzippedFile = True
    try:
        line = fileHandle.readline()
        fileHandle.close()
        return TextIOWrapper(gzip.GzipFile(fileName))
    except:
        return open(fileName)

def fileLen(path):
    """Simple function to get an exact file length."""
    with open(path) as f:
        for i, l in enumerate(f): pass
    return i+1

def filePrefix(f):
    """Returns the filename prefix for a file.  For example:
       /my/dir/myfile.ext --> myfile"""
    head,tail     = os.path.split(f)
    prefix,suffix = os.path.splitext(tail)
    return prefix

def findFile(name, path, delim=':'):
    """Finds the first instance of a file name in the given path string."""
    paths = path.split(delim)
    for p in paths:
        filePath = os.path.join(p, name)
        if os.path.exists(filePath) and os.path.isfile(filePath):
            return filePath

def getAttribute(key, default, **args):
    """Returns the value for the given key in the arguments dict, if found;
    otherwise returns default value."""
    return default if key not in args else args[key]

def getEnvironmentValue(name, default=None):
    """Returns the value for the given environment variable name, if found;
    otherwise returns default value."""
    try:
        return os.environ[name]
    except KeyError:
        return default

def idFactory(pfx='', initial=1):
    """Generates unique ids using the given prefix.  For example,
    idFactory('ev_') will generate 'ev_1', 'ev_2', ..."""
    prefix  = pfx if pfx else ''
    counter = initial
    while 1:
        yield '%s%d' % (prefix, counter)
        counter += 1

def listString(vals, delim=','):
    """Returns a simple string representation of a list of values."""
    return delim.join([str(x) for x in vals])

def logMessage(s, logstream=None):
    """Allows log messages to be output both to stderr and a log file."""
    sys.stderr.write(s)
    if logstream:
        logstream.write(s)

def makeGraphListFile(spliceGraphDir):
    """Given a path to a top-level directory of splice graphs, returns
    a file that contains the paths to all graphs under that directory.
    Follows the standard SpliceGrapher directory structure:
        top-level-dir/chromosome-dir/splice-graph-file"""
    subdirs   = os.path.join(spliceGraphDir, '*')
    target    = os.path.join(subdirs, '*.gff')
    graphList = glob(target)
    if not graphList:
        raise ValueError('No splice graphs found in %s\n' % spliceGraphDir)

    graphList.sort()
    result      = '%s.lis' % spliceGraphDir
    graphStream = open(result,'w')
    graphStream.write('\n'.join(graphList))
    graphStream.close()
    return result

def process_fastq_header(header):
    """Find the chromosome identifier within a FASTQ header."""
    return header.split()[0]

def process_fasta_header(header):
    """Find the chromosome identifier within a FASTA header."""
    return header.split()[0]

def process_labeled_fasta_header(header):
    """ 
    Extract a sequence ID and its label from a fasta file.
    Used with training data FASTA files, so it assumes
    headers have the form "ID label=#".
    """
    id,labelToken = header.split()
    return id.strip(),labelToken.split('=')[1]

def runCommand(s, **args):
    """Announces a command runs it.."""
    logstream   = getAttribute('logstream', None, **args)
    debug       = getAttribute('debug', False, **args)
    exitOnError = getAttribute('exitOnError', True, **args)
    stderr      = getAttribute('stderr', None, **args)
    stdout      = getAttribute('stdout', None, **args)
    message     = '    ' + timeString('%s\n' % s)
    sys.stderr.write(message)
    if logstream:
        logstream.write(message)

    retcode = 0
    if not debug:
        retcode = subprocess.call(s, shell=True, stderr=stderr, stdout=stdout)

    if exitOnError and retcode < 0:
        raise Exception('Error running command: returned %d signal\n%s' % (retcode, s))

def stringToRange(s):
    """Converts a string of the form '1,3,7-10' to a list [1,3,7,8,9,10]."""
    ranges = s.split(',')
    result = set()
    for r in ranges:
        if r.isdigit():
            result.add(int(r))
        else:
            parts = [int(x) for x in r.split('-')]
            assert(len(parts) == 2)
            result.update(list(range(parts[0], parts[1]+1)))
    return sorted(result)

def substringAfter(s, tag):
    """Returns the substring after the given tag, if found; None otherwise."""
    pos = s.find(tag)
    if pos >= 0:
        return s[pos+len(tag):]

def substringBefore(s, tag):
    """Returns the substring before the given tag, if found; None otherwise."""
    pos = s.find(tag)
    if pos >= 0:
        return s[:pos]

def substringBetween(s, tag1, tag2):
    """Returns the substring between the two tags, if found; None otherwise."""
    p1 = s.find(tag1)
    if p1 < 0:
        return None
    p1 += len(tag1)
    p2  = s.find(tag2, p1)
    if p2 > 0:
        return s[p1:p2]

def timestamp(format_string='%Y%m%d%H%M%S'):
    """Returns a timestamp unique to the current second."""
    return time.strftime(format_string, time.localtime())

def timeString(s, format_string='%X', LF=False):
    """Returns the input string with user-readable a timestamp prefix."""
    timestamp = time.strftime(format_string, time.localtime())
    result    =  '%s %s' % (timestamp, s)
    if LF:
        result += '\n'
    return result

def to_numeric(s):
    """Attempts to return the given value as a float or an int, else returns the original string."""
    try:
        return int(s)
    except ValueError:
        # (fails on numbers with decimals like '12.34')
        pass

    try:
        return float(s)
    except ValueError:
        pass

    return s

def validateDir(path):
    """Standard method for validating directory paths."""
    validateFile(path)
    if not os.path.isdir(path):
        raise Exception("'%s' is not a directory; exiting." % path)

def validateFile(path):
    """Standard method for validating file paths."""
    if not path:
        raise Exception("'%s' is not a valid file path; exiting." % path)

    if not os.path.exists(path):
        raise Exception("File '%s' not found; exiting." % path)

def writeStartupMessage():
    """Standardized startup message for all scripts."""
    base = os.path.basename(sys.argv[0])
    sys.stderr.write(timeString('%s Started\n' % base))

class ProgressIndicator(object):
    """A simple progress indicator."""
    def __init__(self, increment, description='', verbose=True):
        self.limit   = increment
        self.barlim  = int(self.limit/10)
        self.dotlim  = int(self.barlim/5)
        self.descr   = description
        self.started = False
        self.verbose = verbose
        self.ctr     = 0

    def count(self):
        """Returns the current count."""
        return self.ctr

    def finish(self):
        """Finishes the progress output by appending a newline,
        if anything has been written."""
        if self.started:
            sys.stderr.write('\n')
        self.started = False

    def reset(self):
        """Resets the indicator to be used again."""
        self.ctr = 0
        self.finish()

    def update(self):
        """Updates the indicator."""
        self.ctr += 1
        if not self.verbose:
            return
        if self.ctr % self.limit == 0:
            sys.stderr.write('%s %s\n' % (commaFormat(self.ctr), self.descr))
        elif self.ctr % self.barlim == 0:
            sys.stderr.write('|')
        elif self.ctr % self.dotlim == 0:
            self.started = True
            sys.stderr.write('.')

class RandomListIterator(object):
    """
    Returns items at random, with replacement, from a list of values.
    For lists with more than ~1000 elements, this is about 30% more
    efficient than using random.sample() repeatedly.
    """
    def __init__(self, values, **args):
        self.values = values
        self.rand   = random.Random()
        self.limit  = len(self.values)-1
        if 'seed' in args:
            self.rand.seed = args['seed']

    def __iter__(self):
        return self

    def __next__(self):
        """Iterator implementation that returns a random value, with
        replacement, from a list."""
        i = self.rand.randint(0,self.limit)
        return self.values[i]
