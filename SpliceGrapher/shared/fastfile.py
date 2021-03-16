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
Permits fast access to large files by using a binary search
on keys found within each record.  Also permits reading part
of a line to avoid reading long records into memory.
"""
from SpliceGrapher.shared.utils import getAttribute
import os

NEWLINE        = '\n'
INT_TYPE       = 'int'
CHAR_TYPE      = 'char'
VALID_KEYTYPES = [INT_TYPE, CHAR_TYPE]

class fastfile(object) :
    """
    Class that provides fast access to large files.  If the file is sorted
    by numeric keys, it can perform a binary search on keys in the file.

    Options:
        delim   - delimiter to use for file records (default=',')
        keypos  - 0-based position where key appears in each record (default=0)
        keytype - type of key: int/char (default=fastfile.INT_TYPE)
    """
    def __init__(self, path, **args) :
        self.stream = open(path, 'r')
        # Find the last position in the file
        self.stream.seek(0, os.SEEK_END)
        self.end    = self.stream.tell()

        # Establish first record location
        self.firstRecPos = 0
        self.setpos(self.firstRecPos)

        self.delim  = getAttribute('delim', ',', **args)
        self.keypos = getAttribute('keypos', 0, **args)
        self.keytype = getAttribute('keytype', INT_TYPE, **args)

        if self.keytype not in VALID_KEYTYPES :
            raise Exception("Invalid key type requested: %s" % self.keytype)

    def close(self) :
        """Convenience wrapper for file stream's close()"""
        self.stream.close()

    def getkey(self, line) :
        """Returns the key for the given line."""
        try :
            result = line.split(self.delim)[self.keypos]
            if self.keytype == INT_TYPE :
                return int(result)
            else :
                return result
        except IndexError :
            raise IndexError("No key in position %d for '%s'" % (self.keypos, line))
        except ValueError :
            raise ValueError("Invalid key in position %d for '%s'" % (self.keypos, line))

    def getline(self, pos) :
        """
        Returns the complete line that ends after the given file position.
        Stored file position is updated to the start of the next line.
        """
        self.stream.seek(pos)
        c = None
        while c != NEWLINE and pos > 0 :
            pos -= 1
            self.stream.seek(pos)
            c = self.stream.read(1)

        result         = self.readline().rstrip()
        self.lastKnown = self.stream.tell()
        return result

    def getlines(self, minkey, maxkey) :
        """Returns all lines with keys between the min and max key values,
        or the end of the file, whichever comes first."""
        result = []
        line   = self.search(minkey)
        key    = self.getkey(line)
        while key <= maxkey :
            if key >= minkey : result.append(line)
            line = self.readline()
            if not line : break
            key  = self.getkey(line)
        return result

    def linelen(self) :
        """Returns the length of the current line."""
        start = self.lastKnown
        self.skipline()
        end   = self.lastKnown
        self.setpos(start)
        return (end-start-1)

    def locateFirstRecord(self) :
        """Based on the key, performs a linear search through the file to find
        the first legitimate record.  Useful for files with headers (e.g., SAM format).
        Note: this is a linear search and thus slow."""
        # Find the first non-header record in the SAM file:
        looking = True
        self.firstRecPos = 0
        while looking :
            self.firstRecPos = self.lastKnown
            line = self.getline(self.firstRecPos)
            try :
                key = self.getkey(line)
                break
            except Exception :
                continue

    def readline(self) :
        """Convenience wrapper for file stream's readline()"""
        return self.stream.readline()

    def search(self, key, lo=None, hi=None) :
        """
        Performs a binary search on a file, looking for the given key.  If the key given
        comes before the first key in the file, it will return the first key in the file.
        Otherwise if the key is not in the file, it returns the highest key below the one given.
        Note: if you wish to grab a set of records within a range, the 'lastKnown' class attribute
        will provide the starting position for the most recent record.
        """
        if lo is None : lo = self.firstRecPos
        if hi is None : hi = self.end-1

        if lo < 0 or hi < 0 :
            raise IndexError("Negative file position specified: %d" % min(lo,hi))
        elif lo > self.end or hi > self.end :
            raise IndexError("File position beyond EOF: %d" % max(lo,hi))

        if lo >= hi :
            return self.getline(hi)
        else :
            midpt  = int(0.5*(lo+hi))
            lokey  = self.getkey(self.getline(lo))
            hikey  = self.getkey(self.getline(hi))
            midkey = self.getkey(self.getline(midpt))

            # Note 1: midpoint file position may fall outside key range 
            # Note 2: calling getline more than once costs little
            #         and keeps the stored file position current
            if midkey < lokey :
                return self.getline(lo)
            elif midkey > hikey :
                return self.getline(hi)
            elif midkey > key :
                return self.search(key, lo, midpt)
            elif midkey < key :
                return self.search(key, midpt, hi-1)
            else :
                return self.getline(midpt)

    def setpos(self, pos) :
        """
        Sets the file pointer to the given position and updates the
        last known position at the same time.
        """
        self.lastKnown = pos
        self.stream.seek(pos)
        
    def skipline(self) :
        """Skips a line in the file."""
        c   = None
        pos = self.stream.tell()
        while c != NEWLINE :
            c = self.stream.read(1)
        self.setpos(self.stream.tell())

    def substring(self, start, end) :
        """Reads part of the current line from the start to the end position."""
        if start > end :
            raise IndexError("Invalid substring specified (%d,%d)" % (start,end))
        elif start < 0 or end < 0 :
            raise IndexError("Negative string position specified (%d,%d)" % (start,end))

        # skip to start position in current line
        self.stream.seek(start, os.SEEK_CUR)
        # read the number of characters requested
        result = self.readline(end-start)
        # reset file pointer back to line start
        self.stream.seek(self.lastKnown)
        return result
