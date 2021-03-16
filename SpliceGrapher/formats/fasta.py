# Copyright (C) 2003, 2004 by BiRC -- Bioinformatics Research Center
#                                     University of Aarhus, Denmark
#                                     Contact: Thomas Mailund <mailund@birc.dk>
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
A parser for FASTA files.

Copyright (C) 2003, 2004 by BiRC -- Bioinformatics Research Center
                                    University of Aarhus, Denmark
                                    Contact: Thomas Mailund <mailund@birc.dk>
with changes by Asa Ben-Hur and Mark Rogers
"""
from SpliceGrapher.shared.utils import ezopen
import os

class MalformedInput :
    "Exception raised when the input file does not look like a fasta file."
    def __init__(self, value) :
        self.param = value
    def __repr__(self) :
        return repr(self.param)

class FastaRecord :
    "a fasta record."
    def __eq__(self, other) :
        return self.header == other.header and self.sequence == other.sequence

    def __init__(self, header, sequence):
        "Create a record with the given header and sequence."
        self.header   = header
        self.sequence = sequence

    def __str__(self) :
        return '>' + self.header + '\n' + self.sequence + '\n'

def _fasta_itr_from_file(file) :
    "Provide an iteration through the fasta records in file."
    h = file.readline().strip()
    if h[0] != '>':
        raise MalformedInput("First header in FASTA file must start with '>' (found %s)" % h)
    h = h[1:]

    seq     = []
    linectr = 0
    for line in file:
        linectr += 1
        line     = line.strip()
        if len(line) < 1 :
            raise MalformedInput('Blank line in FASTA file (line %d)' % linectr)

        if line[0] == '>':
            yield FastaRecord(h,''.join(seq))
            h   = line[1:]
            seq = []
            continue
        seq.append(line)
    yield FastaRecord(h,''.join(seq))

def _fasta_itr_from_name(fname):
    "Provide an iteration through the fasta records in the file named fname. "
    f = ezopen(fname)
    for rec in _fasta_itr_from_file(f) :
        yield rec

def _fasta_itr(src):
    """Provide an iteration through the fasta records in file `src'.
    
    Here `src' can be either a file object or the name of a file.
    """
    if type(src) == str :
        return _fasta_itr_from_name(src)
    elif type(src) == file :
        return _fasta_itr_from_file(src)
    else:
        raise TypeError

def fasta_get_by_name(itr,name, byLength=False):
    "Return the record in itr with the given name."
    x = name.strip()
    n = len(x)
    for rec in itr:
        s = rec.header.strip()
        if byLength and n < len(s) :
            s = s[:n]
        if s == x:
            return rec
    return None

class fasta_itr (object) :
    "An iterator through a sequence of fasta records."
    def __init__(self,src) :
        "Create an iterator through the records in src."
        self.__itr = _fasta_itr(src)

    def __iter__(self) :
        return self

    def __next__(self) :
        return next(self.__itr)

    def __getitem__(self,name) :
        return fasta_get_by_name(iter(self),name)

class fasta_slice (object) :
    """Provide an iteration through the fasta records in 'src', from
    'start' to 'stop'.
    """
    def __init__(self, src, first, last = None):
        """
        :Parameters:
        - `src` - the fasta file/file handle. file can be gzipped.
        - `first` - the first record (either its index in the file or
          its identifier
        - `last` - the last record to be output (index in the file or identifier)
        """
        self.__itr = _fasta_itr(src)
        self.__first = first
        self.__last = last
        if type(first) == int :
            self.__current = 0
        elif type(first) == type('') :
            self.__current = None
        else :
            raise ValueError('bad first')

        self.__foundFirst = False
        if self.__first == 0 or self.__first == '' :
            self.__foundFirst = True

    def __iter__(self) :
        return self

    def __next__(self) :
        """
        Implementation of the iterator interface.
        """
        if not self.__foundFirst :
            for rec in self.__itr :
                if type(self.__first) == int :
                    if self.__first == self.__current :
                        self.__foundFirst = True
                        break
                    self.__current += 1
                else :
                    if rec.header == self.__first :
                        self.__foundFirst = True
                        break
                    self.__current = rec.header
            if not self.__foundFirst :
                raise ValueError('did not find first record')
            return rec
        rec = next(self.__itr)

        if self.__last is not None :
            if type(self.__first) == int :
                self.__current += 1
                if self.__current == self.__last :
                    raise StopIteration
            else :
                if rec.header == self.__last :
                    raise StopIteration
                self.__current = rec.header
        return rec

    def __getitem__(self, name):
        return fasta_get_by_name(iter(self),name)

    def save(self, fileName) :
        outfile = open(fileName, 'w')
        for record in self :
            outfile.write(str(record))

def get_sequence(src, name):
    "Return the record in src with the given name."
    return fasta_itr(src)[name]

def fasta_count(src) :
    """
    count the number of records in a fasta file
    """
    return sum([1 for rec in fasta_itr(src)])

def fasta_split(fileName, num_files, directory = None) :
    """
    split a fasta file into a given number of files
    the resulting files are named by adding a number
    to the provided file name.

    :Parameters:
    - `fileName` - the fasta file to split
    - `num_files` - the number of files to split into
    - `directory` - the directory into which to write the files
    """
    num_records = fasta_count(fileName)

    if directory is None :
        base, ext = os.path.splitext(fileName)
    else :
        dir, name = os.path.split(fileName)
        base, ext = os.path.splitext(name)
        base = os.path.join(directory, base)

    rec_num       = 1
    file_num      = 0
    recs_per_file = float(num_records) / float(num_files)
    rec_limit     = 0

    for rec in fasta_itr(fileName) :
        if rec_num > round(rec_limit) :
            file_num  += 1
            outfile    = open(base + '.' + str(file_num) + ext, 'w')
            rec_limit += recs_per_file
        outfile.write(str(rec))
        rec_num += 1

    ## assert(num_files == file_num)

class FastaRandomizer(object) :
    """
    Class that loads all FASTA sequences in a file and then returns
    random subsets of the sequences.
    """
    def __init__(self, fastaFile) :
        self.records  = [rec for rec in fasta_itr(fastaFile)]
        self.recRange = list(range(len(self.records)))

    def randomRecords(self, n) :
        """
        Returns a list of records from the file in random order.
        Raises a ValueError exception if the number requested is
        greater than the number stored.
        """
        if n > len(self.records) :
            raise ValueError('%d FASTA records requested; only %d stored.' % (n, len(self.records)))

        import random
        idList = random.sample(self.recRange, n)
        return [self.records[i] for i in idList]

    def __str__(self) :
        return 'FastaRandomizer - %d records' % len(self.records)

def truncateSequences(fastaFile, exonSize, intronSize, acceptor=False, outFile=None, verbose=False) :
    """
    Method that truncates FASTA sequences based on given intron and exon sizes.
    """
    if verbose : sys.stderr.write('Loading sequences from %s\n' % fastaFile)
    fiter      = fasta_itr(fastaFile)

    outStream = sys.stdout
    if outFile :
        if verbose : sys.stderr.write('Writing output to %s\n' % outFile)
        outStream = file(outFile, 'w')

    for rec in fiter :
        midpt = (len(rec.sequence)/2)
        if acceptor :
            newRec  = FastaRecord(rec.header, rec.sequence[midpt-intronSize:midpt+exonSize])
        else :
            newRec  = FastaRecord(rec.header, rec.sequence[midpt-exonSize:midpt+intronSize])
        outStream.write(str(newRec))
