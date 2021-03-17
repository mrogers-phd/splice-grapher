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
Module for manipulating SAM-formatted files.
"""
from SpliceGrapher.shared.utils import ProgressIndicator, ezopen, getAttribute
from SpliceGrapher.shared.ShortRead import *
from sys import maxsize as MAXINT
from io import TextIOWrapper
import pysam
import re
import sys

# Header tags
HEADER_HD_TAG = 'HD'
HEADER_LN_TAG = 'LN'
HEADER_SN_TAG = 'SN'
HEADER_SO_TAG = 'SO'
HEADER_SQ_TAG = 'SQ'
HEADER_VN_TAG = 'VN'
HEADER_SQ_LINE = '@SQ'

# SAM files have the following tab-delimited columns:
# Column  Value
# 0       QNAME -- Query pair NAME if paired; or Query NAME if unpaired
# 1       FLAG  -- bitwise FLAG (Sam spec., section 2.2.2)
# 2       RNAME -- Reference sequence NAME 3
# 3       POS   -- 1-based leftmost POSition/coordinate of the clipped sequence
# 4       MAPQ  -- MAPping Quality (phred-scaled posterior probability that the mapping position of this read is incorrect)
# 5       CIGAR -- extended CIGAR string
# 6       MRNM  -- Mate Reference sequence NaMe
# 7       MPOS  -- 1-based leftmost Mate POSition of the clipped sequence
# 8       ISIZE -- inferred Insert SIZE
# 9       SEQ   -- query SEQuence; = for match to the reference; n/N/. for ambiguity; cases are not maintained
# 10      QUAL  -- query QUALity; ASCII-33 gives the Phred base quality
# 11-??   Optional TAG:VTYPE:VALUE triplets
#
# A -- Printable character
# i -- Signed 32-bit integer
# f -- Single-precision floating number
# Z -- Printable string, including space
# H -- Byte array in the Hex format
# B -- Integer or numeric array
VALID_VTYPES = ['A', 'i', 'f', 'Z', 'H', 'B']

ALL_COLUMNS      = [QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, MRNM, MPOS, ISIZE, SEQ, QUAL, TAGS] = list(range(12))
REQUIRED_COLUMNS = ALL_COLUMNS[:-1]
INT_COLUMNS      = [FLAG, POS, MPOS, MAPQ]

# Tag used by Cufflinks and SpliceGrapher to identify alignment strand
STRAND_TAG       = 'XS'

# Tag used by SpliceGrapher SAM files to identify junction type
JCT_CODE_TAG       = 'YC'
KNOWN_JCT_TAG      = 'YC:A:K'
RECOMBINED_JCT_TAG = 'YC:A:U'
PREDICTED_JCT_TAG  = 'YC:A:P'

READ_LEN_TAG     = 'NM'

# Match patterns of valid CIGAR symbol sequences
MATCH_RE    = re.compile('^M(NM)*$')
CIGAR_MATCH = re.compile('^([0-9]+[SH])?[0-9]+[MSH]([0-9]+[NID][0-9]+[MSH])*([0-9]+[SH])?$')
CIGAR_TOKEN = re.compile('[0-9]+[A-Z]')

PYSAM_CIGAR = 'MIDNSHP=X'

# Defined in pysam documentation, but evidently
# not in the Python code.
BAM_CMATCH     = 0
BAM_CINS       = 1
BAM_CDEL       = 2
BAM_CREF_SKIP  = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD       = 6
BAM_CEQUAL     = 7
BAM_CDIFF      = 8

# Our own lists for convenience
PYSAM_MERGE_OPS  = [BAM_CINS, BAM_CSOFT_CLIP, BAM_CHARD_CLIP]
PYSAM_IGNORE_OPS = [BAM_CDEL, BAM_CPAD]
PYSAM_SAVE_OPS   = [BAM_CMATCH, BAM_CREF_SKIP]

NULL_CIGAR      = '*'
NULL_CHROMOSOME = '*'
EXACT_CIGAR     = '='

# Flag values used to assess/interpret records
UNMAPPED_FLAG = 4
REVERSE_FLAG  = 16

class ChromosomeTracker(object):
    def __init__(self):
        self.chromosome = None
        self.finished   = {}

    def allFinished(self, cSet):
        """Returns true if all chromosomes in the given set or list are finished; false otherwise."""
        for c in cSet:
            if not self.isFinished(c) : return False
        return True

    def isFinished(self, c):
        """Returns true if the given chromosome has already been seen; false otherwise."""
        try:
            return self.finished[c.lower()]
        except KeyError:
            return False

    def getFinished(self):
        """Returns a list of all chromosomes that have been seen."""
        return [c for c in self.finished if self.finished[c]]

    def finish(self):
        """Tells the tracker to finish off the last chromosome."""
        if self.chromosome:
            self.finished[self.chromosome] = True

    def update(self, c):
        """Main feature of this class: when it detects a change from one
        chromosome to another, marks the old one as finished and starts a new one."""
        chrom = c.lower()
        if self.chromosome == chrom : return False
        # Mismatch: 
        if self.chromosome:
            self.finished[self.chromosome] = True
        self.finished[chrom] = False
        self.chromosome  = chrom
        return True

def acceptSAMRecord(s, counter, **args):
    """Converts a string to a SAM record and returns the record if it's valid;
    returns None for unmapped SAM records; throws an exception for illegal records.
    NB: assumes record string has been strip()-ed."""
    if not s : return None,None # blank lines
    if s[0] == '@' : return None,None # SAM comments

    try:
        rec = SAMRecord(s)
    except ValueError as ve:
        sys.stderr.write('\n*** Error in SAM file at line %d ***\n' % counter)
        sys.stderr.write('Line: %s\n' % s)
        parts = s.split('\t')
        for i in range(len(parts)):
            sys.stderr.write('Column %d:\t%s\n' % ((i+1), parts[i]))
        sys.stderr.write('\n')
        raise ve

    if rec.attrs[RNAME] == NULL_CHROMOSOME : return None,None
    if bool(rec.attrs[FLAG] & UNMAPPED_FLAG) : return None,None
    if not validCigarString(rec.attrs[CIGAR]) : return None,None

    ungap = '%dM' % len(rec.attrs[SEQ])
    if rec.attrs[CIGAR] in [ungap, EXACT_CIGAR]:
        tokens = [ungap]
    else:
        tokens = cigarTokens(rec.attrs[CIGAR])
        # Convert S,H,I and D tokens into M or N:
        try:
            newCigar = convertCigarString(tokens, len(rec.attrs[SEQ]))
        except ValueError as ve:
            raise ve
        rec.attrs[CIGAR] = newCigar
        tokens = cigarTokens(newCigar)
        newLen = sum([int(s[:-1]) for s in tokens if s[-1] == 'M'])
        if newLen < len(rec.attrs[SEQ]):
            rec.attrs[SEQ] = rec.attrs[SEQ][:newLen]
            rec.attrs[QUAL] = rec.attrs[QUAL][:newLen]
        elif newLen > len(rec.attrs[SEQ]):
            delta = newLen - len(rec.attrs[SEQ])
            rec.attrs[SEQ]  += rec.attrs[SEQ][0]*delta
            rec.attrs[QUAL] += rec.attrs[QUAL][0]*delta

    return rec, tokens

def bamIterator(path):
    """Return an iterator over SAM/BAM records."""
    bamStream = pysam.Samfile(path, 'rb')
    chrMap    = pysamChromosomeMap(bamStream)
    headers   = pysamHeaders(bamStream)
    # As in a SAM file: first the headers, then the alignments.
    for h in headers:
        yield h

    for r in bamStream:
        result = pysamReadToString(r,chrMap)
        yield result

def cigarTokens(cigar):
    """Takes a CIGAR string and returns a list of cigar tokens.
    For example, '12M304N20M' will return ['12M','304N','20M']"""
    return [cigar[m.start():m.end()] for m in CIGAR_TOKEN.finditer(cigar)]

def convertCigarString(tokens, requiredSize):
    """Converts a CIGAR string that contains S, H, I and D tokens
    into one that consists only of M and N tokens.  For read-mapping
    purposes, these are the only aspects that really matter.  Assumes
    that the string is a valid CIGAR string."""
    if not tokens:
        raise ValueError('Unrecognized CIGAR string (no tokens): "%s"' % tokens)
    elif len(tokens) == 1:
        if (tokens[0][-1] == 'M') or (tokens[0] == EXACT_CIGAR):
            return tokens[0]
        else:
            raise ValueError('Unrecognized CIGAR string (single token): "%s"' % tokens)

    # NB: to get here, len(tokens) > 1
    result = []
    prev   = tokens[0]
    if prev[-1] == 'M':
        result.append(prev)

    for curr in tokens[1:]:
        symbol = curr[-1]
        try:
            # Tokens must be an integer followed by a single character
            size = int(curr[:-1])
        except ValueError as ve:
            raise ValueError('Unrecognized CIGAR string: "%s"' % s)

        if symbol == 'M':
            # If a match follows a hard/soft clip, just
            # extend the match to include clipped region.
            if prev[-1] in ['H','S']:
                size += int(prev[:-1])
                result.append('%dM'%size)
            elif prev[-1] == 'M':
                # After an insertion/deletion, the
                # previous token might also be a match
                size      += int(prev[:-1])
                result[-1] = '%dM'%size
            else:
                # Ignore an I or N
                result.append('%dM'%size)
        elif symbol == 'I':
            # Extra characters not in reference, so skip them
            prev = result[-1]
            continue # do not update prev
        elif symbol == 'D':
            # A deletion must follow a match and
            # the lengths should be added together.
            #assert(prev[-1] == 'M')
            size      += int(prev[:-1])
            result[-1] = '%dM'%size
            # Pretend this token was a match
            # instead of a delete
            curr       = result[-1]
        elif symbol == 'N':
            # Splice junction added as is
            result.append(curr)
        elif symbol in ['H','S']:
            # Soft and hard clips should only
            # come at the end, after a match
            #assert(prev[-1] == 'M')
            size      += int(prev[:-1])
            result[-1] = '%dM'%size
        prev = curr

    return ''.join(result)

# A -- Printable character
# i -- Signed 32-bit integer
# f -- Single-precision floating number
# Z -- Printable string, including space
# H -- Byte array in the Hex format
# B -- Integer or numeric array
def pysamAttributesToStrings(pysamAttributes):
    """Converts a pysam list of attribute values (as duples) into a
    canonical SAM attribute string."""
    # NB: pysam does not handle chr/str attributes correctly: it
    #     converts chr types to str types, resulting in ':Z:'
    #     attributes in place of ':A:' attributes in output files.
    result = []
    for name,value in pysamAttributes:
        if type(value) == int:
            newString = '%s:i:%d' % (name, value)
        elif type(value) == float:
            newString = '%s:i:%f' % (name, value)
        elif type(value) == list:
            newString = '%s:B:%s' % (name, value)
        elif type(value) == chr:
            newString = '%s:A:%c' % (name, value)
        elif type(value) == str:
            newString = '%s:Z:%s' % (name, value)
        result.append(newString)
    return result

def pysamChromosomeMap(pysamFile):
    """Builds a map of pysam indexes to their corresponding chromosome names."""
    result = {}
    if HEADER_SQ_TAG in pysamFile.header:
        chromList = pysamFile.header[HEADER_SQ_TAG]
        for i in range(len(chromList)):
            result[str(i)] = chromList[i][HEADER_SN_TAG]
    return result

def pysamCigarToString(pysamCigar):
    """Converts a pysam list of CIGAR values (as duples) into a
    canonical CIGAR string.  (The pysam cigarstring attribute does
    not use the canonical form.)"""
    return ''.join(['%d%s' % (d[1], PYSAM_CIGAR[d[0]]) for d in pysamCigar])

def pysamHeaders(pysamFile):
    """Returns a list of header strings associated with a pysam file."""
    result = []
    hdict  = pysamFile.header
    if HEADER_HD_TAG in hdict:
        valdict = hdict[HEADER_HD_TAG]
        hstring = '@%s' % HEADER_HD_TAG
        for k in [HEADER_VN_TAG, HEADER_SO_TAG]:
            if k in valdict : hstring += '\t%s:%s' % (k, valdict[k])
        result.append(hstring + '\n')

    if HEADER_SQ_TAG in hdict:
        seqList = hdict[HEADER_SQ_TAG]
        for seq in seqList:
            # @SQ SN:Chr4 LN:18585056
            result.append('@%s\t%s:%s\t%s:%s\n' % (HEADER_SQ_TAG, HEADER_SN_TAG, seq[HEADER_SN_TAG], HEADER_LN_TAG, seq[HEADER_LN_TAG]))
    return result

def pysamMergeCigar(pysamCigar):
    """Accepts pysam CIGAR duples that include adjustments (soft/hard
    clipping, inserts, deletes) and returns a set of just match/skip duples
    for corresponding locations in the reference sequence."""
    if len(pysamCigar) == 1 : return pysamCigar
    i = 0
    # skip over initial ignored operations
    while i < len(pysamCigar) and pysamCigar[i][0] in PYSAM_IGNORE_OPS : i += 1
    # starts with an intron; something's wrong
    prev = pysamCigar[i]
    if prev[0] == BAM_CREF_SKIP : raise ValueError('CIGAR error: starts with splice junction')

    result = [prev] if prev[0] in PYSAM_SAVE_OPS else []

    for curr in pysamCigar[i+1:]:
        oper   = curr[0]
        length = curr[1]
        if prev[0] in PYSAM_MERGE_OPS:
            if oper == BAM_CREF_SKIP : raise ValueError('CIGAR error: splice junction preceded by clip/insert/delete')
            curr = (oper, length+prev[1])
        elif oper == prev[0]:
            result[-1] = (prev[0], length+prev[1])
            continue
        elif oper in PYSAM_MERGE_OPS and prev[0] == BAM_CMATCH:
            result[-1] = (prev[0], length+prev[1])

        if oper in PYSAM_SAVE_OPS:
            if result and oper == result[-1][0]:
                result[-1] = (oper, length+result[-1][1])
            else:
                result.append(curr)
        prev = curr
    return result

def pysamReadDepths(bamFile, chromosome, gene, **args):
    """Returns a relative start position and an array of read depths
    for the given gene based on input from a BAM file."""
    margin    = getAttribute('margin', 0, **args)
    verbose   = getAttribute('verbose', False, **args)

    loBound   = max(0, gene.minpos-margin)
    upBound   = gene.maxpos+margin
    allreads  = set([r for r in bamFile.fetch(chromosome, loBound, upBound)])
    nSpliced  = 0
    nUngapped = 0
    result    = [0]*(upBound - loBound + 1)
    for r in allreads:
        if len(r.cigar) > 1 : # Spliced alignment
            if pysamStrand(r) != gene.strand : continue
            nSpliced += 1
            cigar     = pysamMergeCigar(r.cigar)
            pos       = r.pos+1
            for tok in cigar:
                if tok[0] == BAM_CMATCH:
                    start = max(pos, loBound) - loBound
                    end   = min(upBound, pos+tok[1]) - loBound
                    for i in range(start, end):
                        result[i] += 1
                pos += tok[1]
        else : # ungapped alignment
            #assert(len(r.cigar) == 1)
            nUngapped += 1
            start = max(r.pos+1, loBound) - loBound
            end   = min(r.pos+r.qlen+1, upBound) - loBound
            for i in range(start, end+1):
                result[i] += 1

    if verbose : sys.stderr.write('Loaded %d ungapped and %d spliced reads for %s\n' % (nUngapped, nSpliced, gene.id))
    return loBound, result

def pysamReadToString(pysamRecord, chromMap):
    """Converts a pysam record into a normal SAM record.  The chromosome
    map links pysam chromosome indexes to their corresponding names
    (see method pysamChromosomeMap)."""
    parts = str(pysamRecord).split('\t')
    try:
        parts[2] = chromMap[parts[2]]
    except KeyError:
        raise ValueError('Chromosome map is missing pysam index %s:\n%s' % (parts[2], chromMap))

    # pysam records adjust the position by -1
    parts[3] = str(int(parts[3])+1)
    parts[5] = pysamCigarToString(pysamRecord.cigar)
    # pysam stores single-end data 
    try:
        if pysamRecord.next_reference_id < 0 : parts[6] = '*'
        if pysamRecord.next_reference_start < 0 : parts[7] = '0'
    except AttributeError:
        sys.stderr.write('ERROR IN PYSAM RECORD: MISSING ATTRIBUTE\n')
        sys.stderr.write('%s\n' % pysamRecord)
        for k in pysamRecord.__dict__:
            sys.stderr.write('  %s = %s\n' % (k, pysamRecord.__dict__[k]))
        sys.exit(1)

    parts[8] = str(pysamRecord.template_length)

    # pysam adjusts the original chars by subtracting 33, so we must add it back in:
    parts[10] = ''.join([chr(c+33) for c in pysamRecord.query_qualities])

    attrs    = pysamAttributesToStrings(pysamRecord.tags)
    parts    = parts[:-1] + attrs
    return '\t'.join(parts)

def pysamSpliceJunctions(pysamRecord, chrMap):
    """Converts a pysam record to a list of splice junction records.
    Returns a list of SpliceJunction objects if the SAM record
    represents a spliced alignment.  Otherwise returns None.
    Note: assumes the CIGAR string has already been validated
    and that there are at least 2 matches in the list."""
    # Example matches: ['16M', '204N', '20M', '355N', '40M']
    cigarDuples = pysamMergeCigar(pysamRecord.cigar)
    sizes = [d[1] for d in cigarDuples]
    pos   = [pysamRecord.pos+1 + sizes[0]-1]
    for b in sizes[1:-1]:
        nextPos = pos[-1] + b
        pos.append(nextPos)

    # Look for junction type: known='YC:A:K' recombined='YC:A:U' predicted='YC:A:P'
    # (Given by SpliceGrapher, but not other programs.)
    tagDict = dict(pysamRecord.tags)
    jctCode = ''
    try:
        jctCode = tagDict[JCT_CODE_TAG]
    except KeyError:
        pass

    result = []
    exons  = [d[1] for d in cigarDuples if d[0] == BAM_CMATCH]
    k      = 0
    chrom  = chrMap[str(pysamRecord.tid)]
    strand = pysamStrand(pysamRecord,tagDict)
    for i in range(0,len(pos),2):
        jct = SpliceJunction(chrom, pos[i], pos[i+1]+1, [exons[k],exons[k+1]], jctCode, strand)
        result.append(jct)
        k  += 1
    return result

def pysamStrand(pysamRecord, tagDict=None):
    """Convenience method returns the strand given by the pysam record."""
    if not tagDict : tagDict = dict(pysamRecord.tags)
    try:
        return tagDict[STRAND_TAG]
    except KeyError:
        return '-' if pysamRecord.is_reverse else '+'

def getNextSamChromosome(samStream, seed=None, **args):
    """Fetches all records for the next chromosome in the
    stream, along with the first record from the next chromosome."""
    verbose = getAttribute('verbose', False, **args)
    result  = []
    if seed:
        result.append(seed)
        targetChrom = seed.split('\t')[2]

    seedLine    = None
    targetChrom = None
    indicator   = ProgressIndicator(10000000, verbose=verbose)
    for line in samStream:
        indicator.update()
        if line.startswith('@') : continue
        s     = line.strip()
        chrom = s.split('\t')[2]
        if not targetChrom:
            targetChrom = chrom

        if chrom == targetChrom:
            result.append(s)
        else:
            seedLine = s
            break
    indicator.finish()
    return result, seedLine

def getSamAlignments(samRecords, **args):
    """Reads a SAM file and fills a dictionary with the start
    position and length of every read in the file.  Spliced
    alignments are stored as separate short reads."""
    verbose     = getAttribute('verbose', False, **args)
    chromosomes = getAttribute('chromosomes', None, **args)
    chromSet    = makeChromosomeSet(chromosomes)
    indicator   = ProgressIndicator(1000000, verbose=verbose)
    tracker     = ChromosomeTracker()
    result      = {}
    for line in samInput(samRecords):
        indicator.update()
        if line.startswith('@') : continue
        rec, matches = acceptSAMRecord(line.strip(), indicator.ctr)
        if not rec : continue

        c   = rec.chromosome()
        if chromSet:
            if tracker.update(c) and tracker.allFinished(chromSet) : break
            if (c not in chromSet) : continue
        if c not in result:
            result[c] = []

        pos     = rec.attrs[POS]
        for m in matches:
            try:
                code   = m[-1]
                delta  = int(m[:-1])
            except Exception as e:
                sys.stderr.write('\n*** Error in SAM file at line %d ***\n' % indicator.ctr)
                sys.stderr.write('Line: %s\n' % line)
                raise e

            # Update reads for matches only
            if code == 'M' : result[c].append((pos,delta))
            pos += delta

    indicator.finish()
    return result

def getSamDepths(samRecords, **args):
    """Reads a SAM file and fills a read-depth array for each
    chromosome found in the file.  Returns a dictionary of
    read-depth lists indexed by chromosome."""
    verbose     = getAttribute('verbose', False, **args)
    maxpos      = getAttribute('maxpos', MAXINT, **args)
    chromosomes = getAttribute('chromosomes', None, **args)
    chromSet    = makeChromosomeSet(chromosomes)

    if isDepthsFile(samRecords):
        depths,jcts = readDepths(samRecords, junctions=False, **args)
        return depths

    depths = {}
    limit  = {}
    if verbose and maxpos < MAXINT : sys.stderr.write('Loading SAM records up to position %d\n' % maxpos)
    indicator = ProgressIndicator(1000000, verbose=verbose)

    tracker = ChromosomeTracker()
    omitted = 0
    total   = 0
    for line in samInput(samRecords):
        indicator.update()
        if line.startswith('@') : continue
        rec, matches = acceptSAMRecord(line.strip(), indicator.ctr)
        if not rec:
            omitted += 1
            continue

        c   = rec.chromosome()

        # If chromosome changes, check to see if we have loaded all of them
        if chromSet:
            if tracker.update(c) and tracker.allFinished(chromSet) : break
            if (c not in chromSet) : continue

        total += 1

        if c not in limit:
            limit[c] = 0
        if not c in depths:
            depths[c] = [0]*(maxpos+1) if maxpos < MAXINT else [0]

        prvPos  = rec.attrs[POS]
        if prvPos > maxpos : break

        for m in matches:
            try:
                code   = m[-1]
                delta  = int(m[:-1])
            except Exception as e:
                sys.stderr.write('\n*** Error in SAM file at line %d ***\n' % indicator.ctr)
                sys.stderr.write('Line: %s\n' % line)
                raise e
            curPos = prvPos + delta

            # Expand lists as necessary:
            while curPos+1 > len(depths[c]):
                depths[c] += [0]*len(depths[c])

            # Update depth for matches only
            if code == 'M':
                for i in range(prvPos, curPos):
                    depths[c][i] += 1
            prvPos = curPos

        # SAM records should be sorted, but a spliced read may
        # have a lower start position than subsequent reads
        # yet still have a higher end position:
        limit[c] = max(limit[c], curPos+1)

    indicator.finish()

    for c in depths:
        maxDepth = max(maxpos,limit[c]) if maxpos < MAXINT else limit[c]
        if maxDepth < len(depths[c]):
            depths[c] = depths[c][:maxDepth]

    if verbose:
        loaded = total-omitted
        sys.stderr.write('Loaded %d records, omitted %d out of %d total\n' % (loaded, omitted, total))

    return depths

def getSamHeaders(samRecords, **args):
    """Reads a SAM file and returns just the header strings as a list."""
    result = []
    instream = samInput(samRecords)
    for line in instream:
        if not line.startswith('@') : break
        result.append(line.strip())
    if type(instream) == file:
        instream.close()
    return result

def getSamHeaderInfo(samStream, **args):
    """Parses header records in a SAM file and returns a list of
    possible chromosomes along with the first alignment record."""
    verbose  = getAttribute('verbose', False, **args)
    seedLine = None
    result   = set()
    for line in samStream:
        if line.startswith('@'):
            if line.startswith('@SQ'):
                parts = line.strip().split('\t')
                #assert(len(parts) == 3)
                #assert(parts[1].startswith('SN'))
                result.add(parts[1].split(':')[1])
        else:
            seedLine = line
            break
    if verbose:
        sys.stderr.write('Found %d chromosomes in SAM header:\n' % len(result))
        sys.stderr.write('  %s\n' % ','.join(sorted(result)))
    return result, seedLine

def getSamJunctions(samRecords, **args):
    """Reads a SAM file and builds lists of splice-junctions
    for each chromosome found in the file.  Returns a dictionary of
    junction lists indexed by chromosome."""
    verbose     = getAttribute('verbose', False, **args)
    maxpos      = getAttribute('maxpos', MAXINT, **args)
    minanchor   = getAttribute('minanchor', 0, **args)
    minjct      = getAttribute('minjct', 1, **args)
    chromosomes = getAttribute('chromosomes', None, **args)

    if isDepthsFile(samRecords):
        depths,jcts = readDepths(samRecords, depths=False, **args)
        return jcts

    # Restrict records to the given chromosomes
    chromSet = makeChromosomeSet(chromosomes)

    if verbose : sys.stderr.write('Loading splice junctions with minimum anchor %d and minimum junction support %d\n' % (minanchor,minjct))

    indicator  = ProgressIndicator(1000000, verbose=verbose)
    junctions  = {}
    tracker = ChromosomeTracker()
    for line in samInput(samRecords):
        indicator.update()
        if line.startswith('@') : continue
        rec, matches = acceptSAMRecord(line.strip(), indicator.ctr)
        if not rec or len(matches) < 2 : continue

        if rec.attrs[POS] > maxpos : break
        c   = rec.chromosome()

        # If chromosome changes, check to see if we have loaded all of them
        if chromSet:
            if tracker.update(c) and tracker.allFinished(chromSet) : break
            if (c not in chromSet) : continue

        jctList = recordToSpliceJunction(rec, matches)
        if not jctList : continue

        strand = rec.attrs[STRAND_TAG]
        if c not in junctions : junctions[c] = {}
        if strand not in junctions[c] : junctions[c][strand] = {}
        for newJct in jctList:
            junctions[c][strand].setdefault(newJct.p1,{})
            junctions[c][strand][newJct.p1].setdefault(newJct.p2,None)

            if junctions[c][strand][newJct.p1][newJct.p2] is None:
                junctions[c][strand][newJct.p1][newJct.p2] = newJct
            else:
                junctions[c][strand][newJct.p1][newJct.p2].update(newJct)

    result = {}
    for c in sorted(junctions.keys()):
        result[c] = []
        for s in sorted(junctions[c].keys()):
            for a in sorted(junctions[c][s].keys()):
                for b in sorted(junctions[c][s][a].keys()):
                    jct = junctions[c][s][a][b]
                    if jct.count >= minjct and jct.minAnchor() >= minanchor:
                        result[c].append(jct)

    indicator.finish()
    return result

def getSamReadData(samRecords, **args):
    """Reads a SAM file and fills a read-depth array for each
    chromosome found in the file along with a junction array
    for each chromosome.  By default, returns two dictionaries:
    one with read-depth lists indexed by chromosome and another
    with junction lists indexed by chromosome.  If the 'alignments'
    option is set, returns a third dictionary of (pos,len) tuples
    representing every read alignment."""
    alignments  = getAttribute('alignments', False, **args)
    chromosomes = getAttribute('chromosomes', None, **args)
    chromSet    = makeChromosomeSet(chromosomes)
    getJct      = getAttribute('junctions', True, **args)
    maxpos      = getAttribute('maxpos', MAXINT, **args)
    minanchor   = getAttribute('minanchor', 0, **args)
    minjct      = getAttribute('minjct', 1, **args)
    verbose     = getAttribute('verbose', False, **args)

    if isDepthsFile(samRecords):
        if alignments : raise ValueError('Cannot load alignment counts from a depths file')
        depths,jcts = readDepths(samRecords, **args)
        return depths,jcts

    align  = {}
    depths = {}
    jctTmp = {}
    limit  = {}
    if verbose:
        if maxpos < MAXINT:
            sys.stderr.write('   loading SAM records for positions up to %d\n' % maxpos)
        sys.stderr.write('   splice junctions with minimum anchor %d and minimum junction support %d\n' % (minanchor,minjct))

    indicator = ProgressIndicator(1000000, verbose=verbose)
    tracker   = ChromosomeTracker()
    curPos    = 0
    for line in samInput(samRecords):
        indicator.update()
        if line.startswith('@') : continue
        try:
            rec, matches = acceptSAMRecord(line.strip(), indicator.ctr)
        except ValueError:
            continue
        if not rec : continue

        c = rec.chromosome()

        # If chromosome changes, check to see if we have loaded all of them
        if chromSet:
            if tracker.update(c) and tracker.allFinished(chromSet) : break
            if (c not in chromSet) : continue

        if not c in depths:
            depths[c] = [0]*(maxpos+1) if maxpos < MAXINT else [0]

        if c not in limit:
            limit[c] = 0

        if alignments and not c in align:
            align[c] = []

        prvPos  = rec.attrs[POS]
        if prvPos > maxpos : break

        for m in matches:
            try:
                code   = m[-1]
                delta  = int(m[:-1])
            except Exception as e:
                sys.stderr.write('\n*** Error in SAM file: invalid CIGAR string (column 5) at line %d ***\n' % indicator.ctr)
                sys.stderr.write('Invalid record: %s\n' % line)
                sys.stderr.write('CIGAR string:   %s\n' % rec.cigar())
                raise e

            curPos = prvPos + delta

            # Expand lists as necessary:
            while curPos+1 > len(depths[c]):
                depths[c] += [0]*len(depths[c])

            # Update depth for matches only
            if code == 'M':
                for i in range(prvPos, curPos):
                    depths[c][i] += 1
                if alignments:
                    align[c].append((prvPos,delta))
            prvPos   = curPos

        # SAM records should be sorted, but a spliced read may
        # have a lower start position than subsequent reads
        # yet still have a higher end position:
        limit[c] = max(limit[c],curPos+1)

        # Now look for junctions
        if getJct and len(matches) > 1:
            jctList  = recordToSpliceJunction(rec, matches)
            if not jctList : continue

            strand = rec.attrs[STRAND_TAG]
            if not c in jctTmp : jctTmp[c] = {}
            if not strand in jctTmp[c] : jctTmp[c][strand] = {}
            for newJct in jctList:
                jctTmp[c][strand].setdefault(newJct.p1,{})
                jctTmp[c][strand][newJct.p1].setdefault(newJct.p2,None)

                if jctTmp[c][strand][newJct.p1][newJct.p2] is None:
                    jctTmp[c][strand][newJct.p1][newJct.p2] = newJct
                else:
                    jctTmp[c][strand][newJct.p1][newJct.p2].update(newJct)

    indicator.finish()

    for c in depths:
        maxDepth = max(maxpos,limit[c]) if maxpos < MAXINT else limit[c]
        if maxDepth < len(depths[c]):
            depths[c] = depths[c][:maxDepth]

    junctions = {}
    if getJct:
        for c in sorted(jctTmp.keys()):
            junctions[c] = []
            for s in sorted(jctTmp[c].keys()):
                for a in sorted(jctTmp[c][s].keys()):
                    for b in sorted(jctTmp[c][s][a].keys()):
                        jct = jctTmp[c][s][a][b]
                        if jct.count >= minjct and jct.minAnchor() >= minanchor:
                            junctions[c].append(jctTmp[c][s][a][b])

    # Really need a better way to do this
    if alignments:
        return depths, junctions, align
    else:
        return depths, junctions

def getSamSequences(samRecords, **args):
    """Reads a SAM file and returns a dictionary of sequence names and lengths."""
    result  = {}
    headers = getSamHeaders(samRecords, **args)
    for line in headers:
        if not line.startswith(HEADER_SQ_LINE) : continue
        parts = line.strip().split('\t')
        sdict = dict([tuple(s.split(':')) for s in parts[1:]])
        try:
            result[sdict[HEADER_SN_TAG]] = sdict[HEADER_LN_TAG]
        except KeyError:
            continue
    return result

def isBamFile(filePath):
    """Simple heuristic returns True if path is a BAM file; false otherwise."""
    return filePath.lower().endswith('.bam')

def loadSAMRecords(f):
    """Loads SAM records from a file and returns them in a list."""
    hstring = None
    result  = []
    for s in ezopen(f):
        if not s.startswith('@'):
            result.append(SAMRecord(s))
    return result

def makeChromosomeSet(chromList):
    """Converts a string or a list of chromosomes into a set of unique values."""
    chromSet = None
    if chromList:
        if type(chromList) == type(' '):
            chromList = [chromList]
        chromSet = set([c.lower() for c in chromList])
    return chromSet

def recordToSpliceJunction(samRec, matches):
    """Converts a SAM record to a list of splice junction records.
    Returns a list of SpliceJunction objects if the SAM record
    represents a spliced alignment.  Otherwise returns None.
    Note: assumes the CIGAR string has already been validated
    and that there are at least 2 matches in the list."""

    # Example matches: ['16M', '204N', '20M', '355N', '40M']
    sizes = [int(m[:-1]) for m in matches]
    pos   = [samRec.attrs[POS] + sizes[0]-1]
    for b in sizes[1:-1]:
        nextPos = pos[-1] + b
        pos.append(nextPos)

    # Look for junction type: known='YC:A:K' recombined='YC:A:U' predicted='YC:A:P'
    # (Given by SpliceGrapher, but not other programs.)
    try:
        jctCode = samRec[JCT_CODE_TAG]
    except KeyError:
        jctCode = ''

    result = []
    exons  = [int(m[:-1]) for m in matches if m[-1]=='M']
    k      = 0
    for i in range(0,len(pos),2):
        jct = SpliceJunction(samRec.chromosome(), pos[i], pos[i+1]+1, [exons[k],exons[k+1]], jctCode, samRec.attrs[STRAND_TAG])
        result.append(jct)
        k  += 1
    return result

def samInput(source):
    """Convenience method that returns an iterator over a set of SAM records.
    If source is a list, set or file stream it returns the object.  A string
    is interpreted as a file path to be opened for input and the stream returned."""
    if type(source) == list:
        return source

    if type(source) == set:
        return source

    if type(source) == TextIOWrapper:
        return source

    if type(source) == str:
        return samIterator(source)

    raise ValueError('Unrecognized type %s for SAM input' % type(source))

def samIterator(path):
    """Return an iterator over SAM/BAM records."""
    if isBamFile(path):
        return bamIterator(path)
    else:
        return ezopen(path)

def validCigarString(s):
    """Returns true if the given CIGAR string is one that we can use; false otherwise."""
    if s == NULL_CIGAR:
        return False
    elif s == EXACT_CIGAR:
        return True
    else:
        return (CIGAR_MATCH.match(s) is not None)

class SAMRecord(object):
    """Encapsulates all the information in a SAM record."""
    def __init__(self,s):
        self.attrs  = {}
        self.vtypes = {}
        parts       = s.split('\t')
        col         = 0
        while col in REQUIRED_COLUMNS:
            try:
                if col in INT_COLUMNS:
                    self.attrs[col]  = int(parts[col])
                    self.vtypes[col] = 'i'
                else:
                    self.attrs[col]  = parts[col]
                    self.vtypes[col] = 'A'
            except IndexError as ie:
                raise ValueError('Too few columns (%d<%d) in input record:\n%s' % (len(parts), len(REQUIRED_COLUMNS), s))
            col += 1

        while col < len(parts):
            triplets = parts[col].split()
            for triplet in triplets:
                try:
                    (tag,vtype,val)  = triplet.split(':')
                    if vtype not in VALID_VTYPES : raise ValueError()
                    self.vtypes[tag] = vtype
                    self.attrs[tag]  = val
                except ValueError:
                    raise ValueError('Illegal tag:type:value SAM attribute "%s"; record contains %d columns' % (triplet,len(parts)))
            col += 1

        # Only do these once:
        if STRAND_TAG not in self.attrs:
            self.attrs[STRAND_TAG]  = '-' if bool(self.attrs[FLAG] & REVERSE_FLAG) else '+'
            self.vtypes[STRAND_TAG] = 'A'

    def chromosome(self):
        """Convenience method returns the chromosome."""
        return self.attrs[RNAME].lower()

    def cigar(self):
        """Convenience method returns the cigar string."""
        return self.attrs[CIGAR]

    def __cmp__(self, other):
        return self.attrs[POS] - other.attrs[POS]

    def flag(self):
        """Convenience method returns the bitwise flag value as an int."""
        return self.attrs[FLAG]

    def __getitem__(self, key):
        if not key in self.vtypes:
            return self.attrs[key]
        elif self.vtypes[key] in ['A','Z']:
            return self.attrs[key]
        elif self.vtypes[key] == 'i':
            return int(self.attrs[key])
        elif self.vtypes[key] == 'f':
            return float(self.attrs[key])
        elif self.vtypes[key] == 'H':
            return int(self.attrs[key],16)
        else:
            return self.attrs[key]

    # Required for creating dicts/sets
    def __hash__(self):
        return str(self).__hash__()

    def matchpos(self):
        """Convenience method returns the mate position value as an int."""
        return self.attrs[MPOS]

    def pos(self):
        """Convenience method returns the position value as an int."""
        return self.attrs[POS]

    def quality(self):
        """Convenience method returns the phred quality value as an int."""
        return self.attrs[MAPQ]

    def query(self):
        """Convenience method returns the original query sequence."""
        return self.attrs[SEQ]

    def read(self):
        """Convenience method returns the read id."""
        return self.attrs[QNAME]

    def setAttributeString(self, valueString):
        """Sets an attribute in SAM key:vtype:value format.
        Accordingly, the valueString must have the form 'key:vtype:value'."""
        self.setAttribute(valueString.split(':'))

    def setAttribute(self, values):
        """Sets an attribute based on [key, type, value] values using SAM key:type:value format."""
        if type(values) not in [set,list] : raise ValueError('Value must be a set or a list; received %s' % type(values))
        if not (2 <= len(values) <= 3) : raise ValueError('Value list must have 2 or 3 elements; received %d' % len(values))
        key,vtype,value = values
        if vtype not in VALID_VTYPES : raise ValueError('Invalid vtype: %s' % vtype)
        self.vtypes[key] = vtype
        self.attrs[key]  = value

    def __str__(self):
        reqd = [str(self.attrs[x]) for x in REQUIRED_COLUMNS]
        tags = ['%s:%s:%s' % (x,self.vtypes[x],self.attrs[x]) for x in self.attrs if x not in REQUIRED_COLUMNS]
        vals = reqd + tags
        return '\t'.join([x for x in vals if x is not None])

    def strand(self):
        """Convenience method returns the strand given by the bitwise flag."""
        return self.attrs[STRAND_TAG]
