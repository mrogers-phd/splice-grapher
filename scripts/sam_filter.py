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
from SpliceGrapher.predict.SpliceSiteValidator import *
from SpliceGrapher.predict.SiteClassifier      import *
from SpliceGrapher.formats.loader              import *
from SpliceGrapher.formats.FastaLoader         import *
from SpliceGrapher.shared.config               import *
from SpliceGrapher.shared.utils                import *
from SpliceGrapher.shared.streams              import *
from SpliceGrapher.formats.sam                 import *

import sys, gzip, zipfile
from optparse import OptionParser

try :
    from PyML import *
    from PyML.containers import SequenceData
except Exception :
    sys.stderr.write('\n** Unable to load PyML modules required for this script.\n')
    sys.exit(1)

try :
    from pysam import *
    PYSAM_LOADED = True
except Exception :
    PYSAM_LOADED = False

# Match patterns of valid CIGAR symbol sequences
EXACT_CIGAR = '='
NULL_CIGAR  = '*'
NOMINAL_POSITIVE = 1.0

def cmpQuintuple(a,b) :
    # quintuple: chromosome, strand, position, type, score
    # sort by chromosome, then position, then strand
    if a[0] == b[0] :
        return cmp(a[1],b[1]) if a[2]==b[2] else a[2]-b[2]
    else :
        return cmp(a[0],b[0])

def jctString(c,d,a,s) :
    return '%s;%d;%d;%s' % (c,d,a,s)

def mapSitesToGenes(geneModel) :
    result = {}
    for chrom in geneModel.getChromosomes() :
        result[chrom] = {'-':{}, '+':{}}
        genes         = geneModel.getGeneRecords(chrom, gene_type_filter)
        for g in genes :
            sites = g.donorList() + g.acceptorList()
            for x in sites :
                result[chrom][g.strand].setdefault(x,set([]))
                result[chrom][g.strand][x].add(g)
    return result

def bamIterator(path) :
    """Return an iterator over SAM/BAM records."""
    bamStream = Samfile(path, 'rb')
    chrMap    = pysamChromosomeMap(bamStream)
    headers   = pysamHeaders(bamStream)
    for h in headers :
        yield h
    for r in bamStream :
        yield pysamReadToString(r,chrMap)

def samIterator(path, isBam=False) :
    """Return an iterator over SAM/BAM records."""
    if isBam :
        return bamIterator(path)
    else :
        return ezopen(path)

def siteString(c,p,s,stype='d') :
    return '%s;%d;%s;%s' % (c,p,s,stype)

def validSpliceSite(chrom, strand, pos, classifiers, storedValues, **args) :
    """
    Classifies a single splice site given by a chromosome, strand and
    position along the strand, using the set of classifiers provided.
    To improve speed, it uses a 'storedValues' dictionary to avoid
    duplicating predictions for the same location.  Note: the classifiers
    and stored values for donor and acceptor sites should be distinct.
    """
    verbose  = getAttribute('verbose', False, **args)
    siteType = getAttribute('siteType', DONOR_SITE, **args)

    try :
        return (storedValues[chrom][strand][pos] > 0)
    except KeyError :
        storedValues.setdefault(chrom,{})
        storedValues[chrom].setdefault(strand,{})
        ## Maybe add this as an option, but gene models should
        ## validate all non-canonical sites that don't have classifiers
        ## # Assume valid unless proved otherwise.  Allows
        ## # non-canonical sites without classifiers to pass through:
        ## storedValues[chrom][strand][pos] = NOMINAL_POSITIVE
        storedValues[chrom][strand][pos] = 0.0

    for svm in classifiers.values() :
        # The intron/exon size both depend on which classifier is used;
        # unless there are more than 2 classifiers for a site type,
        # this approach is faster in the long run than getting the
        # dimer first and then finding the appropriate classifier.
        (exonSeq, intronSeq, dimer) = getSpliceSiteParts(chrom, pos, siteType, strand, getSequence,
                                                         exonWindow=svm.config.exonSize(),
                                                         intronWindow=svm.config.intronSize())
        if dimer.lower() != svm.config.dimer().lower() : continue

        newSeq = intronSeq + exonSeq if siteType == ACCEPTOR_SITE else exonSeq + intronSeq

        hideStdout() # Hide PyML C stream
        if svm.config.mismatchProfile() :
            ssData = positionalKmerData(svm.config, [newSeq])
        else :
            ssData = SequenceData([newSeq], mink=svm.config.mink(), maxk=svm.config.maxk(), maxShift=svm.config.maxShift())

        if svm.config.normalize() :
            ssData.attachKernel('cosine')
            
        (ssClass,score) = svm.classify(ssData,0)
        showStdout() # Reinstate stream

        storedValues[chrom][strand][pos] = score
        break

    return (storedValues[chrom][strand][pos] > 0)

def getSequence(chrom, pos1, pos2, strand) :
    return seqDict.subsequence(chrom, pos1, pos2, reverse=(strand=='-'))

USAGE = """%prog SAM/BAM-file classifiers [options]

Removes false-positive spliced alignments from a SAM or BAM file.
The classifers may be given as a zip-file or as a list of classifier
configuration (.cfg) files."""

# Establish command-line options:
parser = OptionParser(usage=USAGE)
parser.add_option('-C', dest='cigar',   default=None,          help='File for storing unrecognized CIGAR strings [default: %default]')
parser.add_option('-f', dest='fasta',   default=SG_FASTA_REF,  help='Reference genome FASTA [default: %default]')
parser.add_option('-F', dest='fpsites', default=None,          help='File for storing false-positive SAM alignment records [default: %default]')
parser.add_option('-m', dest='model',   default=SG_GENE_MODEL, help='Gene model file (GTF/GFF3 format) [default: %default]')
parser.add_option('-o', dest='output',  default=None,          help='Output file [default: stdout]')
parser.add_option('-r', dest='report',  default=None,          help='Write classifier scores to file [default: %default]')
parser.add_option('-v', dest='verbose', default=False,         help='Verbose mode [default: %default]', action='store_true')
parser.add_option('-z', dest='gzip',    default=False,         help='Apply gzip compression to output [default: %default]', action='store_true')
opts, args = parser.parse_args(sys.argv[1:])

if len(args) < 2 :
    parser.print_help()
    sys.exit(1)

if not opts.fasta :
    parser.print_help()
    sys.stderr.write('You must specify a reference genome (-f option).\n')
    sys.exit(1)

if opts.gzip and not opts.output :
    parser.print_help()
    sys.stderr.write('You must specify an output file if you want to use GZIP compression.\n')
    sys.exit(1)

for f in args : validateFile(f)

# Check SAM/BAM formats
samFile   = args[0]
bamFormat = samFile.lower().endswith('.bam')
if bamFormat : sys.stderr.write('\nUsing BAM format for input file\n')
if bamFormat and not PYSAM_LOADED :
    parser.print_help()
    sys.stderr.write('Pysam (required for BAM files) not found on your system.\n')
    sys.exit(1)

if opts.model :
    validateFile(opts.model)
else :
    sys.stderr.write('No gene models provided; relying entirely on classifiers.\n')

validateFile(opts.fasta)

writeStartupMessage()

cigarStream = open(opts.cigar,'w') if opts.cigar else None
badStream   = None if not opts.fpsites else open(opts.fpsites,'w')
seqDict     = FastaLoader(opts.fasta, verbose=opts.verbose)

# Accept either regular config files or zip archives
duplicates = set()
cfgFiles   = set()
for f in args[1:] :
    if zipfile.is_zipfile(f) :
        newFiles = set(unzipClassifiers(f))
    else :
        newFiles = set([f])
    duplicates.update(newFiles & cfgFiles)
    cfgFiles |= newFiles

if duplicates :
    sys.stderr.write('** Warning: found duplicate configurations for %s\n' % ', '.join(duplicates))

# Load splice site classifiers for each kind of site
donModels = {}
accModels = {}
hideStdout()
for cfg in sorted(cfgFiles) :
    model = SiteClassifier(cfg, verbose=opts.verbose)
    if model.config.acceptor() :
        accModels[cfg] = model
    else :
        donModels[cfg] = model

    if opts.verbose :
        stype = 'acceptor' if model.config.acceptor() else 'donor'
        sys.stderr.write("loaded classifier for '%s' %s sites\n" % (model.config.dimer(),stype))
showStdout()

# Dictionaries for storing splice site classification results on the fly:
donSites = {}
accSites = {}

# Load gene models
knownDon = {}
knownAcc = {}
if opts.model :
    geneModel  = loadGeneModels(opts.model, verbose=opts.verbose)
    knownAcc   = geneModel.getAllAcceptors(gene_type_filter)
    knownDon   = geneModel.getAllDonors(gene_type_filter)
    siteMap    = mapSitesToGenes(geneModel)
    knownChrom = set(geneModel.getChromosomes())

# Validate junctions
if opts.verbose : sys.stderr.write('Validating junctions in %s\n' % samFile)

novelJct       = set([])
validJct       = set([])
invalidJct     = set([])
intergenicJct  = set([])
transgenicJct  = set([])
knownDonors    = set([])
knownAcceptors = set([])
predDonors     = set([])
predAcceptors  = set([])
ungapped       = 0
spliced        = 0
badChrom       = set([])

if opts.gzip :
    outStream = gzip.open(opts.output,'w')
elif opts.output :
    outStream = open(opts.output,'w')
else :
    outStream = sys.stdout

invalidCigar = 0
unaligned    = 0
unprocessed  = set()
indicator    = ProgressIndicator(1000000, verbose=opts.verbose)
for line in samIterator(samFile, isBam=bamFormat) :
    indicator.update()
    s = line.strip()
    if not s : continue

    if s[0] == '@' :
        outStream.write(line)
        continue

    try :
        rec, matches = acceptSAMRecord(s, indicator.ctr)
    except ValueError,ve :
        invalidCigar += 1
        continue

    if not rec or rec.attrs[CIGAR] == NULL_CIGAR :
        unaligned += 1
        continue

    if len(matches) == 1 :
        outStream.write('%s\n' % rec)
        ungapped += 1
        continue

    jctList = recordToSpliceJunction(rec, matches)
    if not jctList : continue

    # Assume a junction is valid until proved otherwise
    # (necessary to handle long reads with 2+ junctions).
    validJunctions = True
    for jct in jctList :
        spliced += 1
        chrom = jct.chromosome
        don   = jct.donor()      if jct.strand == '+' else jct.donor()-2
        acc   = jct.acceptor()-2 if jct.strand == '+' else jct.acceptor()

        donGene = None
        accGene = None
        if opts.model :
            try :
                donGene = geneModel.getGeneFromLocations(chrom, don, don, jct.strand)
                accGene = geneModel.getGeneFromLocations(chrom, acc, acc, jct.strand)
            except KeyError, ke :
                #sys.stderr.write('\nError looking for genes on chromosome %s (strand %s) at positions %d or %d' \
                #        % (chrom, jct.strand, don, acc))
                #sys.stderr.write('Offending SAM record:\n%s\n' % line)
                #raise ke
                unprocessed.add(jct)
                continue

        isNovel = False
        okDonor = True
        if (chrom in knownDon) and (don in knownDon[chrom][jct.strand]) :
            knownDonors.add(siteString(chrom, don, jct.strand))
        elif validSpliceSite(chrom, jct.strand, don, donModels, donSites, siteType=DONOR_SITE) :
            isNovel = True
            predDonors.add(siteString(chrom, don, jct.strand))
        else :
            isNovel = True
            okDonor = False

        okAcceptor = True
        if (chrom in knownAcc) and (acc in knownAcc[chrom][jct.strand]) :
            knownAcceptors.add(siteString(chrom, acc, jct.strand,'a'))
        elif validSpliceSite(chrom, jct.strand, acc, accModels, accSites, siteType=ACCEPTOR_SITE) :
            predAcceptors.add(siteString(chrom, acc, jct.strand, 'a'))
            isNovel = True
        else :
            isNovel    = True
            okAcceptor = False

        if donGene and accGene and (donGene == accGene) :
            if isNovel :
                novelJct.add(jctString(chrom,don,acc,jct.strand))
                rec.setAttributeString(PREDICTED_JCT_TAG)
            else :
                introns = donGene.getIntrons()
                duple   = (jct.donor(),jct.acceptor())
                if duple in introns :
                    rec.setAttributeString(KNOWN_JCT_TAG)
                else :
                    rec.setAttributeString(RECOMBINED_JCT_TAG)

            if okDonor and okAcceptor :
                try :
                    donGenes    = siteMap[chrom][jct.strand][don]
                    accGenes    = siteMap[chrom][jct.strand][acc]
                    commonGenes = accGenes.intersection(donGenes)
                except KeyError :
                    commonGenes = set()
        
                if commonGenes :
                    validJct.add(jctString(chrom,don,acc,jct.strand))
                else :
                    transgenicJct.add(jctString(chrom,don,acc,jct.strand))
            else :
                validJunctions = False
                # Use line.strip() in case 
                if badStream : badStream.write('%s\n' % s)
                invalidJct.add(jctString(chrom,don,acc,jct.strand))
        else :
            validJunctions = (okDonor and okAcceptor)
            intergenicJct.add(jctString(chrom,don,acc,jct.strand))

    # Don't write the record unless all junctions are valid
    if validJunctions :
        outStream.write('%s\n' % rec)
indicator.finish()

if opts.report :
    if opts.verbose : sys.stderr.write('Writing classifier scores to %s\n' % opts.report)
    quintuples = []
    for c in donSites.keys() :
        for s in donSites[c].keys() :
            for p in donSites[c][s].keys() :
                quintuples.append((c,s,p,'d',donSites[c][s][p]))
    for c in accSites.keys() :
        for s in accSites[c].keys() :
            for p in accSites[c][s].keys() :
                quintuples.append((c,s,p,'a',accSites[c][s][p]))
    quintuples.sort(cmp=cmpQuintuple)
    rstream = open(opts.report,'w')
    rstream.write('%s\n' % '\n'.join(['%s\t%s\t%d\t%s\t%s' % q for q in quintuples]))
    rstream.close()

if opts.verbose :
    sys.stderr.write('Found %s ungapped and %s spliced alignments in %s\n' % (commaFormat(ungapped), commaFormat(spliced), args[0]))

    # First count all valid junctions
    validCount   = len(validJct)
    invalidCount = len(invalidJct)
    total        = validCount + invalidCount
    validPct     = float(100.0*validCount)/total if total > 0 else 0
    invalidPct   = float(100.0*invalidCount)/total if total > 0 else 0
    sys.stderr.write('%s (%.1f%%) TP and %s (%.1f%%) FP junction alignments out of %s total\n' % \
            (commaFormat(validCount), validPct, commaFormat(invalidCount), invalidPct, commaFormat(total)))

    # Next count novel junctions (invalid junctions are novel implicitly)
    novelCount   = len(novelJct)
    total        = novelCount
    validCount   = total - invalidCount
    validPct     = float(100.0*validCount)/total if total > 0 else 0
    invalidPct   = float(100.0*invalidCount)/total if total > 0 else 0
    sys.stderr.write('%s (%.1f%%) TP and %s (%.1f%%) FP novel junction alignments out of %s total\n' % \
            (commaFormat(validCount), validPct, commaFormat(invalidCount), invalidPct, commaFormat(total)))
    sys.stderr.write('\nAlso found %d intergenic and %d transgenic junctions\n' % (len(intergenicJct),len(transgenicJct)))

    if invalidCigar > 0 :
        sys.stderr.write('** Omitted %s records with unrecognized CIGAR patterns.\n' % commaFormat(invalidCigar))

    if unaligned > 0 :
        sys.stderr.write("** Found %s records reported as unaligned CIGAR ('%s') patterns.\n" % (commaFormat(unaligned), NULL_CIGAR))

    if len(unprocessed) > 0 :
        sys.stderr.write('** There were unprocessed records in the following chromosomes:\n')
        chromSet = set([r.chromosome for r in unprocessed])
        sys.stderr.write('  %s\n' % '\n  '.join(sorted(chromSet)))
