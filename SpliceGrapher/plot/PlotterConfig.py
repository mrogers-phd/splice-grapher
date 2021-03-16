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
Module that encapsulates a configuration for a meta-plotter.
"""
from SpliceGrapher.shared.utils import getAttribute
from SpliceGrapher.shared.adjust import MIN_INTRON_SIZE
from ConfigParser import ConfigParser
import re, os

ENV_MATCH = re.compile('\${\w\w*}')

OUTPUT_FORMATS = ['emf', 'eps', 'pdf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz']

# Main section names
SINGLE_SECTION = 'SinglePlotConfig'
MULTI_SECTION  = 'MultiplotConfig'

# Main section tags
VALID_SINGLE_TAGS = [FONT_TAG, LEGEND_TAG, LOGFILE_TAG, HEIGHT_TAG, OUTPUT_TAG, SHRINK_TAG, SFACTOR_TAG, WIDTH_TAG] \
        = ['fontsize', 'legend', 'log_file', 'height', 'output_file', 'shrink_introns', 'shrink_factor', 'width']

VALID_MULTI_TAGS  = [CHROM_TAG, EXT_TAG, GENELIST_TAG] = ['chromosome', 'output_format', 'gene_list']
VALID_MULTI_TAGS += VALID_SINGLE_TAGS

# Plot value names
VALID_PLOT_TAGS = [ACCEPTOR_TAG, ANNOTATE_TAG, BACKGROUND_TAG, CODONS_TAG, COVERAGE_TAG, CLUSTER_TAG, DONOR_TAG,
                   EDGE_TAG, FEATURE_TAG, FORMAT_TAG, GENE_TAG, HIDE_TAG, HIGHLIGHT_TAG,
                   ISOLABEL_TAG, ISOWEIGHT_TAG, LABEL_TAG, LINE_TAG, LOGSCALE_TAG, NEW_JCT_TAG, SIZE_TAG,
                   PLOT_TAG, SOURCE_TAG, TITLE_TAG, XLABEL_TAG, YLIMIT_TAG] \
                = ['acceptors', 'annotate', 'background', 'codons', 'min_coverage', 'clusters', 'donors',
                   'edge_weight', 'features', 'file_format', 'gene_name', 'hide', 'highlight',
                   'iso_labels', 'iso_weights', 'labels', 'line', 'log_threshold', 'novel_jct', 'relative_size',
                   'plot_type', 'source_file', 'title_string', 'x_labels', 'y_limit']

# Recognized plot types:
VALID_PLOT_TYPES = [GENE_MODEL_PLOT, ISOFORM_PLOT, JUNCTION_PLOT, READ_DEPTH_PLOT, SPLICEGRAPH_PLOT, XY_PLOT] \
        = ['gene', 'isoforms', 'junctions', 'read_depth', 'splice_graph', 'xy_plot']

VALID_FILE_FORMATS = [CSV_FORMAT, GENE_MODEL_FORMAT, SAM_FORMAT, SPLICEGRAPH_FORMAT] = ['CSV', 'gene_model', 'SAM', 'splice_graph']
MULTI_FILE_FORMATS = [SPLICEGRAPH_FORMAT]

# Gene models may be specified either using a gene model (GTF/GFF)
# representation or a splice graph representation:
VALID_GENE_FORMATS = [GENE_MODEL_FORMAT, SPLICEGRAPH_FORMAT]

# Plots that provide gene min/max range
MINMAX_PLOT_TYPES  = [GENE_MODEL_PLOT, ISOFORM_PLOT, SPLICEGRAPH_PLOT]

# Parameter type lists
BOOLEAN_TAGS       = [ANNOTATE_TAG, BACKGROUND_TAG, CODONS_TAG, HIDE_TAG, ISOLABEL_TAG, LABEL_TAG, LINE_TAG, LEGEND_TAG, NEW_JCT_TAG,
                      SHRINK_TAG, XLABEL_TAG]
FLOAT_TAGS         = [FONT_TAG, LOGSCALE_TAG, HEIGHT_TAG, WIDTH_TAG, SIZE_TAG, YLIMIT_TAG]
INTEGER_TAGS       = [EDGE_TAG, COVERAGE_TAG, SFACTOR_TAG]
STRING_TAGS        = [ACCEPTOR_TAG, DONOR_TAG, FEATURE_TAG, FORMAT_TAG, GENE_TAG, GENELIST_TAG,
                      HIGHLIGHT_TAG, ISOWEIGHT_TAG, LOGFILE_TAG, OUTPUT_TAG, PLOT_TAG, SOURCE_TAG, TITLE_TAG]

TAG_TYPE = {}
for t in BOOLEAN_TAGS : TAG_TYPE[t] = 'boolean'
for t in FLOAT_TAGS   : TAG_TYPE[t] = 'float'
for t in INTEGER_TAGS : TAG_TYPE[t] = 'int'
for t in STRING_TAGS  : TAG_TYPE[t] = 'string'

# Map plot types to valid file types
VALID_FILE_FORMAT  = { GENE_MODEL_PLOT  : [GENE_MODEL_FORMAT],
                       ISOFORM_PLOT     : [SPLICEGRAPH_FORMAT],
                       JUNCTION_PLOT    : [SAM_FORMAT],
                       READ_DEPTH_PLOT  : [SAM_FORMAT],
                       SPLICEGRAPH_PLOT : [SPLICEGRAPH_FORMAT],
                       XY_PLOT          : [CSV_FORMAT]
                       }

# Recognized file extensions:
VALID_FILE_EXTS    = {CSV_FORMAT         : ['.csv'],
                      GENE_MODEL_FORMAT  : ['.gff', '.gff3', '.gtf'],
                      SAM_FORMAT         : ['.sam'],
                      SPLICEGRAPH_FORMAT : ['.gff']
                      }
# Hack to add gzipped files:
for k in VALID_FILE_EXTS :
    gzList   = ['%s.gz'%s for s in VALID_FILE_EXTS[k]]
    gzipList = ['%s.gzip'%s for s in VALID_FILE_EXTS[k]]
    VALID_FILE_EXTS[k] += gzList + gzipList

def parseEnvironmentVars(value) :
    """Users may use environment variables in the config file.  These
    must be surrounded by ${} as in many shells."""
    envList = ENV_MATCH.findall(value)
    result  = str(value)
    for env in envList :
        name = env[2:-1] # omit '${' and '}'
        try :
            fullValue = os.environ[name]
        except KeyError :
            raise ValueError('Unrecognized environment variable %s found in configuration file.' % env)
        result = result.replace(env,fullValue)
    return result

class SingleConfig(object) :
    """Encapsulates the meta-information provided for a set of plots."""
    def __init__(self) :
        self.legend         = False
        self.height         = None
        self.width          = None
        self.fontsize       = None
        self.log_file       = None
        self.output_file    = None
        self.shrink_introns = False
        self.shrink_factor  = MIN_INTRON_SIZE

    def __setattr__(self, k, v) :
        if k not in VALID_SINGLE_TAGS :
            raise ValueError("Unrecognized %s tag %s" % (SINGLE_SECTION,k))
        self.__dict__[k] = v

    def __eq__(self, o) :
        return self.__str__() == o.__str__()

    def __hash__(self) :
        return self.__str__().__hash__()

    def __str__(self) :
        keys = sorted(self.__dict__.keys())
        return ','.join(['%s=%s'%(x,str(self.__dict__[x])) for x in keys])

class MultiConfig(SingleConfig) :
    """Encapsulates the meta-information provided for a set of multi-plots."""
    def __init__(self) :
        SingleConfig.__init__(self)
        self.chromosome    = None
        self.gene_list     = None
        self.output_format = None

    def __setattr__(self, k, v) :
        if k not in VALID_MULTI_TAGS :
            raise ValueError("Unrecognized %s tag %s" % (MULTI_SECTION,k))
        self.__dict__[k] = v

class PlotInfo(object) :
    """Encapsulates the information required for a given subplot type."""
    def __init__(self, name) :
        self.name           = name
        self.acceptors      = ''
        self.annotate       = False
        self.background     = True
        self.clusters       = False
        self.codons         = False
        self.donors         = ''
        self.edge_weight    = 1
        self.features       = None 
        self.file_format    = None 
        self.gene_name      = None
        self.hide           = False
        self.highlight      = ''
        self.iso_labels     = False
        self.iso_weights    = ''
        self.labels         = False
        self.line           = False
        self.log_threshold  = 100.0
        self.min_coverage   = 1
        self.plot_type      = None
        self.relative_size  = 1.0
        self.novel_jct      = True
        self.source_file    = None
        self.title_string   = None
        self.x_labels       = False
        self.y_limit        = 0.0

    def __setattr__(self, k, v) :
        if k != 'name' and k not in VALID_PLOT_TAGS :
            raise ValueError("Unrecognized %s tag %s" % (self.name,k))
        self.__dict__[k] = v

    def __str__(self) :
        keys = sorted(self.__dict__.keys())
        return ','.join(['%s=%s'%(x,str(self.__dict__[x])) for x in keys])

class SinglePlotConfig(object) :
    """
    Class for reading plot configurations that contain the
    parameters for a meta-plot.
    """
    globalParams = SingleConfig()
    mainSection  = SINGLE_SECTION
    mainTags     = VALID_SINGLE_TAGS
    def __init__(self, cfgFile=None, **args) :
        self.comments = []
        self.verbose  = getAttribute('verbose', False, **args)
        self.config   = ConfigParser()
        if cfgFile :
            self.config.read(cfgFile)
            self.validate()
            self.instantiate()

    def getConfiguration(self) :
        """Returns the global configuration object."""
        return self.globalParams

    def getPlotInfo(self, p) :
        """Returns the given plot information object."""
        return self.plotObjects[p]

    def getPlotIds(self) :
        """Returns a list of all plot sections found in the config file."""
        return [s for s in self.config.sections() if s != self.mainSection]

    def getPlotList(self) :
        """Returns a list of all plot objects found in the config file."""
        return [self.plotObjects[p] for p in self.getPlotIds()]

    def getValue(self, section, name, default=None) :
        """Returns a value from the configuration."""
        tag = name.lower()
        if (section == self.mainSection) and (name not in self.mainTags) :
            raise ValueError("Unrecognized %s tag %s" % (section,name))
        if (section != self.mainSection) and (name not in VALID_PLOT_TAGS) :
            raise ValueError("Unrecognized %s tag %s" % (section,name))

        try :
            # ConfigParser type-based methods can fail in nasty ways,
            # so instead we use a single access point and recast the
            # values by brute force.
            value = self.config.get(section,name)
        except Exception :
            # When ConfigParser fails, it usually throws a TypeError
            # from 'way down in its bowels, but we catch them all.
            return default

        try :
            if tag in FLOAT_TAGS :
                return float(value)
            elif tag in INTEGER_TAGS :
                return int(value)
            elif tag in BOOLEAN_TAGS :
                return value.lower() in ['t','true','0']
            elif tag == SOURCE_TAG :
                return parseEnvironmentVars(value)
            else : # STRING_TAGS
                return value
        except ValueError, ve :
            raise ValueError('Invalid assignment in config file: %s is %s, should be %s\n' % (tag, value, TAG_TYPE[tag]))

    def instantiate(self) :
        """Instantiates all objects related to the configuration."""
        for o in self.config.options(self.mainSection) :
            self.globalParams.__dict__[o] = self.getValue(self.mainSection,o)

        self.plotObjects = {}
        for p in self.getPlotIds() :
            plot = PlotInfo(p)
            for o in self.config.options(p) :
                plot.__dict__[o] = self.getValue(p,o)
            if plot.file_format is None and plot.plot_type != GENE_MODEL_PLOT :
                plot.file_format = VALID_FILE_FORMAT[plot.plot_type][0]
            self.plotObjects[p] = plot

    def validate(self) :
        """Enforces strict adherence to plot format and value types."""
        sections = self.config.sections()
        if self.mainSection not in sections :
            raise ValueError('Plot configuration file is missing [%s] section' % self.mainSection)

        configSet  = set(self.config.options(self.mainSection))
        badOptions = configSet - set(self.mainTags)
        if badOptions :
            raise ValueError('Unrecognized options found in [%s] section: %s\nValid options are: %s' % (self.mainSection, ', '.join(badOptions), ', '.join(self.mainTags)))

        # Can we determine plot boundaries from the plot types provided?
        okMinMax = False
        validSet = set(VALID_PLOT_TAGS)
        for plot in self.getPlotIds() :
            configSet  = set(self.config.options(plot))

            # Validate options generally
            badOptions = configSet - validSet
            if badOptions :
                raise ValueError('Unrecognized options found in %s plot: %s\nValid options are: %s' % (plot, ', '.join(badOptions), ', '.join(validSet)))

            # Validate plot type
            plotType = self.getValue(plot,PLOT_TAG)
            if plotType is None :
                raise ValueError('Missing plot type for %s plot' % plot)
            if plotType not in VALID_PLOT_TYPES :
                raise ValueError('Unrecognized plot type found in %s plot: %s\nValid types are: %s' % (plot, plotType, ', '.join(VALID_PLOT_TYPES)))

            okMinMax |= (plotType in MINMAX_PLOT_TYPES)

            # Validate source file format (must be specified for gene model graphs)
            sourceFormat = self.getValue(plot,FORMAT_TAG)
            if sourceFormat is None and plotType == GENE_MODEL_PLOT :
                raise ValueError('Source file format required for %s plot' % plot)
            elif sourceFormat is None :
                sourceFormat = VALID_FILE_FORMAT[plotType][0]

            if sourceFormat not in VALID_FILE_FORMAT[plotType] :
                raise ValueError("Unrecognized %s file format: %s\nValid formats are: %s" % (plotType,sourceFormat,', '.join(VALID_FILE_FORMAT[plotType])))

            # Validate source file
            sourceFile  = self.getValue(plot,SOURCE_TAG)
            if sourceFile is None :
                raise ValueError('Missing source data for %s plot' % plot)

        if not okMinMax :
            raise ValueError("Insufficient data to determine plot boundaries; add one of the following plot types: %s" % ', '.join(MINMAX_PLOT_TYPES))

class MultiplotConfig(SinglePlotConfig) :
    """
    Class for reading plot configurations that contain the parameters for a multi-plot.
    """
    globalParams = MultiConfig()
    mainSection  = MULTI_SECTION
    mainTags     = VALID_MULTI_TAGS

    def validate(self) :
        """Enforces parent class constraints along with specific directory
        constraints for multi-plot graphs."""
        SinglePlotConfig.validate(self)

        # Check output file format
        outputExt = self.getValue(self.mainSection,EXT_TAG)
        if not outputExt :
            raise ValueError('No output file format specified')
        if outputExt.lower() not in OUTPUT_FORMATS :
            raise ValueError("Output file format '%s' not recognized; must be one of %s" % (outputExt, ', '.join(OUTPUT_FORMATS)))

        for plot in self.getPlotIds() :
            plotType     = self.getValue(plot,PLOT_TAG)
            sourceFormat = self.getValue(plot,FORMAT_TAG)
            sourceFile   = self.getValue(plot,SOURCE_TAG)
            # Splice graphs must use directories; gene models and read depths must use single files
            sourceFormat = self.getValue(plot,FORMAT_TAG)
            if sourceFormat is None :
                sourceFormat = VALID_FILE_FORMAT[plotType][0]

            if sourceFormat in MULTI_FILE_FORMATS :
                if os.path.isfile(sourceFile) :
                    raise ValueError('%s source for %s is a file, %s format requires a directory.' % (plotType,plot,sourceFormat))
            elif os.path.isdir(sourceFile) :
                    raise ValueError('%s source for %s is a directory, not a file.' % (plotType,plot))
