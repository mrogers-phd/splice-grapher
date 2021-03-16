#!/usr/bin/env python
#----------------------------------------------------------------------
# Setup script for the CSU SpliceGrapher package
import sys,os
import textwrap
from distutils.core import setup,Extension

MATPLOTLIB_REQUIRED = '2.1.0'
PYML_REQUIRED       = '0.7.14'

#functions for showing errors and warnings (lifted from the pytables setup script)
def _print_admonition(kind, head, body):
    tw = textwrap.TextWrapper(initial_indent='   ', subsequent_indent='   ')
    print ".. %s:: %s" % (kind.upper(), head)
    for line in tw.wrap(body):
        print line

def print_warning(head, body=''):
    _print_admonition('warning', head, body)

def bad_version(required, found) :
    try :
        a = [int(x) for x in required.split('.')]
        b = [int(x) for x in found.split('.')]
        for i in range(len(a)) :
            if i >= len(b) : return True
            if b[i] < a[i] : return True
        return False
    except Exception :
        return True

# Check existing software versions:
try :
    import matplotlib
    if bad_version(MATPLOTLIB_REQUIRED, matplotlib.__version__) :
        sys.stderr.write('** Warning: SpliceGrapher was written using matplotlib %s and may not work with %s\n' % (MATPLOTLIB_REQUIRED,matplotlib.__version__))
except :
    sys.stderr.write("Could not find Matplotlib in your PYTHONPATH.\nSpliceGrapher visualization modules will not work.\n")

try :
    import PyML
    if bad_version(PYML_REQUIRED, PyML.__version__) :
        sys.stderr.write("** Warning: SpliceGrapher uses PyML %s; found %s\n" % (PYML_REQUIRED,PyML.__version__))
except :
    sys.stderr.write("Could not find PyML in your PYTHONPATH.\nSpliceGrapher splice site classification modules will not work.\n")

PACKAGES     = ['SpliceGrapher',
                'SpliceGrapher.formats',
                'SpliceGrapher.plot',
                'SpliceGrapher.predict',
                'SpliceGrapher.shared',
                'SpliceGrapher.statistics',
                'SpliceGrapher.view',
                ]

SCRIPT_PATHS = ['scripts/build_classifiers.py',
                'scripts/classify_sites.py',
                'scripts/convert_models.py',
                'scripts/ests_to_splicegraph.py',
                'scripts/find_splice_forms.py',
                'scripts/fix_unresolved.py',
                'scripts/genePredToBed',
                'scripts/gene_model_to_splicegraph.py',
                'scripts/generate_known_junctions.py',
                'scripts/generate_predicted_junctions.py',
                'scripts/generate_putative_sequences.py',
                'scripts/generate_roc.py',
                'scripts/generate_splice_site_data.py',
                'scripts/genewise_statistics.py',
                'scripts/get_good_pairs.py',
                'scripts/gtf2gff.py',
                'scripts/isolasso_pipeline.py',
                'scripts/isolasso_update_graphs.py',
                'scripts/plotter.py',
                'scripts/predict_graphs.py',
                'scripts/predict_splicegraph.py',
                'scripts/psginfer_pipeline.py',
                'scripts/psginfer_update_graphs.py',
                'scripts/realignment_pipeline.py',
                'scripts/sam_collate.py',
                'scripts/sam_filter.py',
                'scripts/sam_split.py',
                'scripts/sam_to_depths.py',
                'scripts/select_model_parameters.py',
                'scripts/splicegraph_statistics.py',
                'scripts/splice_junction_pipeline.py',
                'scripts/view_splicegraph_multiplot.py',
                'scripts/view_splicegraphs.py',
                ]

setup(name         = 'SpliceGrapher',
      version      = '0.2.7',
      description  = "predicting splice graphs from diverse evidence",
      author       = "Mark F. Rogers",
      url          = "http://combi.cs.colostate.edu/SpliceGrapher",
      author_email = "rogersma@cs.colostate.edu",
      packages     = PACKAGES,
      scripts      = SCRIPT_PATHS
      )
