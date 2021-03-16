#! /usr/bin/python
from SpliceGrapher.shared.config             import *
from SpliceGrapher.shared.adjust             import *
from SpliceGrapher.shared.GeneModelConverter import *
from SpliceGrapher                           import SpliceGraph
from SpliceGrapher.shared.utils              import *

g = SpliceGraph.getFirstGraph('AT2G04700.gff')
