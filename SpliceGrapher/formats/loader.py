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
from SpliceGrapher.formats.GTFLoader import *
from SpliceGrapher.formats.GeneModel import *

import os

def loadGeneModels(path, **args) :
    """Loads a set of gene models in either GTF or GFF3 format.
    Uses a simple heuristic that looks for '.gtf' anywhere in the
    filename (case insensitive)."""
    base = os.path.basename(path)
    base = base.lower()
    if base.find('.gtf') >= 0 :
        return GTFLoader(path, **args).model
    else :
        return GeneModel(path, **args)
