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
Establish default values used by all modules.
"""
from SpliceGrapher.shared.utils import configMap, findFile, getEnvironmentValue
import os

def getMapValue(cfgMap, section, name, default=None) :
    if not cfgMap : return default
    try :
        return cfgMap[section][name]
    except KeyError :
        return default

SG_CONFIG     = 'SpliceGrapher.cfg'

# Main SpliceGrapher configuration sections
SG_SECTION    = 'SpliceGrapher'
VIEW_SECTION  = 'ViewerOptions'
TRAIN_ACC_SEC = 'TrainingAccept'
TRAIN_REJ_SEC = 'TrainingReject'
SEQ_MAP_SEC   = 'SequenceMap'

# Main configuration option identifiers
SG_FASTA_REF_NAME  = 'FASTA_REFERENCE'
SG_GENE_MODEL_NAME = 'GENE_MODEL'
SG_TEMP_DIR_NAME   = 'TEMP_DIR'
SG_SS_PRED_NAME    = 'SS_PREDICTIONS'

# External software versions
PSGINFER_VERSION   = '1.1.3'

# Commonly-used configuration values:
SG_FASTA_REF  = None
SG_GENE_MODEL = None
TRAIN_INCLUDE = None
TRAIN_REJECT  = None

PATH       = getEnvironmentValue('PATH')
configFile = findFile(SG_CONFIG, PATH) if PATH else None

# Read configuration file into a map:
cfgMap     = configMap(configFile) if configFile else {}

# Load main SpliceGrapher values:
SG_FASTA_REF  = getMapValue(cfgMap, SG_SECTION, SG_FASTA_REF_NAME)
SG_GENE_MODEL = getMapValue(cfgMap, SG_SECTION, SG_GENE_MODEL_NAME)
SG_TEMP_DIR   = getMapValue(cfgMap, SG_SECTION, SG_TEMP_DIR_NAME)
SG_SS_PRED    = getMapValue(cfgMap, SG_SECTION, SG_SS_PRED_NAME)

if not SG_FASTA_REF  : SG_FASTA_REF  = getEnvironmentValue('SG_FASTA_REF')
if not SG_GENE_MODEL : SG_GENE_MODEL = getEnvironmentValue('SG_GENE_MODEL')

# Load training key/value pairs directly from the configuration file
TRAIN_ACCEPT  = cfgMap[TRAIN_ACC_SEC] if TRAIN_ACC_SEC in cfgMap else {}
TRAIN_REJECT  = cfgMap[TRAIN_REJ_SEC] if TRAIN_REJ_SEC in cfgMap else {}

# Load mapping between GFF3 chromosome ids and FASTA ids:
SG_FASTA_ID = {}
try :
    SG_FASTA_ID = cfgMap[SEQ_MAP_SEC]
except Exception :
    pass
# Map IDs both ways:
SG_GFF3_ID  = dict([(SG_FASTA_ID[k],k) for k in SG_FASTA_ID])
