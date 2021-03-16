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
class FeatureGroup(object) :
    """
    Class that represents a group of overlapping features (exons, nodes, junctions, introns).
    A feature must provide 'id', 'minpos' and 'maxpos' attributes.
    """
    def __init__(self, feature) :
        self.group  = {}
        self.group[feature.id] = feature
        self.minpos = feature.minpos
        self.maxpos = feature.maxpos

    def addFeature(self, feature) :
        self.group[feature.id] = feature
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)

    def __contains__(self, feature) :
        return feature.id in self.group

    def __len__(self) :
        return len(self.group)

    def overlaps(self, feature) :
        return self.minpos < feature.maxpos and feature.minpos < self.maxpos

    def values(self) :
        return self.group.values()

class FeatureBundle(object) :
    """
    Class that holds several groups of overlapping features
    """
    def __init__(self, initialGroup=[]) :
        self.bundle = []
        self.addFeatures(initialGroup)

    def addFeature(self, feature) :
        for group in self.bundle :
            if group.overlaps(feature) :
                group.addFeature(feature)
                return
        # no bundle found; create a new one
        self.bundle.append(FeatureGroup(feature))

    def addFeatures(self, featureList) :
        for f in featureList :
            self.addFeature(f)

    def getMultiGroups(self) :
        return [g for g in self.bundle if len(g) > 1]
