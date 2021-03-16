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
from SpliceGrapher.shared.utils import getAttribute, commaFormat

class PathError(Exception):
    pass

class SpliceGraphPath(object):
    """Encapsulates a single path through a splice graph."""
    def __init__(self, node):
        self.nodes   = [node]
        self.lengths = len(node)
        self.minpos  = node.minpos
        self.maxpos  = node.maxpos

    def append(self, n):
        """Append a node to a path."""
        if n in self.nodes:
            raise PathError('%s is already in path %s\n' % (n.id, self))

        self.nodes.append(n)
        self.lengths += len(n)
        self.minpos   = min(self.minpos, n.minpos)
        self.maxpos   = max(self.maxpos, n.maxpos)

    def __cmp__(self, o):
        """Sorting method that prioritizes paths by the number of nodes
        followed by longest exons and total path span."""
        if len(self) != len(o):
            return len(self) - len(o)

        # Second priority: path with longest exons
        if self.lengths != o.lengths:
            return self.lengths - o.lengths

        # Finally try path with longest span
        return self.span() - o.span()

    def __hash__(self):
        return self.__str__().__hash__()

    def __lt__(self, o):
        if self.minpos == o.minpos:
            return self.maxpos < o.maxpos
        return self.minpos < o.maxpos

    def __len__(self):
        """Returns the length given as the number of nodes."""
        return len(self.nodes)

    def span(self):
        """Returns the total number of nucleotides spanned by the path."""
        return self.maxpos-self.minpos+1

    def update(self, o):
        """Update path by appending nodes from another path."""
        self.nodes.extend(o.nodes)
        self.minpos = min(self.minpos, o.minpos)
        self.maxpos = max(self.maxpos, o.maxpos)

    def __str__(self):
        """String representation of a Path instance."""
        return ','.join([n.id for n in self.nodes])

def getAllPaths(G, **args):
    """Returns a list of all paths through the graph, given as lists of node pointers."""
    # We want only resolved roots plus rare cases where a resolved
    # node has an unresolved parent which we also accept as a root:
    verbose = getAttribute('verbose', False, **args)
    limit   = getAttribute('limit', None, **args)

    roots = set([r for r in G.getRoots() if not r.isUnresolved()])
    for n in list(G.nodeDict.values()):
        for p in n.parents:
            if p.isRoot():
                roots.add(p)

    roots = sorted(roots)
    # Memoization improves efficiency
    memo  = {}
    paths = []
    if verbose:
        print("Traversing %d roots in %s:" % (len(roots), G.getName()))
    for r in roots :
        if verbose:
            print("  ", str(r))
        newpaths = getPaths(r, memo, **args)
        paths.extend(newpaths)
        if limit and len(paths) > limit:
            raise ValueError('Too many paths for %s: %s > %s\n' % (G.getName(), commaFormat(len(paths)), commaFormat(limit)))
    return paths

def getPaths(node, memo, **args):
    """Recursive method that returns a list of all paths that start from the given node."""
    limit = getAttribute('limit', None, **args)
    if node in memo:
        return memo[node]
    elif node.isLeaf():
        result = [SpliceGraphPath(node)]
    else:
        result = []
        for c in node.children:
            for p in getPaths(c, memo, **args):
                path = SpliceGraphPath(node)
                path.update(p)
                result.append(path)
        # If children are unresolved, treat node as a leaf:
        if not result:
            result = [SpliceGraphPath(node)]
    memo[node] = result
    if limit and len(result) > limit:
        raise ValueError('Too many paths: %s > %s\n' % (commaFormat(len(result)), commaFormat(limit)))
    return result

def getLongestPath(G, verbose=False):
    """Returns the longest path through the graph, given as a list of node pointers."""
    roots = set([r for r in G.getRoots() if not r.isUnresolved()])
    for n in list(G.nodeDict.values()):
        for p in n.parents:
            if p.isRoot():
                roots.add(p)
    memo   = {}
    result = None
    if verbose:
        print("getLongestPath traversing %d roots in %s:" % (len(roots), G.getName()))
    for r in roots :
        path = getLongest(r, memo)
        if not result or (len(path) > len(result)):
            result = path
    return result

def getLongest(node, memo, verbose=False):
    """Recursive method that returns the longest path from the given node."""
    if node in memo:
        if verbose:
            print("getLongest returning memoized path for", str(node))
        return memo[node]
    elif node.isLeaf():
        if verbose:
            print("getLongest returning leaf node for", str(node))
        result = SpliceGraphPath(node)
    else:
        result = None
        for c in node.children:
            path = SpliceGraphPath(node)
            path.update(getLongest(c,memo,verbose=verbose))
            if not result or (len(path) > len(result)):
                result = path
        if not result:
            result = SpliceGraphPath(node)
        if verbose:
            print("getLongest returning path of length", len(result))
    memo[node] = result
    return result
