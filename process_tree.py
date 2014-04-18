#!/opt/local/bin/python2.6

'''Code for Newick tree analysis. See USAGE below for details.

This code uses the Python Environment for (phylogenetic) Tree
Exploration (ETE), published in:

Jaime Huerta-Cepas, Joaquin Dopazo and Toni Gabaldon. ETE: a python
Environment for Tree Exploration. BMC Bioinformatics 2010, 11:24.

Documentation is at http://packages.python.org/ete2.
'''

import argparse
import sys
from ete2 import Tree

USAGE = '''
Analyzes a Newick tree to find all clusters with the
following properties:
* support > SUPPORT
* mean inter-leaf genetic distance < DISTANCE

Use -f FILE to read the tree from FILE.
Use -r N to generate a random tree with N random leaves.

Writes results to five files, <OUTPUT>.tree, <OUTPUT>.clusters,
<OUTPUT>.sequences, <OUTPUT>.dist, and <OUTPUT>.distt.

<OUTPUT>.tree contains an ascii diagram of the Newick tree with all
internal nodes labeled. The labels correspond to cluster IDs in the
other two files.

<OUTPUT>.clusters contains one line per qualifying cluster:
* <cluster ID> <support> <mean genetic distance>

<OUTPUT>.sequences contains one line per sequence (for sequences
within qualifying clusters):
* <cluster ID> <sequence name>

<OUTPUT>.dist and .distt contain inter-leaf distances. The .dist file
contains inter-leaf distances for all non-leaf nodes. The .distt file
contains inter-leaf distances only for those non-leaf nodes whose
support value exceeds the user-supplied SUPPORT threshold.
'''

VERBOSE = 0

# Newick tree node feature names.
MEAN_LEAF_DIST = 'mean_leaf_dist'
NUM_LEAVES = 'num_leaves'
TOT_DIST_TO_LEAVES = 'tot_dist_to_leaves'

def Log(level, msg):
  '''Prints message if verbosity level is appropriate.'''
  if VERBOSE >= level:
    print msg


class VerboseTree(Tree):        # pylint: disable=R0904
  '''A version of ete2.Tree that prints the support and inter-leaf
  distance of every node.
  '''

  # pylint: disable=C0103
  # pylint: disable=R0914
  def _asciiArt(self, char1='-', show_internal=True, compact=False):
    '''Returns the ASCII representation of the tree. Adapted from
    ete2.Tree.
    '''
    LEN = 5
    PAD = ' ' * LEN
    PA = ' ' * (LEN-1)
    descr = getattr(self, 'name')
    if show_internal:
      descr += (':S%0.2f' % self.support)
      descr += (':D%0.2f' % self.dist)
      if hasattr(self, MEAN_LEAF_DIST):
        descr += (':M%0.2f' % getattr(self, MEAN_LEAF_DIST))
    if not self.is_leaf():
      mids = []
      result = []
      for c in self.children:
        if c is self.children[0]:
          char2 = '/'
        elif c is self.children[-1]:
          char2 = '\\'
        else:
          char2 = '-'
        # pylint: disable=W0212
        (clines, mid) = c._asciiArt(char2, show_internal, compact)
        mids.append(mid+len(result))
        result.extend(clines)
        if not compact:
          result.append('')
      if not compact:
        result.pop()
      (lo, hi, end) = (mids[0], mids[-1], len(result))
      prefixes = [PAD] * (lo+1) + [PA+'|'] * (hi-lo-1) + [PAD] * (end-hi)
      mid = (lo + hi) / 2
      prefixes[mid] = char1 + '-'*(LEN-2) + prefixes[mid][-1]
      result = [p+l for (p, l) in zip(prefixes, result)]
      if show_internal:
        stem = result[mid]
        result[mid] = stem[0] + descr + stem[len(descr)+1:]
      return (result, mid)
    else:
      return ([char1 + '-' + descr], 0)


def Fatal(msg):
  '''Prints msg and exits.'''
  print msg
  sys.exit(1)


def LabelInternalNodes(tree):
  '''Assigns sequential numeric ids to unlabeled internal nodes of a
  tree.
  '''
  node_num = 0
  for node in tree.traverse(strategy='levelorder'):
    if node.name == 'NoName':
      node.name = str(node_num)
      node_num += 1    


def PrintTree(tree, outfile):
  '''Prints <tree> to <outfile>.tree.'''
  fname = outfile + '.tree'
  try:
    f = open(fname, 'w')
    f.write(tree.get_ascii(show_internal=True, compact=False))
    f.write('\n')
    f.close()
  except IOError:
    Fatal('Could not write to ' + fname)
    

def PrintDistances(tree, support_thresh, outfile):
  '''Prints mean inter-leaf distances to two files, <outfile>.dist and
  <outfile>.distt. The former contains inter-leaf distances of all
  non-leaf nodes; the latter contains inter-leaf distances only for
  nodes that meet the threshold criterion given by support_thresh.

  Requires that tree nodes have been annotated with mean_leaf_dist
  attributes.
  '''
  def printDistInternal(node, dist, distt):
    if node.is_leaf():
      return
    if not hasattr(node, MEAN_LEAF_DIST):
      return
    dist.write(str(node.mean_leaf_dist) + '\n')
    if node.support > support_thresh:
      distt.write(str(node.mean_leaf_dist) + '\n')
    for c in node.get_children():
      printDistInternal(c, dist, distt)

  fname1, fname2 = outfile + '.dist', outfile + '.distt'
  try:
    f1, f2 = open(fname1, 'w'), open(fname2, 'w')
    printDistInternal(tree, f1, f2)
    f1.close()
    f2.close()
  except IOError:
    Fatal('Could not write to ' + fname1 + ' or ' + fname2)


def PrintClusters(nodes, outfile):
  '''Given nodes, a list of nodes (cluster roots), prints
  <outfile>.clusters as specified in USAGE.
  '''
  fname = outfile + '.clusters'
  try:
    f = open(fname, 'w')
    for node in sorted(nodes, key=lambda node: int(node.name)):
      f.write(node.name + ' ' + str(node.support) + ' ' +
              str(node.mean_leaf_dist) + '\n')
    f.close()
  except IOError:
    Fatal('Could not write to ' + fname)


def PrintSequences(nodes, outfile):
  '''Given nodes, a list of nodes (cluster roots), prints
  <outfile>.clusters as specified in USAGE.
  '''
  fname = outfile + '.sequences'
  try:
    f = open(fname, 'w')
    for node in sorted(nodes, key=lambda node: int(node.name)):
      for leaf in node.iter_leaves():
        f.write(node.name + ' ' + leaf.name + '\n')
    f.close()
  except IOError:
    Fatal('Could not write to ' + fname)


LEAF_DIST_MEMO = {}
def ComputeMeanInterleafDistance(node):
  '''Computes the mean inter-leaf distance of node. Sets the
  'mean_leaf_dist' attribute to the distance.

  Returns the distance.
  '''
  # Memoize pairwise leaf distance in the LEAF_DIST_MEMO dict to
  # improve performance.
  tot_dist = 0
  leaves = sorted(node.get_leaves(), key=lambda node: node.name)
  for i in xrange(0, len(leaves)-1):
    for j in xrange(i+1, len(leaves)):
      try:
        dist_from_i = LEAF_DIST_MEMO[leaves[i]]
        Log(2, "Found memo row")
      except KeyError:
        dist_from_i = LEAF_DIST_MEMO[leaves[i]] = dict()
      try:
        dist_from_i_to_j = dist_from_i[leaves[j]]
        Log(2, "Found memo item")
      except KeyError:
        dist_from_i_to_j = dist_from_i[leaves[j]] = (
          leaves[i].get_distance(leaves[j]))
        Log(2, "Computing from scratch")
      tot_dist += dist_from_i_to_j
  node.add_feature(MEAN_LEAF_DIST,
                   (tot_dist / (len(leaves) * (len(leaves) - 1) / 2)))
  return node.mean_leaf_dist

def NodesWithSupportAndDistance(node, support_thresh, distance_thresh):
  '''If node's support exceeds support_thresh, and if node's mean
  inter-leaf distance is below distance_thresh, returns a list
  containing only node.

  Otherwise, runs recursively on its children and returns a list of
  the closest (least-deep) descendants that match both criteria.
  '''
  Log(2, "Computing NodesWithSupportAndDistance")
  if node.support > support_thresh and not node.is_leaf():
    if hasattr(node, MEAN_LEAF_DIST):
      d = getattr(node, MEAN_LEAF_DIST)
    else:
      d = ComputeMeanInterleafDistance(node)
    if d < distance_thresh:
      # If a node meets both criteria, we need not examine its children. 
      return [node]
  nodes = []
  for c in node.get_children():
    nodes += NodesWithSupportAndDistance(c, support_thresh, distance_thresh)
  return nodes


def ComputeInterleafDistance(tree):
  '''Annotates each node in tree with the mean inter-leaf distance of the
  cluster rooted at that node.'''

  for node in tree.traverse(strategy='postorder'):
    num_leaves = 1 if node.is_leaf() else 0
    tot_dist_to_leaves = 0
    mean_leaf_dist = 0
    for c in node.children:     # Leaves have no children.
      num_leaves += getattr(c, NUM_LEAVES)
      tot_dist_to_leaves += (c.dist * getattr(c, NUM_LEAVES) +
                             getattr(c, TOT_DIST_TO_LEAVES))
    mean_leaf_dist = (tot_dist_to_leaves / num_leaves) * 2
    node.add_feature(NUM_LEAVES, num_leaves)
    node.add_feature(TOT_DIST_TO_LEAVES, tot_dist_to_leaves)
    node.add_feature(MEAN_LEAF_DIST, mean_leaf_dist)


def ProcessTree(tree, quick, support_thresh, distance_thresh, outfile):
  '''Given tree, find all clusters with the following properties:
  * support > support_thresh
  * mean inter-leaf genetic distance < distance_thresh

  If quick is true, computes mean inter-leaf distances using an
  algorithm linear in the number of leaves. Otherwise, uses the naive
  brute-force quadratic algorithm.

  Writes results to three files, <outfile>.tree, <outfile>.clusters and
  <outfile>.sequences, as specified in USAGE.
  '''
  LabelInternalNodes(tree)
  if quick:
    ComputeInterleafDistance(tree)
  nodes = NodesWithSupportAndDistance(tree, support_thresh, distance_thresh)
  PrintTree(tree, outfile)
  PrintDistances(tree, support_thresh, outfile)
  PrintClusters(nodes, outfile)
  PrintSequences(nodes, outfile)


def main():                     # pylint: disable=C0103
  '''Parses arguments and calls analysis.'''
  parser = argparse.ArgumentParser(
    description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('-d', '--distance', default=0.3, type=float,
                      help='Maximum per-cluster mean inter-leaf distance')
  parser.add_argument('-f', '--file', default='',
                      help='Input file containing Newick tree')
  parser.add_argument('-o', '--output', default='nwk_analysis',
                      help='Basename of output files')
  parser.add_argument('-q', '--quick', default=True,
                      help='Use linear analysis algorithm')
  parser.add_argument('-r', '--random', default=0, type=int,
                      help='Generate a random tree of the given size')
  parser.add_argument('-s', '--support', default=0.9, type=float,
                      help='Minimum support of cluster root')
  parser.add_argument('-v', '--verbose', action='count',
                      help='Show debug info (repeat to increase verbosity)')
  args = parser.parse_args()
  global VERBOSE                # pylint: disable=W0603
  VERBOSE = args.verbose
  if args.file != '':
    tree = VerboseTree(args.file)
  elif args.random > 0:
    tree = VerboseTree(format=2)
    tree.populate(args.random, random_branches=True)
  else:
    Fatal(parser.format_help())
  ProcessTree(tree, args.quick, args.support, args.distance, args.output)


if __name__ == '__main__':
  main()
