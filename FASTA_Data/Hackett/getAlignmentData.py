import treeswift
import sys
import os


newick_file = sys.argv[1]
tree = treeswift.read_tree_newick(newick_file)
print(tree.diameter())
