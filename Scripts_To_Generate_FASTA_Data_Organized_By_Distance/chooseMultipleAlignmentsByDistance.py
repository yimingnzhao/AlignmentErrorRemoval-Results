import treeswift
import sys
import random
import os
import numpy as np



"""
Checks if an input string is a float

Args:
    num (str): the string to check

Return:
    bool: whether the string is a float or not
"""
def isFloat( num ):
    try:
        float(num)
        return True
    except ValueError:
        return False



"""
Generates an alignment file of a subtree

Args:
    leaf_labels (str list): list of labels of alignments
    aln_file(str): the path to the alignment file
    file_label(str): the label for the output file

Return:
    None
"""
def generateAlignmentFile( leaf_labels, aln_file, file_label ):
    chosen_tree_labels = leaf_labels
    output_f = open("ChosenAlignments/chosen_alignment_" + str(file_label) + ".fasta", "a")
    f = open(aln_file, "r")
    line = f.readline().strip('\n')
    while line:
        if line in chosen_tree_labels:
            output_f.write(line + '\n')
            output_f.write(f.readline())
        else:
            f.readline()
        line = f.readline().strip('\n')
    output_f.close()
    f.close()


USAGE = "python chooseMultipleAlignmentsByDistance.py [newick tree file] [alignment file] [minimum diameter] [maximum diameter]"
DESCRIPTION = "Chooses the top 10 trees (or all trees if the treeset is less than 10) within the diameter range and outputs to a directory 'ChosenAlignments'"
if not len(sys.argv) == 5:
    print()
    print("\tError: Incorrect number of parameters\n")
    print("\tUSAGE: " + USAGE)
    print()
    print("\tDESCRIPTION: " + DESCRIPTION)
    print()
    sys.exit()

if not isFloat(sys.argv[3]):
    print()
    print("\tError: [minimum diameter] is not a float")
    print("\tUSAGE: " + USAGE)
    sys.exit()

if not isFloat(sys.argv[4]):
    print()
    print("\tError: [maximum diameter] is not a float")
    print("\tUSAGE: " + USAGE)
    sys.exit()

newick_file = sys.argv[1]
align_file = sys.argv[2]
min_diameter = float(sys.argv[3])
max_diameter = float(sys.argv[4])


print("Traversing through the tree to find subtrees of the right diameters...", flush=True)
tree = treeswift.read_tree_newick(newick_file)
print("Total Tree Diamter: " + str(tree.diameter()), flush=True)
trees = []
num_nodes_in_tree = []
max_alignments = 0
for node in tree.traverse_postorder(internal=True, leaves=False):
    subtree = tree.extract_subtree(node)
    diameter = subtree.diameter()
    if diameter <= max_diameter and diameter > min_diameter:
            num_nodes_in_tree.append(subtree.num_nodes(internal=False))
            trees.append(subtree)

if len(trees) == 0:
    print("Error: No clades found with the given diameter parameters");
    sys.exit()


# Creates a directory to put all output files
print("Creating directory to store subtrees...", flush=True);
try:
    os.mkdir("ChosenAlignments")
except OSError:
    print("Error: Failed to create a new directory")
    sys.exit()


# Sorts the subtrees by number of leaves
print("Sorting subtrees by leaf counts...", flush=True);
np_trees = np.array(trees)
np_num_nodes_in_tree = np.array(num_nodes_in_tree)
inds = np_num_nodes_in_tree.argsort()
sorted_trees =  np_trees[inds]


# Chooses the 10 trees with the most number of leaves
i = len(sorted_trees) - 1
while i >= 0 and i >= len(sorted_trees) - 10:
    print("Generating alignment file, repitition: " + str(len(sorted_trees)-1-i), flush=True);

    subtree_diameter = sorted_trees[i].diameter();
    max_pairwise_distance = 0;
    mean_pairwise_distance = 0;
    subtree_leaves = 0;
    leaf_set = []
    leaf_set_labels = []
    for leaf in sorted_trees[i].traverse_postorder(internal=False):
        leaf_set.append(leaf);
        leaf_set_labels.append(">" + leaf.get_label());

    for leaf1 in leaf_set:
        for leaf2 in leaf_set:
            if leaf1 == leaf2:
                continue
            distance = sorted_trees[i].distance_between(leaf1, leaf2);
            mean_pairwise_distance += distance;
            max_pairwise_distance = max( max_pairwise_distance, distance )
            subtree_leaves += 1

    subtree_diameter = round(float(subtree_diameter), 6);
    max_pairwise_distance = round(float(max_pairwise_distance), 6);
    mean_pairwise_distance = round(float(float(mean_pairwise_distance)/subtree_leaves), 6);
    file_label = "repitition" + str(len(sorted_trees)-1-i) + "_diameter" + str(subtree_diameter) + "_meanPairDist" + str(mean_pairwise_distance) + "_leaves" + str(len(leaf_set));
    generateAlignmentFile( leaf_set_labels, align_file, file_label )
    i = i - 1

