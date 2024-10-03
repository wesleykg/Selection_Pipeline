'''Usage: 0b_prune-tree.py <tree> <taxa_to_keep> <clade>'''

from ete3 import Tree  # Manipulating trees
import os  # Manipulating filenames

# Check if running interactively in an iPython console, or in a script from the
# command line
def in_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False
# Run in a script from the command line
if in_ipython() is False:
    from docopt import docopt  # Command line argument handler
    cmdln_args = docopt(__doc__)
    tree_file = cmdln_args.get('<tree>')
    taxa_to_keep_file = cmdln_args.get('<taxa_to_keep>')
    clade = cmdln_args.get('<clade>')
# Run interactively in an iPython console
if in_ipython() is True:
    tree_file = 'RAxML_bestTree.DNA_part_Jan2018.tre'
    taxa_to_keep_file = 'taxa_burmPartials.txt'
    clade = 'burmPartials'

# Retrieve name of the original tree and append with the clade being pruned
# down to
tree_name = os.path.splitext(tree_file)[0]
out_tree_filename = tree_name + '_' + clade + '.tre'

# Read in by ete3
full_tree = Tree(tree_file)

# Retrieve the list of taxa to be pruned down to
taxa_to_keep = []
with open(taxa_to_keep_file, 'r') as taxa_list:
    for taxon in taxa_list:
        taxon = taxon.rstrip('\n')    
        taxa_to_keep.append(taxon)

# Prune the full tree down to just the taxa listed
full_tree.prune(taxa_to_keep, preserve_branch_length=True)

# PAML requires an unrooted tree
full_tree.unroot()

# Write the pruned tree to a new file with the above filename
full_tree.write(outfile=out_tree_filename)
