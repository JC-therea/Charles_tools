from Bio import Phylo
import matplotlib
import matplotlib.pyplot as plt
from io import StringIO
import sys
from Bio.Seq import Seq
import glob

Path_to_files = sys.argv[1]

print("Used species tree:")
print(Path_to_files)

species_tree_file = Path_to_files
sp_tree = Phylo.read(species_tree_file, "newick")
Phylo.draw(sp_tree)
