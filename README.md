plot_LST
========

# Basic use

    python plot_LST_vXX.py id f l

where

- f is int, first line in .lst file
- l is int, last line in .lst file
- id being the name of a replicon in the gembase format : ggee.000.c00.00

If an id.config file exist, it will color the genes as many colors as there are classes in the .config file which is formatted as follow :

    GENE001c01 classe1
    GENE002c01 classe2
    GENE003c01 classe1
  
  
# with annotations 
    
   python plot_LST_vXX.py id f l -a

# circular representation

   python plot_LST_vXX.py id f l -c

# On a tree
   python plot_LST_vXX.py id_tree -t

with -t option, an id_tree.treeconfig file is needed, and is formatted as follow :

    id1 f1 l1
    id2 f2 l2
    id3 f3 l3

where idx are names of leaves in the necessary i_tree.tree file, in Newick format, like :

    ((id1,id2),id3)

id_tree is the name of the tree.

To colorize genes, one .config file is needed where all desired genes from all classes are concatenated in it.
