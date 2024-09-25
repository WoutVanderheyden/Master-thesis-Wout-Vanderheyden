## ANALYZE THE E2 DATA OF FLAVIVIRIDAE ##

#comparePhylo is part of the ape package 
library(ape)

#list.files()
tree_3di = read.tree("/data/leuven/358/vsc35887/master_thesis/iqtree_output/refolded_fullglyco_E2_famsa_iqtree/refolded_fullglyco_E2_3di_famsa.treefile")
tree_AA = read.tree("/data/leuven/358/vsc35887/master_thesis/iqtree_output/refolded_fullglyco_E2_famsa_iqtree/refolded_fullglyco_E2_AA_famsa.treefile")

rf_distance <- dist.topo(tree_3di, tree_AA)
print(rf_distance) #230 splits that differ

result <- comparePhylo(tree_3di, tree_AA)
print(result)

# Both trees have the same number of tips: 189.
# Both trees have the same tip labels.
# Both trees have the same number of nodes: 187.
# Both trees are unrooted.
# Both trees are not ultrametric.
# 71 splits in common.

#############################################################################
####### Now I want and try to compare my tree with the Author's tree ########
#############################################################################

tree_3di_author = read.tree("/vsc-hard-mounts/leuven-data/358/vsc35887/master_thesis/data_flaviviridae/glycoprotein_structural_alignments_and_trees/3di/fn_3di_trees/refolded_fullglyco_E2_3di_famsa.fas.treefile")
tree_AA_author = read.tree("/vsc-hard-mounts/leuven-data/358/vsc35887/master_thesis/data_flaviviridae/glycoprotein_structural_alignments_and_trees/3di/fn_3di_trees/refolded_fullglyco_E2_AA_famsa.fas.treefile")


#####  COMPARE THE 3DI TREES ######
rf_distance <- dist.topo(tree_3di, tree_3di_author)
print(rf_distance) #114 splits in which the trees differ

result <- comparePhylo(tree_3di, tree_3di_author)
print(result)

# Both trees have the same number of tips: 189.
# Both trees have the same tip labels.
# Both trees have the same number of nodes: 187.
# Both trees are unrooted.
# Both trees are not ultrametric.
# 129 splits in common.



###### COMPARE THE AA TREES #######
rf_distance <- dist.topo(tree_AA, tree_AA_author)
print(rf_distance) #50 splits in which the trees differ

result <- comparePhylo(tree_AA, tree_AA_author)
print(result)

# Both trees have the same number of tips: 189.
# Both trees have the same tip labels.
# Both trees have the same number of nodes: 187.
# Both trees are unrooted.
# Both trees are not ultrametric.
# 161 splits in common.