################################################################################
### TAXONOMIC CLASSIFICATION AND PHYLOGENETIC TREE VISUALIZATION

# Obtain classification for species list, create phylogenetic trees.
# Trees are created using the same features as ggplot2, so 
# are highly adaptable/customizable.


################################################################################
# SETUP

install.packages(c("ggtree","dplyr","ggplot2","tidyr"))
install.packages(c("bold"))
install.packages("remotes")
install.packages("phytools")
devtools::install_github("nishanmudalige/BOLD.R", force = TRUE)
remotes::install_github("ropensci/bold")
remotes::install_github("ropensci/taxize")
remotes::install_github("ropensci/ggtree")


library(dplyr); library(taxize); library(ggtree); library(ggplot2); library(tidyr)

taxon.names = c("Myxicola infundibulum", "Peprilus triacanthus", "Lacuna vincta", "Calanoida",
                "Pseudopleuronectes americanus", "Centropages typicus", "Osmerus mordax",               
                "Fundulus heteroclitus", "Metridia lucens", "Priapulus caudatus",
                "Meganyctiphanes norvegica", "Urophycis tenuis", "Polydora cornuta",
                "Nanomia cara","Calanus finmarchicus","Zygochlamys delicatula",
                "Brevoortia tyrannus", "Bougainvillia muscus", "Scalibregma inflatum",
                "Alosa pseudoharengus", "Acartia tonsa", "Scomber scombrus",
                "Flabellina verrucosa", "Chaetogaster diastrophus", "Hydra vulgaris",
                "Chone infundibuliformis", "Salvelinus fontinalis", "Mytilus edulis",
                "Aurelia aurita", "Ascidia ahodori", "Leptodiaptomus minutus" )



#Gulf 2020 CMs 5000 Class Size
taxon.names = c('Acartia', 'Ascidiacea', 'Bivalvia', 'Brachyura', 'Bryozoa', 'Calanoida', 'Calanus', 'Centropages', 'Cirripedia', 'Copepoda', 'Decapoda', 'Echinodermata', 'Eurytemora', 'Evadne', 'Fritillaria', 'Gastropoda', 'Harpacticoida', 'Hydrozoa', 'Labidocera', 'Obelia', 'Oithona', 'Osteichthyes', 'Paracalanus', 'Podon', 'Polychaeta', 'Pseudocalanus', 'Pseudodiaptomus', 'Temora', 'Tortanus')

#Gulf 2020 CMs Max Class Size
taxon.names = c('Acartia', 'Ascidiacea', 'Bivalvia', 'Brachyura', 'Bryozoa', 'Calanoida', 'Calanus', 'Centropages', 'Cirripedia', 'Copepoda', 'Decapoda', 'Echinodermata', 'Eurytemora', 'Evadne', 'Fritillaria', 'Gastropoda', 'Harpacticoida', 'Hydrozoa', 'Labidocera', 'Obelia', 'Oithona', 'Osteichthyes', 'Paracalanus', 'Podon', 'Polychaeta', 'Pseudocalanus', 'Pseudodiaptomus', 'Temora')

# PA 2021 CMs
taxon.names = c('Polychaeta', 'Cirripedia', 'Evadne', 'Podon', 'Acartia', 'Centropages', 'Corycaeidae', 'Calanoida', 'Paracalanus', 'Pseudocalanus', 'Tortanus', 'Copepoda', 'Oithona', 'Brachyura', 'Decapoda', 'Bryozoa', 'Oikopleura', 'Ascidiacea', 'Calycophorae', 'Hydrozoa', 'Echinodermata', 'Bivalvia', 'Gastropoda', 'Platyhelminthes', 'Fritillaria')

# NL 2020 CMs
taxon.names = c('Acartia', 'Bivalvia', 'Bryozoa', 'Calanoida', 'Calanus', 'Centropages', 'Chaetognatha', 'Cirripedia', 'Copepoda', 'Echinodermata', 'Euphausia', 'Eurytemora', 'Evadne', 'Fritillaria', 'Gastropoda', 'Harpacticoida', 'Microcalanus', 'Obelia', 'Oikopleura', 'Oithona', 'Podon', 'Polychaeta', 'Pseudocalanus', 'Temora')

# NL 2021 CMs
taxon.names = c('Acartia', 'Bivalvia', 'Bryozoa', 'Calanoida', 'Calanus', 'Centropages', 'Chaetognatha', 'Cirripedia', 'Copepoda', 'Echinodermata', 'Euphausiacea', 'Eurytemora', 'Evadne', 'Fritillaria', 'Gastropoda', 'Harpacticoida', 'Microcalanus', 'Obelia', 'Oikopleura', 'Oithona', 'Podon', 'Polychaeta', 'Pseudocalanus', 'Temora')




# Gulf Ind Bar Graphs
taxon.names = c('Acartia', 'Ascidiacea', 'Bivalvia', 'Brachyura', 'Bryozoa', 'Calanoida', 'Calanus', 'Centropages', 'Chaetognatha', 'Cirripedia', 'Cnidaria', 'Copepoda', 'Decapoda', 'Echinodermata', 'Eurytemora', 'Evadne', 'Fritillaria', 'Gastropoda', 'Harpacticoida', 'Hydrozoa', 'Labidocera', 'Microsetella', 'Monstrillidae', 'Obelia', 'Oikopleura', 'Oithona', 'Osteichthyes', 'Ostracoda', 'Paracalanus', 'Podon', 'Polychaeta', 'Pseudocalanus', 'Pseudodiaptomus', 'Temora', 'Tortanus')

# Pacific Ind Bar Graphs
taxon.names = c('Acartia', 'Alciopidae', 'Amphipoda', 'Ascidiacea', 'Bivalvia', 'Brachyura', 'Bryozoa', 'Calanoida', 'Calanus', 'Calycophorae', 'Centropages', 'Chaetognatha', 'Cirripedia', 'Copepoda', 'Corycaeidae', 'Ctenophora', 'Cubozoa', 'Decapoda', 'Echinodermata', 'Epilabidocera', 'Euphysa', 'Evadne', 'Fritillaria', 'Gastropoda', 'Harpacticoida', 'Hydrozoa', 'Isopoda', 'Leuckartiara', 'Nemertea', 'Neoturris', 'Obelia', 'Oikopleura', 'Oithona', 'Osteichthyes', 'Ostracoda', 'Paracalanus', 'Platyhelminthes', 'Podon', 'Polychaeta', 'Pseudocalanus', 'Sarsia', 'Scyphozoa', 'Siphonophorae', 'Tentaculata', 'Tortanus')

# NL 2020 Ind Bar Graphs
taxon.names = c('Acartia', 'Aglantha', 'Amphipoda', 'Anthoathecata', 'Bivalvia', 'Bryozoa', 'Calanoida', 'Calanus', 'Centropages', 'Chaetognatha', 'Chiridius', 'Cirripedia', 'Cnidaria', 'Copepoda', 'Decapoda', 'Echinodermata', 'Euphausia', 'Euphausiacea', 'Eurytemora', 'Evadne', 'Fritillaria', 'Gastropoda', 'Harpacticoida', 'Hydrozoa', 'Metridia', 'Microcalanus', 'Obelia', 'Oikopleura', 'Oithona', 'Podon', 'Polychaeta', 'Pseudocalanus', 'Temora', 'Tomopteris')

# NL 2021 Ind Bar Graphs
taxon.names = c('Acartia', 'Aglantha', 'Bivalvia', 'Bryozoa', 'Calanoida', 'Calanus', 'Centropages', 'Chaetognatha', 'Cirripedia', 'Copepoda', 'Decapoda', 'Echinodermata', 'Euphausia', 'Euphausiacea', 'Eurytemora', 'Evadne', 'Foraminifera', 'Fritillaria', 'Gastropoda', 'Harpacticoida', 'Hydrozoa', 'Microcalanus', 'Monstrillidae', 'Obelia', 'Oikopleura', 'Oithona', 'Ostracoda', 'Podon', 'Polychaeta', 'Pseudocalanus', 'Sarsia', 'Temora', 'Tortanus')


################################################################################
# CLASSIFICATION 
print(taxon.names)
# taxon.names <- unique(taxon.names)
# Finding the WoRMS ID first will speed up classification 
taxa.id = taxize::get_wormsid(taxon.names, rows=1, marine_only = FALSE)
# Creates a list of ranks for each species
highertaxa = taxize::classification(taxa.id[!is.na(taxa.id)], db="worms") %>%
  lapply(as.data.frame) # Each list element must be data.frame for tree building later

new_taxa = dplyr::bind_rows(
  lapply(highertaxa, function(x) 
    select(x, -id) %>%
    pivot_wider(
      names_from = "rank",
      values_from = "name")), 
  .id = "id") %>%
  mutate(taxonID = paste0("urn:lsid:marinespecies.org:taxname:",id)) %>%
  select(c("Kingdom","Phylum","Class","Order","Family","Genus","taxonID"))

print(highertaxa)
################################################################################
# TREE BUILDING


# Build phylogenetic tree with ggtree
taxize_tree <- class2tree(highertaxa, check = TRUE)



# If needed, extract the tree tip labels in the order they are listed on the tree 
# for use in other plots
is_tip <- taxize_tree$phylo$edge[,2] <= length(taxize_tree$phylo$tip.label)
ordered_tips <- taxize_tree$phylo$edge[is_tip, 2]
spp.ord <- rev(taxize_tree$phylo$tip.label[ordered_tips]) 
print(spp.ord)

ordered_leaf_tips <- taxize_tree$phylo$tip.label

# Print the ordered leaf tips
print(ordered_leaf_tips)

taxize_tree$edge.length[1] <- 200  # Change the length of the first branch
taxize_tree$edge.length[2] <- 0.5  # Change the length of the second branch

# Create tree
ggtree(taxize_tree$phylo) +
  geom_tiplab() +
  xlim(0,200)  # the max value may change, so play around with it 
               # if the tip labels get cut off 

print(taxize_tree)
ggtree(taxize_tree$phylo, branch.length = "none") # depending on how you want the visualization

taxize_tree$phylo$edge.length <- taxize_tree$phylo$edge.length * 0.5  # Shrink branch lengths by half
ggtree(taxize_tree$phylo, branch.length = "none") +
  #geom_nodelab(geom = "label") + # if you want higher classification labels
  #geom_tiplab() + # if you want species labels
  xlim(-.5,10)
  ylim(0, 40)
  
  
ggtree(taxize_tree$phylo, branch.length = "none") +
  geom_nodelab(geom = "label", size = 3) +  # Adjust size of node labels
  geom_tiplab(size = 3) +  # Adjust size of tip labels
  xlim(-.5, 12) +
  ylim(0, 20)
  
tree_dict <- as.list(taxize_tree)
print(tree_dict)

library(ggtree)
library(ape)

# Assuming `tree` is your ggtree object (you can convert it to an ape tree if needed)
# If your tree is already an ape object, you can skip the next step
tree_ape <- as.phylo(taxize_tree)

# Export the tree to Newick format
newick_tree <- write.tree(tree_ape)
cat(newick_tree, file = "tree.nwk")



library(phytools)

# Export to JSON format (assuming your tree is in `tree_ape` form)
json_tree <- phyloToJSON(taxize_tree)
writeLines(json_tree, "tree.json")


# Convert `classtree` to a nested list (for easy JSON conversion)
tree_list <- function(tree) {
  if (is.null(tree)) return(NULL)
  
  # Recursively extract the tree structure into a list
  lapply(tree, function(x) {
    list(
      name = x$name,
      rank = x$rank,
      children = tree_list(x$children)
    )
  })
}

# Apply the function to your tree
nested_tree <- tree_list(taxize_tree)

# Convert the nested list to a JSON object
tree_json <- toJSON(nested_tree, pretty = TRUE)

# Save the JSON output to a file
write(tree_json, file = "tree.json")
