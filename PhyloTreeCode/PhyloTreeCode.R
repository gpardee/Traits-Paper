################################################################
## Creating the phylogentic tree... This is used in the fourth corner test to be able to permute by phylogeny
## Code was originally created by Tina Harrison for use in her GCB paper
## Code modified on 21 June 2019 by GLP for the RMBL dataset

#data organization
library(reshape2)
library(ggplot2)

#packages for permutation
library(permute)

#packages for phylogenies
require(rncl)
require(ape)
require(phytools)
require(Hmisc)
library(devtools)
install_github("liamrevell/phytools")
require(phytools)
##create phylogenetic tree
get.bee.phylo <- function(com, which.tree=1, tip.length=0.005)####
  
  which.tree=1
  tip.length=0.005
  
  require(rncl)
  require(ape)
  require(phytools)
  require(Hmisc)
  library(devtools)
  #install_github("liamrevell/phytools")
  require(phytools)
  #make sure phanghorn and phytools are loaded and vegan is unloaded
  
  tree0 <- read.newick(file = "~/Documents/R_analysis/Fourth corner/bee_tree_Hedtke2013BMC.rtf.phy")
  # look at a tree
  plot(tree0, edge.width = 1)
  
  # make sure tree number (which.tree) is valid, use tree number 1, which was aggined at the top
  #if(!which.tree %in% c(1:length(tree0))) print("invalid tree number") don't need to use this
  
  #assign which tree to 1, when you troubleshoot a function you have to first assign variables that show up in bracket
  which.tree=1
  tree1 <- collapse.singles(tree0[[which.tree]]) 
  
  plot.phylo(tree1, show.tip.label=T, cex=0.2)
  
  # tree1$edge.length
  #GLP comment: need to read in file that has all the species in study system 
  com <- read.csv (file = "~/Documents/R_analysis/Fourth corner/PhyloTreeR_27May2019.csv", header=FALSE)
  # ah. make a species-genus membership matrix (binary)
  yy <-as.data.frame(cbind((com$V1), unlist(lapply((com$V1), first.word))))
  #yy <- as.data.frame(rbind(yy, c("Melitta_sp", "Melitta")))
  yy
  # prune tree to genus.list (phytools drop.tip command)
  ##GLP comment: prune tree allows you to only use the genera you're interested in for your data
  pruned.tree1 <- drop.tip(tree1, tree1$tip.label[-match(unique(yy[,2]), tree1$tip.label)])
  plot.phylo(pruned.tree1, align.tip.label=T) # wow so easy!!!
  #GLP comment: this should create a tree with only the genera you want
  
  pt1 <- root(pruned.tree1, outgroup = "Hylaeus")
  plot.phylo(pt1, align.tip.label=T)
  # note could also use setdiff
  PatristicDistMatrix <- cophenetic(pruned.tree1)
  new.pdm1 <- PatristicDistMatrix[order(rownames(PatristicDistMatrix)),]
  new.pdm2 <- new.pdm1[,order(colnames(new.pdm1))]
  
  genmember <- as.matrix(table(com$V1, yy$V2))
  expand.PDM <- genmember %*% t(genmember %*% new.pdm2)
  # use this for species-level analyses:
  PDM <- as.dist(expand.PDM)
  as.matrix(dist(PDM))
  
  # expand genus tips into polytomies
  source("~/Documents/R_analysis/Fourth corner/genera_to_polytomies.R")
  #run genera_to_polytomies
  tree <- make.composite.with.polytomies(pruned.tree1, yy$V2, com$V1, max.genus.age=tip.length)
  
  # make.composite.with.polytomies
  plot(tree, cex=0.4) # tada!!!
  
  tree
  writeNexus(tree, file="FinalTree")
  library(ade4)
  ################################################################