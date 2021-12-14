################################################################
########## Fourth Corner analysis- 22 June 2019  ###############

# Code orginally created by Tina Harrison for her GCB paper, but modified by G Pardee for the RMBL paper


#First load in required packages

### Load packages  ###

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

################################################################



#Step 2: create phylogenetic tree

get.bee.phylo <- function(com, which.tree=1, tip.length=0.005)####
  
  which.tree=1
tip.length=0.005

tree0 <- read.newick(file = "~/Dropbox/fourth corner paper/4th corner codes/updated analysis/R codes/TinaCodesPhyloTree/bee_tree_Hedtke2013BMC.rtf.phy")
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
#all bees, including bees only ID'ed to genus
com <- read.csv (file = "~/Dropbox/fourth corner paper/4th corner codes/updated analysis/Emergence/PhyloTreeR_27May2019_Emergence.csv", header=FALSE)
#no bombus
#com <- read.csv (file = "~/Dropbox/fourth corner paper/4th corner codes/updated analysis/Emergence/PhyloTreeR_27May2019_EmergenceNoBombus.csv", header=FALSE)
#excluding bees only Id'ed to genus
#com <- read.csv (file = "~/Dropbox/fourth corner paper/4th corner codes/updated analysis/Emergence/PhyloTreeR_27May2019_Emergence_no_sp.csv", header=FALSE)

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
source("~/Dropbox/fourth corner paper/4th corner codes/updated analysis/R codes/TinaCodesPhyloTree/genera_to_polytomies.R")
#run genera_to_polytomies
tree <- make.composite.with.polytomies(pruned.tree1, yy$V2, com$V1, max.genus.age=tip.length)

# make.composite.with.polytomies
plot(tree, cex=0.4) # tada!!!

tree
writeNexus(tree, file="FinalTree")

################################################################
#run fourth corner analysis

#download permute package
library(ade4)
#Com is community, with site/years as rows, species as columns, abundance data

#all bees, including bees only Id'ed to genus
Com=read.csv("~/Dropbox/fourth corner paper/4th corner codes/updated analysis/Emergence/Abundance_matrix_27May2019_Emergence.csv",header=TRUE)
#no bombus
#Com=read.csv("~/Dropbox/fourth corner paper/4th corner codes/updated analysis/Emergence/Abundance_matrix_27May2019_EmergenceNoBombus.csv",header=TRUE)
#excluding bees only ID'ed to genus
#Com=read.csv("~/Dropbox/fourth corner paper/4th corner codes/updated analysis/Emergence/Abundance_matrix_27May2019_Emergence_no_sp.csv",header=TRUE)

#Env is environment, with site/years as rows, environmental variables as columns
#original environmental variables
#Env=read.csv("~/Dropbox/fourth corner paper/4th corner codes/updated analysis/FinalMatrices/Site_matrix_27May2019.csv",header=TRUE)

#Drought index
Env=read.csv("~/Dropbox/fourth corner paper/4th corner codes/updated analysis/FinalMatrices/Site_matrix_Drought_11January2021.csv",header=TRUE)

#Trait is traits, with species as rows, traits as columns. We permuted one trait at a time 
#including all bees, even ones only ID'ed to genus
Trait=read.csv("~/Dropbox/fourth corner paper/4th corner codes/updated analysis/Emergence/TraitsMatrix_27May2019_Emergence.csv",header=TRUE)
#no bombus
#Trait=read.csv("~/Dropbox/fourth corner paper/4th corner codes/updated analysis/Emergence/TraitsMatrix_27May2019_EmergenceNoBombus.csv",header=TRUE)
#excluding bees only Id'ed to genus
#Trait=read.csv("~/Dropbox/fourth corner paper/4th corner codes/updated analysis/Emergence/TraitsMatrix_27May2019_Emergence_no_sp.csv",header=TRUE)

#first, set 1st column of each dataset as the rownames. 
Trait <- data.frame(Trait[,-1], row.names=Trait[,1])
Com <- data.frame(Com[,-1], row.names=Com[,1])
Env <- data.frame(Env[,-1], row.names=Env[,1])


#I manually set the Its to 9999. 

Fourth.corner <- function(Trait, Its = 9999) {
  #this currently 
  Com <- Com[, match(rownames(Trait), colnames(Com))]
  tree
  #data should not have any NAs or else it won't run, but all the traits have been found for each species
  
  # shuffle rows (breaks environment x site links)
  f.row <- fourthcorner(Env, Com, Trait, modeltype = 2, p.adjust.method.G="none", p.adjust.method.D = "none", nrepet=9999)
  # shuffle columns (breaks trait x species links)
  #phylogenetic shuffling here
  f.col <- fourthcorner(Env, Com, Trait, modeltype = 4, p.adjust.method.G="none", p.adjust.method.D = "none", nrepet=9999)
  
  # shuffle columns phylogenetically
  # source Harmon & Glor 2010 code
  source("~/Dropbox/fourth corner paper/4th corner codes/updated analysis/R codes/TinaCodesPhyloTree/phylogenetic_permutation_functions.R")
  #source("/Users/Gabriellapardee/Dropbox/fourth corner paper/data/phylogenetic_permutation_functions.R")
  
  # empty matrices for null values
  pp.D2 <- matrix(nrow=9999, ncol=length(f.col$tabD2$obs))
  pp.D <- matrix(nrow=9999, ncol=length(f.col$tabD$obs))
  pp.G <- matrix(nrow=9999, ncol=length(f.col$tabG$obs))
  
  # get fourthcorner observervation for Its PP shuffles. PP=phylogenetic permutation
  for(i in 1:9999) {
    
    x <- tryCatch(phyloPermute(tree, k=1), error=function(e) e)
    if(!inherits(x, "error")){ 
      Ldat.new <- Com[, match(x, colnames(Com))]
      
      pfc <- fourthcorner(Env, Ldat.new, Trait, nrepet=0)
      
      pp.D2[i, ] <- pfc$tabD2$obs
      pp.D[i, ] <- pfc$tabD$obs
      pp.G[i, ] <- pfc$tabG$obs			}
    
    ### test: should produce same results as regular Model 4
    # PP <- colnames(Ldat)[sample(1:ncol(Ldat), size = ncol(Ldat), replace = F)]
    
  }
  # we added this permutation in because we want the rows to be permuted within block. therefore, with the shuffle function we can break the rownames into two columns to permute within block. 
  # set this permutation up the same as above where we built empty matrices for null values
  
  # first read in dataset to shuffle 
  Blocks<-read.csv("~/Dropbox/Fourth corner paper/Fourth corner analysis/Block_year.csv") #made a new csv file that only contains the block_years
  y=as.vector(Blocks$Block_year)  #created a vector 
  CTRL=how(within=Within(type="free"), plots=Plots(strata=gl(3,8, labels=paste(y))))  #Plots(strata=gl(3,8) means that we have 3 blocks and 8 years
  CTRL 
  
  # empty matrices for null values, sitep= site permutations
  sitep.D2 <- matrix(nrow=9999, ncol=length(f.col$tabD2$obs))
  sitep.D <- matrix(nrow=9999, ncol=length(f.col$tabD$obs))
  sitep.G <- matrix(nrow=9999, ncol=length(f.col$tabG$obs))
  # get fourthcorner observervation for Its sitep shuffles
  
  
  #for(i in 1:9999) {
  # x <-shuffle(24, CTRL)#this is the function that shuffles within sites
  #x_names <- y[x] #want to convert numbers to names
  #if(!inherits(x_names, "error")){ 
  # Ldat.new1 <- Com[match(x_names, rownames(Com)),]
  
  #  sfc <- fourthcorner(Env, Ldat.new1, Trait, nrepet=0)
  
  #  sitep.D2[i, ] <- sfc$tabD2$obs
  # sitep.D[i, ] <- sfc$tabD$obs
  #sitep.G[i, ] <- sfc$tabG$obs			}
  
  
  #}
  # p-values
  # sim is matrix of simulated values (reps as rows, n tests as columns)
  #these are p values for PP
  ppd2 <- with(f.col, as.krandtest(sim = pp.D2, obs = tabD2$obs, alter = tabD2$alter))
  ppd <- with(f.col, as.krandtest(sim = pp.D, obs = tabD$obs, alter = tabD$alter))
  ppG <- with(f.col, as.krandtest(sim = pp.G, obs = tabG$obs, alter = tabG$alter))
  
  #get p values for within site shuffling (sitep)
  sitepd2 <- with(f.col, as.krandtest(sim = sitep.D2, obs = tabD2$obs, alter = tabD2$alter))
  sitepd <- with(f.col, as.krandtest(sim = sitep.D, obs = tabD$obs, alter = tabD$alter))
  sitepG <- with(f.col, as.krandtest(sim = sitep.G, obs = tabG$obs, alter = tabG$alter))
  
  # final output tables; look at all 3 tests together for simplicity
  TabD2 <- with(f.col$tabD2, data.frame(
    Test = gsub(" ", "", names),
    Obs = round(obs, 4),
    SD = round(sqrt(expvar$Variance), 4),
    Alter = alter,
    Pvalue.sites = f.row$tabD2$pvalue,
    Pvalue.sp = f.col$tabD2$pvalue,
    Pvalue.sp.PP = round(ppd2$pvalue,5),
    Pvalue.sp.site=round(sitepd2$pvalue, 5)
  ))
  
  TabD <- with(f.col$tabD, data.frame(
    Test = gsub(" ", ".", names),
    Obs = round(obs, 4),
    Alter = alter,
    Pvalue.sites = f.row$tabD$pvalue,
    Pvalue.sp = f.col$tabD2$pvalue,
    Pvalue.sp.PP = round(ppd$pvalue,5),
    Pvalue.sp.site=round(sitepd$pvalue,5)
  ))
  
  TabG <- with(f.col$tabG, data.frame(
    Test = gsub(" ", "", names),
    Obs = round(obs, 4),
    Alter = alter,
    Pvalue.sites = f.row$tabG$pvalue,
    Pvalue.sp = f.col$tabD2$pvalue,
    Pvalue.sp.PP = round(ppG$pvalue,5),
    Pvalue.sp.site=round(sitepG$pvalue,5)
  ))
  
  list(TabG = TabG, TabD2 = TabD2, TabD = TabD, n.sp = nrow(Trait), n.spec = sum(Com))	
  
}

Fourthcorner1=Fourth.corner(Trait, 9999)
#changed the maxprint options to print out more

print(Fourthcorner1$TabD2)
write.csv(Fourthcorner1$TabD2, file="EmergenceOutputD2_DroughtCon.csv")

print(Fourthcorner1$TabG)
write.csv(Fourthcorner1$TabG, file="EmergenceTabG_DroughtCon.csv")
