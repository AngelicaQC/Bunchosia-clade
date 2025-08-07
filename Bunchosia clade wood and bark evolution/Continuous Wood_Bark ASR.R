##Ancestral state reconstruction for Continuous Traits in the Bunchosia Clade###

# Set working directory
setwd()

##1. Charge libraries####
library(ape)
library(phytools)


##2. Load and prune the phylogenetic tree####
tree <- read.nexus("MalpighiaceaeTree_301121.tree")

# Define the species to retain
speciesTree<- c("Bunchosia_glandulifera","Bunchosia_maritima_AL"
                ,"Bunchosia_montana_WC","Bunchosia_nitida"
                ,"Bunchosia_polystachia","Dicella_nucifera"
                ,"Echinopterys_eglandulosa","Echinopterys_setosa_WC"
                ,"Heladena_multiflora","Henleophytum_echinatum"
                ,"Hiptage_benghalensis","Malpighia_glabra"
                ,"Niedenzuella_stannea","Thryallis_longifolia"
                ,"Tetrapterys_schiedeana","Tristellateia_australasiae")
#prune tree to include only selected species
pruned_tree<-drop.tip(tree,tree$tip.label[!(tree$tip.label %in% speciesTree)])
pruned_tree<-ladderize(pruned_tree, right=FALSE)
class(pruned_tree)<-"phylo"
# Plot the pruned tree
plotTree(pruned_tree, fsize= 0.8, offset=0.3, lwd=2, ftype="i")

# Check if the tree is binary
is.binary(pruned_tree) #true

##3.Load trait data####
data<- read.csv("BunchosiaCladeData2.csv", header=TRUE, row.names = 1)
colnames(data) <- gsub(" ", ".", colnames(data))


#Check for name matching between tree and trait dataset
name.check(pruned_tree, data) #Should return "OK"

#3.1. Select continuous trait: Vessel frequency
VesselFreq<-setNames(data$Vessel.frequency, rownames(data))
str(VesselFreq)

##4. Ancestral state reconstruction (ASR) for continuous trait####
#Estimate ancestral states using maximum likelihood under a Brownian motion model of trait evolution
#We will also compute variances & 95% confidence intervals for each node
fitVF<-fastAnc(pruned_tree, VesselFreq,vars=TRUE,CI=TRUE)

#print
print(fitVF)
range(fitVF) ## compare to root node


##5. Visualization of continuous trait on the phylogeny ####
# Use contMap to visualize trait evolution
VesselF_tree<-contMap(pruned_tree, VesselFreq, plot = FALSE)
VesselF_tree

# Customize color gradient
VesselF_tree<-setMap(VesselF_tree,colors=c("yellow","orange2","violetred2","magenta4","midnightblue"))

# Plot the trait map with customized color scheme
plot(VesselF_tree,fsize=c(1,1), lwd=5, outline=FALSE,
     leg.txt="Vessel frequency", legend=0.5*max(nodeHeights(pruned_tree)))


##6. Phylogenetic Signal Analysis####
# Estimate Pagel's lambda and test significance using 1000 simulations
VesselF_lambda<-phylosig(pruned_tree, VesselFreq, method = "lambda", test = TRUE, nsim = 1000)
print(VesselF_lambda)

# Plot log-likelihood profile for lambda
plot(VesselF_lambda)
title(main = "Vessel frequency/ mm^2- Pagel's lambda")

#########################################################

#3.1. Select continuous trait: Vessel width
VesselW<-setNames(data$Vessel.width, rownames(data))
str(VesselW)

##4. Ancestral state reconstruction 
fitVw<-fastAnc(pruned_tree, VesselW,vars=TRUE,CI=TRUE)
print(fitVw)
range(VesselW)

##5. Visualization of continuous trait
VesselW_tree<-contMap(pruned_tree, VesselW, plot = FALSE)
VesselW_tree<-setMap(VesselW_tree,colors=c("yellow","orange2","violetred2","magenta4","midnightblue"))
plot(VesselW_tree,fsize=c(1,1), lwd=5, outline=FALSE,
     leg.txt="Vessel width", legend=0.5*max(nodeHeights(pruned_tree)))

##6. Phylogenetic Signal Analysis####
VesselW_lambda<-phylosig(pruned_tree, VesselW, method = "lambda", test = TRUE, nsim = 1000)
print(VesselW_lambda)
plot(VesselW_lambda)
title(main = "Vessel width (µm)- Pagel's lambda")

#########################################################

#3.1. Select continuous trait: Intervessel pits size
Pits<-setNames(data$Intervessel.pits, rownames(data))
str(Pits)

##4. Ancestral state reconstruction 
fitPits<-fastAnc(pruned_tree,Pits,vars=TRUE,CI=TRUE)
print(fitPits)
range(Pits)

##5. Visualization of continuous trait
Pits_tree<-contMap(pruned_tree, Pits, plot = FALSE)
Pits_tree<-setMap(Pits_tree,colors=c("yellow","orange2","violetred2","magenta4","midnightblue"))
plot(Pits_tree,fsize=c(1,1), lwd=5, outline=FALSE,
     leg.txt="Intervessel pits size", legend=0.5*max(nodeHeights(pruned_tree)))

##6. Phylogenetic Signal Analysis####
Pits_lambda<-phylosig(pruned_tree, Pits, method = "lambda", test = TRUE, nsim = 1000)
print(Pits_lambda)
plot(Pits_lambda)
title(main = "Intervessel pits size (µm)- Pagel's lambda")

#########################################################

#3.1. Select continuous trait: Parenchyma proportion
Parenchyma<-setNames(data$Parenchyma.proportion, rownames(data))
str(Parenchyma)

##4. Ancestral state reconstruction 
fitParenchyma<-fastAnc(pruned_tree,Parenchyma,vars=TRUE,CI=TRUE)
print(fitParenchyma)
range(Parenchyma)

##5. Visualization of continuous trait
Parenchyma_tree<-contMap(pruned_tree, Parenchyma, plot = FALSE)
Parenchyma_tree<-setMap(Parenchyma_tree,colors=c("yellow","orange2","violetred2","magenta4","midnightblue"))
plot(Parenchyma_tree,fsize=c(1,1), lwd=5, outline=FALSE,
     leg.txt="Parenchyma proportion", legend=0.5*max(nodeHeights(pruned_tree)))

##6.Phylogenetic Signal Analysis####
Parenchyma_l<-phylosig(pruned_tree, Parenchyma, method = "lambda", test = TRUE, nsim = 1000)
print(Parenchyma_l)
plot(Parenchyma_l)
title(main = "Parenchyma proportion (mm^2)- Pagel's lambda")

#########################################################

# Define the species to retain
speciesTree2<- c("Bunchosia_maritima_AL","Bunchosia_montana_WC"
                ,"Bunchosia_polystachia","Dicella_nucifera"
                ,"Echinopterys_eglandulosa","Echinopterys_setosa_WC"
                ,"Heladena_multiflora","Henleophytum_echinatum"
                ,"Hiptage_benghalensis","Malpighia_glabra"
                ,"Niedenzuella_stannea","Thryallis_longifolia"
                ,"Tetrapterys_schiedeana","Tristellateia_australasiae")
#prune tree to include only selected species
pruned_tree2<-drop.tip(tree,tree$tip.label[!(tree$tip.label %in% speciesTree2)])
pruned_tree2<-ladderize(pruned_tree2, right=FALSE)

#3.1. Select continuous trait: STE area
STEArea<-setNames(data$STE.area, rownames(data))
STEArea <- na.omit(STEArea)
str(STEArea)

STEArea <- STEArea[names(STEArea) %in% pruned_tree2$tip.label]

##4. Ancestral state reconstruction 
fitSTEArea<-fastAnc(pruned_tree2,STEArea,vars=TRUE,CI=TRUE)
print(fitSTEArea)
range(STEArea)

##5. Visualization of continuous trait
STEArea_tree<-contMap(pruned_tree2, STEArea, plot = FALSE)
STEArea_tree<-setMap(STEArea_tree,colors=c("yellow","orange2","violetred2","magenta4","midnightblue"))
plot(STEArea_tree,fsize=c(1,1), lwd=5, outline=FALSE,
     leg.txt="STE Area", legend=0.5*max(nodeHeights(pruned_tree2)))

##6. Phylogenetic Signal Analysis####
STEArea_lambda<-phylosig(pruned_tree2, STEArea, method = "lambda", test = TRUE, nsim = 1000)
print(STEArea_lambda)
plot(STEArea_lambda)
title(main = "STE Area (mm^2)- Pagel's lambda")



