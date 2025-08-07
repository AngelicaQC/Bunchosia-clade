### Correlation Analyses for Bunchosia Clade ###

## 1. Load requiered libraries####
library("phytools")
library("ape")
library("geiger")
library("maps")

##2. Load and prune the tree####
tree <- read.nexus("MalpighiaceaeTree_301121.tree")
###Prune the tree
speciesTree<- c("Bunchosia_glandulifera"
                ,"Bunchosia_maritima_AL"
                ,"Bunchosia_montana_WC"
                ,"Bunchosia_nitida"
                ,"Bunchosia_polystachia"
                ,"Dicella_nucifera"
                ,"Echinopterys_eglandulosa"
                ,"Echinopterys_setosa_WC"
                ,"Heladena_multiflora"
                ,"Henleophytum_echinatum"
                ,"Hiptage_benghalensis"
                ,"Malpighia_glabra"
                ,"Niedenzuella_stannea"
                ,"Thryallis_longifolia"
                ,"Tetrapterys_schiedeana"
                ,"Tristellateia_australasiae")
pruned_tree<-drop.tip(tree,tree$tip.label[!(tree$tip.label %in% speciesTree)])
pruned_tree<-ladderize(pruned_tree, F)

# Check if the tree is binary
is.binary(pruned_tree) ## Should return TRUE

##3. Load trait data####
data<- read.csv("BunchosiaCladeData.csv", header=TRUE, row.names = 1)

#Verify that tree tip labels match dataset row names
name.check(pruned_tree, data) #should return "Ok"

#Select traits for correlation analysis
Habit<- setNames(data$Habit.binary, rownames(data))
Biphasic<- setNames(data$Biphasic.development, rownames(data))
VesselDimorphism<- setNames(data$Vessel.dimorphism, rownames(data))

#4. Pagel's Correlation Test ####

#Habit vs. Biphasic Development
Habit_Biphasic<-fitPagel(pruned_tree,x=Habit, y=Biphasic, pi="fitzjohn")
print(Habit_Biphasic)

#Plot Pagel's model
plot(Habit_Biphasic, cex.main=1,cex.sub=0.8,
     cex.traits=0.8,cex.rates=0.7,
     lwd.by.rate=TRUE,max.lwd=4.5, spacer=0.4)

#Habit vs. Vessel dimorphism
Habit_Dimorphism<-fitPagel(pruned_tree,x=Habit, y=VesselDimorphism, pi="fitzjohn")
print(Habit_Dimorphism)

#Plot Pagel's model
plot(Habit_Dimorphism, cex.main=1,cex.sub=0.8,
     cex.traits=0.8,cex.rates=0.7,
     lwd.by.rate=TRUE,max.lwd=4.5, spacer=0.4)

