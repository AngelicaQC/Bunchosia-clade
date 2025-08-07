##Ancestral state reconstruction for the Bunchosia clade##

####
##nalysis of the evolution of biphasic development in wood##

# Set working directory
setwd()

##1.Charge libraries####
library(pacman)
p_load(phytools, ape, geiger, maps)


##2.Load and prune the phylogenetic tree####
tree <- read.nexus("MalpighiaceaeTree_301121.tree")
tree #167 tips and 166 internal nodes
# Define the species to retain
speciesTree<- c("Bunchosia_glandulifera","Bunchosia_maritima_AL"
                ,"Bunchosia_montana_WC","Bunchosia_nitida"
                ,"Bunchosia_polystachia","Dicella_nucifera"
                ,"Echinopterys_eglandulosa","Echinopterys_setosa_WC"
                ,"Heladena_multiflora","Henleophytum_echinatum"
                ,"Hiptage_benghalensis","Malpighia_glabra"
                ,"Niedenzuella_stannea","Thryallis_longifolia"
                ,"Tetrapterys_schiedeana","Tristellateia_australasiae")

#prune the tree to retain only the selected taxa
pruned_tree<-drop.tip(tree,tree$tip.label[!(tree$tip.label %in% speciesTree)])
pruned_tree<-ladderize(pruned_tree, right=FALSE)
class(pruned_tree)<-"phylo"
# Plot the pruned tree
plotTree(pruned_tree, fsize= 0.8, offset=0.3, lwd=2, ftype="i")

# Check if the tree is binary
is.binary(pruned_tree) #true

##3.Load trait data####
# Trait data must have tip labels as row names
data<- read.csv("BunchosiaCladeData2.csv", header=TRUE, row.names = 1)

#Verify if tree tip labels match dataset row names
name.check(pruned_tree, data) #Should return "OK"

#3.1. Select the trait of interest: Biphasic development ####
BD<-setNames(data$Biphasic.development, rownames(data))

##4.Model fitting for ASR####
# Fit Mk models (Equal Rates (ER), All Rates Different (ARD))
fitER<-fitMk(pruned_tree, BD,model="ER", pi="fitzjohn")
fitARD<-fitMk(pruned_tree, BD,model="ARD", pi="fitzjohn")
#Compare model fits using AIC
BD_aov<-anova(fitER, fitARD)

#Plot transition models
par(mfrow=c(1,2))
plot(fitER, cex.traits=0.65,spacer=0.4,lwd=2, cex.rates=0.8)
legend("topleft",legend=paste("AIC =",round(AIC(fitER),1)),bty="n", cex=0.8)
mtext(cex=0.9,"a)Biphasic development transition rates -\"ER\" model", font.main=1)

plot(fitARD, width=TRUE, color=TRUE, cex.traits=0.65,spacer=0.4, lwd=2, cex.rates=0.8)
legend("topleft",legend=paste("AIC =",round(AIC(fitARD),1)),bty="n", cex=0.8)
mtext(cex=0.9, "b)Biphasic development transition rates -\"ARD\" model", font.main=1)

dev.off()

##5.Stochastic Character Mapping (500 Simulations)####
#Generate stochastic maps
BD_simmap<- make.simmap(pruned_tree, BD, model="ER",nsim=500)
BD_summ<-summary(BD_simmap)
BD_summ

#Define colors for different habits
cols<-setNames(c("midnightblue","orange2")
               [1:length(unique(BD))],sort(unique(BD)))
#Plot summary of stochastic maps
plot(BD_summ,fsize=1,ftype="i",cex=0.6,lwd=3,offset=0.5, colors=cols)
add.simmap.legend(x=2,y=16,prompt=FALSE, shape="circle", fsize=1, 
                  colors= cols, c("Absent","Present"))
cladelabels(text="Bunchosia Clade", node=c(18), offset=15, wing.length=0.5, cex=1.2,
            orientation="vertical")

#Plot all 500 simulations (optional, very dense output)
par(mfrow=c(20,25))
sapply(BD_simmap, plotSimmap, colors=cols, lwd=1, ftype="off")
dev.off()

#Plot density of character state changes
BD_density<-density(BD_simmap)
pdf(file = "BD_density.pdf", height = 8, width = 11)
plot(BD_density, Cex.axis=0.5, cex.lab=0.6)
dev.off()

##6.Phylogenetic Signal Analysis####
# Estimate Pagel's lambda using the ER model
BD_ER<-fitDiscrete(pruned_tree, BD, model ="ER", transform = "lambda")
print(BD_ER)
#lambda = 1.000000



# Note: Steps 3.1 to 6 can be repeated for the next discrete trait:
#"Vessel.dimorphism", "Porosity", "Vessel.arrangement", 
#"Vessel.grouping", "Apotracheal.parenchyma", "Paratracheal.parenchyma",
#"Banded.parenchyma", "Septate.fibers", "Gelatinous.fibers", 
#"Interconnected.rays", "Ray.width", "Rays.composition",
#"Prismatic.crystals.in.axial.parenchyma", "Prismatic.crystals.in.ray.parenchyma",
#"STE.distribution", "Druses.in.phloem", "Sclerenchyma.arrangement"


