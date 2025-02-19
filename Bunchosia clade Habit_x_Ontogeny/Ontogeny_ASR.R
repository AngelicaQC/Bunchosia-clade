##Ancestral state reconstruction for the Bunchosia clade##
##Ontogeny##

# Set working directory
setwd()

##1.Charge libraries####
library("phytools")
library("ape")
library("geiger")
library("maps")

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
#prune tree
pruned_tree<-drop.tip(tree,tree$tip.label[!(tree$tip.label %in% speciesTree)])
pruned_tree<-ladderize(pruned_tree, right=FALSE)
class(pruned_tree)<-"phylo"
# Plot the pruned tree
plotTree(pruned_tree, fsize= 0.8, offset=0.3, lwd=2, ftype="i")

# Check if the tree is binary
is.binary(pruned_tree)  #Should return TRUE

##3.load trait data####
data<- read.csv("BunchosiaCladeData.csv", header=TRUE, row.names = 1)

#Verify if tree tip labels match dataset row names
name.check(pruned_tree, data) #Should return "OK"

#Select Ontogeny trait
Ontogeny<-setNames(data$Cambial.Variant.type, rownames(data))

##4.Model fitting for ASR####
# Fit Mk models (Equal Rates (ER), All Rates Different (ARD), Symmetric (SYM))
fitER<-fitMk(pruned_tree,Ontogeny,model="ER", pi="fitzjohn")
fitARD<-fitMk(pruned_tree,Ontogeny,model="ARD", pi="fitzjohn")
fitSYM<-fitMk(pruned_tree,Ontogeny, model="SYM", pi="fitzjohn")
#Compare models
Ontogeny_aov<-anova(fitER, fitARD, fitSYM)

#Plot transition models
layout(matrix(c(1,1,2,2,4,3,3,4),2,4,byrow=TRUE))
plot(fitER, cex.traits=0.7,spacer=0.35,lwd=2, cex.rates=0.8)
legend("topleft",legend=paste("AIC =",round(AIC(fitER),1)),bty="n", cex=1)
mtext(cex=0.9,"a)Ontogeny transition rates -\"ER\" model", font.main=1)

plot(fitARD, width=TRUE, color=TRUE, cex.traits=0.7,spacer=0.35,lwd=2, cex.rates=0.8)
legend("topleft",legend=paste("AIC =",round(AIC(fitARD),1)),bty="n", cex=1)
mtext(cex=0.9, "b)Ontogeny transition rates -\"ARD\" model", font.main=1)

plot(fitSYM, width=TRUE, color=TRUE, cex.traits=0.7,spacer=0.35,lwd=2, cex.rates=0.8)
legend("topleft",legend=paste("AIC =",round(AIC(fitSYM),1)),bty="n", cex=1)
mtext(cex=0.9, "c)Ontogeny transition rates -\"SYM\" model", font.main=1)
dev.off()

##5.Stochastic Character Mapping (500 Simulations)####
#Generate stochastic maps
Ontogeny_simmap<- make.simmap(pruned_tree, Ontogeny, model="ER",nsim=500)
Ontogeny_summ<-summary(Ontogeny_simmap)
Ontogeny_summ
#Define colors for different habits
cols<-setNames(c("midnightblue","magenta4","violetred3", "orange2","yellow")
               [1:length(unique(Ontogeny))],sort(unique(Ontogeny)))
#Plot summary of stochastic maps
plot(Ontogeny_summ,fsize=1,ftype="i",cex=0.5,lwd=3,offset=0.3, colors=cols)
add.simmap.legend(x=2,y=16,prompt=FALSE, shape="circle", fsize=1.2, 
                  colors= cols, c("Conspicuous rays", "Interxylary phloem", "Phloem wedges+Interxylary phloem",
                                  "Regular stem", "Wavy cambium"))
cladelabels(text="Bunchosia Clade", node=c(22), offset=13, wing.length=0.5, cex=1,
            orientation="vertical")

#Plot all 500 simulations
par(mfrow=c(20,25))
sapply(Ontogeny_simmap, plotSimmap, colors=cols, lwd=1, ftype="off")
dev.off()

#Plot density of character state changes
O_density<-density(Ontogeny_simmap)
pdf(file = "Ontogeny_density.pdf", height = 8, width = 11)
plot(O_density, cex.axis=0.5, cex.lab=0.6)

dev.off()

##6.Phylogenetic Signal Analysis####
Onto_ER<-fitDiscrete(pruned_tree, Ontogeny, model ="ER", transform = "lambda")
print(Onto_ER)

