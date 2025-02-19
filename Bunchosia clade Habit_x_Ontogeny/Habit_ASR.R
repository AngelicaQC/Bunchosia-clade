##Ancestral state reconstruction for the Bunchosia clade##
##Habit##

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
is.binary(pruned_tree) #true

##3.load trait data####
data<- read.csv("BunchosiaCladeData.csv", header=TRUE, row.names = 1)

#Verify if tree tip labels match dataset row names
name.check(pruned_tree, data) #Should return "OK"

#Select Habit trait
Habit<-setNames(data$Habit, rownames(data))

##4.Model fitting for ASR####
# Fit Mk models (Equal Rates (ER), All Rates Different (ARD), Symmetric (SYM))
fitER<-fitMk(pruned_tree,Habit,model="ER", pi="fitzjohn")
fitARD<-fitMk(pruned_tree,Habit,model="ARD", pi="fitzjohn")
fitSYM<-fitMk(pruned_tree,Habit, model="SYM", pi="fitzjohn")
#Compare models
Habit_aov<-anova(fitER, fitARD, fitSYM)

#Plot transition models
layout(matrix(c(1,1,2,2,4,3,3,4),2,4,byrow=TRUE))
plot(fitER, cex.traits=0.7,spacer=0.35,lwd=2, cex.rates=0.8)
legend("topleft",legend=paste("AIC =",round(AIC(fitER),1)),bty="n", cex=1)
mtext(cex=0.9,"a)Habit transition rates -\"ER\" model", font.main=1)

plot(fitARD, width=TRUE, color=TRUE, cex.traits=0.7,spacer=0.35,lwd=2, cex.rates=0.8)
legend("topleft",legend=paste("AIC =",round(AIC(fitARD),1)),bty="n", cex=1)
mtext(cex=0.9, "b)Habit transition rates -\"ARD\" model", font.main=1)

plot(fitSYM, width=TRUE, color=TRUE, cex.traits=0.7,spacer=0.35,lwd=2, cex.rates=0.8)
legend("topleft",legend=paste("AIC =",round(AIC(fitSYM),1)),bty="n", cex=1)
mtext(cex=0.9, "c)Habit transition rates -\"SYM\" model", font.main=1)

dev.off()

##5.Stochastic Character Mapping (500 Simulations)####
#Generate stochastic maps
Habit_simmap<- make.simmap(pruned_tree, Habit, model="ER",nsim=500)
Habit_summ<-summary(Habit_simmap)
Habit_summ

#Define colors for different habits
cols<-setNames(c("magenta4", "orange2","yellow", "midnightblue")
               [1:length(unique(Habit))],sort(unique(Habit)))
#Plot summary of stochastic maps
plot(Habit_summ,fsize=1,ftype="i",cex=0.5,lwd=3,offset=0.3, colors=cols)
add.simmap.legend(x=2,y=16,prompt=FALSE, shape="circle", fsize=1.2, 
                  colors= cols, 
                  c("Climbing shrub","Liana", "Shrub","Tree"))
cladelabels(text="Bunchosia Clade", node=c(22), offset=13, wing.length=0.5, cex=1.2,
           orientation="vertical")

#Plot all 500 simulations
par(mfrow=c(20,25))
sapply(Habit_simmap, plotSimmap, colors=cols, lwd=1, ftype="off")

dev.off()

#Plot density of character state changes
H_density<-density(Habit_simmap)
pdf(file = "Habit_density.pdf", height = 8, width = 11)
plot(H_density, Cex.axis=0.5, cex.lab=0.6)
dev.off()

##6.Phylogenetic Signal Analysis####
SF_ER<-fitDiscrete(pruned_tree, Habit, model ="ER", transform = "lambda")
print(SF_ER)

