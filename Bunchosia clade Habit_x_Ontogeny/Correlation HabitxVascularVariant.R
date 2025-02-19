### Correlation Analyses for Bunchosia Clade ###

##1.Load libraries####
library("phytools")
library("ape")
library("geiger")
library("maps")

##2.Load and prune the tree####
tree <- read.nexus("MalpighiaceaeTree_301121.tree")
###pruned tree
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
class(pruned_tree)<-"phylo"

# Check if the tree is binary
is.binary(pruned_tree) ## Should return TRUE

##3.Load trait data####
data<- read.csv("BunchosiaCladeData.csv", header=TRUE, row.names = 1)

#Verify that tree tip labels match dataset row names
name.check(pruned_tree, data) #should return "Ok"

#Select traits for correlation analysis
Habit1<-setNames(data$Habit.binary, rownames(data))
Variant<-setNames(data$Cambial.Variant.presence, rownames(data))

#4.Pagel's Correlation Test####
#Fit Pagel's (1994) model to test dependency between traits
Habit_Variant<-fitPagel(pruned_tree,x=Habit1, y=Variant, pi="fitzjohn")
Habit_Variant
anova(Habit_Variant)

#Plot Pagel's model
plot(Habit_Variant, cex.main=1,cex.sub=0.8,
     cex.traits=0.8,cex.rates=0.7,
     lwd.by.rate=TRUE,max.lwd=4.5, spacer=0.4)

##5.Plot both traits####
#Model fitting for Cambial Variants
fitER<-fitMk(pruned_tree,Variant,model="ER")
fitARD<-fitMk(pruned_tree,Variant,model="ARD", pi="estimated")
Variant_aov<-anova(fitER, fitARD) #ER is best
#Stochastic Character Mapping for Trait
Variant_simmap<- make.simmap(pruned_tree, Variant, model="ER",nsim=500)
Variant_summ<-summary(Variant_simmap)

#Model fitting for Habit
fitER<-fitMk(pruned_tree,Habit1,model="ER")
fitARD<-fitMk(pruned_tree,Habit1,model="ARD", pi="estimated")
Variant_aov<-anova(fitER, fitARD) #ER is best
#Stochastic Character Mapping for Trait
Habit1_simmap<- make.simmap(pruned_tree, Habit1, model="ER",nsim=500)
Habit1_summ<-summary(Habit1_simmap)

#Define colors for traits
cols_h<-setNames(c("magenta4", "orange2")[1:length(unique(Habit1))],sort(unique(Habit1)))
cols_v<-setNames(c("orange2","magenta4")[1:length(unique(Variant))],sort(unique(Variant)))

#Generate PDF for visualization
pdf(file="Pagel_HabitxVariant.pdf", height=11,width=15)
par(mfrow=c(1,2))
plot(Habit1_summ,fsize=0.8,ftype="i",cex=0.8,
     lwd=3,colors=cols_h, offset=0.95)
add.simmap.legend(x=1,y=15,prompt=FALSE, shape="circle", fsize=1.2, 
                  colors= cols_h, c("Climbing","Self-supporting"))
title(main="Habit", line=-2, cex.main=1.5)

plot(Variant_summ, invert=TRUE,lwd=3, direction="leftwards",ftype="off",
     colors=cols_v, cex=0.8)
add.simmap.legend(x=27,y=15,prompt=FALSE, shape="circle", fsize=1.2, 
                  colors= cols_v, c("Vascular variant absent","Vascular variant present"))
title(main="Vascular Variant", line=-2, cex.main=1.5)

dev.off()

#Density plots of character state changes
par(mfrow=c(2,1))
H<-density(Habit1_simmap)
plot(H, Cex.axis=0.4, cex.lab=0.5)
V<-density(Variant_simmap)
plot(V, Cex.axis=0.4, cex.lab=0.5)
