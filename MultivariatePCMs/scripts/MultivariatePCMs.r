# Multivariate PCMs

#Read Tree and data; prune
library(geomorph)
library(geiger)
library(phytools)
tree.best<-read.nexus("TutorialData/Consensus of 1000 salamander trees.nex") #Maximum Credible Tree
plethdata<-read.csv("TutorialData/meandata-CinGlutOnly.csv",row.names=1, header=TRUE )
 Y<-plethdata[,c(3:8)]
     Y<-apply(Y,2,as.numeric);row.names(Y)<-row.names(plethdata)
 size<-plethdata[,2];names(size)<-row.names(plethdata) 
 gp<-as.factor(plethdata[,1]); names(gp)<-row.names(plethdata)
plethtree<-treedata(phy = tree.best,data = plethdata, warnings = FALSE)$phy
plot(plethtree)
axisPhylo(1)

gdf <- geomorph.data.frame(Y=Y, size=size,gp=gp,phy = plethtree) #needed for geomorph

Y.pc<-prcomp(Y)$x   #principal components of Y
plot(Y.pc,asp=1, pch=21,bg="black")
text(x=2+Y.pc[,1],y=.5+Y.pc[,2],row.names(Y.pc),cex=.5)

#### Analyses
#0: Phylomorphospace (ie, plot the data)
phylomorphospace(tree = plethtree,X = Y.pc[,1:2])

#Note: in next version of geomorph, this will be deprecated. gm.prcomp will be used instead
plotGMPhyloMorphoSpace(phy = plethtree,A = arrayspecs(Y,p=1,k=6), ancStates = FALSE, tip.labels = FALSE) 
plotGMPhyloMorphoSpace(phy = plethtree,A = arrayspecs(Y,p=1,k=6),ancStates = FALSE, zaxis = "time") 

#1: Phylogenetic regression via phylo-transform and RRPP: Regression of multivariate data~size
PGLS.reg<-procD.pgls(Y ~ size, phy = phy, data = gdf, iter = 999, print.progress = FALSE) # randomize residuals
summary(PGLS.reg)
plot(PGLS.reg, type = "regression", predictor=size,reg.type = "RegScore")

#2: Phylogenetic ANOVA via phylo-transform and RRPP: ANOVA of multivariate data~group
PGLS.aov<-procD.pgls(Y ~ gp, phy = phy, data = gdf, iter = 999, print.progress = FALSE) # randomize residuals
summary(PGLS.aov)  #no difference once phylogeny considered

#3: Phylogenetic PLS
PLS.Y <- phylo.integration(A = Y[,1:3], A2 = Y[,4:6], phy= plethtree, print.progress = FALSE)
summary(PLS.Y)
plot(PLS.Y)

#4: Phylogenetic signal
PS.shape <- physignal(A=Y,phy=plethtree,iter=999, print.progress = FALSE)
summary(PS.shape)
plot(PS.shape)

#5a: Compare Evolutionary Rates Among Clades
ER<-compare.evol.rates(A=Y, phy=plethtree,gp=gp,iter=999, print.progress = FALSE)
summary(ER)   #significantly higher rate of morphological evolution 'large' Plethodon
plot(ER)

#5b: Compare Evolutionary Rates Among Traits
var.gp <- c("B","A","A","B","B","B")  #head (A) vs. body (B)
EMR<-compare.multi.evol.rates(A=Y,gp=var.gp, Subset=TRUE, phy= plethtree,iter=999, print.progress = FALSE)
summary(EMR) #Limb traits evolve faster than head traits
plot(EMR)

