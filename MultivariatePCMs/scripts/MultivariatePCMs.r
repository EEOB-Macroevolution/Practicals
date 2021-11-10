
## ----read_data ----
library(geomorph)
library(geiger)
library(phytools)
tree.best<-read.nexus("../data/Consensus of 1000 salamander trees.nex") #Maximum Credible Tree
plethdata<-read.csv("../data/meandata-CinGlutOnly.csv",row.names=1, header=TRUE )
Y<-plethdata[,c(3:8)] ##This is our multivariate data
Y<-apply(Y,2,as.numeric);row.names(Y)<-row.names(plethdata) ##This is our formated data
svl<-plethdata[,2];names(svl)<-row.names(plethdata) 
groups<-as.factor(plethdata[,1]); names(groups)<-row.names(plethdata)
plethtree<-treedata(phy = tree.best,data = plethdata, warnings = FALSE)$phy
plot(plethtree)
axisPhylo(1)


PCA<-gm.prcomp(Y,phy=plethtree)   #principal components of Y
plot(PCA, pch=21,bg="black") ##we can put the PCA directly into the plot function

## ----Phylomorphospace ----
plot(PCA,phylo=TRUE)

## ----phylo_pca ----

phylo_PCA<-gm.prcomp(Y,phy=plethtree)   #principal components of Y
plot(phylo_PCA, pch=21,bg="black") ##we can put the PCA directly into the plot function

## ----phylo_reg----
geo_df <- geomorph.data.frame(Y=Y, svl=svl,groups=groups,phy = plethtree) #needed for geomorph
PGLS.reg<-procD.pgls(Y ~ svl, phy = phy, data = geo_df, iter = 999, print.progress = FALSE) # randomize residuals
summary(PGLS.reg)
plot(PGLS.reg, type = "regression", predictor=svl,reg.type = "RegScore")


## ---- phylo_ANOVA----
PGLS.aov<-procD.pgls(Y ~ groups, phy = phy, data = geo_df, iter = 999, print.progress = FALSE) # randomize residuals
summary(PGLS.aov)  #no difference once phylogeny considered

## ---- phylo_PLS ----
PLS.Y <- phylo.integration(A = Y[,1:3], A2 = Y[,4:6], phy= plethtree, print.progress = FALSE)
summary(PLS.Y)
plot(PLS.Y)

## ----phylo_signal----
PS.shape <- physignal(A=Y,phy=plethtree,iter=999, print.progress = FALSE)
summary(PS.shape)
plot(PS.shape)

## ----evol_rates----
ER<-compare.evol.rates(A=Y, phy=plethtree,gp=groups,iter=999, print.progress = FALSE)
summary(ER)   #significantly higher rate of morphological evolution 'large' Plethodon
plot(ER)

var.gp <- c("B","A","A","B","B","B")  #head (A) vs. body (B)
EMR<-compare.multi.evol.rates(A=Y,gp=var.gp, Subset=TRUE, phy= plethtree,iter=999, print.progress = FALSE)
summary(EMR) #Limb traits evolve faster than head traits
plot(EMR)

## ----pleth_read ----
data("plethspecies")
pleth_tree<-plethspecies$phy
landmarks<-plethspecies$land
procD_landmarks <- gpagen(landmarks)



