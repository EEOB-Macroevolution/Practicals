## ----Ordination_prep ----
library(geomorph)
data(plethspecies) 
Y.gpa <- gpagen(plethspecies$land, print.progress = F)    #GPA-alignment

## ----PCA_example ----
PCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
plot(PCA, main = "PCA", pch = 21, bg = "red")

## ----Phylomorph_example ----
plot(PCA, phylo = TRUE, pch = 21, bg = "red",
     main = "Phylomorphospace")

## ----phyloPCA_example ----
phylo.PCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, GLS = TRUE)
plot(phylo.PCA, phylo = TRUE, pch = 21, bg="red",
     main = "phylogenetic PCA")

## ----PaCA_example ----
PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
                      align.to.phy = TRUE, GLS = TRUE)
plot(PaCA.gls, phylo = TRUE, pch = 21, bg = "red",
     main = "PaCA: Phylogenetically-Aligned Component Analysis")

## ----read_data ----
library(RRPP)
library(geiger)
library(phytools)
tree.best<-read.nexus("../data/Consensus of 1000 salamander trees.nex") #Maximum Credible Tree
plethdata<-read.csv("../data/meandata-CinGlutOnly.csv",row.names=1, header=TRUE )
svl<-plethdata[,2];names(svl)<-row.names(plethdata) 
groups<-as.factor(plethdata[,1]); names(groups)<-row.names(plethdata)
Y<-as.matrix(plethdata[,c(3:8)])  ##This is our multivariate data
plethtree<-treedata(phy = tree.best,data = plethdata, warnings = FALSE)$phy
plot(plethtree)
axisPhylo(1)
PCA<-gm.prcomp(Y,phy=plethtree)   #principal components of Y
plot(PCA, pch=21,bg="black") ##we can put the PCA directly into the plot function

## ----phylo_reg----
PhyCov <- vcv.phylo(plethtree) #generate phylogenetic covariance matrix
rdf <- rrpp.data.frame(Y=Y, svl=svl,groups=groups,Cov = PhyCov)  #RRPP data frame
PGLS.reg <- lm.rrpp(Y ~ svl, Cov = PhyCov, data = rdf, iter = 999, 
                   print.progress = FALSE,) # randomize residuals  
anova(PGLS.reg)

## ---- phylo_ANOVA----
PGLS.aov<-lm.rrpp(Y ~ groups, Cov = PhyCov, data = rdf, iter = 999, print.progress = FALSE) # randomize residuals
anova(PGLS.aov)  #no difference once phylogeny considered

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



