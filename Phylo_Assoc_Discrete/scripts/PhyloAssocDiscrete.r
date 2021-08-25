
## ----read_data----
library(geiger, warn.conflicts=F, quietly=T) #also loads ape
library(phytools, warn.conflicts=F, quietly=T)


tree<-read.tree("TutorialData/tree.64.tre",tree.names=T)
mydata<-read.csv('TutorialData/DiscreteData.csv', row.names=1, header=TRUE)
head(mydata) ##Head can be used to show the first few rows of a data frame 

#Match data with tree
data.pruned<-treedata(phy=tree,data = mydata, warnings=FALSE)
tree<-data.pruned$phy
mydata<-data.pruned$data


## ----plot_data----

plot.phylo(tree,show.tip.label = F) ## setting show.tip.label = F prevents the species names from appearing


## We make a colour key to denote which state gets which color
colorkey1 <- c('red','black') ##trait V1
colorkey2 <- c('green','orange') ##trait V3 
names(colorkey1) <-c(0,1) ##Here we are saying that 0s get red and 1s get black
names(colorkey2) <-c(0,1) ##Here we are saying that 0s get green and 1s get orange
V1cols<-colorkey1[as.character(mydata[,1])] ##putting the trait values in as.character() is important!
V3cols<-colorkey2[as.character(mydata[,3])]


tiplabels(pch = 17,col = V1cols,cex=0.6)
tiplabels(pch = 19,col = V3cols,offset = 0.2,cex=0.6)

## ----single_trait----
library(corHMM)

plot.phylo(tree,show.tip.label = F)
tiplabels(pch = 17,col = V1cols)

#Set up data. corHMM is picky about how it wants the data formatted
trt1<-cbind(row.names(mydata),mydata[,1])

#Set up initial rate matrices
trait1_model_er <- getStateMat4Dat(trt1,"ER")  #equal transition rates
trait1_model_ard <- getStateMat4Dat(trt1,"ARD")  #All rates different


## ----single_trait_plot
trait1_model_er
trait1_model_ard

plotMKmodel(trait1_model_er$rate.mat,rate.cat = 1)
plotMKmodel(trait1_model_ard$rate.mat,rate.cat = 1)


## ----single_trait_fit----

#Fit models: TRAIT 1
trait1_fit_er <-corHMM(phy = tree, data = trt1, rate.cat = 1, rate.mat = trait1_model_er$rate.mat)
trait1_fit_ard <-corHMM(phy = tree, data = trt1, rate.cat = 1, rate.mat = trait1_model_er$rate.mat)

##The estimated rates
trait1_fit_er$solution
trait1_fit_ard$solution

#compare models: logL and AIC
trait1_fit_er$loglik
trait1_fit_ard$loglik

##We can formally perform a likelihood ratio test of nested models by doing the following
summary_stat <- -2 * (trait1_fit_er$loglik - trait1_fit_ard$loglik) ## -2 * (model_constrained - model_full)
##This is our p-value
1-pchisq(summary_stat,df = 1 ) ##degrees of freedom is the difference in the number of parameters between the two models. In our case 2-1=1

trait1_fit_er$AIC
trait1_fit_ard$AIC




## ---- single_trait2_and_trait3 ----

##This code is outdated but the conclusions stand
# 
# #Fit models: TRAIT 2
# plot.phylo(tree,show.tip.label = F)
# tiplabels(pie = to.matrix(mydata[,3],sort(unique(mydata[,3]))),piecol=c("red", "black"),cex=.3, offset=0)
# res2.er<-corHMM(tree,trt2,rate.cat=1,rate.mat=rmat.er,node.states="marginal") 
# res2.ard<-corHMM(tree,trt2,rate.cat=1,rate.mat=rmat.ard,node.states="marginal") 
# #compare models: logL and AIC
# res2.er$loglik 
# res2.ard$loglik  #identical in this case: no need for formal LRT 
# res2.er$AICc
# res2.ard$AICc    #Choose simpler model (lower AICc b/c fewer parameters and same logL)
# res2.er$solution  #rate transition parameters
# 
# #Fit models: TRAIT 3
# plot.phylo(tree,show.tip.label = F)
# tiplabels(pie = to.matrix(mydata[,5],sort(unique(mydata[,5]))),piecol=c("red", "black"),cex=.3, offset=0)
# res3.er<-corHMM(tree,trt3,rate.cat=1,rate.mat=rmat.er,node.states="marginal") 
# res3.ard<-corHMM(tree,trt3,rate.cat=1,rate.mat=rmat.ard,node.states="marginal") 
# #compare models: logL and AIC
# res3.er$loglik 
# res3.ard$loglik
# LRT<- -2*(res3.er$loglik - res3.ard$loglik)  #LRT test
# LRT
# 1-pchisq( LRT,df = 1)  #probability of LRT
# res3.er$AICc
# res3.ard$AICc    #ARD preferred
# res3.ard$solution  #rate transition parameters
#    #transition from 

## ---- two_trait ----


##first plot the data
V1cols<-colorkey1[as.character(mydata[,1])] ##putting the trait values in as.character() is important!
V2cols<-colorkey2[as.character(mydata[,2])]
plot(tree,show.tip.label = F)
tiplabels(pch = 17,col = V1cols)
tiplabels(pch = 19,col = V2cols,offset = 0.2)

##we set up the same format of data frame as before but now we include two traits: trait 1 and trait 2
trtset12<-cbind(row.names(mydata),mydata[,1:2])

trait12_model_er <- getStateMat4Dat(trtset12,"ER")  #equal transition rates
trait12_model_ard <- getStateMat4Dat(trtset12,"ARD")  #All rates different

trait12_model_er
trait12_model_ard


## ----two_trait_fixed ----

##coax the full transition matrix by having a data frame with all trait combinations
## First get all the unique states for the traits in question
states1<-unique(mydata[,1])
states2<-unique(mydata[,2])

traits12_expanded<-expand.grid('species',states1,states2) ##Having some string at the beginning is important! Our data needs the first column to be populated with names, this achieves that

##try making matrices again
trait12_model_er <- getStateMat4Dat(traits12_expanded,"ER")  #equal transition rates
trait12_model_ard <- getStateMat4Dat(traits12_expanded,"ARD")  #All rates different

##The full transition matrix. Hot dog!
trait12_model_er
trait12_model_ard

plotMKmodel(trait12_model_er$rate.mat,rate.cat = 1)
plotMKmodel(trait12_model_ard$rate.mat,rate.cat = 1)


## ----two_trait_pagel-----
tr1<-mydata[,1]; names(tr1)<-row.names(mydata)
tr2<-mydata[,2]; names(tr2)<-row.names(mydata)
pagel_fit <- fitPagel(tree,x=tr1,y=tr2)
pagel_fit

## ---- two_trait_fit ----

##Code is out of date and won't run
# LRT<- -2*(disc.res12.er$loglik - disc.res12.ard$loglik)  #LRT test
# 1-pchisq( LRT,df = 4)  #probability of LRT
# disc.res12.er$AIC
# disc.res12.ard$AIC   #STRONG SUPPORT for correlated (dependent) evolution
#    #Implemented in phytools
# 
# disc.res12.b
#   #NOTE: LOOK AT PHYLOGENY AND DATA HERE!  WE have strong support for correlated evolution, but both traits are co-distributed and clustered
#     #essentially, 1 evolutionary shift could explain the pattern: the 'BiSSE' problem (need ancestral state estimation: next week)
#  ##### ALWAYS PLOT YOUR DATA! DON'T JUST RUN STATISTICS (statistical analysis alone here leads to mis-interpretation)
# 
# #Traits 3&4
# plot.phylo(tree,show.tip.label = F)
# tiplabels(pie = to.matrix(mydata[,3],sort(unique(mydata[,3]))),piecol=c("red", "black"),cex=.3, offset=0)
# tiplabels(pie = to.matrix(mydata[,4],sort(unique(mydata[,4]))),piecol=c("green", "orange"),cex=.3, offset=0.2)
# disc.res34.er<-corDISC(tree,trtset34,ntraits=2,model="ER",node.states="marginal") 
# disc.res34.ard<-corDISC(tree,trtset34,ntraits=2,model="ARD",node.states="marginal") 
# disc.res34.er$loglik
# disc.res34.ard$loglik 
# disc.res34.er$AIC
# disc.res34.ard$AIC # VERY STRONG SUPPORT FOR CORRELATED EVOLUTION


## ----two_trait_order----

trait_model_er <- getStateMat4Dat(traits12_expanded,"ER")  #Transition matrix with equal transition rates for to two binary traits


##We will use the ER transition matrix as a starting point for our directional model
trait_model_dir<-trait_model_er

##Give some of the transitions their own rate
trait_model_dir$rate.mat[1,2]<-2
trait_model_dir$rate.mat[3,4]<-3

trait_model_dir

## ----two_trait_plot----


#Two trait analysis: Order to transitions matters (changes in trait 2 DEPEND on values of trait 1: Maddison 1990)
#Must define the two models for comparison
plot.phylo(tree,show.tip.label = F)
tiplabels(pie = to.matrix(mydata[,3],sort(unique(mydata[,3]))),piecol=c("red", "black"),cex=.3, offset=0)
tiplabels(pie = to.matrix(mydata[,4],sort(unique(mydata[,4]))),piecol=c("green", "orange"),cex=.3, offset=0.2)

##----two_trait_dir_fit----

traits34<-cbind(row.names(mydata),mydata[,3:4])

trait34_fit_er  <-corHMM(phy = tree, data = traits34, rate.cat = 1, rate.mat = trait_model_er$rate.mat)
trait34_fit_dir <-corHMM(phy = tree, data = traits34, rate.cat = 1, rate.mat = trait_model_dir$rate.mat)

trait34_fit_er$loglik
trait34_fit_dir$loglik

summart_stat <- -2*(trait34_fit_er$loglik-trait34_fit_dir$loglik)
1-pchisq(summary_stat,df=2) ##2 degrees of freedom b/c 3-1=2, the difference in the number of parameters 
##high p value means we fail to reject the null model that allows for different transitions of one trait based on being in the other.

trait34_fit_er$AIC
trait34_fit_dir$AIC ##Strong preference for trait34_fit_dir 

