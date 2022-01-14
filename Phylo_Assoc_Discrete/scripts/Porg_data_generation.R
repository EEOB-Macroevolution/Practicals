library(ape)
library(phytools)
library(corHMM)

set.seed(50014)
tree <- pbtree(n=30)

species<-c('Green','Red','Blue-footed','Wrinkley-necked',"Dean's", "Stubby","Yellow-crested",'Striped', "Sapsucker",'White-bellied',
           'Hoth',"Yavin 4", "Emperor", "Brush-tailed","True","Common","Red-chested", "Lesser",'Soft-billed',"Bearded","Sun-chested",
           'Gregarious','Horned',"Purple-throated","Oasis",'Short-tailed',"Long-tailed","Buff","Marbled","Whistling")
species<-paste(species,'porg', sep=' ')
tree$tip.label<-species

temp<-cbind(species,0:1,c(1,1,0)) ## get a dummy data frame to create the rate matrix
rate_mat<-getStateMat4Dat(temp,'ER')$rate.mat
rate_mat[2,1]<- rate_mat[4,2]<-0.3 ##seldom have have an unwebbed beakless porg -> they have nothing to eat on land
rate_mat[1,3]<- rate_mat[2,4] <-2 ##make it easy to go from unwebbed to webbed -> more niches to exploit
rate_mat[3,4] <-1.5 ##Lots of fish to be speared -> good strategy for survival

#Make rows sum to zero
for(rw in 1:4){
  rate_mat[rw,rw]<- -sum(rate_mat[rw,]) ##make the diagonal the negative of the rowsum
}
rate_mat
trait_sim<-sim.char(tree,rate_mat,model='discrete') ##simulate characters
state1<- as.numeric(trait_sim %in% c(3,4))
state2<- as.numeric(trait_sim %in% c(2,4))

my_data<-as.data.frame(cbind(species,state1,state2))

colnames(my_data)<-c('species','webbed feet','beak')

write.csv(x=my_data,file='../data/porg_data.csv')
write.tree(tree,file='../data/porg.tre')
