library(ape)
library(phytools)
library(corHMM)

set.seed(50011)
tree <- pbtree(n=30)

species<-c('Burley','Black','Spiny','Winkly',"Soft", "Stubby","Yellow-bellied",'Striped', "Jen's",'Two-toed',
           'Puffy',"Light", "Coarse", "Star-nosed","Mountain","Common","Pygmy", "Pocket",'Polar',"Bearded","Giant",
           'Shy','Brown',"Flying","Cave",'Spotted',"Tree","Stout","Blue-eyed","Red")
species<-paste(species,'fizzgig', sep=' ')
tree$tip.label<-species

temp<-cbind(species,0:1) ## get a dummy data frame to create the rate matrix
rate_mat<-getStateMat4Dat(temp,'ER')$rate.mat


#Make rows sum to zero
for(rw in 1:2){
  rate_mat[rw,rw]<- -sum(rate_mat[rw,]) ##make the diagonal the negative of the rowsum
}
rate_mat
disc_trait_sim<-sim.char(tree,rate_mat,model='discrete')[,,1] ##simulate characters
disc_trait_sim[disc_trait_sim==1]<-'insectivore'
disc_trait_sim[disc_trait_sim==2]<-'aviovore'
cont_trait_sim<-sim.char(tree,par=0.2)[,,1]


my_data<-as.data.frame(cbind(species,disc_trait_sim,cont_trait_sim))

colnames(my_data)<-c('species','type','bite_radius')

write.csv(x=my_data,file='../data/fizzgig_data.csv',row.names = F)
write.tree(tree,file='../data/fizzgig.tre')
