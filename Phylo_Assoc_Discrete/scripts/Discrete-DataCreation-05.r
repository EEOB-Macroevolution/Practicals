### CREate data


n<-64
tree<-compute.brlen(stree(n,type="balanced"),method=rep(1,n))
plot(tree, show.tip.label = FALSE)
write.tree(tree,'tree.64.tre')

#create some data
traits<-cbind(rep(c("0","1"), each=n/2), rep(c(0,1), each=n/2), rep(c("0","1"), 32), rep(c(0,1), 32),c(rep(0:1, 4), rep(0,(n-9)), 1),sample(0:1, n, replace = T))

row.names(traits)<-tree$tip.label
write.csv(traits,"DiscreteData.csv")
