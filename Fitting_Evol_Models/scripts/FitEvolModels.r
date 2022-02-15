#Evolutionary models of continuous trait evolution in R

## ----read_data----
library(phytools)
library(geiger)
library(OUwie)

tree<-read.tree("../data/anole.gp.tre",tree.names=T)
group<-read.csv('../data/anole.gp.csv', row.names=1, header=TRUE,colClasses=c('factor'))
gp<-as.factor(t(group)); names(gp)<-row.names(group)
svl<-read.csv('../data/anole.svl2.csv', row.names=1, header=TRUE)
svl<-as.matrix(treedata(phy = tree,data = svl, warnings=FALSE)$data)[,1]  #match data to tree

####----plot_data----

##BM plot data
tree.col<-contMap(tree,svl,plot=FALSE)  #runs Anc. St. Est. on branches of tree
plot(tree.col,legend=0.7*max(nodeHeights(tree)),
     fsize=c(0.7,0.9))

#PLOT data
tree.col<-contMap(tree,svl,plot=FALSE)  #runs Anc. St. Est. on branches of tree
plot(tree.col,type="fan",legend=0.7*max(nodeHeights(tree)),
     fsize=c(0.7,0.9))
cols<-setNames(palette()[1:length(unique(gp))],sort(unique(gp)))
tiplabels(pie=model.matrix(~gp-1),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.8*par()$usr[1],
                  y=-max(nodeHeights(tree)-.6),fsize=0.8)

## ----BM_sim----
nsim <- 100
t <- 0:100  # time
sig2 <- 0.005
X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2)), nsim, length(t) - 1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))
plot(t, X[1, ], xlab = "time", ylab = "phenotype", ylim = c(-2, 2), type = "l")
apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)

## ----fit_BM----

fit.BM1<-fitContinuous(tree, svl, model="BM")  #Brownian motion model
fit.BM1

##----OU_sim----

#from Lars Shmitz tutorial: http://schmitzlab.info/BMandOU.html
OU.sim <- function(n, theta, alpha, sigma,x0){
  dw  <- rnorm(n, 0)
  x <- c(x0)
  for (i in 2:(n+1)) {
    x[i]  <-  x[i-1] + alpha*(theta-x[i-1]) + sigma*dw[i-1]
  }
  return(x);
}
OU1.sim1 <- replicate(100, OU.sim(n=100, theta=0.5, alpha=0.5, sigma=0.03, x0=0.5), simplify=FALSE) #No change
OU1.sim2 <- replicate(100, OU.sim(n=100, theta=0.75, alpha=0.5, sigma=0.03, x0=0.25), simplify=FALSE) #start and theta differ
par(mfcol = c(1, 2))
plot(OU1.sim1[[1]],xlab = "time", ylab = "phenotype", ylim = c(0,1), type = "l")
for(i in 2:100){lines(OU1.sim1[[i]])}
plot(OU1.sim2[[1]],xlab = "time", ylab = "phenotype", ylim = c(0,1), type = "l")
for(i in 2:100){lines(OU1.sim2[[i]])}
par(mfcol = c(1,1))

## ----OU_fit----
options(warn=-1)
fit.OU1<-fitContinuous(tree, svl,model="OU")    #OU1 model
fit.OU1


## ----BM_OU_compare----
LRT<- -(2*(fit.BM1$opt$lnL-fit.OU1$opt$lnL))
prob<-pchisq(LRT, 1, lower.tail=FALSE)
LRT
prob

fit.BM1$opt$aic
fit.OU1$opt$aic


## ----evol_models_geiger----
fit.BMtrend<-fitContinuous(tree, svl, model="trend")   #Brownian motion with a trend
fit.EB<-fitContinuous(tree, svl, model="EB")   #Early-burst model
fit.lambda<-fitContinuous(tree, svl,
        bounds = list(lambda = c(min = exp(-5), max = 2)), model="lambda")  #Lambda model
options(warn=-1)
fit.K<-fitContinuous(tree, svl, model="kappa")   #Early-burst model

#Examine AIC
c(fit.BM1$opt$aic,fit.BMtrend$opt$aic,fit.EB$opt$aic,fit.lambda$opt$aic,fit.OU1$opt$aic)
  #NOTE: none of these generate dAIC > 4.  So go with simplest model (BM1)



## ----OUwie_model_basic----

data<-data.frame(Genus_species=names(svl),Reg=gp,X=svl)  #input data.frame for OUwie
tree.simmap<-make.simmap(tree,gp)  # perform & plot stochastic maps (we would normally do this x100)
plot(tree.simmap,type='fan')
fitBM1<-OUwie(tree.simmap,data,model="BM1",simmap.tree=TRUE)
fitOU1<-OUwie(tree.simmap,data,model="OU1",simmap.tree=TRUE)
fitOUM<-OUwie(tree.simmap,data,model="OUM",simmap.tree=TRUE)

fitBM1
fitOU1
fitOUM  #OUM is strongly preferred (examine AIC)

#How it SHOULD be run
#trees.simmap<-make.simmap(tree = tree,x = gp,nsim = 100)  # 100 simmaps
#fitOUM.100<-lapply(1:100, function(j) OUwie(trees.simmap[[j]],data,model="OUM",simmap.tree=TRUE))
#OUM.AIC<-unlist(lapply(1:100, function(j) fitOUM.100[[j]]$AIC))
#hist(OUM.AIC)
#abline(v=fitBM1$AIC, lwd=2) #add value for BM1


##----BM_many_rates----

# 3A: BM1 vs BMM: Comparing evolutionary rates
tree.simmap<-make.simmap(tree,gp)  # perform & plot stochastic maps (we would normally do this x100)
BMM.res<-OUwie(tree.simmap,data,model="BMS",simmap.tree = TRUE)
BMM.res


## ----Rate_shifts----
# 3B: Identifying rate shifts on phylogeny
  #Bayesian MCMC: single rate shift (Revell et al. 2012: Evol.)
BM.MCMC<-evol.rate.mcmc(tree=tree,x=svl, quiet=TRUE)
post.splits<-minSplit(tree,BM.MCMC$mcmc)  #summarize
MCMC.post<-posterior.evolrate(tree=tree,mcmc=BM.MCMC$mcmc,tips = BM.MCMC$tips,ave.shift = post.splits)
   #Plot rescaled to rates
ave.rates(tree,post.splits,extract.clade(tree, node=post.splits$node)$tip.label,
    colMeans(MCMC.post)["sig1"], colMeans(MCMC.post)["sig2"],post.splits)

   #Reversible-jump MCMC: multiple rate shifts (Eastman et al. 2011: Evol.)
BM.RJMC<-rjmcmc.bm(phy = tree,dat = svl)
#Run again for plotting
r <- paste(sample(letters,9,replace=TRUE),collapse="")
rjmcmc.bm(phy=tree, dat=svl, prop.width=1.5, ngen=20000, samp=500, filebase=r, simple.start=TRUE, type="rbm")
outdir <- paste("relaxedBM", r, sep=".")
ps <- load.rjmcmc(outdir)
plot(x=ps, par="shifts", burnin=0.25, legend=TRUE, show.tip=FALSE, edge.width=2,type='fan')



## ----lecture_lambda----

#1: what does lambda do?
TreeLambda0 <- rescale(tree, model = "lambda", 0)
TreeLambda5 <- rescale(tree, model = "lambda", 0.5)

par(mfcol = c(1, 3))
plot(tree, show.tip.label = FALSE)
plot(TreeLambda5,show.tip.label = FALSE)
plot(TreeLambda0,show.tip.label = FALSE)
par(mfcol=c(1,1))


