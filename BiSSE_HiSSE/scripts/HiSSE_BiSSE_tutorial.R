# BiSSE and HiSSE in R
# Macroevolution -- EEOB 565


##set the seed for reproducibility
set.seed(9177207) ##The exponent of the largest known repunit prime. ((10^9177207)-1)/9 is the prime

# load packages
library(diversitree)
library(hisse)

# Simulate tree
sim_parameters <- c(1.0, 4, 0.7, 0.3, 0.2, 0.8)
names(sim_parameters)<-c('lambda0','lambda1','mu0','mu1','q01','q10')
tree_bisse <- tree.bisse(sim_parameters, max.taxa = 100)

tree_bisse

names(tree_bisse) ##See which components exist for our simulated phylogeny.

tree_bisse$tip.state


treehistory = history.from.sim.discrete(tree_bisse, 0:1)
plot(treehistory, tree_bisse, cols = c("red", "blue"))

sim_parameters

bisse_model = make.bisse(tree_bisse, tree_bisse$tip.state)
bisse_model

bisse_model(c(1.0, 4.0, 0.5, 0.5, 0.2, 0.8))
bisse_model(c(4.0, 1.0, 0  , 1  , 0.5, 0.5))

initial_pars<-starting.point.bisse(tree_bisse)
fit_model<-find.mle(bisse_model,initial_pars)

fit_model$lnLik
round(fit_model$par,digits=2) ##we round to two digits for easier reading

sim_parameters

constrained_bisse_model <-constrain(bisse_model,mu0 ~ mu1)

constrained_initial_pars<-initial_pars[-3]
fit_constrained_model <- find.mle(constrained_bisse_model,constrained_initial_pars)
round(fit_constrained_model$par,digits = 2) ##Round estimates to 2 digits for simplicity

anova(fit_model,constrained=fit_constrained_model)

##simulate binary trait data on the tree
transition_rates<-c(0.5,1)
secondary_trait <- sim.character(tree = tree_bisse, pars =transition_rates , x0 = 1, model='mk2')


incorrect_bisse_model <- make.bisse(tree_bisse, secondary_trait)
incorrect_bisse_mle <- find.mle(incorrect_bisse_model, initial_pars)


incorrect_null_bisse_model <- constrain(incorrect_bisse_model, lambda0 ~ lambda1)
incorrect_null_bisse_model <- constrain(incorrect_null_bisse_model, mu0 ~ mu1)
constrained_initial_pars<-initial_pars[c(-1,-3)]
incorrect_null_bisse_mle <- find.mle(incorrect_null_bisse_model, constrained_initial_pars)

anova(incorrect_bisse_mle, constrained = incorrect_null_bisse_mle)

null_rate_matrix <- TransMatMakerHiSSE(hidden.traits = 1,make.null=T)
null_rate_matrix


null_net_turnover <- c(1,1,2,2)
null_extinction_fraction <- c(1,1,2,2)


hisse_states <- cbind(names(secondary_trait), secondary_trait)

null_hisse <- hisse(phy = tree_bisse, data = hisse_states, hidden.states=TRUE,
                   turnover = null_net_turnover, eps = null_extinction_fraction,
                   trans.rate = null_rate_matrix)
null_hisse

null_hisse$AIC
AIC(incorrect_bisse_mle)
AIC(incorrect_null_bisse_mle)

hisse_states2 <- cbind(names(tree_bisse$tip.state), tree_bisse$tip.state)
null_hisse2 <- hisse(phy = tree_bisse, data = hisse_states2, hidden.states=TRUE, 
                   turnover = null_net_turnover, eps = null_extinction_fraction,
                   trans.rate = null_rate_matrix)

null_hisse2$AIC
AIC(fit_model)
AIC(fit_constrained_model)

