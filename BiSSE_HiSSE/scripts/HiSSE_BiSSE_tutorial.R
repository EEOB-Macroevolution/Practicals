# BiSSE and HiSSE in R
# Macroevolution -- EEOB 565


##set the seed for reproducibility
set.seed(9177207) ##The exponent of the largest known repunit prime. ((10^9177207)-1)/9 is the prime

# load packages
library(diversitree)
library(hisse)

# Simulate tree
sim_parameters <- c(1.0, 4, 0.7, 0.3, 0.2, 0.8) # True parameter values
names(sim_parameters)<-c('lambda0','lambda1','mu0','mu1','q01','q10')
tree_bisse <- tree.bisse(sim_parameters, max.taxa = 100)

# View details
tree_bisse

names(tree_bisse) ##See which components exist for our simulated phylogeny.

# Display the states for each tip
tree_bisse$tip.state

# Plot the character state history 
treehistory = history.from.sim.discrete(tree_bisse, 0:1)
plot(treehistory, tree_bisse, cols = c("red", "blue"))

# View the true parameter values
sim_parameters

# Create a bisse-model variable with our tree and tip states
bisse_model = make.bisse(tree_bisse, tree_bisse$tip.state)
bisse_model

# Compute the likelihood under different parameter values
## Parameter order: lambda0, lambda1, mu0, mu1, q01, q10 
bisse_model(c(1.0, 4.0, 0.5, 0.5, 0.2, 0.8))
bisse_model(c(4.0, 1.0, 0  , 1  , 0.5, 0.5))

# Use greedy search algorithm to find optimum
initial_pars<-starting.point.bisse(tree_bisse) # provide starting values
fit_model<-find.mle(bisse_model,initial_pars) # run search algorithm

# View the log likelihood of the optimal model
fit_model$lnLik

# View the optimal parameters
round(fit_model$par,digits=2) ##we round to two digits for easier reading

# Compare to true parameters
sim_parameters

# Create an alternative model w/ equal trait-associated extinction rates  
constrained_bisse_model <-constrain(bisse_model,mu0 ~ mu1)
# create a set of starting parameters that matches constrained model
constrained_initial_pars<-initial_pars[-3] 
# search for optimum
fit_constrained_model <- find.mle(constrained_bisse_model,constrained_initial_pars)
round(fit_constrained_model$par,digits = 2) ##Round estimates to 2 digits for simplicity

# Compare the fit of the constrained and unconstrained models
anova(fit_model,constrained=fit_constrained_model)

# Simulate an independent binary trait on our tree 
transition_rates<-c(0.5,1)
secondary_trait <- sim.character(tree = tree_bisse, pars =transition_rates , x0 = 1, model='mk2')

# Initialize a bisse model for the independent trait
# in this case, the bisse model is not the true model of evolution for this trait
incorrect_bisse_model <- make.bisse(tree_bisse, secondary_trait)
incorrect_bisse_mle <- find.mle(incorrect_bisse_model, initial_pars)

# Set up a null model where the speciation and extinction rates are the same for each trait
incorrect_null_bisse_model <- constrain(incorrect_bisse_model, lambda0 ~ lambda1)
incorrect_null_bisse_model <- constrain(incorrect_null_bisse_model, mu0 ~ mu1)
constrained_initial_pars<-initial_pars[c(-1,-3)]
incorrect_null_bisse_mle <- find.mle(incorrect_null_bisse_model, constrained_initial_pars)

# Compare the two models
anova(incorrect_bisse_mle, constrained = incorrect_null_bisse_mle)

#### Apply HiSSE model ####

# Construct a transition rate matrix 
null_rate_matrix <- TransMatMakerHiSSE(hidden.traits = 1,make.null=T) #  hidden.traits=1 specifies the hisse model
null_rate_matrix

# the hisse package paramterizes the model using transformed parameters
null_net_turnover <- c(1,1,2,2) # lambda + mu
null_extinction_fraction <- c(1,1,2,2) # mu/lambda

# Transform the tip state data to the hisse format
hisse_states <- cbind(names(secondary_trait), secondary_trait)

# perform ML estimation under the hisse null model (note that it takes some time to run this)
null_hisse <- hisse(phy = tree_bisse, data = hisse_states, hidden.states=TRUE,
                   turnover = null_net_turnover, eps = null_extinction_fraction,
                   trans.rate = null_rate_matrix)
null_hisse

# compare models
null_hisse$AIC
AIC(incorrect_bisse_mle)
AIC(incorrect_null_bisse_mle)

# specify a hisse model based on the DEPENDENT trait (our first trait)
hisse_states2 <- cbind(names(tree_bisse$tip.state), tree_bisse$tip.state)
null_hisse2 <- hisse(phy = tree_bisse, data = hisse_states2, hidden.states=TRUE, 
                   turnover = null_net_turnover, eps = null_extinction_fraction,
                   trans.rate = null_rate_matrix)

# Compare models
null_hisse2$AIC
AIC(fit_model)
AIC(fit_constrained_model)

