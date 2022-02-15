#Load Needed packages to knit rmarkdown
library(rmarkdown)
library(knitr)

##Set working directory to this file's location!

print('building Onboarding Tutorial')
render(input = "../Onboarding_Rintro/Setting_up_R_Rstudio.Rmd",quiet=T)

print('building Intro R Tutorial')
render(input = "../Onboarding_Rintro/Intro_to_R.Rmd",quiet=T)

print('building Ancestral State Estimation Tutorial')
render(input = "../Ancestral_State_Estimation/scripts/AncStateEstimation_Tutorial.Rmd",
       output_dir = "/Ancestral_State_Estimation",quiet = T)

print('building BiSSE/HiSSE Tutorial')
render(input = "../BiSSE_HiSSE/scripts/HiSSE_BiSSE_Tutorial.Rmd",
       output_dir = "../BiSSE_HiSSE",quiet = T)

print('building Fitting Evolutionary Models Tutorial')
render(input = "../Fitting_Evol_Models/scripts/Fit_Evol_Models_Tutorial.Rmd",
       output_dir = "../Fitting_Evol_Models",quiet = T)

print('building Introduction to Phylogenetics Tutorial')
render(input = "../Intro_to_Phylo/scripts/intro_to_phylo.Rmd",
       output_dir = "../Intro_to_Phylo",quiet = T)

print('building Multivariate PCM Tutorial')
render(input = "../MultivariatePCMs/scripts/MultivariatePCM_Tutorial.Rmd",
       output_dir = "../MultivariatePCMs",quiet = T)

print('building Discrete Phylogenetic Association Tutorial')
render(input = "../Phylo_Assoc_Discrete/scripts/Phylo_Assoc_Discrete_Tutorial.Rmd",
       output_dir = "../Phylo_Assoc_Discrete",quiet = T)

print('building Phylogenetic Regression Tutorial')
render(input = "../Phylo_Regression/scripts/Phylo_Regression_Tutorial.Rmd",
       output_dir = "../Phylo_Regression",quiet = T)

