# EEOB565X- Practical Course Files Repository

This repository contains the tutorials for Macroevolution (BIOL 465 / EEOB 565), a graduate and advanced undergraduate course at Iowa State University. 

## Topics Covered

In this repo you will find all of the materials related to practicals for EEOB565X. The materials are divided by macroevolutionary topic. These topics include:
* [**Onboarding**](https://eeob-macroevolution.github.io/Practicals/Onboarding_RIntro/Setting_up_R_Rstudio.html) to R, Rsutdio, and the tutorials themselves
* [**An Introduction to R**](https://eeob-macroevolution.github.io/Practicals/Onboarding_RIntro/intro_to_R.html)
* [**An Introduction to Phylogenetics in R**](https://eeob-macroevolution.github.io/Practicals/Intro_to_Phylo/intro_to_phylo.html)
* [**Phylogenetic Association of Discrete Characters**](https://eeob-macroevolution.github.io/Practicals/Phylo_Assoc_Discrete/Phylo_Assoc_Discrete_Tutorial.html)
* [**Phylogenetic Regression**](https://eeob-macroevolution.github.io/Practicals/Phylo_Regression/Phylo_Regression_Tutorial.html)
* [**Ancestral State Estimation**](https://eeob-macroevolution.github.io/Practicals/Ancestral_State_Estimation/AncStateEstimation_Tutorial.html)
* [**State Dependent Diversification**](https://eeob-macroevolution.github.io/Practicals/BiSSE_HiSSE/HiSSE_BiSSE_tutorial.html)
* [**Fitting Evolutionary Models**](https://eeob-macroevolution.github.io/Practicals/Fitting_Evol_Models/Fit_Evol_Models_Tutorial.html)
* [**Multivariate Phylogenetic Comparative Methods**](https://eeob-macroevolution.github.io/Practicals/MultivariatePCMs/MultivariatePCM_Tutorial.html)

### Folder Structure

In each folder you should find the following files:

* **Macroevolutionary_Topic_Tutorial.html** This file will contain the tutorial for the given unit, this is what you'll be using the most.
* **readme.md** This is a text document that details the folders contents and any specific notes about the unit. It will also list the required R package needed for the tutorial, so this might be good to check before looking the tutorial.

You will also see a few different subfolders: 

* **scripts** This folder contains all the scripts used to generate data, run analyses, and create the tutorial. All of the code you see in the html tutorial can be found in the file that ends in `_code.r`. You will not find any code in the file ending in `.rmd`, that is the file used to create the html.
* **data** This folder will contain all the needed data files, usually trait data matrices and phylogenetic trees. We will cover in the first few tutorials on how to read these data into R.
* **figs** This folder contains any figures used in the tutorial. You likely won't need to access this folder unless you want to reference some of the figures. 

If new tutorials are added you can copy `Tutorial_Skeleton` for a blank copy of the folder structure

----

## Other Notes

Be sure to check `maintenance_notes.text` at the beginning of each semester for code maintenance, package updating, etc. 

If new tutorials are made, be sure to check `Onboarding/Tutorial_Onboarding.rmd` for conventions and example code for the conventions 
