# EEOB565X-Spring2020 Course Files Repository

GitHub pages are turned on for this repository. So to link to any file, you just need to append the path to the repo URL [https://eeob-macroevolution.github.io/EEOB565X-Spring2020/](https://eeob-macroevolution.github.io/EEOB565X-Spring2020/). For example: [https://eeob-macroevolution.github.io/EEOB565X-Spring2020/practicals/01-intro-to-phylo/01-intro-to-phylo.html](https://eeob-macroevolution.github.io/EEOB565X-Spring2020/practicals/01-intro-to-phylo/01-intro-to-phylo.html).

## Topics Covered

In this repo you will find all of the materials related to practicals for EEOB565X. The materials are divided by macroevolutionary topic. These topics include:
* [**Onboarding**](https://eeob-macroevolution.github.io/Practicals/Onboarding/Setting_up_R_Rstudio.html) to R, Rsutdio, and the tutorials themselves
* **An Introduction to R and Phylogenetics in R**
* **Phylogenetic Inference**
* **State Dependent Diversification**
* **Phylogenetic Regression**
* **Phylogenetic Association of Discrete Characters**
* **Ancestral State Estimation**
* **Fitting Evolutionary Models**
* **Multivariate Phylogenetic Comparative Methods**

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
