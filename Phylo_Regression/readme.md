

###Contents
This file contains the contents for the Phylogenetic Regression

You will find the tutorial labeled PhyloRegressionhtml here. 


You will find the following subdirectories here:
+ **scripts** This folder contains all the scripts used to generate data, code presented in the tutorial, and markdown for the tutorial itself
+ **data**  This folder contains all the data used in the tutorial for analysis
+ **figs** This folder contains the figures used in the tutorial



## Going thru the tutorial

If following along with the code in this tutorial, we assume the working directory is the same as the tutorial itself (i.e. `/Tutorial_Skeleton`)
You will need the following R packages to run the tutorial:


ape
geiger
nlme
RRPP
geomorph
MASS

**Note:** If when plotting you get an error regarding not being able to find the function `getYMult()`, then you may need to additionally install the `Plotrix` package.


## Year to year maintenance 

1. Update the packages to their latest versions
2. try knitting the tutorial again
    + If things break try fixing them 
    + If no reasonable fix exists, find package versions where the tutorial is stable and make note of the version above 