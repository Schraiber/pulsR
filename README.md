# Fitting puncutated and gradual evolutionary models in R

## Installation

For now, just download all the R files into a directory. We will soon implement a devtools install and ultimately submit the package to cran.

## Data format

Phylogenies are expected to be in ape format. Data is expected to be in a named vector, in which names correspond to tips of the phylogeny.

## Example data

Vertebrate data used in Landis and Schraiber (2017) can be found in vertebrate.body\_size.clean.rds

## Fitting data

Data can be fit using the function `fit_reml_levy` in the `levy_pruning_optim.r` script. The main parameters are `phy`, which is an ape format phylogeny; `dat`, which is a named vector of measurements; and `model`, which indicates which model to fit. Possible models are `BM`, `OU`, `JN`, `VG`, `NIG`, `BMJN`, `BMVG`, and `BMNIG`. Additional paremeters are described in the script file. 

This function returns a list with an element `params` indicating the maximum likelihood parameters. Other returned elements are explained in the script file.

## Simulating data

Data can be simulated using the function `rlevy` in `levy_prurning_sim.r`. Simulation requires an ape format phylogeny and a model (same choices as above). Models require parameters; for explanation of parameters, see Landis and Schraiber (2017). 
