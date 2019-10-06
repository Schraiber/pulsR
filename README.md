# Fitting puncutated and gradual evolutionary models in R

## Installation

The `pulsR` package is installed using `devtools`:
```
# install the devtools package (if needed)
install.packages("devtools")

# load the devtools library
library(devtools)

# install the pulsR package
install_github("https://github.com/Schraiber/pulsR")
```

## Library

Once installed, `pulsR` is loaded like other `R` packages:
```
# load the pulsR library
library(pulsR)
```

## Data format

Phylogenies are expected to be in ape format. Data are expected to be in a named vector, in which names correspond to tips of the phylogeny.

## Example data

Vertebrate data used in Landis and Schraiber (2017) can be found in vertebrate.body\_size.clean.rds

## Fitting data

Data can be fit using the function `fit_reml_levy` in the `levy_pruning_optim.r` script. The main parameters are `phy`, which is an ape format phylogeny; `dat`, which is a named vector of measurements; and `model`, which indicates which model to fit. Possible models are `BM`, `OU`, `JN`, `VG`, `NIG`, `BMJN`, `BMVG`, and `BMNIG`. Additional paremeters are described in the script file. 

This function returns a list with an element `params` indicating the maximum likelihood parameters. Other returned elements are explained in the script file.

## Simulating data

Data can be simulated using the function `rlevy` in `levy_prurning_sim.r`. Simulation requires an ape format phylogeny and a model (same choices as above). Models require parameters; for explanation of parameters, see Landis and Schraiber (2017). 
