# Evo-Demo Hyperstate Matrix Model Code Repository

This repository contains code associated with the manuscript "Toward a unified approach to modeling adaptation among demographers and evolutionary ecologists" submitted to Methods in Ecology and Evolution (MEE). The repository was created by S. Jenouvrier and includes code for computing the IBM and Evo-Demo Hyperstate Matrix (EvoDemo-Hyper MPM).

## Overview
 
The EvoDemo-Hyper MPM model introduces a new approach that incorporates the genetic inheritance of quantitative traits, enabling rapid computation of evolutionary and demographic dynamics. The model is evaluated against individual-based simulations (IBM) and provides analytical approximations for adaptation rates across six distinct scales in response to selection. The study demonstrates that, while different modeling approaches yield equivalent results for the same biological scenario, semantic distinctions between fields may sometimes obscure these similarities.

## Contents

R Code (IBM and Theoretical Derivations): Developed by Timothée Bonnet, the R scripts include code to run IBM simulations and compute theoretical derivations.
MATLAB Code (EvoDemo-Hyper MPM): Developed by Stephanie Jenouvrier, Jimmy Garnier, and Joanie Van de Walle, the MATLAB scripts implement the EvoDemo-Hyper MPM model and calculate adaptation rates from theoretical formulas.

## File Descriptions

**MATLAB Code Files**

* `CompRA_THEO_definition_methods_GITHUB.m`: Main MATLAB script to calculate the rate of adaptation using analytical theoretical formulas.
* `Emat.m` and `BD_proj_mat.m`: Functions for matrix operations used within the EvoDemo-Hyper MPM framework.
* `Main_EvoDemo_Hyper_MPM_GITHUB.m`: Main function for executing the EvoDemo-Hyper MPM simulations.

**R Code Files**

* `IBM_Function.R`: R functions to run the individual-based model (IBM) simulations and do post treatment of simulation output.
* `running_IBM.Rmd` and compiled version `running_IBM.html`: Run the individual-based model (IBM) simulations across vital rates and species to reproduce the results presented in the manuscript.
* `MPM_Functions.R`: R functions for executing the EvoDemo-Hyper MPM simulations (translated from MATLAB).
* `running_MPM.Rmd` and compiled version `running_MPM.html`: Demonstration of running EvoDemo-Hyper MPM simulations. 
* `TheoreticalDerivations.R`: R script to compute theoretical derivations (translated from MATLAB) as described in the manuscript.

**Miscellaneous files**

* `Eco-EvoHyperstateModel.Rproj`: R-Studio project properties (not essential)
* `.gitignore`: file extensions to be ignored by git (not essential)

## How to Use

**IBM Simulations (R)**: Use the R script `running_IBM.Rmd` to run IBM simulations. These scripts also include theoretical derivations `TheoreticalDerivations.R` for computing analytical adaptation rates.

**EvoDemo-Hyper MPM (MATLAB)**: Run `CompRA_THEO_definition_methods_GITHUB.m` in MATLAB to compute analytical adaptation rates. 
Run `Main_EvoDemo_Hyper_MPM_GITHUB.m` to compute adaptation rates using the EvoDemo-Hyper MPM.  Ensure that `Emat.m` and `BD_proj_mat.m` are in the same directory or accessible via your MATLAB path.

## Authors

Timothée Bonnet (R code for IBM simulations and theoretical derivations); Stéphanie Jenouvrier, Jimmy Garnier, Joanie Van de Walle (MATLAB code for EvoDemo-Hyper MPM)

## Citation

If you use this code, please cite the corresponding manuscript once it is published.

