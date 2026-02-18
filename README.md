# BayesianAMRPrevalence

This repository houses the code for the academic paper "Bayesian estimation of the prevalence of antimicrobial resistance: a mathematical modelling study" in *Journal of Antimicrobial Chemotherapy* 2024 Sep 3;79(9):2317-2326.  doi: 10.1093/jac/dkae230.

The code contains randomisation elements, and a seed was not set when the original results were generated - however, comparable results can be reproduced by downloading the *microbiologyevents.csv* file in MIMIC-IV version 2.2 from PhysioNet, installing the packages listed in the *Functions* script, replacing anywhere #FILEPATH# is stated with your working directory, and running the scripts in this repository in the following order:

1. Functions.R
2. BEAR precision validation.R
3. OP precision validation.R
4. BEAR resource validation.R
5. OP resource validation.R
6. Plots.R
