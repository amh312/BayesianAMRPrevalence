# BayesianAMRPrevalence

This repository houses the code for the academic paper "Bayesian estimation of the prevalence of antimicrobial resistance: a mathematical modelling study" in *Journal of Antimicrobial Chemotherapy* 2024 Sep 3;79(9):2317-2326.  doi: 10.1093/jac/dkae230.

The code contains randomisation elements, and a seed was not set when the original results were generated - however, comparable results can be reproduced by downloading the *microbiologyevents.csv* file in MIMIC-IV version 2.2 from *PhysioNet*, installing the packages listed in the *Functions.R* script, replacing anywhere #FILEPATH# is stated with your working directory, and running the scripts in this repository in the following order:

1. Functions.R
2. BEAR precision validation.R
3. OP precision validation.R
4. BEAR resource validation.R
5. OP resource validation.R
6. Plots.R



## Notes

1. To reduce computation time in the resource analyses (*BEAR resource validation.R* and *OP resource validation.R*) when reproducing the code, starting numbers (argument *startingnum* in function *number_validate()*) of samples informing the likelihood have been started at the values from the manuscript results. To reproduce the time taken to complete the full analysis, these numbers should be changed to 2.
2. The resource validation scripts may sporadically throw non-reproducible errors. If such an error occurs, re-running the code from the block that generated the error will usually resolve the problem. Examples of such errors that have been observed are: "Error in 'filter()': Can't specify an argument named 'by' in this verb. Did you mean to use '.by' instead?;"Error in 'arg_match()': ! 'arg' must be a symbol, not a function."; "Error in UseMethod('depth') :  no applicable method for 'depth' applied to an object of class 'NULL'").
