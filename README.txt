Overview - A function to extract fixed effects for use in a causal forest and example files.

There are four (4) files included.

FILE NAME 		-> 	FILE PURPOSE
fe_extract.R 	->	includes the fe.extract function that extracts group-level fixed effect estimates using LASSO

avg_tes.R 		->	includes the avg.tes function that estimates ATE, ATC, ATT, and ATO both with and without one-way 
						clustered standard errors
							
simulate_data.R ->	includes the simulate.data function that simulates an unbalanced panel of data with one group 
						identifier along with helper functions for the simulation
							
example.R 		->	includes sample code that shows how to use the fixed effect estimates from the fe.extract function 
						as part of a causal forest and estimate clustered standard errors using the avg.tes function on 
						simulated data

To use the above files, you will need to have the following libraries installed in R:
Matrix, glmnet, dplyr, grf, and sandwich