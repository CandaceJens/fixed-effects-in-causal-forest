Overview - These files take Monte Carlo results and plot Figures 1, 2, and 3 from Jens, Page, and Reeder (2021).

There are seven (7) files included.

FILE NAME 		-> 	FILE PURPOSE
complex.continuous.R	-> 	makes right panel of Figure 2

complex.discrete.R	->	makes middle panel of Figure 2

complex.fixed.R		->	makes left panel of Figure 2	

fig3.R			-> 	makes Figure 3

linear.continuous.R	-> 	makes right panel of Figure 1

linear.discrete.R	->	makes middle panel of Figure 1

linear.fixed.R		->	makes left panel of Figure 1

To use the above files, you will need to have the following libraries installed in R:
tidyverse, ggridges

Monte Carlo results needed - should be placed in the "MCresults" folder
	complex.continuous.all.csv
	complex.discrete.all.csv
	complex.fixed.all.csv
	linear.continuous.all.csv
	linear.discrete.all.csv
	linear.fixed.all.csv

Single simulation results needed - should be placed in the "singlesim" folder
	hte_sample.csv

Figures subfolder is empty but included because it is the name of the output folder for each file.
Change the project directory at the beginning of each file to run and save locally.

