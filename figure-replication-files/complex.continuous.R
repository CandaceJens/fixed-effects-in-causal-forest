rm(list=ls())


##Loads required packages 

require(tidyverse)
require(ggridges)


##Setting directory
project.dir <- "C:/Users/cjens/Desktop/ML FE/MCplots"
##Set project directory here


setwd(project.dir)



###### Load Data ######

data.first <- read.csv(paste(project.dir, "MCresults/complex.continuous.all.csv", sep="/"), 
	header = TRUE)

###What we want:
##at.reg.no.fe
##at.reg.fe
##base.ate
##cluster.ate
##ols.fe.ate
##lasso.fe.ate
##mean.ate
##demeaned.ate
##part.demeaned.ate



data.complex.continuous.1 <- data.first %>% select(at.reg.no.fe, at.reg.fe,
		base.ate, cluster.ate, ols.fe.ate, lasso.fe.ate, mean.ate, demeaned.ate, 
		part.demeaned.ate) 

data.complex.continuous.1 <- (data.complex.continuous.1 - data.first$mean_true_effect) / data.first$mean_true_effect

		
data.complex.continuous.1.melted <- data.complex.continuous.1 %>% 
				gather(data.complex.continuous.1)



model.names <- c("OLS without FE", "OLS with FE", "Causal forest",
	"Cluster-robust causal forest", "Causal forest with OLS FE",
	"Causal forest with LASSO FE", "Causal forest with mean y",
	"Causal forest with demeaned data", "Causal forest with demeaned y-hat")

###Plot:

Fig2PanelR <- data.complex.continuous.1.melted %>%
 	mutate(data.complex.continuous.1 = fct_relevel(data.complex.continuous.1, c("at.reg.no.fe", "at.reg.fe",
		"base.ate", "cluster.ate", "ols.fe.ate", "lasso.fe.ate", "mean.ate", "demeaned.ate", 
		"part.demeaned.ate"))) %>%
	ggplot(aes(y = data.complex.continuous.1, x = value)) + 
	stat_density_ridges(quantile_lines = TRUE,
                      quantile_fun=function(x,...)mean(x), alpha = 0.5) +
	xlab("bias") +
	xlim(-0.5, 0.5) + 
	scale_y_discrete(labels = model.names) + 
	theme_minimal() +
	theme(legend.position = "bottom",
        legend.text=element_text(margin=margin(r=0.5,unit="inch"))) + 
	ggtitle("Densities of bias for complex functional form \n
			Continuous treatment effects") +
	theme(plot.title = element_text(hjust = 0.5, size = 12, face="bold"),
		axis.title.y = element_blank(),
		axis.text.y = element_blank(),
		axis.text.x = element_text(size = 12),
		axis.title.x = element_text(size = 12)) +	
	geom_vline(xintercept = 0, linetype = 2) 


Fig2PanelR


###Save:
ggsave(Fig2PanelR, filename = paste(project.dir, "figures/Fig2PanelR.pdf", sep="/"),
       width = 3.6, height = 6, units = "in")








