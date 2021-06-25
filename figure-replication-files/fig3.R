rm(list=ls())


##Loads required packages 

require(tidyverse)
require(ggridges)



##Setting directory
project.dir <- "C:/Users/cjens/Desktop/ML FE/MCplots"
####Set a project directory here


setwd(project.dir)




multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



###### Load Data ######

data.subset <- read.csv(paste(project.dir, "singlesim/hte_sample.csv", sep="/"), 
	header = TRUE)

###List of estimators:
names(data.subset)

##base.hte  (drop this one for a 2x3)
##cluster.hte
##ols.fe.hte
##lasse.fe.hte
##mean.hte
##demeaned.hte
##part.demeaned.hte 

base.hte.plot <- ggplot(data.subset, aes(x = base.hte, y = true_effect)) +
		geom_point(shape = 4, colour = "#56B4E9") +
		geom_abline(intercept = 0, slope = 1, color="gray50", 
                 size=1) +
		theme_minimal() +
		geom_smooth(method = lm, linetype="dashed", se = FALSE, colour = "black") +
		xlab("estimated treatment effect") +
		xlim(0, 1.5)  +
		ylab("true treatment effect") +
		ylim(0, 1.5)  +
		ggtitle("Causal forest") +
		theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
		axis.text.y = element_text(size = 10),
		axis.title.y = element_text(size = 10),
		axis.text.x = element_text(size = 10),
		axis.title.x = element_text(size = 10))
	
base.hte.plot


cluster.hte.plot <- ggplot(data.subset, aes(x = cluster.hte, y = true_effect)) +
		geom_point(shape = 4, colour = "#56B4E9") +
		geom_abline(intercept = 0, slope = 1, color="gray50", 
                 size=1) +
		theme_minimal() +
		geom_smooth(method = lm, linetype="dashed", se = FALSE, colour = "black") +
		xlab("estimated treatment effect") +
		xlim(0, 1.5)  +
		ylab("true treatment effect") +
		ylim(0, 1.5)  +
		ggtitle("Cluster-robust causal forest") +
		theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
		axis.text.y = element_text(size = 10),
		axis.title.y = element_text(size = 10),
		axis.text.x = element_text(size = 10),
		axis.title.x = element_text(size = 10))
	
cluster.hte.plot



ols.fe.hte.plot <- ggplot(data.subset, aes(x = ols.fe.hte, y = true_effect)) +
		geom_point(shape = 4, colour = "#56B4E9") +
		geom_abline(intercept = 0, slope = 1, color="gray50", 
                 size=1) +
		theme_minimal() +
		geom_smooth(method = lm, linetype="dashed", se = FALSE, colour = "black") +
		xlab("estimated treatment effect") +
		xlim(0, 1.5)  +
		ylab("true treatment effect") +
		ylim(0, 1.5)  +
		ggtitle("Causal forest with OLS fixed effect") +
		theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
		axis.text.y = element_text(size = 10),
		axis.title.y = element_text(size = 10),
		axis.text.x = element_text(size = 10),
		axis.title.x = element_text(size = 10))
	
ols.fe.hte.plot


lasso.fe.hte.plot <- ggplot(data.subset, aes(x = lasso.fe.hte, y = true_effect)) +
		geom_point(shape = 4, colour = "#56B4E9") +
		geom_abline(intercept = 0, slope = 1, color="gray50", 
                 size=1) +
		theme_minimal() +
		geom_smooth(method = lm, linetype="dashed", se = FALSE, colour = "black") +
		xlab("estimated treatment effect") +
		xlim(0, 1.5)  +
		ylab("true treatment effect") +
		ylim(0, 1.5)  +
		ggtitle("Causal forest with LASSO fixed effect") +
		theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
		axis.text.y = element_text(size = 10),
		axis.title.y = element_text(size = 10),
		axis.text.x = element_text(size = 10),
		axis.title.x = element_text(size = 10))
	
lasso.fe.hte.plot


mean.hte.plot <- ggplot(data.subset, aes(x = mean.hte, y = true_effect)) +
		geom_point(shape = 4, colour = "#56B4E9") +
		geom_abline(intercept = 0, slope = 1, color="gray50", 
                 size=1) +
		theme_minimal() +
		geom_smooth(method = lm, linetype="dashed", se = FALSE, colour = "black") +
		xlab("estimated treatment effect") +
		xlim(0, 1.5)  +
		ylab("true treatment effect") +
		ylim(0, 1.5)  +
		ggtitle("Causal forest with mean y") +
		theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
		axis.text.y = element_text(size = 10),
		axis.title.y = element_text(size = 10),
		axis.text.x = element_text(size = 10),
		axis.title.x = element_text(size = 10))
	
mean.hte.plot

demeaned.hte.plot <- ggplot(data.subset, aes(x = demeaned.hte, y = true_effect)) +
		geom_point(shape = 4, colour = "#56B4E9") +
		geom_abline(intercept = 0, slope = 1, color="gray50", 
                 size=1) +
		theme_minimal() +
		geom_smooth(method = lm, linetype="dashed", se = FALSE, colour = "black") +
		xlab("estimated treatment effect") +
		xlim(0, 1.5)  +
		ylab("true treatment effect") +
		ylim(0, 1.5)  +
		ggtitle("Causal forest with demeaned data") +
		theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
		axis.text.y = element_text(size = 10),
		axis.title.y = element_text(size = 10),
		axis.text.x = element_text(size = 10),
		axis.title.x = element_text(size = 10))
	
demeaned.hte.plot


part.demeaned.hte.plot <- ggplot(data.subset, aes(x = part.demeaned.hte, y = true_effect)) +
		geom_point(shape = 4, colour = "#56B4E9") +
		geom_abline(intercept = 0, slope = 1, color="gray50", 
                 size=1) +
		theme_minimal() +
		geom_smooth(method = lm, linetype="dashed", se = FALSE, colour = "black") +
		xlab("estimated treatment effect") +
		xlim(0, 1.5)  +
		ylab("true treatment effect") +
		ylim(0, 1.5)  +
		ggtitle("Causal forest with demeaned y-hat") +
		theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
		axis.text.y = element_text(size = 10),
		axis.title.y = element_text(size = 10),
		axis.text.x = element_text(size = 10),
		axis.title.x = element_text(size = 10))
	
part.demeaned.hte.plot


##base.hte  (drop this one for a 2x3)
##cluster.hte
##ols.fe.hte
##lasse.fe.hte
##mean.hte
##demeaned.hte
##part.demeaned.hte 


####Stack 'em up:
multiplot(cluster.hte.plot, ols.fe.hte.plot,  demeaned.hte.plot, lasso.fe.hte.plot,  mean.hte.plot, 
		part.demeaned.hte.plot, cols = 2)

ggsave(multiplot(cluster.hte.plot, ols.fe.hte.plot,  demeaned.hte.plot, lasso.fe.hte.plot,  mean.hte.plot, 
		part.demeaned.hte.plot, cols = 2), filename = paste(project.dir, "figures/figure3.pdf", sep="/"),
       width = 7, height = 7, units = "in")








