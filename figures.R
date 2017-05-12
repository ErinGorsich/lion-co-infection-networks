###########################################################
###########################################################
###########################################################
# Figures for Lion Co-infection manuscript
# Erin E. Gorsich
# 12 - May - 2017
###########################################################
###########################################################
###########################################################

###########################################################
# read in HB's data, packages and scripts
###########################################################
rm(list = ls())
setwd("~/GitHub/lion-co-infection-networks")
alldata<-read.csv("~/Documents/postdoc_buffology/HB-lion_coinfection_network/COINCLIQ.csv")
drop <- c('THEILERIABABESIA', 'OTHER', 'BABESA') # take out some
alldata <- alldata[, c(!names(alldata) %in% drop)]
dataneg<-alldata[alldata$FIV=="Negative",]
datapos<-alldata[alldata$FIV=="Positive",]
plotdata<-alldata[,c(11:length(colnames(alldata)))]

library(lattice)  
library(bipartite)
library(igraph)  # tried bigraph but not avaliable
library(gplots)
library(HiveR)
library(grid)

# functions that groom raw data into network matricies
source("prep_data_to_matrix.R")
# functions that calculate network statistics   
source("calculate_net_and_node_stats.R")
# functions to make hive plots  
source("make_bipart_hive_plot.R")   
source("mod.edge2HPD.R")
#####################


###########################################################
###########################################################
# Prep analyses for figures
###########################################################
###########################################################

# make networks
# all animals, FIV+ animals, FIV- animals
##################################################
# turn data into matricies for monopartite and bipartite nets
all2<-make_all_matrix(alldata)  
all3<-make_all_lionmatrix(alldata) 
allb<-make_all_bipartite_matrix(alldata)
FIVneg2<-make_matrix(dataneg); colnames(FIVneg2) <- colnames(plotdata[-1])
FIVneg3<-make_lionmatrix(dataneg); rownames(FIVneg3) <- colnames(plotdata[-1])
FIVpos2<-make_matrix(datapos); colnames(FIVpos2) <- colnames(plotdata[-1])
FIVpos3<-make_lionmatrix(datapos); rownames(FIVpos3) <- colnames(plotdata[-1])

# parnet: nodes=parasites; edges= # of shared hosts, symmetric
all_parnet<- tcrossprod(t(all2)) 
FIVneg_parnet<- tcrossprod(t(FIVneg2))
FIVpos_parnet<- tcrossprod(t(FIVpos2))
# all_lionnet: nodes=lions; edges=number of shared parasites, symmetric
all_lionnet<-tcrossprod(t(all3))
FIVneg_lionnet<-tcrossprod(t(FIVneg3))  
FIVpos_lionnet<-tcrossprod(t(FIVpos3))

# Bipartite networks: nodes = parasites and lions, linkes = infection
FIVpos_binet<- make_bipartite_matrix(datapos)
FIVneg_binet<- make_bipartite_matrix(dataneg)

# summarize and write results to csv: 
calculate_net_and_node_stats(all_parnet, name="allparasitenetwork")
calculate_net_and_node_stats(all_lionnet, name="alllionnetwork")
calculate_net_and_node_stats(FIVpos_parnet, name="FIVposparasitenetwork")
calculate_net_and_node_stats(FIVneg_parnet, name="FIVnegparasitenetwork")
calculate_net_and_node_stats(FIVpos_lionnet, name="FIVposlionnetwork")
calculate_net_and_node_stats(FIVneg_lionnet, name="FIVneglionnetwork")


###########################################################
###########################################################
# Figure 2: Richness by FIV status
###########################################################
###########################################################
# read only necessary results back in ...
bineg<- read.csv("~/Documents/postdoc_buffology/HB-lion_coinfection_network/lion_coinfection_networks_files/nodedfbinet_FIVneg.csv")
binegl<- bineg[bineg$axis==1,]
bipos<- read.csv("~/Documents/postdoc_buffology/HB-lion_coinfection_network/lion_coinfection_networks_files/nodedfbinet_FIVpos.csv")
biposl<- bipos[bipos$axis==1,]

breaks <- seq(0, 10, 1)
# Note: Value divided by maximum value in the network...
poswt <- hist(biposl$degree_nonorm, breaks = breaks, plot = FALSE)
negwt <- hist(binegl$degree_nonorm, breaks = breaks, plot = FALSE)

# weighted
png("Figure2_richness.png", height = 480, width = 700, units = "px")
par(mar = c(6, 7, 4, 2))  # increase bottom and left by one and two
plot(negwt$mids, negwt$density, xlim = c(0, 10), 
	ylim= c(0, 0.36), xlab = "Parasite richness",
	las= 1, ylab= "",
	type = "b", cex = 1.5, pch= 19, lty = 1, lwd = 1, 
	tck = -0.025, cex.axis = 1.5, cex.lab = 1.8,
	bty = "n")
lines(poswt$mids, poswt$density, lwd = 1, lty = 3, cex = 2)
points(poswt$mids, poswt$density, pch= 21, lty = 3, cex = 2)
legend("topright", c("FIV -", "FIV +"), lty=c(1, 3), 
	pch=c(19, 21), lwd=1, bty="n", cex = 2)
mtext(text = "Density (proportion of lions)", 
	cex = 1.8, side = 2, line = 4.5)
dev.off()


###########################################################
###########################################################
# Figure 1: Richness by FIV status
###########################################################
###########################################################
# Hive plots
make_bipart_hive_plot(allb, datatype="all", demogdata=alldata, save_nodedf=FALSE, name="nosave")
