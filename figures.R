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
# Read in HB's data, packages and scripts
###########################################################
rm(list = ls())
setwd("~/GitHub/lion-co-infection-networks")
alldata <- read.csv("~/Documents/postdoc_buffology/HB-lion_coinfection_network/COINCLIQ.csv")
drop <- c('THEILERIABABESIA', 'OTHER', 'BABESA') # take out some
alldata <- alldata[ , c(!names(alldata) %in% drop)]
dataneg <- alldata[alldata$FIV == "Negative",]
datapos <- alldata[alldata$FIV == "Positive",]
plotdata <- alldata[ , c(11:length(colnames(alldata)))]

library(lattice)  
library(bipartite)
library(igraph)  # tried bigraph but not avaliable
library(gplots)
library(HiveR)
library(grid)
library('circlize')
library('RColorBrewer')
library('dplyr')

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
all2 <- make_all_matrix(alldata)  
all3 <- make_all_lionmatrix(alldata) 
allb <- make_all_bipartite_matrix(alldata)
FIVneg2 <- make_matrix(dataneg); colnames(FIVneg2) <- colnames(plotdata[-1])
FIVneg3 <- make_lionmatrix(dataneg); rownames(FIVneg3) <- colnames(plotdata[-1])
FIVpos2 <- make_matrix(datapos); colnames(FIVpos2) <- colnames(plotdata[-1])
FIVpos3 <- make_lionmatrix(datapos); rownames(FIVpos3) <- colnames(plotdata[-1])

# parnet: nodes=parasites; edges= # of shared hosts, symmetric
all_parnet <- tcrossprod(t(all2))
FIVneg_parnet <- tcrossprod(t(FIVneg2))
FIVpos_parnet <- tcrossprod(t(FIVpos2))
# all_lionnet: nodes=lions; edges=number of shared parasites, symmetric
all_lionnet <- tcrossprod(t(all3))
FIVneg_lionnet <- tcrossprod(t(FIVneg3))  
FIVpos_lionnet <- tcrossprod(t(FIVpos3))

# Bipartite networks: nodes = parasites and lions, linkes = infection
FIVpos_binet <- make_bipartite_matrix(datapos)
FIVneg_binet <- make_bipartite_matrix(dataneg)

# summarize and write results to csv: 
calculate_net_and_node_stats(all_parnet, name = "allparasitenetwork")
calculate_net_and_node_stats(all_lionnet, name ="alllionnetwork")
calculate_net_and_node_stats(FIVpos_parnet, name ="FIVposparasitenetwork")
calculate_net_and_node_stats(FIVneg_parnet, name ="FIVnegparasitenetwork")
calculate_net_and_node_stats(FIVpos_lionnet, name ="FIVposlionnetwork")
calculate_net_and_node_stats(FIVneg_lionnet, name ="FIVneglionnetwork")


###########################################################
###########################################################
# Figure 0: Ugly network diagram & Overall Hive Plot with FIV highlihgted
###########################################################
###########################################################
# specify color choices by type
gut = c("ASCARIDS", "TAPES", "HOOKS", "COCCIDIA", "TOXO", "WHIPS")
bloodborne = c("ERLICHIAANAPLASMA", "BFELIS", "BLEO", 
	"BMICROTI", "BLENGUA", "TBICORNIS", "HEPATOZOON", 
	"TANNAE", "BCANIS", "BROSSI", "BVOGELI")
virus = c("FIV", "CDV", "FPV", "CALICI", "CORONA")

newdf<-data.frame(parasite=c(virus, gut, bloodborne), 
    type=c(rep("virus", length(virus)), rep("gut", length(gut)), 
    		rep("bloodborne", length(bloodborne)) ), 
    prev=NA, 
    col= c("#000033", rep("#663399", length(virus)-1), 
    		rep("#006666", length(gut)), 
		rep("#FF6600", length(bloodborne))) )

newdf<-data.frame(parasite=c(virus, gut, bloodborne), 
    type=c(rep("virus", length(virus)), rep("gut", length(gut)), 
    		rep("bloodborne", length(bloodborne)) ), 
    prev=NA, 
    col= c("red4", rep("#663399", length(virus)-1), 
    		rep("#66CCCC", length(gut)), 
		rep("#FF9933", length(bloodborne))) )

# Ugly network diagram
# all parasites: 
set.seed(1)
ga<-graph.adjacency(all_parnet, mode="upper", weighted=TRUE) # undirected/named/weighted/self
ga_all<- simplify(ga, remove.multiple=FALSE, remove.loops=TRUE)
V(ga_all)$color <- as.character(newdf$col[match( V(ga_all)$name, newdf$parasite)])

ga<-graph.adjacency(FIVneg_parnet, mode="upper", weighted=TRUE) # undirected/named/weighted/self
ga_neg <- simplify(ga, remove.multiple=FALSE, remove.loops=TRUE)
V(ga_neg)$color <- as.character(newdf$col[match( V(ga_neg)$name, newdf$parasite)])

ga <- graph.adjacency(FIVpos_parnet, mode="upper", weighted=TRUE) # undirected/named/weighted/self
ga_pos <- simplify(ga, remove.multiple=FALSE, remove.loops=TRUE)
V(ga_pos)$color <- as.character(newdf$col[match(V(ga_pos)$name, newdf$parasite)])
# test1
par(mfrow = c(1,3))
plot(ga_all, main = "All parasites", vertex.label = NA)
plot(ga_neg, main = "FIV -", vertex.label = NA)
plot(ga_pos, main = "FIV +", vertex.label = NA)

# test 2: ring layout
ga_pos$layout <- layout.circle
ga_neg$layout <- layout.circle
par(mfrow = c(1,2), mar = c(1, 0.5, 1, 0.1))
plot(ga_neg, main = "FIV -", vertex.label = NA)
plot(ga_pos, main = "FIV +", vertex.label = NA)

# test3
# size of circle proportional to abundance...
n <- data.frame(degree(ga_neg, mode = "all"))[,1] 
p <- data.frame(degree(ga_pos, mode = "all"))[,1]
V(ga_neg)$size <- (n + 1)# / (6*max(c(n, p))) # default size is 15
V(ga_pos)$size <- (p + 1)# / (6*max(c(n, p)))
par(mfrow = c(1,2), mar = c(1, 0.5, 1, 0.1))
plot(ga_neg, main = "FIV -", vertex.label = NA)
plot(ga_pos, main = "FIV +", vertex.label = NA)

# test3 with random layout
n <- data.frame(degree(ga_neg, mode = "all"))[,1] 
p <- data.frame(degree(ga_pos, mode = "all"))[,1]
V(ga_neg)$size <- (n + 1)# / (6*max(c(n, p))) # default size is 15
V(ga_pos)$size <- (p + 1)# / (6*max(c(n, p)))
par(mfrow = c(1,2), mar = c(1, 0.5, 1, 0.1))
plot(ga_neg, main = "FIV -", vertex.label = NA, layout = layout.random)
plot(ga_pos, main = "FIV +", vertex.label = NA, layout = layout.random)

# test4- width of line proporitonal to number of shared hosts
#E(ga_neg)$width <-  # default is 1
#E(ga_pos)$width <-

numlions <- length(rownames(FIVneg2))
par(mfrow = c(4, 5), mar = c(1, 0.5, 1, 0.1))
plot(ga_neg, main = "FIV -", vertex.label = NA, layout = layout.random)
for(i in 1:19){
	ss_FIVpos2 <- sample_n(data.frame(FIVpos2), numlions, replace = FALSE)
	FIVpos_parnet<- tcrossprod(t(ss_FIVpos2))

	# make networks and plot	
	ga<-graph.adjacency(FIVpos_parnet, mode="upper", weighted=TRUE) 
	ga_pos<- simplify(ga, remove.multiple=FALSE, remove.loops=TRUE)
	V(ga_pos)$color <- as.character(newdf$col[match(V(ga_pos)$name, newdf$parasite)])
	V(ga_pos)$size <- 1 + data.frame(degree(ga_pos, mode = "all"))[,1]

	plot(ga_pos, main = "Subsample FIV +", vertex.label = NA, layout = layout.random)
	rm(ga, gapos, ss_FIVpos2)
}

plot(ga_neg, main = "FIV -", vertex.label = NA, layout = layout.circle)
for(i in 1:19){
	ss_FIVpos2 <- sample_n(data.frame(FIVpos2), numlions, replace = FALSE)
	FIVpos_parnet<- tcrossprod(t(ss_FIVpos2))

	# make networks and plot	
	ga<-graph.adjacency(FIVpos_parnet, mode="upper", weighted=TRUE) 
	ga_pos<- simplify(ga, remove.multiple=FALSE, remove.loops=TRUE)
	V(ga_pos)$color <- as.character(newdf$col[match(V(ga_pos)$name, newdf$parasite)])
	V(ga_pos)$size <- 1 + data.frame(degree(ga_pos, mode = "all"))[,1]

	plot(ga_pos, main = "Subsample FIV +", vertex.label = NA, layout = layout.circle)
	rm(ga, gapos, ss_FIVpos2)
}

# function makes hive plot
png("Figure1_hive_alldata.png", height = 600, width = 1200, units = "px")
make_bipart_hive_plot(allb, datatype="all", demogdata=alldata, save_nodedf=FALSE, name="nosave")
dev.off()

###########################################################
###########################################################
# Figure 2: Richness by FIV status
###########################################################
###########################################################
# read only necessary results back in ...
bineg <- read.csv("~/Documents/postdoc_buffology/HB-lion_coinfection_network/lion_coinfection_networks_files/nodedfbinet_FIVneg.csv")
binegl <- bineg[bineg$axis == 1,]
bipos <- read.csv("~/Documents/postdoc_buffology/HB-lion_coinfection_network/lion_coinfection_networks_files/nodedfbinet_FIVpos.csv")
biposl <- bipos[bipos$axis == 1,]

breaks <- seq(0, 10, 1)
# Note: Value divided by maximum value in the network...
poswt <- hist(biposl$degree_nonorm, breaks = breaks, plot = FALSE)
negwt <- hist(binegl$degree_nonorm, breaks = breaks, plot = FALSE)

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

# Still need to do resampling

###########################################################
###########################################################
# Figure 3: Cluster Figures
###########################################################
###########################################################
par <- read.csv("~/Documents/postdoc_buffology/HB-lion_coinfection_network/lion_coinfection_networks_files/allparnet_com.csv")
lion <- read.csv("~/Documents/MATLAB/lionparasites.csv")

# Make a list of all parasite pairs and their associated community
data <- data.frame(
	opar = rep(par$Parasite, length(par$Parasite)), 
	Ocommunity = rep(par$Community, length(par$Parasite)),
	dpar = rep(par$Parasite, each = length(par$Parasite)), 
	Dcommunity = rep(par$Community, each = length(par$Parasite)), 
	oname = rep(par$Name, length(par$Parasite)), 
	dname = rep(par$Name, each = length(par$Parasite)),
	wt = NA)
for(i in 1:length(data[,1])){
	data$wt[i] <- lion[which(as.character(lion$names) == data$opar[i]), colnames(lion) == data$dpar[i]]
}
# remove duplicates
#data$wt[data$opar == data$dpar] <- 0

EL = cbind(as.character(data$oname), as.character(data$dname))   # needs to be in matrix form for igraph
graph=graph.edgelist(EL, directed=FALSE)
E(graph)$weight <- data$wt
adjmat=get.adjacency(graph, attr="weight")  # symmetric, weighted, no self loops.

# make colors
colmatch <- data.frame(
	com = seq(1, 4), 
	col = brewer.pal(8, "Paired")[c(1,3,5,7)], 
	coldark = brewer.pal(8, "Paired")[c(2,4,6,8)] )
comcol <- data.frame(counties = par$Name, com = par$Community)
comcol$coldark= as.character(colmatch$coldark[match(comcol$com, colmatch$com)])
comcol$col= as.character(colmatch$col[match(comcol$com, colmatch$com)])

# make df1: holds data to set spacing
df1 <- data.frame(order=c(1:length(unique(data$oname))),
	region=sort(unique(data$oname)), xmin=rep(0, length(unique(data$oname))), 
	xmax=NA, xmax_mid=NA, neworder = NA)
df1$rcol=as.character(comcol$col[match(df1$region, comcol$counties)])
df1$lcol=as.character(comcol$coldark[match(df1$region, comcol$counties)])
df1$com <- comcol$com[match(df1$region, comcol$counties)]

# match order of matrix names and df1 names
ord1 = numeric()
for(i in 1:dim(adjmat)[1]){
	ord1[i]=which(df1$region==colnames(adjmat)[i]) 
}  # order to match that in comcol
adjmat.ord=adjmat[order(ord1),order(ord1)] 
mat=data.matrix(adjmat.ord)
cn=colnames(mat)
rn=rownames(mat)

#Table Marginal values
col_sum = colSums(mat)
row_sum = rowSums(mat)

#Limits for table
df1$xmax = row_sum+col_sum  #requires correct ordering...
df1$xmax_mid= col_sum

# Really want order on plot to be by community type, then xmax. 
orderdf = data.frame( 
	new = c("TANNAE", "BROSSI", "BCANIS", "BVOGELI", "FPV",
	"HEPATO", "BLENGUA", "FIV", "BFELIS", "CALICI", "ASCARIDS", 
	"BLEO", "ERLI-ANAPL", "BMICROTI", "TBICORNIS", "CORONA", 
	"WHIPS", "TOXO", "CDV", "COCC", "TAPES", "HOOKS"), 
	order = seq(1:length(df1$region)) )

for (i in 1:length(df1$neworder)){
	df1$neworder[i] <- orderdf$order[orderdf$new == df1$region[i]]
}
df1<-df1[order(df1$neworder, decreasing=FALSE),]  # reorder df1 by xmax, will plot by this order

#Number of unique columns
factors = colnames(mat)[df1$order]
factors = factor(factors, levels = factors) # ordered alphabetically to match df1
colnames(mat)[df1$order] == df1$region  # test, should be the same



####################
png(filename="Figure3_parasitecluster.png", height=800, width=800, units="px")
par(mar = c(2,2,2,2))  #0000 initially
circos.clear()
#Set how Circos will draw external ring, direction and starting point 
circos.par(cell.padding = c(0, 0, 0, 0), clock.wise = TRUE, start.degree = 90, gap.degree=0.6)  # first round figs at 0.8
#Initialize circos, limits of table and values (factors)
circos.initialize(factors = factors, xlim = cbind(df1$xmin, df1$xmax))
##########################
# Plot outside

#circos.trackPlotRegion(
	#Unique row values
	#factors = as.factor(df1$order),  #neworder to go by largest size
 	ylim = c(0, 1),
	# bg.border = NA,
 
	#Colors for outer ring by category
 	#bg.col =df1$rcol,
 	#track.height = 0.05,
# 	panel.fun = function(x, y){
#		name = get.cell.meta.data("sector.index")
#  		i = get.cell.meta.data("sector.numeric.index")
#  		xlim = get.cell.meta.data("xlim")
#  		ylim = get.cell.meta.data("ylim")
  		#Details for text style in outer ring
  		#circos.text(mean(xlim), 1, name, adj = c(0.5, 0),facing="bending",cex=0.4)  # previous label as sector.name didn't work?
  		#circos.text(mean(xlim), 1, name, adj = c(0.5, 0), facing = "clockwise", cex = 0.6)
#  		circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], col = df1$rcol[i], border=df1$lcol[i], lwd=0.4)
#  		circos.lines(x=rep(df1$xmax_mid[i],2), y=c(0,1), col=df1$lcol[i], lwd=0.4)
#  		print(name)
#	}
#)

circos.trackPlotRegion(
	#Unique row values
	#factors = as.factor(df1$order),  #neworder to go by largest size
 	ylim = c(0, 1.4),
	# bg.border = NA,
	#Colors for outer ring by category
 	#bg.col =df1$rcol,
 	#track.height = 0.05,
 	panel.fun = function(x, y){
		name = get.cell.meta.data("sector.index")
  		i = get.cell.meta.data("sector.numeric.index")
  		xlim = get.cell.meta.data("xlim")
  		ylim = get.cell.meta.data("ylim")
  		theta = mean(get.cell.meta.data("xplot"))
  		#Details for text style in outer ring
  		#circos.text(mean(xlim), 1, name, adj = c(0.5, 0),facing="bending",cex=0.4)  # previous label as sector.name didn't work?
  		circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], col = df1$rcol[i], border=df1$lcol[i], lwd=0.4)
  		#circos.lines(x=rep(df1$xmax_mid[i],2), y=c(0,1), col=df1$lcol[i], lwd=0.4)
  		if (theta < 90 || theta > 270){
  			circos.text(mean(xlim), 0.3, name, adj = c(0, 0), facing = "clockwise", cex = 0.6)
  		} else {
  			circos.text(mean(xlim), 1.2, name, adj = c(0, 0), facing = "reverse.clockwise", cex = 0.6)
  		}
  		print(name)
	}
)


# make df2, holds flows
n<-nrow(mat)
df1<-df1[order(df1$order),]  # to match order of m
m<-mat
df1$sum1<-df1$xmax_mid
df1$sum2<-numeric(n)
df1$sum3<-df1$xmax_mid
df1<-df1[order(df1$neworder),]  # to match order of m

#df2 <- cbind(as.data.frame(mat),orig=rownames(mat),  stringsAsFactors=FALSE)
#df2 <- cbind(as.data.frame(mat), orig = seq(1:length(rownames(mat))))
#df2 <- reshape(df2, idvar="orig", varying=list(1:n), direction="long",
#	timevar="dest", time=seq(1:length(rownames(mat))),  v.names = "m")

# added by me to df2, adds neworig and dest columns to df2
#df2$neworig<-NA
#df2$newdest<-NA
#for (i in 1: length(df2$orig)){
#	for (j in 1: length(df1$order)){
#		if (df2$orig[i]==df1$order[j]) df2$neworig[i]<-df1$neworder[j]
#		if (df2$dest[i]==df1$order[j]) df2$newdest[i]<-df1$neworder[j]
#	}
#}

#df2 <- arrange(df2, neworig, desc(m))  # here m is the column
#loose zero links
#df2 <- subset(df2, m>0)

#for(k in 1:nrow(df2)){
#  #i,j reference of flow matrix
#  i<-match(df2$neworig[k], df1$neworder)
#  j<-match(df2$newdest[k], df1$neworder)
#  q<-match(df2$orig[k], df1$neworder)  # at k=1 want q=26
#  p<-match(df2$dest[k], df1$neworder)

# new, non-directed
df2 <- cbind(as.data.frame(mat), orig = seq(1:length(rownames(mat))))
df2 <- reshape(df2, idvar="orig", varying=list(1:n), direction="long",
	timevar="dest", time=seq(1:length(rownames(mat))),  v.names = "m")
df2 <- df2[df2$orig > df2$dest,]

# added by me to df2, adds neworig and dest columns to df2
df2$neworig<-NA
df2$newdest<-NA
for (i in 1: length(df2$orig)){
	for (j in 1: length(df1$order)){
		if (df2$orig[i]==df1$order[j]) df2$neworig[i]<-df1$neworder[j]
		if (df2$dest[i]==df1$order[j]) df2$newdest[i]<-df1$neworder[j]
	}
}

#loose zero links
df2 <- subset(df2, m>0)
df2$ocol <- df1$rcol[match(df2$neworig, df1$neworder)]
df2$dcol <- df1$rcol[match(df2$newdest, df1$neworder)]
df2$within <- 1
df2$within[!(df2$ocol == df2$dcol)] <- 0
df1$sum1 <- 0
df2 <- arrange(df2, within, neworig, desc(m))  # here m is the column

for (k in 1:nrow(df2)){
	#i,j reference of flow matrix
	i<-match(df2$neworig[k], df1$neworder)
	j<-match(df2$newdest[k], df1$neworder)
	q<-match(df2$orig[k], df1$neworder)  # at k=1 want q=26
	p<-match(df2$dest[k], df1$neworder)
 
#plot link  
# this way makes it so you plot incomming shipments followed by outgoing shipments. 
  	if(df1$rcol[i] == df1$rcol[j]){
		circos.link(sector.index1 = df1$region[i], point1=c(df1$sum1[i], df1$sum1[i] + abs(2*m[q, p])), # specify origin link (half way )
			sector.index2=df1$region[j], point2=c(df1$sum2[j], df1$sum2[j] + abs(2*m[q, p])),       # specify destination link
			col = df1$rcol[i],
			border= df1$lcol[i], lwd=0.4, rou=0.75)
	} else {
		circos.link(sector.index1 = df1$region[i], point1=c(df1$sum1[i], df1$sum1[i] + abs(2*m[q, p])), # specify origin link (half way )
			sector.index2=df1$region[j], point2=c(df1$sum2[j], df1$sum2[j] + abs(2*m[q, p])),       # specify destination link
			col = "lightgray",
			border= "darkgray", lwd=0.4, rou=0.75)
	}        
#update sum1 and sum2 for use when plotting the next link
	df1$sum1[i] = df1$sum1[i] + abs(2*m[q, p])
	df1$sum1[j] = df1$sum1[j] + abs(2*m[q, p])
	df1$sum2[i] = df1$sum2[i] + abs(2*m[q, p])
	df1$sum2[j] = df1$sum2[j] + abs(2*m[q, p])
}
dev.off()

###########################################################
###########################################################
# Figure 3b: Cluster Figures with FIV elements highlighted
###########################################################
###########################################################
par <- read.csv("~/Documents/postdoc_buffology/HB-lion_coinfection_network/lion_coinfection_networks_files/allparnet_com.csv")
lion<- read.csv("~/Documents/MATLAB/lionparasites.csv")

# Make a list of all parasite pairs and their associated community
data <- data.frame(
	opar = rep(par$Parasite, length(par$Parasite)), 
	Ocommunity = rep(par$Community, length(par$Parasite)),
	dpar = rep(par$Parasite, each = length(par$Parasite)), 
	Dcommunity = rep(par$Community, each = length(par$Parasite)), 
	oname = rep(par$Name, length(par$Parasite)), 
	dname = rep(par$Name, each = length(par$Parasite)),
	wt = NA)
for(i in 1:length(data[,1])){
	data$wt[i] <- lion[which(as.character(lion$names) == data$opar[i]), colnames(lion) == data$dpar[i]]
}
# remove duplicates
#data$wt[data$opar == data$dpar] <- 0

EL=cbind(as.character(data$oname), as.character(data$dname))   # needs to be in matrix form for igraph
graph=graph.edgelist(EL, directed=FALSE)
E(graph)$weight <- data$wt
adjmat=get.adjacency(graph, attr="weight")  # symmetric, weighted, no self loops.

# make colors
colmatch <- data.frame(
	com = seq(1, 4), 
	col = brewer.pal(8, "Paired")[c(1,3,5,7)], 
	coldark = brewer.pal(8, "Paired")[c(2,4,6,8)] )
comcol <- data.frame(counties = par$Name, com = par$Community)
comcol$coldark= as.character(colmatch$coldark[match(comcol$com, colmatch$com)])
comcol$col= as.character(colmatch$col[match(comcol$com, colmatch$com)])

	
# make df1: holds data to set spacing
df1 <- data.frame(order=c(1:length(unique(data$oname))),
	region=sort(unique(data$oname)), xmin=rep(0, length(unique(data$oname))), 
	xmax=NA, xmax_mid=NA, neworder = NA)
df1$rcol=as.character(comcol$col[match(df1$region, comcol$counties)])
df1$lcol=as.character(comcol$coldark[match(df1$region, comcol$counties)])
df1$com <- comcol$com[match(df1$region, comcol$counties)]

# match order of matrix names and df1 names
ord1 = numeric()
for(i in 1:dim(adjmat)[1]){
	ord1[i]=which(df1$region==colnames(adjmat)[i]) 
}  # order to match that in comcol
adjmat.ord=adjmat[order(ord1),order(ord1)] 
mat=data.matrix(adjmat.ord)
cn=colnames(mat)
rn=rownames(mat)

#Table Marginal values
col_sum = colSums(mat)
row_sum = rowSums(mat)

#Limits for table
df1$xmax = row_sum+col_sum  #requires correct ordering...
df1$xmax_mid= col_sum

# Really want order on plot to be by community type, then xmax. 
orderdf = data.frame( 
	new = c("TANNAE", "BROSSI", "BCANIS", "BVOGELI", "FPV",
	"HEPATO", "BLENGUA", "FIV", "BFELIS", "CALICI", "ASCARIDS", 
	"BLEO", "ERLI-ANAPL", "BMICROTI", "TBICORNIS", "CORONA", 
	"WHIPS", "TOXO", "CDV", "COCC", "TAPES", "HOOKS"), 
	order = seq(1:length(df1$region)) )

for (i in 1:length(df1$neworder)){
	df1$neworder[i] <- orderdf$order[orderdf$new == df1$region[i]]
}
df1<-df1[order(df1$neworder, decreasing=FALSE),]  # reorder df1 by xmax, will plot by this order

#Number of unique columns
factors = colnames(mat)[df1$order]
factors = factor(factors, levels = factors) # ordered alphabetically to match df1
colnames(mat)[df1$order] == df1$region  # test, should be the same



####################
png(filename="Figure3b_parasitecluster_FIV.png", height=800, width=800, units="px")
par(mar = c(2,2,2,2))  #0000 initially
circos.clear()
#Set how Circos will draw external ring, direction and starting point 
circos.par(cell.padding = c(0, 0, 0, 0), clock.wise = TRUE, start.degree = 90, gap.degree=0.6)  # first round figs at 0.8
#Initialize circos, limits of table and values (factors)
circos.initialize(factors = factors, xlim = cbind(df1$xmin, df1$xmax))
##########################
# Plot outside

#circos.trackPlotRegion(
	#Unique row values
	#factors = as.factor(df1$order),  #neworder to go by largest size
 	ylim = c(0, 1),
	# bg.border = NA,
 
	#Colors for outer ring by category
 	#bg.col =df1$rcol,
 	#track.height = 0.05,
# 	panel.fun = function(x, y){
#		name = get.cell.meta.data("sector.index")
#  		i = get.cell.meta.data("sector.numeric.index")
#  		xlim = get.cell.meta.data("xlim")
#  		ylim = get.cell.meta.data("ylim")
  		#Details for text style in outer ring
  		#circos.text(mean(xlim), 1, name, adj = c(0.5, 0),facing="bending",cex=0.4)  # previous label as sector.name didn't work?
  		#circos.text(mean(xlim), 1, name, adj = c(0.5, 0), facing = "clockwise", cex = 0.6)
#  		circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], col = df1$rcol[i], border=df1$lcol[i], lwd=0.4)
#  		circos.lines(x=rep(df1$xmax_mid[i],2), y=c(0,1), col=df1$lcol[i], lwd=0.4)
#  		print(name)
#	}
#)

circos.trackPlotRegion(
	#Unique row values
	#factors = as.factor(df1$order),  #neworder to go by largest size
 	ylim = c(0, 1.4),
	# bg.border = NA,
	#Colors for outer ring by category
 	#bg.col =df1$rcol,
 	#track.height = 0.05,
 	panel.fun = function(x, y){
		name = get.cell.meta.data("sector.index")
  		i = get.cell.meta.data("sector.numeric.index")
  		xlim = get.cell.meta.data("xlim")
  		ylim = get.cell.meta.data("ylim")
  		theta = mean(get.cell.meta.data("xplot"))
  		#Details for text style in outer ring
  		#circos.text(mean(xlim), 1, name, adj = c(0.5, 0),facing="bending",cex=0.4)  # previous label as sector.name didn't work?
  		circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], col = df1$rcol[i], border=df1$lcol[i], lwd=0.4)
  		#circos.lines(x=rep(df1$xmax_mid[i],2), y=c(0,1), col=df1$lcol[i], lwd=0.4)
  		if (theta < 90 || theta > 270){
  			circos.text(mean(xlim), 0.3, name, adj = c(0, 0), facing = "clockwise", cex = 0.6)
  		} else {
  			circos.text(mean(xlim), 1.2, name, adj = c(0, 0), facing = "reverse.clockwise", cex = 0.6)
  		}
  		print(name)
	}
)


# make df2, holds flows
n<-nrow(mat)
df1<-df1[order(df1$order),]  # to match order of m
m<-mat
df1$sum1<-df1$xmax_mid
df1$sum2<-numeric(n)
df1$sum3<-df1$xmax_mid
df1<-df1[order(df1$neworder),]  # to match order of m

#df2 <- cbind(as.data.frame(mat),orig=rownames(mat),  stringsAsFactors=FALSE)
#df2 <- cbind(as.data.frame(mat), orig = seq(1:length(rownames(mat))))
#df2 <- reshape(df2, idvar="orig", varying=list(1:n), direction="long",
#	timevar="dest", time=seq(1:length(rownames(mat))),  v.names = "m")

# added by me to df2, adds neworig and dest columns to df2
#df2$neworig<-NA
#df2$newdest<-NA
#for (i in 1: length(df2$orig)){
#	for (j in 1: length(df1$order)){
#		if (df2$orig[i]==df1$order[j]) df2$neworig[i]<-df1$neworder[j]
#		if (df2$dest[i]==df1$order[j]) df2$newdest[i]<-df1$neworder[j]
#	}
#}

#df2 <- arrange(df2, neworig, desc(m))  # here m is the column
#loose zero links
#df2 <- subset(df2, m>0)

#for(k in 1:nrow(df2)){
#  #i,j reference of flow matrix
#  i<-match(df2$neworig[k], df1$neworder)
#  j<-match(df2$newdest[k], df1$neworder)
#  q<-match(df2$orig[k], df1$neworder)  # at k=1 want q=26
#  p<-match(df2$dest[k], df1$neworder)

# new, non-directed
df2 <- cbind(as.data.frame(mat), orig = seq(1:length(rownames(mat))))
df2 <- reshape(df2, idvar="orig", varying=list(1:n), direction="long",
	timevar="dest", time=seq(1:length(rownames(mat))),  v.names = "m")
df2 <- df2[df2$orig > df2$dest,]

# added by me to df2, adds neworig and dest columns to df2
df2$neworig<-NA
df2$newdest<-NA
for (i in 1: length(df2$orig)){
	for (j in 1: length(df1$order)){
		if (df2$orig[i]==df1$order[j]) df2$neworig[i]<-df1$neworder[j]
		if (df2$dest[i]==df1$order[j]) df2$newdest[i]<-df1$neworder[j]
	}
}

#loose zero links
df2 <- subset(df2, m>0)
df2$ocol <- "lightgray"
df2$dcol <- "lightgray"
df2$ocol[df2$orig == 14 | df2$dest == 14] <- "#E31A1C"
df2$dcol[df2$orig == 14 | df2$dest == 14] <- "#E31A1C"
df2$within <- 1
df2$nodeocol <- df1$rcol[match(df2$neworig, df1$neworder)]
df2$nodedcol <- df1$rcol[match(df2$newdest, df1$neworder)]
df2$within[!(df2$nodeocol == df2$nodedcol)] <- 0
df1$sum1 <- 0
df2
# hopefully plot lightgrays followed by colored
df2 <- arrange(df2, desc(ocol), desc(within), desc(m))  # here m is the column

for (k in 1:nrow(df2)){
	#i,j reference of flow matrix
	i<-match(df2$neworig[k], df1$neworder)
	j<-match(df2$newdest[k], df1$neworder)
	q<-match(df2$orig[k], df1$neworder)  # at k=1 want q=26
	p<-match(df2$dest[k], df1$neworder)
 
#plot link  
  	if(df1$rcol[i] == df1$rcol[j]){
		circos.link(sector.index1 = df1$region[i], point1=c(df1$sum1[i], 
			df1$sum1[i] + abs(2*m[q, p])), # specify origin link (half way )
			sector.index2=df1$region[j], point2=c(df1$sum2[j], df1$sum2[j] + abs(2*m[q, p])),# specify destination link
			col = df2$ocol[k],
			border= df1$lcol[i], lwd=0.4, rou=0.75)
	} else {
		circos.link(sector.index1 = df1$region[i], point1 = c(df1$sum1[i],
			df1$sum1[i] + abs(2*m[q, p])), # specify origin link (half way )
			sector.index2 = df1$region[j], point2 = c(df1$sum2[j], df1$sum2[j] + abs(2*m[q, p])),       # specify destination link
			col = df2$ocol[k],
			border = df2$ocol[k], lwd=0.4, rou=0.75)
	}        
#update sum1 and sum2 for use when plotting the next link
	df1$sum1[i] = df1$sum1[i] + abs(2*m[q, p])
	df1$sum1[j] = df1$sum1[j] + abs(2*m[q, p])
	df1$sum2[i] = df1$sum2[i] + abs(2*m[q, p])
	df1$sum2[j] = df1$sum2[j] + abs(2*m[q, p])
}
dev.off()


