library('igraph')
library('circlize')
library('RColorBrewer')
library('dplyr')
setwd("~/Documents/postdoc_buffology/HB-lion_coinfection_network/lion_coinfection_networks_files")
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
#tiff(filename="beef_circle2", height=90, width=90, units="mm", res=600, compression="lzw")
par(mar = c(1,1,1,1))  #0000 initially
circos.clear()
#Set how Circos will draw external ring, direction and starting point 
circos.par(cell.padding = c(0, 0, 0, 0), clock.wise = TRUE, start.degree = 90, gap.degree=0.6)  # first round figs at 0.8
#Initialize circos, limits of table and values (factors)
circos.initialize(factors = factors, xlim = cbind(df1$xmin, df1$xmax))
##########################
# Plot outside

circos.trackPlotRegion(
	#Unique row values
	#factors = as.factor(df1$order),  #neworder to go by largest size
 	ylim = c(0, 1),
	# bg.border = NA,
 
	#Colors for outer ring by category
 	#bg.col =df1$rcol,
 	#track.height = 0.05,
 	panel.fun = function(x, y){
		name = get.cell.meta.data("sector.index")
  		i = get.cell.meta.data("sector.numeric.index")
  		xlim = get.cell.meta.data("xlim")
  		ylim = get.cell.meta.data("ylim")
  		#Details for text style in outer ring
  		circos.text(mean(xlim), 1, name, adj = c(0.5, 0),facing="bending",cex=0.4)  # previous label as sector.name didn't work?
  		circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], col = df1$rcol[i], border=df1$lcol[i], lwd=0.4)
  		circos.lines(x=rep(df1$xmax_mid[i],2), y=c(0,1), col=df1$lcol[i], lwd=0.4)
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
df2 <- cbind(as.data.frame(mat), orig = seq(1:length(rownames(mat))))
df2 <- reshape(df2, idvar="orig", varying=list(1:n), direction="long",
	timevar="dest", time=seq(1:length(rownames(mat))),  v.names = "m")

# added by me to df2, adds neworig and dest columns to df2
df2$neworig<-NA
df2$newdest<-NA
for (i in 1: length(df2$orig)){
	for (j in 1: length(df1$order)){
		if (df2$orig[i]==df1$order[j]) df2$neworig[i]<-df1$neworder[j]
		if (df2$dest[i]==df1$order[j]) df2$newdest[i]<-df1$neworder[j]
	}
}

df2 <- arrange(df2, neworig, desc(m))  # here m is the column
#loose zero links
df2 <- subset(df2, m>0)

for(k in 1:nrow(df2)){
  #i,j reference of flow matrix
  i<-match(df2$neworig[k], df1$neworder)
  j<-match(df2$newdest[k], df1$neworder)
  q<-match(df2$orig[k], df1$neworder)  # at k=1 want q=26
  p<-match(df2$dest[k], df1$neworder)

 #plot link  
 # this way makes it so you plot incomming shipments followed by outgoing shipments. 
  circos.link(sector.index1 = df1$region[i], point1=c(df1$sum1[i], df1$sum1[i] + abs(m[q, p])),  # specify origin link
          sector.index2=df1$region[j], point2=c(df1$sum2[j], df1$sum2[j] + abs(m[q, p])),      # specify destination link
          col = df1$rcol[i],
          border= df1$lcol[i], lwd=0.4, rou=0.75
          )
  #update sum1 and sum2 for use when plotting the next link
  df1$sum1[i] = df1$sum1[i] + abs(m[q, p])
  df1$sum2[j] = df1$sum2[j] + abs(m[q, p])
}
