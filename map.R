###########################################
###########################################
# Figures and Map of Lion Samples for HB
# 17-Jan-2017
# http://www.nickeubank.com/wp-content/uploads/2015/10/RGIS3_MakingMaps_part1_mappingVectorData.html#basemaps
###########################################
###########################################
library(sp)
library(maptools)   # for geospatial services; also loads foreign and sp
library(rgdal)      # for map projection work; also loads sp
library(RgoogleMaps)
library(dismo)		# basemaps
library(RColorBrewer)
# stupid function to make colors transparent- there are way easier ways to do this
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

###########################################
# Load Park boundary & rivers
###########################################
# To plot any shp file, you read them and then project them, this is the common projection
# plot then the name plots them, you can use any of the normal commands to change things (fill, col to specify border and inside colors, lwd to specify line type ect.)
boundary<-readShapePoly("~/Documents/GIS/KNP_data/kruger_boundary.shp")
proj4string(boundary) <- "+proj=longlat +datum=WGS84"
rivers<-readShapeLines("~/Documents/GIS/KNP_data/rivers_main.shp")
proj4string(rivers) <- "+proj=longlat +datum=WGS84"
plot(boundary, axes=TRUE, border="black", col = "gray", lwd = 1.2, main = "") 
plot(rivers, axes=TRUE, col="dark blue", add= TRUE, lwd= 0.8) 

###########################################
# Load base image; reproject boundary so same
###########################################
# Use RgoogleMaps AND the dismo package's basemap
# specify Lat/Long range
x = c(3400000, 3588000); y = c(-2549398, -2940755)
xy = cbind(x, y)
base.map <- gmap(xy, type = "satellite", scale = 2) # type = hybrid, terrain, satellite, roadmap
#base.map <- gmap(boundary, type = "satellite", scale = 2) # type = hybrid, terrain, satellite, roadmap
reprojected.boundary <- spTransform(boundary, base.map@crs)

###########################################
# Load infection data
###########################################
viral <- read.csv("~/Documents/postdoc_buffology/HB-lion_coinfection_network/data_for_map/VIRALRICH.csv")
viral$LONGITUDE <- -viral$LONGITUDE
viral$LATITUDE <- - viral$LATITUDE
hemo <- read.csv("~/Documents/postdoc_buffology/HB-lion_coinfection_network/data_for_map/PRIDERICHHEMO.csv")
hemo$LONGITUDE <- -hemo$LONGITUDE
hemo$LATITUDE <- - hemo$LATITUDE
gi <- read.csv("~/Documents/postdoc_buffology/HB-lion_coinfection_network/data_for_map/PRIDERICHGI.csv")
gi$LONGITUDE <- - gi$LONGITUDE
gi$LATITUDE <- - gi$LATITUDE

# set up colors, max richness = 9.  Using the color pallettes in RColorBrewer
max(gi$RICHNESS, hemo$RICHNESS, viral$RICHNESS)
colmatch = data.frame(richness = seq(0,9,1), col = c("white", brewer.pal(9, "Reds")))
viral$col <- colmatch$col[match(viral$RICHNESS, colmatch$richness)]
hemo$col <- colmatch$col[match(hemo$RICHNESS, colmatch$richness)]
gi$col <- as.character(colmatch$col[match(gi$RICHNESS, colmatch$richness)])

# switch the coordinates of the points to same projection as map
coordinates(viral) <- c("LATITUDE", "LONGITUDE")  # named backwards, x, y
proj4string(viral) <- CRS("+proj=longlat +datum=WGS84")
reprojected.viral.points <- spTransform(viral, base.map@crs)

coordinates(hemo) <- c("LATITUDE", "LONGITUDE")  # named backwards, x, y
proj4string(hemo) <- CRS("+proj=longlat +datum=WGS84")
reprojected.hemo.points <- spTransform(hemo, base.map@crs)

coordinates(gi) <- c("LATITUDE", "LONGITUDE")  # named backwards, x, y
proj4string(gi) <- CRS("+proj=longlat +datum=WGS84")
reprojected.gi.points <- spTransform(gi, base.map@crs)

###########################################
# Plots
###########################################
# viral
plot(base.map, axes = FALSE)
plot(reprojected.boundary, add = TRUE, border = "black", 
	col = makeTransparent("gray", 110), lwd = 1.2, title = "Viral Richness")
points(reprojected.viral.points, pch = 21, col = "black", 
	bg = as.character(reprojected.viral.points@data$col),add = TRUE, 
	cex = (reprojected.viral.points@data$RICHNESS + 1)/2)
#hemo
plot(base.map, axes = FALSE)
plot(reprojected.boundary, add = TRUE, border = "black", col = makeTransparent("gray", 110), lwd = 1.2)
points(reprojected.hemo.points, pch = 21, col = "black", 
	bg = as.character(reprojected.hemo.points@data$col), 
	add = TRUE, cex = (reprojected.hemo.points@data$RICHNESS + 1)/2)
# gi
plot(base.map, axes = FALSE)
plot(reprojected.boundary, add = TRUE, border = "black", 
	col = makeTransparent("gray", 110), lwd = 1.2)
points(reprojected.gi.points, pch = 21, bg = as.character(reprojected.gi.points@data$col),
	col = "black", add = TRUE, cex = (reprojected.gi.points@data$RICHNESS + 1)/2)
# legend
legend_image <- as.raster(matrix(colfunc(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)

legend.col <- function(col, lev){
 
opar <- par
 
n <- length(col)
 
bx <- par("usr")
 
box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
box.cy <- c(bx[3], bx[3])
box.sy <- (bx[4] - bx[3]) / n
 
xx <- rep(box.cx, each = 2)
 
par(xpd = TRUE)
for(i in 1:n){
 
yy <- c(box.cy[1] + (box.sy * (i - 1)),
box.cy[1] + (box.sy * (i)),
box.cy[1] + (box.sy * (i)),
box.cy[1] + (box.sy * (i - 1)))
polygon(xx, yy, col = col[i], border = col[i])
 
}
par(new = TRUE)
plot(0, 0, type = "n",
ylim = c(min(lev), max(lev)),
yaxt = "n", ylab = "",
xaxt = "n", xlab = "",
frame.plot = FALSE)
axis(side = 4, las = 2, tick = FALSE, line = .25)
par <- opar
}
plot.new()
legend.col(col = as.character(colmatch$col), lev = colmatch$richness)

# for additional style information: style="feature:road.local|element:geometry|hue:0x00ff00|saturation:100 &style=feature:landscape|element:geometry|lightness:-100"


