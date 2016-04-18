## FUNCTION TO MAKE HIVE PLOTS for HB's data
make_bipart_hive_plot<-function(matdata, datatype, demogdata=alldata, save_nodedf, name) {
  ########################
  ########################
  ## INPUT: matdata: data matricies for bipartite network
  ## datatype= "all" or "FIV" when add classes of FIV will have to think more
  ## demogdata= dataframe containing Lion ID, age, sex, pride status, bcs (alldata)
  ## save_nodedf= TRUE/FALSE, should the nodedata frame be saved. 
  ## name= what to call output image
  ## OUTPUT: hiveplot image- use R studio for now
  ## need to set up saving and high res photo output sections
  ## colmatch not in function, need saved 
  ########################
  ########################		
  # make is so bloodbornes go on one side of the plot, rest on other!
 bloodborne = c("ERLICHIAANAPLASMA", "THEILERIABABESIA", "BABESA", 
                 "BFELIS", "BLEO", "BMICROTI", "BLENGUA", "TBICORNIS", 
                 "HEPATOZOON", "TANNAE", "BCANIS", "BROSSI", "BVOGELI")
 nonbloodborne<-c("ASCARIDS", "TAPES", "HOOKS", "COCCIDIA", "TOXO", "WHIPS",
                  "FIV", "CDV", "FPV", "CALICI", "CORONA", "OTHER")
 if(datatype=="all"){
  uppermat<-matdata[1:111,]
  lowermat<-matdata[112:136,]
  for (i in 1:length(colnames(uppermat))){
    if (colnames(uppermat)[i] %in% bloodborne){
      uppermat[,i]<-0
    }
  }
  for (j in 1:length(rownames(lowermat))){
    if (rownames(lowermat)[j] %in% nonbloodborne){
      lowermat[j,]<-0
    }
  }
  } else if (datatype=="FIV"){
   uppermat<-matdata[1:111,]
   lowermat<-matdata[112:125,]
   for (i in 1:length(colnames(uppermat))){
     if (colnames(uppermat)[i] %in% bloodborne) {
       uppermat[,i]<-0
     }
   }
   for (j in 1:length(rownames(lowermat))){
     if (rownames(lowermat)[j] %in% nonbloodborne) {
       lowermat[j,]<-0
     }
   } 
  } else {
    print("datatype is not correct")
  }        
 matdata<-rbind(uppermat, lowermat) 
 
 
  # turn adjacency matrix into an i-graph object
  ga<-graph.adjacency(matdata)	
  is.simple(ga)  # should be a simple dataframe
  
  if(datatype=="all"){
  V(ga)$type<-c(rep(TRUE, length(colnames(matdata))-length(colnames(all_parnet))),
                rep(FALSE, length(colnames(all_parnet)))) 
  } else if (datatype=="FIV"){
    V(ga)$type<-c(rep(TRUE, length(colnames(matdata))-length(colnames(FIVneg_parnet))),
                  rep(FALSE, length(colnames(FIVneg_parnet))))                        
  } else {
    print("datatype is not correct")
  }                                               
  is.bipartite(ga)
  
  nodedf<- data.frame(name=V(ga)$name, degree=NA, radius=NA, color=NA, axis= NA, 
                      symbol=19, malefemale=NA, agecat=NA, 	size=NA)
  
  # specify node information (circle size= total degree *1.5; axis by age category)
  nodedf$agecat<-as.character(alldata$LIFESTAGE[match(nodedf$name, alldata$LION)])
  nodedf$agecat[is.na(nodedf$agecat)]<- "parasite"

  nodedf$malefemale<-as.character(alldata$SEX[match(nodedf$name, alldata$LION)])
  nodedf$malefemale[is.na(nodedf$agecat)]<- "parasite"

  nodedf$color<- as.character(newdf$col[match(nodedf$name, newdf$parasite)])
  nodedf$color[is.na(nodedf$color)]<- "thistle4" # label lions

  nodedf$axis<-1
  nodedf$axis[V(ga)$type==FALSE]<-2
  
  nodedf$degree<-bipartite.degree.centrality(ga)$Bipartite.Degree.Centrality
 
  #nodedf$symbol[nodedf$sex=="female"]<-19
  #nodedf$symbol[nodedf$agecat=="male"]<-17
  
  nodedf<-nodedf[order(nodedf$axis, nodedf$degree),]
  nodedf$size[1:length(nodedf$axis[nodedf$axis==1])]<-nodedf$degree[1:length(
    nodedf$axis[nodedf$axis==1])]*1.2
  nodedf$size[length(nodedf$axis[nodedf$axis==1]):length(nodedf$axis)]<-nodedf$degree[length(
    nodedf$axis[nodedf$axis==1]):length(nodedf$axis)]*1.2 
  # specify location of node on axis 1 as the sum of the total size.
  nodedf$radius<-NA
  nodedf$radius[1]<-nodedf$size[1]
  
  for (i in 2:length(nodedf$axis[nodedf$axis==1])){
      nodedf$radius[i]<-nodedf$size[i]+ nodedf$radius[i-1]+0.2
  }
  
  # specify location of node on axis 2 as the sum of the total size.
  nodedf$radius[length(nodedf$axis[nodedf$axis==1])+1]<-nodedf$size[length( 
  nodedf$axis[nodedf$axis==1])+1]+.5
  
    for (i in (length(nodedf$axis[nodedf$axis==1])+2):((length(nodedf$axis[nodedf$axis==1]))+ 
                                                         length(nodedf$axis[nodedf$axis==2])) ){
      nodedf$radius[i]<-nodedf$size[i]+ nodedf$radius[i-1]+0.2
    }
  
  if (save_nodedf){
    write.csv(nodedf, paste("nodedf", name, ".csv", sep=""))
  }
  
  # make an edge-list with weight attributes
  edgelist<-get.edgelist(ga)
  edf<-data.frame(onode=as.character(edgelist[,1]), dnode=as.character(edgelist[,2]))
  edf$wt<-0.6
  edf$col<-NA
  edf$col<-nodedf$color[match(edf$dnode, nodedf$name)]
  for (i in 1:length(edf[,1]))
    if (edf$dnode[i] %in% alldata$LION){
      edf$col[i]<-"lightblue3"      
    }
  #edf$onodenew<-as.character(edf$onode); edf$dnodenew<-as.character(edf$dnode)
  #edf$onodenew[edf$col=="lightblue3"]<-NA
  #for (i in 1:length(edf[,1])){
  #  if (is.na(edf$onodenew[i])){
  #    edf$onodenew[i]<-edf$dnode[i]
  #    edf$dnodenew[i]<-edf$onode[i]
  #  }
  #}
  nodedf$axis<-as.integer(nodedf$axis)
  hive1<-mod.edge2HPD(edge_df= edf[,1:2], edge.weight=edf$wt, edge.color=edf$col, 
                      node.color= nodedf[,c("name", "color")], 	node.size=nodedf[, c("name", "size")], 
                      node.radius= nodedf[,c("name", "radius")], node.axis=nodedf[,c("name", "axis")] )
  
  par(mar=c(0,0,0,0))
  plotHive(hive1, bkgnd = "white", axLabs=c("L", "P"), 
           axLab.pos = 5, axLab.gpar = gpar(col = "black"))
  
  rm(ga, nodedf)
}
