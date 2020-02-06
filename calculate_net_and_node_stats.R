calculate_net_and_node_stats<-function(matdata, name){
    #####################################
    # INPUT: 
    # matdata= matrix, symmetric
    # name= filename for saving; .csv is appended below
    # OUTPUT: 
    # 2 csv files: netstats_filename.csv; nodestats_filename.csv
    #####################################
    net.stats = data.frame(matrix(NA, nrow=1, ncol=11, dimnames= list(NULL, 
        c("NumNodes", "MeanDegree", "MeanDegreeUnWt", "VarDegree", 
        "VarDegreeUnWt", "Diameter", "WeightedDiameter","Reciprocity",
        "Transitivity","VertexConnectivity", "WeightedAssortativity"))))
    node.stats = data.frame(matrix(NA, nrow = length(colnames(matdata)), 
        ncol = 8, dimnames = list(NULL, 
        c("NodeID", "DegreeWt", "DegreeUnWt", "Betweenness", "BetweennessUnWt", 
            "Closeness", "Transitivity", "ANND")) ) )
    # make igraph object
    ga <- graph.adjacency(matdata, mode = "upper", weighted = TRUE) # undir/named/wt/self
    ga2 <- graph.adjacency(matdata, mode = "upper", weighted = NULL) # undir/named/unwt
    ga3 <- simplify(ga, remove.multiple = FALSE, remove.loops = TRUE) # undir/named/wt/noself
  
    # calculate node.stats
    node.stats$NodeID = V(ga)$name
    node.stats$DegreeWt = graph.strength(ga)
    node.stats$DegreeUnWt = degree(ga3, mode = "all") # without self loops
    node.stats$Betweenness = betweenness(ga, directed = FALSE)  # weighted
    node.stats$BetweennessUnWt = betweenness(simplify(ga2), directed = FALSE) # unwt betweenness
    node.stats$Closeness = closeness(ga, mode = "all")
    node.stats$ANND = graph.knn(ga3, weights = E(ga3)$weight)$knn
    node.stats$Transitivity = transitivity(ga, type = c("local"))
    net.stats$NumNodes = length(node.stats$NodeID)                                        
    net.stats$MeanDegree = mean(node.stats$DegreeWt)
    net.stats$VarDegree = var(node.stats$DegreeWt)
    net.stats$MeanDegreeUnWt = mean(node.stats$DegreeUnWt)
    net.stats$VarDegreeUnWt = var(node.stats$DegreeUnWt)
    net.stats$Diameter = diameter(ga,weights=(1/E(ga))) 
    net.stats$WeightedDiameter = diameter(ga2) 
    net.stats$Reciprocity = reciprocity(ga)
    net.stats$Transitivity = transitivity(ga,type = c("global"))
    net.stats$WeightedAssortativity = cor(node.stats$DegreeWt, node.stats$ANND)

    write.csv(node.stats, paste("node_stats_", name, ".csv", sep=""))
    write.csv(net.stats, paste("net_stats_", name, ".csv", sep=""))
}


# NOT UPDATED


##########################################################
# Functions called in calculate_net_and_node_stats
##########################################################
bipartite.betweenness.centrality= function(g){
  # code in BiGraph package but not avaliable for my version of R
  # http://www.inside-r.org/packages/cran/biGraph/docs/bipartite.betweenness.centrality
  # This function is the adaption to bipartite graphs as presented in Borgatti and Everett (1997). 
  # Borgatti, S. P. and Everett, M. G. (1997) Network analysis of 2--mode data. Social Networks 19, 243--269.
  
  if (!is.null(V(g)$type)){
    # determine maximal raw scores for both vertex subsets
    if (length(V(g)[type==FALSE])<length(V(g)[type==TRUE])){
      mrs_TRUE <- 2*(length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-1)
    }
    else{
      mrs_TRUE <- 0.5*length(V(g)[type==FALSE])*(length(V(g)[type==FALSE])-1)+0.5*(length(V(g)[type==TRUE])-1)*(length(V(g)[type==TRUE])-2)+(length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-1)
    }
    if (length(V(g)[type==TRUE])<length(V(g)[type==FALSE])){
      mrs_FALSE <- 0.5*length(V(g)[type==TRUE])*(length(V(g)[type==TRUE])-1)+0.5*(length(V(g)[type==FALSE])-1)*(length(V(g)[type==FALSE])-2)+(length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1)
      
    }
    else{
      mrs_FALSE <- 2*(length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1)
    }
    
    # get raw betweenness centrality scores from igraph
    betweenness_rs <- betweenness(g)
    # "bipartite" normalization of scores
    for (i in V(g)){
      if (V(g)[i]$type==TRUE){
        V(g)[i]$betweenness.centrality <- betweenness_rs[i+1]/mrs_TRUE
      }
      else{
        V(g)[i]$betweenness.centrality <- betweenness_rs[i+1]/mrs_FALSE
      }
    }
    # return value as list
    return(list("Bipartite.Betweenness.Centrality"=V(g)$betweenness.centrality))
  }
  else {
    # boolean vertex attribute 'type' is required
    cat("vertex attribute <type> is missing")
  }
}

# works
bipartite.closeness.centrality = function(g){
  if (!is.null(V(g)$type)){
    # determine maximal raw scores for both vertex subsets
    mrs_TRUE <- length(V(g)[type==FALSE]) + 2*length(V(g)[type==TRUE]) - 2
    mrs_FALSE <- length(V(g)[type==TRUE]) + 2*length(V(g)[type==FALSE]) - 2
    # get sum of all geodesic paths for each vertex
    rowsums_shortest_paths <- rowSums(shortest.paths(g))
    # "bipartite" normalization of scores
    for (i in V(g)){
      if (V(g)[i]$type==TRUE){
        V(g)[i]$closeness.centrality <- mrs_TRUE/rowsums_shortest_paths[i+1]
      }
      else{
        V(g)[i]$closeness.centrality <- mrs_FALSE/rowsums_shortest_paths[i+1]
      }
    }
    # return value as list
    return(list("Bipartite.Closeness.Centrality"=V(g)$closeness.centrality))
  }
  else {
    # boolean vertex attribute 'type' is required
    cat("vertex attribute <type> is missing")
  }
}


# doesn't work: 
bipartite.degree.centrality = function(g,loops = FALSE){
    # if boolean vertex attribute <type> is present, 
    # calc bipartite degree centrality, otherwise monopartite deg centrality
    if (!is.null(V(g)$type)){
        for (i in V(g)){
            V(g)[i]$degree.centrality <- degree(g, v = i) / length(
                V(g)[type == !V(g)[i]$type])
        }
        # return value vector as list item
        return(list("Bipartite.Degree.Centrality"=V(g)$degree.centrality))
    }
    else {
        for (i in V(g)){
            if (!loops){
                V(g)[i]$degree.centrality <- degree(g,v = i,loops = FALSE) / (
                    length(V(g)) - 1)
        }
        else{
          V(g)[i]$degree.centrality <- degree(g, v = i,loops = TRUE) / (
              length(V(g)))
        }
      }
      # return value vector as list item
      return(list("Monopartite.Degree.Centrality"=V(g)$degree.centrality))
    }
  }

#calculate_all_bipartite_net_and_node_stats<-function(edgedata, name){
#####################################
# INPUT: 
# matdata
#datatype= "all" OR "FIV"
# OUTPUT: 
# 2 csv files: netstats_filename.csv; nodestats_filename.csv
#####################################
#	require("bipartite")  
#  net.stats= data.frame()
#  node.stats=data.frame(matrix(NA, nrow= length(c(unique(edgedata$lion), unique(edgedata$parasite))), 
#                               ncol=7, dimnames=list(NULL, 
#                                                     c("NodeID","NodeType","Degree", 
#                                                        "Bipart_Degree", "Bipart_Betweenness", 
#                                                       "Bipart_Closeness", "ANND"

# make bipartite object and define it as bipartite                              
#  temp<-frame2webs(all, varnames=c("parasite", "lion"))       
#  igraph<-graph.data.frame(edgedata)
#  V(igraph)$type <- V(igraph)$name %in% edgedata[,1]
#c(rep(TRUE, length(unique(edgedata$lion))), rep(FALSE, length(unique(edgedata$parasite)) ))  
#  bipartite.projection(igraph)

# calculate node.stats
#  node.stats$NodeID = V(igraph)$name
#  node.stats$NodeType= V(igraph)$type
#  node.stats$Degree = degree(igraph)
#  node.stats$Bipart_Degree = bipartite.degree.centrality(igraph)
#  node.stats$Bipart_Betweenness = bipartite.betweenness.centrality(igraph)  # This function is the adaption to bipartite graphs as presented in Borgatti and Everett (1997).
#  node.stats$Bipart_Closeness = bipartite.closeness.centrality(igraph)
#  node.stats$ANND = graph.knn(igraph)$knn
#  node.stats$modularity

#  net.stats$connectance = 
#  net.stats$nestedness = 
#  net.stats$modularity = 
#  net.stats$meandegree_lion = 
#  net.stats$meandegree_parasite = 
#}


#calculate_bipartite_net_and_node_stats<-function(matdata, datatype=="all", name){
#####################################
# INPUT: 
# matdata= matrix, symmetric
# name= filename for saving; .csv is appended below
# OUTPUT: 
# 2 csv files: netstats_filename.csv; nodestats_filename.csv
# see http://toreopsahl.com/tnet/two-mode-networks/
#####################################
#    
#    net.stats= data.frame(matrix(NA, nrow=4, ncol=8, dimnames= 
#                                   list(NULL, c("NumNodes", "NumParasites", "NumHosts", "MeanDegree"))))
#    
#    node.stats=data.frame(matrix(NA, nrow = length(colnames(matdata)), 
#                                 ncol = 6, dimnames = list(NULL, 
#                                                           c("NodeID", "NodeType" "Degree", "Betweenness",
#                                                             "Closeness", "Transitivity", "ANND")) ) )

#    ga<-graph.adjacency(matdata, mode="upper")  
#    is.simple(ga)   
#    if(datatype=="all"){
# set lion verticies=TRUE, parasite verticies to FALSE
#      V(ga)$type<-c(rep(TRUE, length(colnames(matdata))-length(colnames(all_parnet))),
#                    rep(FALSE, length(colnames(all_parnet))))                                                    
#    } else if (datatype=="FIV"){
#      V(ga)$type<-c(rep(TRUE, length(colnames(matdata))-length(colnames(FIVneg_parnet))),
#                    rep(FALSE, length(colnames(FIVneg_parnet))))                        
#    } else {
#      print("datatype is not correct")
#    }                                               
#    is.bipartite(ga)                                            

#proj_net <- bipartite.projection(ga)
#set.seed(123)  # for reproducible plot
#plot(ga,vertex.color=ifelse(V(ga)$type,"green","red"))

# centrality measures borrowed from bipartite package
#"NodeID", "NodeType" "Degree", "Betweenness",
#"Closeness", "Transitivity", "ANND"
#  node.stats$NodeID<- V(ga)
# node.stats$NodeType<-c(rep("Lion", length(V(ga)$type==TRUE), rep("Parasite", length(V(ga)$type==FALSE)))
#  node.stats$Degree<-bipartite.degree.centrality(ga)
#  node.stats$Betweenness<-  bipartite.betweenness.centrality(ga)
#   node.stats$Closeness<-bipartite.closeness.centrality(ga)

# transitivity needs done in tnet package (http://toreopsahl.com/tnet/two-mode-networks/clustering/)
# clustering: Opsahl (2012) proposed a new coefficient for two-mode networks that measures 
# closure among three nodes from the primary node set instead of only two primary nodes (e.g., Robins and Alexander, 2004).      
#}
