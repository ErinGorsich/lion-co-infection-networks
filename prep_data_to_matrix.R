make_edgelist= function(data){
  temp<- NA
  start<-data[1,]
  edgelist<-data.frame(parasite=colnames(temp[which(start=="Positive")]), 
                       lion=rep(data$LION[1], length(colnames(temp[which(start=="Positive")]))))
  rm(start)
  for (i in 2:length(data[,1])){
    temp<-data[i,]
    df<-data.frame(parasite=colnames(temp[which(temp=="Positive")]), 
                   lion=rep(data$LION[i], length(colnames(temp[which(temp=="Positive")]))))
    edgelist<-rbind(edgelist, df)
    rm(df, temp)
  }
  return(edgelist)
}

make_matrix= function(data){
  tmat<-matrix(NA, nrow=length(data[,1]), ncol=21,  # makes lions columns 
               dimnames=list(data$LION, colnames(colnames(plotdata[-1])))) # all infections minus FIV
  for (i in 1:nrow(tmat)){
    for (j in 1:ncol(tmat)){
      ifelse(data[i,11+j]=="Positive", tmat[i,j] <- 1, tmat[i, j] <- 0)  
    }
  }
  return(tmat)
}

make_lionmatrix= function(data){
  tmat<-matrix(NA, ncol=length(data[,1]), nrow=21,  # 21 rows here bc with FIV 
               dimnames=list(colnames(colnames(plotdata[-1])), data$LION)) # all infections minus FIV
  for (i in 1:nrow(tmat)){
    for (j in 1:ncol(tmat)){
      ifelse(data[j,11+i]=="Positive", tmat[i,j] <- 1, tmat[i, j] <- 0)  
    }
  }
  return(tmat)
}

make_all_matrix= function(data){
  tmat<-matrix(NA, nrow=length(data[,1]), ncol=22,  # 22 rows here bc with FIV 
               dimnames=list(data$LION, colnames(plotdata))) # all infections minus FIV
  for (i in 1:nrow(tmat)){
    for (j in 1:ncol(tmat)){
      ifelse(data[i,10+j]=="Positive", tmat[i,j] <- 1, tmat[i, j] <- 0)  
    }
  }
  return(tmat)
}
make_all_lionmatrix= function(data){
  tmat<-matrix(NA, ncol=length(data[,1]), nrow=22,  # 22 rows here bc with FIV 
               dimnames=list(colnames(plotdata), data$LION)) # all infections minus FIV
  for (i in 1:nrow(tmat)){
    for (j in 1:ncol(tmat)){
      ifelse(data[j,10+i]=="Positive", tmat[i,j] <- 1, tmat[i, j] <- 0)  
    }
  }
  return(tmat)
}

make_all_bipartite_matrix= function(data){
  tmat<-matrix(0, ncol=length(c(as.character(data$LION), colnames(plotdata)) ), 
               nrow=length(c(as.character(data$LION), colnames(plotdata))) , dimnames=list(
                 c(as.character(data$LION), colnames(plotdata)), 
                 c(as.character(data$LION), colnames(plotdata)) ) )
  # lion to lion transitions are 0
  # parastie to parasite connections are 0
  # fill in rows & columns= from lion (1:111) to parasite (112:136)
  for (i in 1:length(data[,1])){   # for each lion... 
      for (j in 1:length(plotdata[1,])){  # for each parasite
            if(data[i, 10+j]=="Positive") {
               tmat[i, 111+j] <- 1
               tmat[111+j, i] <-1}  # symmetric
            else {
              tmat[i, 111+j] <- 0
              tmat[111+j, i] <- 0
            }  
        }                          
  }
 return(tmat)
}

make_bipartite_matrix= function(data){
  tmat<-matrix(0, ncol=length(c(as.character(data$LION), colnames(plotdata)) )-1, 
               nrow=length(c(as.character(data$LION), colnames(plotdata)))-1 , dimnames=list(
                 c(as.character(data$LION), colnames(plotdata)[-1]), 
                 c(as.character(data$LION), colnames(plotdata)[-1]) ) )
  # fill in rows and columns of tmat with presence absence of parasite
  for (i in 1:length(data[,1])){   # for each lion... 
    for (j in 1:c(length( colnames(tmat))-length(data[,1]) ) ){  # for each parasite
      if(data[i, 11+j]=="Positive") {

        tmat[i, c(length(data[,1]))+j] <- 1
        tmat[c(length(data[,1]))+j, i] <-1}  # symmetric
      else {
        tmat[i, c(length(data[,1]))+j] <- 0
        tmat[c(length(data[,1]))+j, i] <- 0
      }  
    }                          
  }
  # keep rows and columns with no contacts... yes
  return(tmat)             
}
