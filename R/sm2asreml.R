"sm2asreml"<-function(A=NULL, rownames=NULL){
  ginv<-data.frame(Row=A@i+1, Column=rep(1:length(A@p[-1]), diff(A@p)), Ainverse=A@x)
  ginv<-ginv[order(ginv$Row),]
  ginv<-ginv[which(ginv$Row>=ginv$Column),]
  attr(ginv, "rowNames")<-rownames
  return(ginv)
}

