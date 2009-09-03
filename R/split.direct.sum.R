split.direct.sum<-function(x){

  if(is.na(x)){
   return(NULL)
  }else{
    openB<-gregexpr("\\(", x)[[1]]
    closeB<-gregexpr("\\)", x)[[1]]
    plus<-gregexpr("\\+", x)[[1]]
    internals<-as.matrix(mapply(function(x,y){(plus>x & plus<y)}, x=openB, y=closeB))
    rterms<-strsplit(x, "")[[1]]
    rterms[plus[which(rowSums(internals)!=0)]]<-"leaveMCMCleave"
    rterms<-paste(rterms, collapse="")
    rterms<-strsplit(rterms, " *\\+ *")[[1]]
    rterms<-gsub("leaveMCMCleave", "+", rterms)
    return(rterms)
  }
}
