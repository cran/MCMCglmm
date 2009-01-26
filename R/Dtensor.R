"DD" <- function(expr, name) {
   order <- length(name)
   if(order < 1) stop("'order' must be >= 1")
   if(order == 1) D(expr,name[1])
   else DD(D(expr, name[1]), name[-1])
}

"Dtensor"<-function(expr, name=NULL, mu=NULL, m=1, evaluate=TRUE){
	
	if(is.null(mu)){
		if(evaluate){stop("can't evaluate tensor: mu not given")}
	}else{	
      if(is.data.frame(mu)==FALSE){
        mu<-as.data.frame(t(mu)) 
      }
   }	 
   if(is.null(name)){
	name<-all.vars(expr)
    rankx<-length(name)
    comb.pos<-expand.grid(lapply(1:m, function(x){1:rankx}))
    if(evaluate==TRUE){
	  names(mu)<- name
      D<-to.tensor(1:(rankx^m), rep(rankx,m))
      for(i in 1:length(comb.pos[,1])){
        D[i]<-eval(DD(expr, name[unlist(comb.pos[i,])]), mu)
      }
    }else{
     D<-as.list(1:length(comb.pos[,1]))
     for(i in 1:length(comb.pos[,1])){
        D[[i]]<-DD(expr, name[unlist(comb.pos[i,])])
      }
    }
  }else{
    if(evaluate==TRUE){
	  names(mu)<- name
      D<-eval(DD(expr, name), mu)
    }else{
      D<-DD(expr, name)
    }
  }
  return(D)
}


"evalDtensor"<-function(x, mu){

   if(is.data.frame(mu)==FALSE){
     mu<-as.data.frame(t(mu))
   }
  eDtensor<-lapply(x, function(x){eval(x, mu)})  
  rankx<-length(mu)
  if(rankx!=1){
    m<-log(length(x))/log(rankx)
  }else{
    m<-1
  }
  if(m%%1!=0){stop("mu wrong length")}
  to.tensor(unlist(eDtensor), rep(rankx,m))
}
