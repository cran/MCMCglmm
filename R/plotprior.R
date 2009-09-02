"plotprior"<-function(prior, cor=FALSE, n=1000, ask = dev.interactive(), ...){

  if(is.null(prior$R)==FALSE){
    nV<-length(prior$R$V)
  }else{
    nV<-0
  }
  if(is.null(prior$G)==FALSE){
    nV<-nV+sum(unlist(lapply(prior$G, function(x){length(x$V)})))
  }

  PD<-matrix(0, n, 0)

  if(is.null(prior$R)==FALSE){
    if(is.matrix(prior$R$V)==FALSE){stop("prior V for R-structure is not a matrix")}
	if(isSymmetric(prior$R$V)==FALSE){stop("prior V for R-structure is not symmetric")}
    if(is.positive.definite(prior$R$V)==FALSE){stop("prior V for R-structure is not positive definite")}
    nd<-dim(as.matrix(prior$R$V))[1]
    if(prior$R$n<nd){
      warning("prior for the R-structure is improper cannot be simulated")
    }else{
      RS<-rIW(prior$R$V, nu=prior$R$n, prior$R$fix, n=n)
      colnames(RS)<-paste("R[", outer(1:nd,1:nd, paste, sep=","), "]", sep="") 
        if(cor){
          v.elements<-diag(matrix(1:(nd^2),nd,nd))
          sd<-sqrt(RS[,v.elements])
          sd<-t(apply(sd, 1, function(x){outer(x, x)}))
          RS<-RS/sd
          RS<-RS[,-v.elements]
        }
      PD<-cbind(PD,RS)
    }
  }
  if(is.null(prior$G)==FALSE){
    for(i in 1:length(prior$G)){
      if(is.matrix(prior$G[[i]]$V)==FALSE){stop(paste("prior V for G-structure ", i, " is not a matrix", sep=""))}
		if(isSymmetric(prior$G[[i]]$V)==FALSE){stop(paste("prior V for G-structure ", i,  " is not symmetric", sep=""))}
		if(is.positive.definite(prior$G[[i]]$V)==FALSE){stop(paste("prior V for G-structure ", i, " is not positive definite", sep=""))}
      nd<-dim(as.matrix(prior$G[[i]]$V))[1]
      if(prior$G[[i]]$n<nd){
        warning(paste("prior for G-structure ", i,  " is improper and cannot be simulated", sep=""))
      }else{
        RS<-rIW(as.matrix(prior$G[[i]]$V), nu=prior$G[[i]]$n, prior$G[[i]]$fix, n=n)
        colnames(RS)<-paste(names(prior$G)[i], "[", outer(1:nd,1:nd, paste, sep=","), "]", sep="")
        if(cor){
          v.elements<-diag(matrix(1:(nd^2),nd,nd))
          sd<-sqrt(RS[,v.elements])
          sd<-t(apply(sd, 1, function(x){outer(x, x)}))
          RS<-RS/sd
          RS<-RS[,-v.elements]
        }
        PD<-cbind(PD,RS)
      }
    }
  }
  plot.dim<-min(ceiling(sqrt(dim(PD)[2])),3)
  par(mfrow=c(plot.dim, plot.dim))
  for(j in 1:dim(PD)[2]){
   xlim<-list(...)$xlim
   if(is.null(xlim)){
     hist(PD[,j],main=colnames(PD)[j], ..., xlab="")
   }else{
     rm.extremes<-which(PD[,j]<xlim[1] | PD[,j]>xlim[2])
     if(length(rm.extremes)>0){
       hist(PD[,j][-rm.extremes],main=colnames(PD)[j], ..., xlab="")
     }
   }
   if(j%%(plot.dim^2)==0 & ask & j!=dim(PD)[2]){
     ask(msg = "Hit <RETURN> to see next plot: ")  
   }
  }
}


