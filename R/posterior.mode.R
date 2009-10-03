"posterior.mode"<-function(x, adjust=0.01, ...){

   if(is.mcmc(x)==FALSE){
     warning("posterior.mode expecting mcmc object")
   }

  find.mode<-function(x,...){
    dx<-density(x, ...)
    dx$x[which.max(dx$y)]
  }

  apply(as.matrix(x), 2, find.mode, ...)
}

