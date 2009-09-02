"posterior.mode"<-function(x, ...){

  find.mode<-function(x,...){
    dx<-density(x, ...)
    dx$x[which.max(dx$y)]
  }

   if(is.mcmc(x)==FALSE){
     stop("posterior.mode expecting mcmc object")
   }else{
     apply(as.matrix(x), 2, find.mode)
   }

 }

