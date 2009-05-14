"posterior.mode"<-function(x){

   if(is.mcmc(x)==FALSE){
     stop("posterior.mode expecting mcmc object")
   }else{
     apply(as.matrix(x), 2, function(x){dx<-density(x)
       dx$x[which.max(dx$y)]})
   }

 }

