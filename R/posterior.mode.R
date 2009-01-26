"posterior.mode"<-function(x){
 dx<-density(x)
 dx$x[which.max(dx$y)]
 }

