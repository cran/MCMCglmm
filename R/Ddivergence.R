Ddivergence<-function(CA=NULL, CB=NULL, n=10000){

  if(dim(CA)[1]!=dim(CB)[1] | dim(CA)[2]!=dim(CB)[2] | dim(CA)[1]!=dim(CA)[2]){
     stop("matrices must be the same dimension and square")
  }

  xi<-mvrnorm(n, rep(0,dim(CA)[1]), cov2cor(CA))
  fx<-dmvnorm(xi, rep(0,dim(CA)[1]), cov2cor(CA))
  gx<-dmvnorm(xi, rep(0,dim(CA)[1]), cov2cor(CB))

  mean(sqrt(0.5*((fx-gx)^2)/(fx+gx)))

}

