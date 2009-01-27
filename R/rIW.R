"rIW"<-function(nu,V,split=NULL, n=1){
 
  if(is.matrix(V)==FALSE | dim(V)[1]!=dim(V)[2]){stop("V must be a square matrix")}
  if(is.positive.definite(V)==FALSE){stop("V must be positive definite")}
  if(nu<=0){stop("nu must be greater than zero")}
  if(is.null(split)){
    split=-998
  }else{
    if(split%%1!=0 | split<1 | split>dim(V)[1]){
      stop(paste("split must be an integer between 1 and ", dim(V)[1]))
    }
  }
  if(split==1){
    if(n==1){
      return(matrix(rep(V,n), dim(V)[1], dim(V)[1]))
    }else{
      return(t(matrix(rep(V,n), dim(V)[1]^2, n)))
    }
  }else{

    output<-.C("rIW",
      as.double(nu),
      as.double(c(V)),
      as.integer(dim(V)[1]),
      as.integer(split-1),
      as.integer(n),
      as.double(matrix(0,n,dim(V)[1]^2))
    )
    if(n==1){
      return(matrix(output[[6]], dim(V)[1], dim(V)[1]))
    }else{
      return(t(matrix(output[[6]], dim(V)[1]^2, n)))
    }
  }
}
