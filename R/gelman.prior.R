gelman.prior<-function(formula, data, scale=1, intercept=scale){
  X1<-model.matrix(formula, data)
  X2<-get_all_vars(formula, data)
  X2<-as.data.frame(lapply(X2, function(x){if(is.numeric(x)){scale(x, scale=sd(x)*2*(length(x)-1)/length(x))}else{x}}))
  X2<-model.matrix(formula, data=X2)
  if(all(X2[,1]==1)){
    X2[,-1]<-apply(X2[,-1,drop=FALSE], 2, function(x){if(any(!x%in%c(0,1))){x}else{scale(x, center=sum(x)/length(x), scale=1)}})
  }else{
    X2<-apply(X2, 2, function(x){if(any(!x%in%c(0,1))){x}else{scale(x, center=sum(x)/length(x), scale=1)}})
  }
  P<-solve(t(X1)%*%X1, t(X1)%*%X2)
  I<-diag(nrow(P))*scale
  I[1,1]<-intercept
  P%*%I%*%t(P)*scale
}

