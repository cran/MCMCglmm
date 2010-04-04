spl<-function(x, k=10, knots=NULL, type="LRTP"){

  if(is.null(knots)){
    knots<-quantile(x, 1:k/(k+1), na.rm=T)
  }

  if(type=="LRTP"){  # low-rank thin plate slpine
    Z<-outer(x, knots, function(x,z){abs(x-z)^3})
    B<-outer(knots, knots, function(x,z){abs(x-z)^3})
  }
  if(type=="cubic"){ # non-cyclic cubic spline
  }
  if(type=="cyc-cubic"){ # cyclic cubic spline
  }
  if(type=="res-cubic"){ # restricted cubic spline
  }

  Bsvd<-svd(B)
  Z<-t(solve(t(Bsvd$v%*%(t(Bsvd$u)*sqrt(Bsvd$d))), t(Z)))

  return(Z)
} 




