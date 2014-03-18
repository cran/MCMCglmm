"predict.MCMCglmm"<-function(object, newdata=NULL, marginal=object$Random$formula, type="response", interval="none", level=0.95, ...){

 # warning("predict.MCMCglmm is still developmental - be careful")

  if(type%in%c("response", "terms")==FALSE){stop("type must be response or terms")}
  if(interval%in%c("none", "confidence", "prediction")==FALSE){stop("interval must be none, confidence or prediction")}
  if(!is.null(marginal)){
    if(class(marginal)!="formula"){stop("marginal should be NULL or a formula")}
  }

  rcomponents<-split.direct.sum(as.character(object$Random$formula)[2])
  mcomponents<-split.direct.sum(as.character(marginal)[2])

  if(any(mcomponents%in%rcomponents==FALSE)){stop("marginal formula does not correspond to model formula")}

  marginal<-rep(as.numeric(rcomponents%in%mcomponents), object$Random$nrt)

  if(is.null(newdata)==FALSE){
   stop("sorry newdata not implemented yet")
  }
  if(any(marginal==0) & dim(object$Sol)[2]==dim(object$X)[2]){
     stop("posterior distribution of random effects not saved: pass pr=TRUE to MCMCglmm")
  }
  if(any(marginal==0) & is.null(object$Z)){
     stop("random effect design matrix not saved: pass saveZ=TRUE to MCMCglmm")
  }

  if(is.null(object$X) & is.null(newdata)){
     stop("either specify saveX=TRUE when model fitting or pass a data frame")
  }

  if(is.null(object$Random$nfl)==FALSE){  # there are random effects
    st<-c(1,cumsum(object$Random$nrl*object$Random$nfl)+1)  # starting column for random effects of each component 
    st<-st[-length(st)]
    end<-cumsum(object$Random$nrl*object$Random$nfl)        # ennding column for random effects of each component 
    keep<-unlist(mapply(st[which(marginal==0)], end[which(marginal==0)], FUN=":"))    # random effects to be kept
    v.terms<-as.factor(rep(1:sum(object$Random$nfl), rep(object$Random$nrl,object$Random$nfl))) # denotes the variance term that the column of Z refers to.
    rm.v<-rep(marginal==0, object$Random$nfl)
  }else{
    keep<-NULL
  }

  object$Sol<-object$Sol[,c(1:object$Fixed$nfl, object$Fixed$nfl+keep),drop=FALSE]

  W<-cBind(object$X, object$Z)
  W<-W[,c(1:object$Fixed$nfl, object$Fixed$nfl+keep), drop=FALSE]
  
  if((type=="response" & any(object$family!="gaussian" & object$family!="cengaussian")) | interval=="prediction"){
     
    # need to calculate variance when
    # i) prediction to be made on data scale for non-gaussian data 
    # ii) obtaining prediction intervals
    pos.error<-sum(object$Random$nfl^2)+cumsum(unlist(sapply(object$Residual$nfl, simplify=FALSE, function(x){c(1,rep(x+1,x-1))})))

    postvar<-t(apply(object$VCV, 1, function(x){x[pos.error[object$error.term]]}))

    if(is.null(object$Random$nfl)==FALSE){  # there are random effects

       same.block<-which(outer(rep(1:length(object$Random$nfl), object$Random$nfl), rep(1:length(object$Random$nfl), object$Random$nfl), "=="))
       # gets positions in V that could have non-zero (co)variances

       if(nlevels(v.terms)>1){

         vpred<-apply(object$Z, 1, function(x){
           nz<-which(x!=0)
           ZZ<-crossprod(t(x[nz]))
           M<-model.matrix(~v.terms[nz]-1)
           M[,which(rm.v)]<-0
           ZZ<-t(M)%*%ZZ%*%M
           ZZ[same.block]
         })

      }else{
        vpred<-apply(object$Z, 1, function(x){
          nz<-which(x!=0)
          ZZ<-crossprod(t(x[nz]))
          M<-matrix(1,1,1)
          M[,which(rm.v)]<-0
          ZZ<-t(M)%*%ZZ%*%M
          ZZ[same.block]
        })
      }

      # sum((zz')*V) this gives the expected variance due to the random effects.
      # zeroing out rows/columns of (zz') that pertian to non-marginalised (co)variances gives the variance due marginalised components
      # removing off-diagonal elements of (zz') that pertian to random effects in different components reduces the storage space.

      postvar<-postvar+object$VCV[,1:sum(object$Random$nfl^2), drop=FALSE]%*%vpred
    }
  }
  post.pred<-t(apply(object$Sol,1, function(x){(W%*%x)@x}))

  if(interval=="prediction"){
   post.pred<-matrix(rnorm(prod(dim(post.pred)), post.pred, sqrt(postvar)), dim(post.pred)[1],  dim(post.pred)[2])
  }

  if(type=="response"){
    if(any(object$family%in%c("poisson","cenpoisson","multinomial","categorical","gaussian","cengaussian", "ordinal", "threshold")==FALSE)){
      stop("sorry - prediction on data scale not implemented for this family")
    }
    if(any(object$family%in%c("poisson","cenpoisson"))){
      keep<-which(object$family%in%c("poisson","cenpoisson"))
      if(interval=="prediction"){
         post.pred[,keep]<-exp(post.pred[,keep])
      }else{
         post.pred[,keep]<-exp(post.pred[,keep]+0.5*postvar[,keep])
      }
    }
    if(any(object$family%in%c("multinomial","categorical"))){
      c2 <- (16 * sqrt(3)/(15 * pi))^2
      keep<-which(object$family%in%c("multinomial","categorical"))
      if(interval=="prediction"){
         post.pred[,keep]<-plogis(post.pred[,keep])
      }else{
         post.pred[,keep]<-plogis(post.pred[,keep]/sqrt(1 + c2 * postvar[,keep]))
      }
    }
    if(any(object$family%in%c("ordinal"))){      
      keep<-which(object$family%in%c("ordinal"))
      CP<-cbind(-Inf, 0, object$CP, Inf)
      q<-matrix(0,dim(post.pred)[1], length(keep))
      if(interval=="prediction"){         
        for(i in 2:(dim(CP)[2]-1)){
          q<-q+(pnorm(CP[,i+1]-post.pred[,keep])-pnorm(CP[,i]-post.pred[,keep]))*(i-1)
        }
      }else{
        for(i in 2:(dim(CP)[2]-1)){
           q<-q+(pnorm(CP[,i+1]-post.pred[,keep],0,sqrt(postvar[,keep]+1))-pnorm(CP[,i]-post.pred[,keep],0,sqrt(postvar[,keep]+1)))*(i-1)
        }
      }
      post.pred[,keep]<-q
      rm(q) 
    }
    if(any(object$family%in%c("threshold"))){      
      keep<-which(object$family%in%c("threshold"))
      if(interval=="prediction"){         
           post.pred[,keep]<-post.pred[,keep]>0
      }else{
           post.pred[,keep]<-pnorm(post.pred[,keep],0,sqrt(postvar[,keep]))
      }
    }
  }
  pred<-matrix(colMeans(post.pred), dim(post.pred)[2],1)
    
  if(interval!="none"){
    pred<-cbind(pred, coda::HPDinterval(mcmc(post.pred), prob=level))   
    colnames(pred)<-c("fit", "lwr", "upr")
  }
  
  rownames(pred)<-1:dim(pred)[1]

  return(pred)
}  
