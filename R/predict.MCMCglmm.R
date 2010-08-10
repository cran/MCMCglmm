"predict.MCMCglmm"<-function(object, newdata=NULL, marginal=object$Random$formula, type="response", interval="none", level=0.95, ...){

  warning("predict.MCMCglmm is still developmental - be careful")

  if(type%in%c("response", "terms")==FALSE){stop("type must be response or terms")}
  if(interval%in%c("none", "confidence", "prediction")==FALSE){stop("interval must be none, confidence or prediction")}

  rcomponents<-split.direct.sum(as.character(object$Random$formula)[2])
  mcomponents<-split.direct.sum(as.character(marginal)[2])

  if(any(mcomponents%in%rcomponents==FALSE)){stop("marginal formula does not correspond to model formula")}

  marginal<-rep(rep(as.numeric(rcomponents%in%mcomponents), object$Random$nrt), object$Random$nfl)
  
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
    st<-c(1,cumsum(rep(object$Random$nrl, object$Random$nfl))+1)
    st<-st[-length(st)]
    end<-cumsum(rep(object$Random$nrl, object$Random$nfl)) 
    comp<-rep(1:length(object$Random$nfl), object$Random$nfl)
    keep<-unlist(mapply(st[which(marginal==0)], end[which(marginal==0)], FUN=":"))
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

    vpred<-matrix(0,dim(object$X)[1], sum(object$Random$nfl[which(marginal==1)]^2, na.rm=T)+sum(object$Residual$nfl^2))  # obtain (co)variances for all marginised effects and error term

    cnt<-0

    if(any(marginal==1)){

      st<-st[which(marginal==1)]
      end<-end[which(marginal==1)]
      comp<-comp[which(marginal==1)]

      for(i in 1:length(st)){   # iterate through G-terms        
        for(j in 1:length(st)){   
          if(comp[i]==comp[j]){
            cnt<-cnt+1
            vpred[,cnt]<-diag(object$Z[,st[i]:end[i]]%*%t(object$Z[,st[j]:end[j]]))
          }
        }
      }
    }

    comp<-rep(1:length(object$Residual$nfl), object$Residual$nfl)
 
    for(i in 1:length(comp)){   # iterate through R-terms        
      for(j in 1:length(comp)){    
        if(comp[i]==comp[j]){
         cnt<-cnt+1
         vpred[,cnt][which(object$error.term==i & object$error.term==j)]<-1 
        }
      }  
    } 
 
    keep<-which(rep(rep(as.numeric(rcomponents%in%mcomponents), object$Random$nrt), object$Random$nfl^2)==1)

    keep<-c(keep, which(rep(rep(rep(1, length(object$Residual$nrt)), object$Residual$nrt), object$Residual$nfl^2)==1)) 

    postvar<-t(apply(object$VCV[,keep,drop=FALSE],1, function(x){(vpred%*%x)}))
  }

  post.pred<-t(apply(object$Sol,1, function(x){(W%*%x)@x}))

  if(interval=="prediction"){
   post.pred<-matrix(rnorm(prod(dim(post.pred)), post.pred, sqrt(postvar)), dim(post.pred)[1],  dim(post.pred)[2])
  }

  if(type=="response"){
    if(any(object$family%in%c("poisson","cenpoisson","multinomial","categorical","gaussian","cengaussian", "ordinal")==FALSE)){
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

  }

  pred<-matrix(colMeans(post.pred), dim(post.pred)[2],1)
    
  if(interval!="none"){
    pred<-cbind(pred, coda::HPDinterval(mcmc(post.pred), prob=level))   
    colnames(pred)<-c("fit", "lwr", "upr")
  }
  
  rownames(pred)<-1:dim(pred)[1]

  return(pred)
}  
