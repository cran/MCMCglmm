priorformat<-function(prior, start, nfl, meta, diagR){

       if(is.null(prior)){
          prior<-list(V=diag(sum(nfl)), nu=0, fix=as.numeric(meta), alpha.mu=rep(1,sum(nfl)), alpha.V=diag(sum(nfl))*0)
       }else{
         if(is.null(prior$V)){
           stop("V not specified for some prior$G/prior$R elements")
         }
         if(is.matrix(prior$V)==FALSE){
            prior$V<-as.matrix(prior$V)
         }
         if(!is.null(prior$alpha.mu) & diagR>0){
            stop("Parameter exapnded priors not implemented for residual structures")
         }
         if(!is.null(prior$alpha.V) & diagR>0){
           stop("Parameter exapnded priors not implemented for residual structures")
         }
         if(is.null(prior$fix)){
           prior$fix<-0
         }
         if(any(names(prior)%in%c("V", "n", "nu", "alpha.mu", "alpha.V", "fix")==FALSE)){
            paste(paste(names(prior)[which(names(prior)%in%c("V", "n", "nu", "alpha.mu", "alpha.V", "fix")==FALSE)], sep=" "), " are not valid prior specifications for G/R-structures")
         }
         if(diagR==3){     # need to expand trait:units prior to us(trait):units prior      
           if(dim(prior$V)[1]!=1){
             stop("V is the wrong dimension for some prior$G/prior$R elements")
           }
           prior$V<-diag(sum(nfl))*as.numeric(prior$V)
         }else{
           if(any(dim(prior$V)!=sum(nfl))){
             stop("V is the wrong dimension for some prior$G/prior$R elements")
           }
         }
         if(diagR==1){
           if(!all(diag(prior$V)>0)){
             stop("V is not positive definite for some prior$G/prior$R elements")
           }
         }else{
           if(is.positive.definite(prior$V)==FALSE){
             stop("V is not positive definite for some prior$G/prior$R elements")
           }
         }
         if(is.null(prior$alpha.V)){
            prior$alpha.V<-prior$V*0
         }else{
           if(is.matrix(prior$alpha.V)==FALSE){
             prior$alpha.V<-as.matrix(prior$alpha.V)
           }
           if(any(dim(prior$alpha.V)!=dim(prior$V))){
             stop("alpha.V is the wrong dimension for some prior$G/prior$R elements")
           }
           if(is.positive.definite(prior$alpha.V)==FALSE & all(prior$alpha.V==0)==FALSE){
             stop("alpha.V is not positive definite for some prior$G/prior$R elements")
           }
         } 
         if(is.null(prior$alpha.mu)){
            prior$alpha.mu<-matrix(1, nrow(prior$V), 1)
         }else{
           if(is.matrix(prior$alpha.mu)==FALSE){
             prior$alpha.mu<-matrix(prior$alpha.mu, length(prior$alpha.mu), 1)
           }
           if(length(prior$alpha.mu)!=nrow(prior$alpha.V)){
             stop("alpha.mu is the wrong length for some prior$G/prior$R elements")
           }
         } 
         if(prior$fix!=0){
           CM<-prior$V[prior$fix:nrow(prior$V),prior$fix:nrow(prior$V)]
           if(sum(CM!=0)>nrow(prior$V) & prior$fix>1){
             stop("sorry - matrices to be conditioned on must be diagonal")
           }           
           if(prior$fix!=1){
             if(is.null(prior$n)){
               stop("nu not specified for some prior$G/prior$R elements")
             }
           }else{
             prior$nu=1
           }
         }else{
           if(is.null(prior$n)){
             stop("nu not specified for some prior$G/prior$R elements")
           }
         }
       }
       if(is.null(start)){
         if(det(prior$V)<1e-8 & prior$fix==0){
           start<-prior$V+diag(nrow(prior$V))
         }else{
           start<-prior$V
         }
       }else{
         if(is.matrix(start)==FALSE){
           start<-as.matrix(start)
         }	
         if(any(dim(start)!=sum(nfl))){
           stop("V is the wrong dimension for some start$G/start$R elements")
         }
         if(is.positive.definite(start)==FALSE){
           stop(paste("some start$G/start$R structures are not positive definite"))
         }
       }

       pfix<-function(x,y){
         if((prior$fix>y | prior$fix<x) & prior$fix!=0){
           if(prior$fix>y){
             match(prior$fix, x:y, nomatch=0)
           }else{
             match(prior$fix, x:y, nomatch=1)
           }
         }else{
             match(prior$fix, x:y,, nomatch=0)
         }
       }

       prior<-mapply(x=cumsum(nfl)-(nfl-1), y=cumsum(nfl),  function(x,y){list(V=as.matrix(prior$V[x:y,x:y, drop=FALSE]), nu=prior$n, fix=pfix(x,y), alpha.mu=prior$alpha.mu[x:y], alpha.V=as.matrix(prior$alpha.V[x:y,x:y, drop=FALSE]))}, SIMPLIFY=FALSE)
       start<-mapply(x=cumsum(nfl)-(nfl-1), y=cumsum(nfl),  function(x,y){list(start=as.matrix(start[x:y,x:y, drop=FALSE]))}, SIMPLIFY=FALSE)

       return(list(prior=prior, start=start))
}

