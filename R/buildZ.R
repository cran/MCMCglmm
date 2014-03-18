buildZ<-function(x, data, formZ=TRUE, nginverse=NULL, augmiss=FALSE){

  vtype="idh"   # form of covariances for the random effets of a random term classified by some factor 
  rtype="iid"   # form of covariances between random effets of random terms

  if(length(grep("^us\\(", x))>0){
    vtype<-"us"
  }
  if(length(grep("^corh\\(", x))>0){
    vtype<-"corh"
    stop("sorry - corh not yet implemented")
  }
  if(length(grep("^cor\\(", x))>0){
    vtype<-"cor"
    stop("sorry - cor not yet implemented. Note that cor as used in versions <2.18 has now been replaced by corg")
  }
  if(length(grep("^corg\\(", x))>0){
    vtype<-"corg"
  }
  if(length(grep("^corgh\\(", x))>0){
    vtype<-"corgh"
  }
  if(length(grep("^idv\\(", x))>0){
    vtype<-"idv"
  }
  if(length(grep("(^|:)mm\\(", x))>0){
    rtype<-"mm"
  }
  if(length(grep("(^|:)str\\(", x))>0){
    rtype<-"str"
  }

  fformula<-gsub("^(us|corg|corgh|corh|cor|idh|idv)\\(", "", x)

  if(grepl("^str\\(|^mm\\(", fformula)){
    fformula<-paste("):", fformula, sep="")
  }

  openB<-gregexpr("\\(", fformula)[[1]]
  closeB<-gregexpr("\\)", fformula)[[1]]

  if(openB[1]!=-1){
    while(openB[1]<closeB[1]){
      openB<-openB[-1]
      closeB<-closeB[-1]
      if(length(closeB)==0 | length(openB)==0){break}
    }
  }

  rterms<-substr(fformula,closeB[1]+2,nchar(fformula))
  rterms<-gsub("(^|:)(str|mm|)\\(", "", rterms)
  rterms<-gsub("\\)$", "", rterms)

  if(rtype!="iid" & any(grepl(":", rterms))){
    stop("interactions not permitted in str and mm structures")
  }
  rterms<-strsplit(rterms, ":| \\+ ")[[1]]

  if(is.null(nginverse)==FALSE){
    Aterm<-na.omit(match(rterms,nginverse))[1]
    if(is.na(Aterm)){
      Aterm<-0
    }
  }else{
    Aterm<-0
  }

  fformula<-substr(fformula,1,closeB-1)
  fformula<-gsub(" ", "", fformula)

  idv.vnames<-paste(fformula, collapse="")

  if(fformula!=""){fformula<-c("+", fformula)}

  fformula<-as.formula(paste("~", "-1",  paste(fformula, collapse=""), sep=""))

  orig.na.action<-options("na.action")[[1]]
  options("na.action"="na.pass")

  X<-model.matrix(fformula, data)
  Xterms<-names(attr(X, "contrasts"))

  X[which(is.na(X))]<-1e-64        # structural zeros needed when us() because may be things like slope terms in leg for animal

  options("na.action"=orig.na.action)

  X<-as(X, "sparseMatrix")

  if(is.null(Xterms)==FALSE){
    colnames(X)<-substr(colnames(X),nchar(Xterms)+1, nchar(colnames(X)))
  }

  nfl<-ncol(X)

  if(formZ){

    ZZ<-list()

    for(k in 1:max(1,(1-(rtype=="iid"))*length(rterms))){  # iterate through Rterms or evalute once if iid

      if(rtype=="iid"){
        select.terms<-1:length(rterms)
      }else{
        select.terms<-k
      } 

      if(length(rterms)==0){
        rfactor<-rep(as.factor(1), dim(data)[1]) 	
      }else{
        if(any(rterms%in%colnames(data)==FALSE)){stop(paste("object", paste(rterms, collapse=" and "), "not found"))}
        rfactor<-interaction(data[,rterms[select.terms]], drop=(Aterm==0 & rtype=="iid"))
        # random effects are dropped if they are not animal or part of a mm or str structure 
      }

      nrl<-nlevels(rfactor)
      nd<-length(rfactor)
              
      data_pos<-as.numeric(rfactor)
      ZZ[[k]]<-Matrix(0,nd,nrl)
      ZZ[[k]][,1][2]<-1              # force it out of vbeing upper triangle!!!!
      ZZ[[k]]@p<-as.integer(c(0,cumsum(table(rfactor))))    
      cnt<-0
      for(j in 1:nrl){
        hit<-which(data_pos==j)
        hit<-hit-nd*(ceiling(hit/nd)-1)
        if(length(hit)>0){
          ZZ[[k]]@i[cnt+1:length(hit)]<-as.integer(hit-1)
          cnt<-cnt+length(hit)
        }
      }
      ZZ[[k]]@x<-rep(1,length(ZZ[[k]]@i))

      colnames(ZZ[[k]])<-paste(paste(rterms[select.terms], collapse=":"),levels(rfactor), sep=".")
      if(any(data$MCMC_dummy==0 & is.na(rfactor))){warning("missing values in random predictors")}

    }

    if(nfl==0){

      nfl<-1  

      for(k in 1:max(1,(1-(rtype=="iid"))*length(rterms))){  # iterate through Rtrems or evalute once if iid

        missing<-which(diff(ZZ[[k]]@p)==0)  # non-represented random effects (or zero covariates for RR)

        #################################################################################################
        # It woule be more efficient to allow complete sparse columns in Z - but need to change C code! #
        #################################################################################################
 
        if(length(missing)>0){        
          if(Aterm==0 & rtype=="iid"){                       # if not associated with ginv drop terms     
            ZZ[[k]]<-ZZ[[k]][,-missing,drop=FALSE]
            nrl[i]<-ncol(ZZ[k])              
          }else{                              # if associated with ginv set arbitrary row to structural zero   
            for(j in 1:length(missing)){
              ZZ[[k]][j,missing[j]]<-1e-64     
            }
          }
        }
        if(k>1){
          if(any(levels(data[,rterms[select.terms]])!=levels(data[,rterms[1]]))){
            stop("terms involved in mm/str structures must have identical levels")
          }
          if(rtype=="mm"){
            colnames(ZZ[[k]])<-colnames(ZZ[[1]])
            ZZ[[1]]<-ZZ[[1]]+ZZ[[k]]
          }
          if(rtype=="str"){
            ZZ[[1]]<-cBind(ZZ[[1]],ZZ[[k]])
          }
        }
      }
      vnames<-paste(rterms, collapse=":")
      if(rtype=="mm"){
         vnames<-paste(rterms, collapse="+")
      }
      if(rtype=="str"){
        vnames<-apply(expand.grid(rterms,rterms), 1, paste, collapse=".")
        nfl<-nfl*length(rterms)
        vtype<-rep("us", length(vtype))
      }
  
      return(list(Z=ZZ[[1]], nfl=nfl, nrl=nrl, Aterm=Aterm, vtype=vtype, vnames=vnames, ordering=NULL, trait.ordering=NULL))

    }else{

      if(vtype=="idh"){
        nrl<-rep(nrl, nfl)
      }
      
      for(i in 1:nfl){
    
        Xtmp<-Matrix(0,nd,nd)
        Xcol<-X[,i,drop=FALSE]
        Xtmp@p<-as.integer(rep(0:c(length(Xcol@i)-1),diff(c(0,Xcol@i+1))))
        Xtmp@p[(length(Xtmp@p)+1):(nd+1)]<-length(Xcol@i)
        Xtmp@i<-Xcol@i
        Xtmp@x<-Xcol@x

        Z<-list()

        for(k in 1:max(1,(1-(rtype=="iid"))*length(rterms))){  # iterate through Rtrems or evalute once if iid

          if(rtype=="iid"){
            select.terms<-1:length(rterms)
          }

          Z[[k]]<-Xtmp%*%ZZ[[k]]
          colnames(Z[[k]])<-paste(paste(rterms[select.terms], collapse=":"), colnames(X)[i], colnames(Z[[k]]),sep=".") 

          missing<-which(diff(Z[[k]]@p)==0)

          #################################################################################################
          # It woule be more efficient to allow complete sparse columns in Z - but need to change C code! #
          #################################################################################################

          if(length(missing)>0){
            if(vtype=="idh" | vtype=="idv"){
              if(Aterm==0 & rtype=="iid"){                     # if univariate G-structure and no Ginv term drop levels
                Z[[k]]<-Z[[k]][,-missing,drop=FALSE]
                nrl[i]<-ncol(Z[[k]])
              }else{                            # if univariate G-structure but Ginv term replace with structural zeros
                for(j in 1:length(missing)){
                  Z[[k]][j,missing[j]]<-1e-64
                }
              }
            }else{                              # if multiivariate G-structure replace with structural zeros
              for(j in 1:length(missing)){
                Z[[k]][j,missing[j]]<-1e-64
              }
            }
          }
          if(k>1){
            if(any(levels(data[,rterms[select.terms]])!=levels(data[,rterms[1]]))){
              stop("terms involved in mm/str structures must have identical levels")
            }
            if(rtype=="mm"){
              colnames(Z[[k]])<-colnames(Z[[1]])
              Z[[1]]<-Z[[1]]+Z[[k]]
            }
            if(rtype=="str"){
              Z[[1]]<-cBind(Z[[1]],Z[[k]])
            }
          }
        }
        if(i==1){
          Zsave<-Z[[1]]
        }else{
          Zsave<-cBind(Zsave,Z[[1]])
        }
      }
      vnames<-colnames(X)  

      if(vtype=="idh"){
        Aterm<-rep(Aterm, nfl)
        vtype<-rep(vtype, nfl)
        nfl<-rep(1, nfl)
      }else{
        if(vtype[1]=="idv"){
          nrl<-nfl*nrl[1]
	  nfl<-1
	  vnames<-idv.vnames
        }
      }
      if(length(rterms)==0){
        rterms<-""
      }
      if(rtype=="iid"){
        if(vtype[1]=="us" | grepl("cor", vtype[1])){
          vnames<-paste(expand.grid(vnames, vnames)[,1],expand.grid(vnames, vnames)[,2], sep=":")
        }
        vnames<-paste(vnames, paste(rterms[select.terms], collapse=":"), sep=".")
      }
      if(rtype=="mm"){
        if(vtype[1]=="us"){
          vnames<-paste(expand.grid(vnames, vnames)[,1],expand.grid(vnames, vnames)[,2], sep=":")
        }
        vnames<-paste(vnames, paste(rterms, collapse="+"), sep=".")
      }
      if(rtype=="str"){
        nfl<-nfl*length(rterms)
        if(vtype[1]!="idh"){
          vnames<-paste(expand.grid(rterms, vnames)[,1],expand.grid(rterms, vnames)[,2], sep=":")
          vnames<-paste(expand.grid(vnames, vnames)[,1], expand.grid(vnames, vnames)[,2], sep=".")
        }else{
          vnamestmp<-vnames
          vnames<-c()
          for(i in 1:length(vnamestmp)){
            vnames<-c(vnames, paste(expand.grid(paste(rterms, vnamestmp[i], sep=":"),paste(rterms, vnamestmp[i], sep=":"))[,1], expand.grid(paste(rterms, vnamestmp[i], sep=":"),paste(rterms, vnamestmp[i], sep=":"))[,2], sep="."))
          }
        }
        vtype<-rep("us", length(vtype))
      }
      return(list(Z=Zsave, nfl=nfl, nrl=nrl, Aterm=Aterm, vtype=vtype, vnames=vnames, ordering=NULL, trait.ordering=NULL))
    }

  }else{

#############################################################################
### this is for residual structuers which currently have simple structure ###
#############################################################################

     if(length(rterms)==0){
       rfactor<-rep(as.factor(1), dim(data)[1]) 	
     }else{
       if(any(rterms%in%colnames(data)==FALSE)){stop(paste("object", paste(rterms, collapse=" and "), "not found"))}
       rfactor<-interaction(data[,rterms], drop=(Aterm==0 & rtype=="iid"))  # random effects are dropped if they are not animal 
     }

     if(nfl==0){

       nfl=1

       if(any(duplicated(rfactor))){stop("R-structure miss-specified: each residual must be unique to a data point")}
       if(length(rfactor)!=dim(data)[1]){stop("R-structure miss-specified: each data point must have a residual")}

       nrl<-nlevels(rfactor)
       ordering<-order(as.numeric(rfactor))
       trait.ordering<-data$trait[1]
       vnames<-paste(rterms, collapse=":")

     }else{
      
       if(any(!X@x%in%c(0,1,1e-64) | any(rowSums(X)>=2))){stop(paste("The formula for", vtype, "R-structures must generate a design matrix where each row is all zero's except for one element equal to one."))}

       # With augmented data (for animal models or us structures in the random effects) the ffac term in idh(ffac):rfac may be missing. 
       # Replace with arbitrary levels sampled from the frequency of observed levels. Only required for the first R-structure.

       if(augmiss){
         miss.ffac<-which(data$MCMC_dummy==1 & apply(X, 1, function(x){all(x==1e-64)}))
         if(length(miss.ffac)>0){
           X[miss.ffac,]<-diag(ncol(X))[sample(1:ncol(X), length(miss.ffac), TRUE, colMeans(X[drop=FALSE,-miss.ffac,])),]
         }
       }

       vnames<-paste(colnames(X), paste(rterms, collapse=":"), sep=".")

       trait.ordering<-as.numeric(data$trait[X@i[X@p[1:ncol(X)]+1]+1])
       ordering<-(length(unique(rfactor))+1)*as.numeric(X%*%Matrix(1:nfl,nfl,1))+as.numeric(rfactor)
       ordering[which(as.integer(as.numeric(X%*%Matrix(1:nfl,nfl,1)))==0)]<-NA  # records not associated with this R-structure in block-diagonals

       ordering<-order(ordering, na.last=NA)

       nrl<-as.integer(colSums(X))

       if(vtype=="idh"){
         nfl<-ncol(X)
         Aterm<-rep(Aterm, nfl)
         vtype<-rep(vtype, nfl)
         nfl<-rep(1, nfl)
       }else{
         if(var(as.integer(nrl))!=0){stop("problem in buildZ.R: expecting identical random levels per ffac")}
         nrl<-nrl[1]
         if(vtype=="idv"){
           nrl<-nfl*nrl[1]
   	   nfl<-1
           vnames<-idv.vnames          
         }else{
           vnames<-paste(paste(expand.grid(colnames(X), colnames(X))[,1],expand.grid(colnames(X), colnames(X))[,2], sep=":"),paste(rterms, collapse=":"), sep=".")
         }
       }
     }

     return(list(Z=NULL, nfl=nfl, nrl=nrl, Aterm=Aterm, vtype=vtype, vnames=vnames, ordering=ordering, trait.ordering=trait.ordering))
  }
}
