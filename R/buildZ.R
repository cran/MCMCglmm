buildZ<-function(x, data, formZ=TRUE){

  vtype="idh"
  if(length(grep("^us\\(", x))>0){
    vtype<-"us"
  }
  if(length(grep("^corh\\(", x))>0){
    vtype<-"corh"
  }
  if(length(grep("^cor\\(", x))>0){
    vtype<-"cor"
  }
  if(length(grep("^idv\\(", x))>0){
    vtype<-"idv"
  }

  fformula<-gsub("^(us|cor|corh|idh|idv)\\(", "", x)

  openB<-gregexpr("\\(", fformula)[[1]]
  closeB<-gregexpr("\\)", fformula)[[1]]

  if(openB!=-1){
    while(openB[1]<closeB[1]){
      openB<-openB[-1]
      closeB<-closeB[-1]
      if(length(closeB)==0 | length(openB)==0){break}
    }
  }

  rterms<-substr(fformula,closeB[1]+2,nchar(fformula))
  rterms<-strsplit(rterms, ":")[[1]]

  if(any(rterms=="animal")){
    Aterm<-1
  }else{
    Aterm<-0
  }

  fformula<-substr(fformula,1,closeB-1)
  if(fformula!=""){fformula<-c("+", fformula)}

  fformula<-as.formula(paste("~", "-1",  paste(fformula, collapse=""), sep=""))

  orig.na.action<-options("na.action")[[1]]
  options("na.action"="na.pass")

  X<-model.matrix(fformula, data)

  X[which(is.na(X))]<-1e-64        # structural zeros needed when us() because may be things like slope terms in leg for animal


  options("na.action"=orig.na.action)

  X<-as(X, "sparseMatrix")

  if(length(rterms)==0){
    rfactor<-rep(as.factor(1), dim(data)[1]) 	
  }else{
    if(any(rterms%in%colnames(data)==FALSE)){stop(paste("object", paste(rterms, collapse=" and "), "not found"))}
    rfactor<-interaction(data[,rterms], drop=(Aterm==0))  # random effects are dropped if they are not animal 
  }

  nfl<-ncol(X)
  nrl<-nlevels(rfactor)
  nd<-length(rfactor)

  if(formZ){
             
    data_pos<-as.numeric(rfactor)

    Ztmp<-Matrix(0,nd,nrl)
    Ztmp[,1][2]<-1              # force it out of vbeing upper triangle!!!!
    Ztmp@p<-as.integer(c(0,cumsum(table(rfactor))))    
    cnt<-0
    for(j in 1:nrl){
      hit<-which(data_pos==j)
      hit<-hit-nd*(ceiling(hit/nd)-1)
      if(length(hit)>0){
        Ztmp@i[cnt+1:length(hit)]<-as.integer(hit-1)
        cnt<-cnt+length(hit)
      }
    }
    Ztmp@x<-rep(1,length(Ztmp@i))

    colnames(Ztmp)<-paste(paste(rterms, collapse=":"),levels(rfactor), sep=".")

    if(nfl==0){

      missing<-which(colSums(Ztmp)==0)  # non-represented random effects (or zero covariates for RR)

      #################################################################################################
      # It woule be more efficient to allow complete sparse columns in Z - but need to change C code! #
      #################################################################################################
 
      if(length(missing)>0){        
        if(Aterm==0){                       # if not associated with ginv drop terms     
          Ztmp<-Ztmp[,-missing,drop=FALSE]
          nrl[i]<-ncol(Ztmp)              
        }else{                              # if associated with ginv set arbitrary row to structural zero   
          for(j in 1:length(missing)){
            Ztmp[j,missing[j]]<-1e-64     
          }
        }
      }

      nfl<-1
      return(list(Z=Ztmp, nfl=nfl, nrl=nrl, Aterm=Aterm, vtype=vtype, vnames=paste(rterms, collapse=":"), ordering=NULL, trait.ordering=NULL))

    }else{
  
       vnames<-paste(colnames(X), paste(rterms, collapse=":"), sep=".")

      if(vtype=="idh"){
        nrl<-rep(nrl, nfl)
      }else{
        vnames<- paste(expand.grid(vnames, vnames)[,1],expand.grid(vnames, vnames)[,2], sep=".")
      }
      for(i in 1:nfl){
        Xtmp<-Matrix(0,nd,nd)
        Xcol<-X[,i,drop=FALSE]
        Xtmp@p<-as.integer(rep(0:c(length(Xcol@i)-1),diff(c(0,Xcol@i+1))))
        Xtmp@p[(length(Xtmp@p)+1):(nd+1)]<-length(Xcol@i)
        Xtmp@i<-Xcol@i
        Xtmp@x<-Xcol@x
        Z<-Xtmp%*%Ztmp
        colnames(Z)<-paste(paste(rterms, collapse=":"), colnames(X)[i], colnames(Z),sep=".") 

        missing<-which(colSums(Z)==0)

        #################################################################################################
        # It woule be more efficient to allow complete sparse columns in Z - but need to change C code! #
        #################################################################################################

        if(length(missing)>0){
          if(vtype=="idh" | vtype=="idv"){
            if(Aterm==0){                     # if univariate G-structure and no Ginv term drop levels
              Z<-Z[,-missing,drop=FALSE]
              nrl[i]<-ncol(Z)
            }else{                            # if univariate G-structure but Ginv term replace with structural zeros
              for(j in 1:length(missing)){
                 Z[j,missing[j]]<-1e-64
              }
            }
          }else{                              # if multiivariate G-structure replace with structural zeros
            for(j in 1:length(missing)){
              Z[j,missing[j]]<-1e-64
            }
          }
        }


        if(i==1){
          Zsave<-Z
        }else{
          Zsave<-cBind(Zsave,Z)
        }
      }
      if(vtype=="idh"){
        Aterm<-rep(Aterm, nfl)
        vtype<-rep(vtype, nfl)
        nfl<-rep(1, nfl)
      }
      if(vtype[1]=="idv"){
        nrl<-nfl*nrl[1]
	nfl<-1
	vnames<-vnames[1]  
      }
      return(list(Z=Zsave, nfl=nfl, nrl=nrl, Aterm=Aterm, vtype=vtype, vnames=vnames, ordering=NULL, trait.ordering=NULL))
    }

  }else{

    if(nfl==0){
      nfl=1
      ordering<-order(as.numeric(rfactor))
      trait.ordering<-data$trait[1]
      vnames<-paste(rterms, collapse=":")
    }else{
      vnames<-paste(colnames(X), paste(rterms, collapse=":"), sep=".")
      trait.ordering<-as.numeric(data$trait[X@i[X@p[1:ncol(X)]+1]+1])
      ordering<-order(100000000*as.numeric(X%*%Matrix(1:nfl,nfl,1))+as.numeric(rfactor))
      if(vtype=="idh"){
        nfl<-ncol(X)
        Aterm<-rep(Aterm, nfl)
        vtype<-rep(vtype, nfl)
        nfl<-rep(1, nfl)
        nrl<-colSums(X)
      }else{
        if(vtype=="idv"){
          nrl<-nfl*nrl[1]
	  nfl<-1
	  vnames<-vnames[1]  
        }else{
          vnames<-paste(expand.grid(vnames, vnames)[,1],expand.grid(vnames, vnames)[,2], sep=".")
        }
      }
    }
    return(list(Z=NULL, nfl=nfl, nrl=nrl, Aterm=Aterm, vtype=vtype, vnames=vnames, ordering=ordering, trait.ordering=trait.ordering))
  }
}
