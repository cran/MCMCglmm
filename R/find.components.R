find.components<-function(x, data){

  notrait=is.null(data$trait)

  vtype="us"
  if(length(grep("^idh\\(", x))>0){
    vtype<-"idh"
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

  fformula<-substr(fformula,1,closeB-1)
  fformula<-gsub(" ", "", fformula)

  interact<-strsplit(fformula, "\\+")[[1]]
  interact<-interact[grep(":|\\*",interact)]
  interact<-strsplit(interact, ":|\\*")

  if(fformula!=""){fformula<-c("+", fformula)}

  fformula<-as.formula(paste("~", "-1", paste(fformula, collapse=""), sep=""))

  ndata=data[1:min(dim(data)[1], 10),]

  if(notrait){
    ndata$trait<-gl(2, 1,min(dim(data)[1], 10))
  }
  X<-model.matrix(fformula, data=ndata)

  fformula=names(attr(X, "contrasts"))

  # IDEALLY - WHEN ANIMAL SPECIFIED at.level/set.level TERMS SHOULD BE RETAINED IN fformula

  if(is.null(fformula)==FALSE){
    if(any(fformula=="trait") & notrait){
      fformula<-fformula[-which(fformula=="trait")]
      if(length(fformula)==0){
        fformula<-NULL
      }
    }
  }
  if(any(rterms=="animal")){
     if(is.null(fformula)){
       fformula<-"MCMC_dummy"
     }
  }else{
    if(is.null(fformula) | vtype=="idh"| vtype=="idv"){
      rterms<-NULL
      fformula=NULL
    }
  }
  # if any terms are interacted they have to be added
  if(length(interact)>0){
    for(i in 1:length(interact)){
      if(sum(interact[[i]]%in%fformula)>1){
        interact[[i]]<-interact[[i]][which(interact[[i]]%in%fformula)]
        fformula<-fformula[-which(fformula%in%interact[[i]])]
      }else{
        interact<-interact[-i]
      }
    }  
  }

  return(list(rterms=rterms,fformula=c(as.list(fformula), interact)))
}
