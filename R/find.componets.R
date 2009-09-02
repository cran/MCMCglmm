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


  fformula<-gsub("^(us|cor|corh|idh)\\(", "", x)

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

  if(fformula!=""){fformula<-c("+", fformula)}

  fformula<-as.formula(paste("~", "-1", paste(fformula, collapse=""), sep=""))

  ndata=data[1:2,]
  if(notrait){
    ndata=data[1:2,]
    ndata$trait<-as.factor(c("A", "B"))
  }
  X<-model.matrix(fformula, data=ndata)
  fformula=names(attr(X, "contrasts"))

  if(is.null(fformula)==FALSE){
    if(fformula=="trait" & notrait){
      fformula=NULL
    }
  }
  if(rterms=="animal"){
     if(is.null(fformula)){
       fformula<-"MCMC_dummy"
     }
  }else{
    if(is.null(fformula) | vtype=="idh"){
      rterms<-NULL
      fformula=NULL
    }
  }
  return(c(rterms,fformula))
}
