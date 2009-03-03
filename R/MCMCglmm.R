"MCMCglmm"<-function(fixed, random=NULL, rcov=~units, family="gaussian", mev=NULL, data=NULL, start=NULL, prior=NULL, tune=NULL, pedigree=NULL, nodes="ALL",scale=TRUE,  nitt=13000, thin=10, burnin=3000, pr=FALSE, pl=FALSE, verbose=TRUE, DIC=TRUE){
   
    options("na.action"="na.pass")
    row.names(data)=NULL	

    if(is.null(data)){stop("data argument is NULL")}
    if(is.null(fixed)){stop("fixed is NULL")}

    if(any(names(data)%in%c("units", "trait", "MCMC_y", "MCMC_y.additional","MCMC_liab","MCMC_meta", "MCMC_mev", "MCMC_family.names"))){
      stop(paste(names(data)[which(names(data)%in%c("units", "trait", "MCMC_y", "MCMC_y.additional","MCMC_liab","MCMC_meta", "MCMC_mev", "MCMC_family.names"))], " is a reserved variable please rename it"))
    }

    family.types<-c("gaussian", "poisson", "multinomial", "notyet_weibull", "exponential", "cengaussian", "cenpoisson", "notyet_cenweibull", "cenexponential",  "notyet_zigaussian", "zipoisson", "notyet_ziweibull", "notyet_ziexponential")

    if(((is.null(start$G) & is.null(random)==FALSE) & is.null(start$R)==FALSE) | (is.null(start$R) & is.null(start$G)==FALSE)){stop("need both or neither starting R and G structures")}
    if(((is.null(prior$G) & is.null(random)==FALSE) & is.null(prior$R)==FALSE) | (is.null(prior$R) & is.null(prior$G)==FALSE)){stop("either both or neither R and G structures need a prior")}
 
    if(is.null(start$QUASI)){
      QUASI=TRUE
    }else{
      QUASI<-start$QUASI
      if(is.logical(QUASI)==FALSE){stop("start$QUASI should be TRUE or FALSE")}
    }

    original.fixed<-fixed                                                                            # original model specification
    original.random<-random 
    original.rcov<-rcov

    response.names<-names(get_all_vars(as.formula(paste(as.character(fixed)[2], "~1")), data))       # response variable names 

    if(length(grep("\\+", rcov))>0){stop("R-structure must currently be defined by a single component")}

######################################################################################
# for phyloegnetic/pedigree analyses form A and augment with missing nodes if needed #
###################################################################################### 	
	
	if(is.null(pedigree)==FALSE){
	  if(is.null(data$animal) |length(grep("animal", random))==0){stop("pedigree/phylogeny has been passed, but no 'animal' variable")}
	  Ai<-inverseA(pedigree, nodes=nodes, scale=scale)
          if(any(data$animal%in%Ai$node.names==FALSE)){stop("individuals in data not appearing in pedigree")}
	  data$animal<-factor(data$animal, levels=Ai$node.names)                           # add factor levels to animal in data for additional nodes
	  if(length(unique(data$animal))!=length(Ai$node.names)){
	    warning(paste("some pedigree/phylogeny nodes not appearing in data: adding",length(Ai$node.names)-length(unique(data$animal)), "missing records"))
            missing.combinations<-setdiff(levels(data$animal), data$animal)
            if(length(missing.combinations)>0){
	      data[dim(data)[1]+1:length(missing.combinations),]<-as.data.frame(NA)
	      data$animal[dim(data)[1]-(length(missing.combinations)-1):0]<-missing.combinations
	    }	
	  }
	}
        
#######################################################################################################
# for interactions involving pedigrees/phylogenies or us structures augment with missing combinations #
#######################################################################################################

	rterms<-strsplit(as.character(random)[2], " *\\+ *")[[1]]
	
	for(i in 1:length(rterms)){
	  components=NULL
	  rterms.split<-lapply(rterms, strsplit,"\\(|\\):|:")
	  if(length(rterms.split[[i]][[1]])==2 & any(rterms.split[[i]][[1]]=="animal")){
	    components<-rterms.split[[i]][[1]]
	  }
	  if(length(rterms.split[[i]][[1]])==3 & rterms.split[[i]][[1]][1]=="us" & rterms.split[[i]][[1]][2]!="trait"){
	    components<-rterms.split[[i]][[1]][2:3]
	  }
	  if(is.null(components)==FALSE){
	    Ztmp<-model.matrix(formula(paste(response.names[1], "~", components[2], ":", components[1], "-1", sep="")), data)
		  
            missing.combinations<-which(colSums(Ztmp, na.rm=TRUE)==0)
		  
            if(length(missing.combinations)>0){
			  
	      warning(paste("some combinations in", rterms[i], "do not exist and", length(missing.combinations), "missing records have been generated"))
			  
              missing.fixed<-ceiling(missing.combinations/nlevels(data[,components[2]]))
	      missing.random<-missing.combinations-nlevels(data[,components[2]])*(missing.fixed-1)
	      missing.fixed<-levels(data[,components[1]])[missing.fixed]
	      missing.random<-levels(data[,components[2]])[missing.random]
	      missing.combinations<-cbind(missing.fixed, missing.random)
			  
	      missing.comb1<-which(is.na(data[,components[1]]) & is.na(data[,components[2]])==FALSE)
	      missing.comb2<-which(is.na(data[,components[1]])==FALSE & is.na(data[,components[2]]))
	      missing.comb12<-which(is.na(data[,components[1]]) & is.na(data[,components[2]]))
			  
	      matching<-match(unique(missing.combinations[,2]), data[,components[2]][missing.comb1])			
	      while(any(is.na(matching)==FALSE)){
	        data[,components[1]][missing.comb1[na.omit(matching)]]<-missing.combinations[,1][match(data[,components[2]][missing.comb1[na.omit(matching)]],missing.combinations[,2])]
       	        missing.comb1<-missing.comb1[-match(data[,components[2]][missing.comb1[na.omit(matching)]],missing.combinations[,2])]	
		missing.combinations<-missing.combinations[-match(data[,components[2]][missing.comb1[na.omit(matching)]],missing.combinations[,2]),, drop=FALSE]
		matching<-match(unique(missing.combinations[,2]), data[,components[2]][missing.comb1]) 
	      }
	      matching<-match(unique(missing.combinations[,1]), data[,components[1]][missing.comb2])
	      while(any(is.na(matching)==FALSE)){
		data[,components[2]][missing.comb2[na.omit(matching)]]<-missing.combinations[,2][match(data[,components[1]][missing.comb2[na.omit(matching)]],missing.combinations[,1])]
		missing.comb2<-missing.comb2[-match(data[,components[1]][missing.comb2[na.omit(matching)]],missing.combinations[,1])]	
		missing.combinations<-missing.combinations[-match(data[,components[1]][missing.comb2[na.omit(matching)]],missing.combinations[,1]),, drop=FALSE]
		matching<-match(unique(missing.combinations[,1]), data[,components[1]][missing.comb2]) 
	      }
	      if(length(missing.comb12)>0 & length(missing.combinations)>0){
		data[missing.comb12,]<-missing.combinations[1:min(length(missing.comb12), length(missing.combinations[,1])),]
		missing.combinations<-missing.combinations[-1:min(length(missing.comb12), length(missing.combinations[,1])),,drop=FALSE,]
	      }			
	      if(length(missing.combinations)>0){                                                      # add dummy records if still needed
		data[dim(data)[1]+1:dim(missing.combinations)[1],]<-NA                  
		data[,components][dim(data)[1]-(dim(missing.combinations)[1]-1):0,]<-missing.combinations
	      }
	    }		  
	  }
	}

##############################################################################################
# if R-structure is of form idh(!=trait) assign new records with arbitrary levels of !-trait #
##############################################################################################
	
	 rterms.split<-lapply(strsplit(as.character(rcov)[2], " *\\+ *")[[1]], strsplit,"\\(|\\):|:")[[1]][[1]]

  	 components=NULL
	 if(length(rterms.split)==2){
	    components<-rterms.split
	 }
	 if(length(rterms.split)==3){
	    components<-rterms.split[2:3]
	 }
	 if(is.null(components)==FALSE){
	   if(any(components!="trait" & components!="units")){
             components<-components[which(components!="trait" & components!="units")]
	     data[,components][which(is.na(data[,components]))]<-sample(levels(data[,components]), sum(is.na(data[,components])), replace=TRUE)
	   }
	 }
      

    MVasUV=FALSE  # Multivaraite as Univariate
    if(is.null(family)){
       if(is.null(data$family)){
         stop("no family specified")
       }else{
         if(length(unique(data$family))==1){
           family.names<-as.character(data$family[1])
         }else{     
           family.names<-as.character(data$family)   
           MVasUV=TRUE
         }
         if(length(grep("cen|multinomial|zi", family.names))>0){ 
           stop("For setting up multi-trait models as univariate the responses cannot come from distributions that require more than one data column (i.e. censored, multinomial, zero-inflated): set it up as multivariate using cbind(...)")
         }
       }
    }else{
       if(is.null(data$family)==FALSE){
          stop("family column exists in data and in family argument to MCMCglmm: specify family=NULL")
       }else{
          family.names<-family 
       }                                                                            # response distribution
    }

    # zero-infalted and multinomial>2 need to be preserved in the same R-structure even if idh 
	
    if(length(grep("zi|multinomial[3:19]", family.names))>0){ 
      rcov=as.formula(paste(gsub("idh\\(", "us\\(", rcov), collapse=""))
      diagR=TRUE	  
    }else{
      diagR=FALSE
    }
	
    nS<-dim(data)[1]                             # number of subjects
    y.additional<-matrix(NA, nS,0)               # matrix equal in dimension to y holding the additional parameters of the distribution (n, upper interval etc.)

    nt<-1                                        # number of traits (to be iterated because y may change dimension with multinomial/categorical/censored distributions)
    mfac<-c()

    if(MVasUV){
      if(any(family.names=="categorical")){
        y.additional<-matrix(NA, nS,1)
        y.additional[which(family.names=="categorical")]<-1
        family.names[which(family.names=="categorical")]<-"multinomial"
        family.names[which(is.na(family.names))]<-"gaussian"
        mfac<-1
      }
    }else{
      for(i in 1:length(family)){

        dist.preffix<-substr(family[i],1,2)                  

        if(any(dist.preffix%in%c("ce", "mu", "ca", "tr", "zi"))){

######################
# categorical traits #
######################

          if(dist.preffix=="ca"){
            cont<-as.matrix(model.matrix(~ as.factor(data[[response.names[nt]]]))[,-1])                            # form new J-1 variable   
            nJ<-dim(cont)[2]                                                                                       # number of J-1 categories 
            mfac<-c(mfac, rep(nJ,nJ))             
            new.names<-paste(response.names[nt], ".", levels(as.factor(data[[response.names[nt]]]))[-1], sep="")   # give new variables names
            colnames(cont)<-new.names                        
            data<-data[,-which(names(data)==response.names[nt])]                                                   # remove original variable
            data<-cbind(data, cont)                                                                                # add new variables to data.frame
            ones<-rep(1, length(response.names))
            ones[which(response.names==response.names[nt])]<-nJ
            response.names<-rep(response.names, ones)
            family.names<-rep(family.names, ones)
            family.names[which(response.names==response.names[nt])]<-"multinomial"
            response.names[which(response.names==response.names[nt])]<-new.names
            nt<-nt+nJ
            y.additional<-cbind(y.additional, matrix(1,nS,nJ))
          }

######################
# multinomial traits #
######################

         if(dist.preffix=="mu"){
           nJ<-as.numeric(substr(family[i],12,nchar(family[i])))-1                                                                            # number of J-1 categories
		   if(nJ<1){stop("Multinomial must have at least 2 categories")}	 
           mfac<-c(mfac, rep(nJ,nJ))  
		   if(all(data[,match(response.names[0:nJ+nt], names(data))]%%1==0, na.rm=T)==FALSE | all(data[,match(response.names[0:nJ+nt], names(data))]>=0, na.rm=T)==FALSE){stop("multinomial data must be positive integers")}
           y.additional<-cbind(y.additional, matrix(rowSums(data[,match(response.names[0:nJ+nt], names(data))]), nS,nJ))             # get n of the multinomial
           data<-data[,-which(names(data)==response.names[nt+nJ]),drop=FALSE]                                                                      # remove first category
           response.names<-response.names[-(nt+nJ)]
           family.names[0:nJ+nt]<-"multinomial"
           family.names<-family.names[-(nt+nJ)]
           nt<-nt+nJ
         }

###################
# censored traits #
###################

         if(dist.preffix=="ce"){         
		   if(any(data[,which(names(data)==response.names[nt+1])]<data[,which(names(data)==response.names[nt])], na.rm=T)){stop("for censored traits left censoring point must be less than right censoring point")}
           y.additional<-cbind(y.additional, data[,which(names(data)==response.names[nt+1])])                               # get upper interval
           if(family.names[nt]=="cenpoisson"){
			  if(all(data[,response.names[0:1+nt]]%%1==0, na.rm=T)==FALSE | all(data[,response.names[0:1+nt]]>=0, na.rm=T)==FALSE){stop("Poisson data must be positive integers")}
             data[[response.names[nt]]][which(data[[response.names[nt]]]!=data[[response.names[nt+1]]])]<-data[[response.names[nt]]][which(data[[response.names[nt]]]!=data[[response.names[nt+1]]])]-1
           }
		   if(family.names[nt]=="cenexponential"){
			   if(any(data[,response.names[0:1+nt]]<0, na.rm=T)){stop("Exponential data must be positive")}  
		   }
		   data<-data[,-which(names(data)==response.names[nt+1]),drop=FALSE]                                                # remove upper interval from the response
           response.names<-response.names[-(nt+1)]
           nt<-nt+1
         }

###################
# truncated traits #
###################
		  
	 if(dist.preffix=="tr"){   
	   y.additional<-cbind(y.additional, data[,which(names(data)==response.names[nt+1])])                     # get upper interval
	   data<-data[,-which(names(data)==response.names[nt+1]),drop=FALSE]                                      # remove upper interval from the response
	   response.names<-response.names[-(nt+1)]
	   nt<-nt+1
	 }

########################
# zero-infalted traits #
########################
		  
	 if(dist.preffix=="zi"){                                                                         
	   y.additional<-cbind(y.additional, rep(1,nS), rep(0,nS))
	   if(all(data[,response.names[nt]]%%1==0, na.rm=T)==FALSE | all(data[,response.names[nt]]>=0, na.rm=T)==FALSE){stop("Poisson data must be integers")}
	   cont<-as.matrix(as.numeric(data[,which(names(data)==response.names[nt])]==0))
	   colnames(cont)<-paste("zi", response.names[nt], sep=".")
	   data<-cbind(data, cont)
	   mfac<-c(mfac, rep(1,1))
	   ones<-rep(1, length(response.names))
	   ones[which(response.names==response.names[nt])]<-2
	   response.names<-rep(response.names, ones)
	   family.names<-rep(family.names, ones)
	   family.names[which(response.names==response.names[nt])]<-family.names[nt]
	   response.names[which(response.names==response.names[nt])]<-c(response.names[nt], paste("zi", response.names[nt], sep="."))
	   nt<-nt+2			  
	 }
		  
#######################################
# gaussian/poisson/exponential traits #
#######################################

        }else{
			if(family.names[nt]=="poisson"){
				if(all(data[,response.names[nt]]%%1==0, na.rm=T)==FALSE | all(data[,response.names[nt]]>=0, na.rm=T)==FALSE){stop("Poisson data must be positive integers")}
			}
			if(family.names[nt]=="exponential"){
				if(any(data[,response.names[nt]]<0, na.rm=T)){stop("Exponential data must be positive")}
			}	
          y.additional<-cbind(y.additional,matrix(NA,nS,1))     
          nt<-nt+1
        }
      }	
      nt<-nt-1
    }

    if(sum((family.names%in%family.types)==FALSE)!=0){stop(paste(unique(family[which((family.names%in%family.types)==FALSE)]), "not a supported distribution"))}

###**************************************########################

    data<-data.frame(data, units=as.factor(1:nS))
    data<-reshape(data, varying=response.names, v.names="MCMC_y", direction="long", timevar="trait")       # reshape the data into long format 
    data$trait<-factor(response.names[data$trait], response.names)
    data$MCMC_y.additional<-c(y.additional) 

    if(MVasUV){
      data$MCMC_family.names<-family.names
    }else{
      if(length(response.names)!=length(family.names)){stop("family must have the same length as the number of responses")}
      data$MCMC_family.names<-rep(family.names, each=nS)       
    }
###**************************************########################

    fixed<-update(fixed,y~.,) 

######################################################
# for (random) meta-analysis add weights/model terms #
###################################################### 	

    if(is.null(mev)==FALSE){
      if(any(dim(mev)!=dim(y.additional))){stop("mev has to be the same dimension as y")}
      data$MCMC_mev<-sqrt(mev)*sqrt(2/3)  # needs to be multiplied by sqrt(2/3) because legendre polynomial slope is sqrt(3/2)
      data$MCMC_meta<-factor(1:dim(data)[1], levels=1:dim(data)[1])
      if(is.null(random)){
        random = ~leg(MCMC_mev, -1):MCMC_meta
        if(is.null(prior$R)==FALSE){
          prior$G<-list(G1=list(V=as.matrix(1), n=1, fix=1))
        }
      }else{
        random<-update(random,~.+leg(MCMC_mev, -1):MCMC_meta)
        if(is.null(start$G)==FALSE){
          start$G[[length(start$G)+1]]<-as.matrix(1)
        }
        if(is.null(prior$G)==FALSE){
          prior$G[[length(prior$G)+1]]<-list(V=as.matrix(1), n=1, fix=1)
        }
      } 
    }

    rmodel.terms<-strsplit(as.character(random)[2], " *\\+ *")[[1]]
    rmodel.terms<-c(na.omit(rmodel.terms), strsplit(as.character(rcov)[2], " *\\+ *")[[1]]) 
 
    nfl<-c()                                                       # number of fixed levels the random term is structured by
    nrl<-c()                                                       # number of random levels
    variance.names<-c()
    GR<-list()
    GRprior<-list()
    Aterm<-c()


    if(is.null(start$R)){
      NOstartG=TRUE
    }else{
      NOstartG=FALSE
      start$G[[length(start$G)+1]]<-start$R
    }
    if(is.null(prior$R)){
      NOpriorG=TRUE
    }else{
      NOpriorG=FALSE
      prior$G[[length(prior$G)+1]]<-prior$R
    }

    nr<-1

    for(r in 1:length(rmodel.terms)){
 
      if(r==length(rmodel.terms)){nG<-nr-1}  # number of G structures

#########################
# normal random effects #
#########################

      if(length(grep("\\(", rmodel.terms[r]))==0){
        
        nfl[nr]<-1

        if(length(grep("\\:", rmodel.terms[r]))>0){
          components<-strsplit(rmodel.terms[r], "\\:")[[1]]
          if(is.factor(data[,components[1]])==FALSE){stop(paste(components[1], " is not a factor"))}
          if(is.factor(data[,components[2]])==FALSE){stop(paste(components[2], " is not a factor"))}
	  if(any(levels(data[,components[1]])%in%unique(data[,components[1]])==FALSE)){stop(paste("some levels not observed in ", components[1]))}
	  if(any(levels(data[,components[2]])%in%unique(data[,components[2]])==FALSE)){stop(paste("some levels not observed in ", components[2]))}
          rmodel.terms[r]<-paste(components, collapse=".")
          data[[rmodel.terms[r]]]<-as.factor(paste(data[,components[1]], data[,components[2]], sep=""))
        }
		  
        if(is.factor(data[,rmodel.terms[r]])==FALSE){stop(paste(rmodel.terms[r], " is not a factor"))}
	if(any(levels(data[,rmodel.terms[r]])%in%unique(data[,rmodel.terms[r]])==FALSE)){stop(paste("some levels not observed in ", rmodel.terms[r]))}		  
		  
        nrl[nr]<-nlevels(data[,rmodel.terms[r]])
        variance.names<-c(variance.names, rmodel.terms[r])
		  
        if(NOpriorG==TRUE){
            GRprior[[nr]]<-list(V=matrix(1), n=0)
        }else{

          if(length(prior$G)<r){stop("priorG/priorR have the wrong number of structures")}
          if(is.null(prior$G[[r]]$V)){stop(paste("V does not exist for G/R structure", r))}
	  if(is.matrix(prior$G[[r]]$V)==FALSE){prior$G[[r]]$V<-as.matrix(prior$G[[r]]$V)}	
          if(dim(prior$G[[r]]$V)[1]!=1  | dim(prior$G[[r]]$V)[2]!=1){stop(paste("V is the wrong dimension for G/R structure", r))}
          if(is.null(prior$G[[r]]$n)){stop(paste("n does not exist for G/R structure", r))}
          GRprior[[nr]]<-prior$G[[r]]
        }

        if(NOstartG==TRUE){

          GR[[nr]]<-GRprior[[nr]]$V

        }else{

          if(length(start$G)<r){stop("starting G/R has the wrong number of structures")}
	  if(is.matrix(start$G[[r]]$V)==FALSE){start$G[[r]]$V<-as.matrix(start$G[[r]]$V)}	
          if(dim(start$G[[r]])[1]!=1  | dim(start$G[[r]])[2]!=1){stop(paste("starting G/R structure", r, "is the wrong dimension"))}
          if(is.positive.definite(start$G[[r]])==FALSE){stop(paste("starting G/R structure", r, " is not positive definite"))}

          GR[[nr]]<-start$G[[r]]
        }


        if(length(grep("animal", rmodel.terms[r]))>0){
          Aterm[nr]<-1
        }else{
          Aterm[nr]<-0
        }

        nr<-nr+1
      }

#######################################################################
# checking all combination of factors exist for us and cor structures #
#######################################################################

      if(length(c(grep("us\\(", rmodel.terms[r]), grep("cor\\(", rmodel.terms[r])))>0){
        components<-strsplit(rmodel.terms[r], "\\(|\\)\\:")[[1]][2:3]                               # find component terms
         if(is.factor(data[,components[1]])==FALSE){stop(paste(components[1], " is not a factor"))}
         if(is.factor(data[,components[2]])==FALSE){stop(paste(components[2], " is not a factor"))}
 	 if(any(levels(data[,components[1]])%in%unique(data[,components[1]])==FALSE)){stop(paste("some levels not observed in ", components[1]))}
	 if(any(levels(data[,components[2]])%in%unique(data[,components[2]])==FALSE)){stop(paste("some levels not observed in ", components[2]))}

		  
         nfl[nr]<-nlevels(data[,components[1]]) 
         nrl[nr]<-nlevels(data[,components[2]]) 
         variance.names<-c(variance.names, paste(components[2],components[1], expand.grid(levels(data[,components[1]]), levels(data[,components[1]]))[,1], expand.grid(levels(data[,components[1]]), levels(data[,components[1]]))[,2], sep="."))
	rmodel.terms[r]<-paste(rev(components), collapse=":")


        if(NOpriorG==TRUE){
          GRprior[[nr]]<-list(V=diag(nfl[nr]), n=0)
        }else{

          if(length(prior$G)<r){stop("priorG/priorR have the wrong number of structures")}
          if(is.null(prior$G[[r]]$V)){stop(paste("V does not exist for G/R structure", r))}
	  if(is.matrix(prior$G[[r]]$V)==FALSE){prior$G[[r]]$V<-as.matrix(prior$G[[r]]$V)}	
          if(dim(prior$G[[r]]$V)[1]!=nfl[nr] | dim(prior$G[[r]]$V)[2]!=nfl[nr]){stop(paste("V is the wrong dimension for G/R structure", r))}
          if(is.null(prior$G[[r]]$n)){stop(paste("n does not exist for G/R structure", r))}

          GRprior[[nr]]<-prior$G[[r]]
        }

        if(NOstartG==TRUE){

          GR[[nr]]<-GRprior[[nr]]$V

        }else{

          if(length(start$G)<r){stop("starting G has the wrong number of structures")}
	  if(is.matrix(start$G[[r]]$V)==FALSE){start$G[[r]]$V<-as.matrix(start$G[[r]]$V)}	
          if(dim(start$G[[r]])[1]!=nfl[nr]  | dim(start$G[[r]])[2]!=nfl[nr]){stop(paste("starting G structure", r, "is the wrong dimension"))}
          if(is.positive.definite(start$G[[r]])==FALSE){stop(paste("starting G/R structure", r, " is not positive definite"))}

          GR[[nr]]<-start$G[[r]]
        }

        if(r==length(rmodel.terms)){                                                                         # rearrange data to match R-structure
         data<-data[order(100000000*as.numeric(data[,components[1]])+as.numeric(data[,components[2]])),]
        }

        if(length(grep("animal", rmodel.terms[r]))>0){
          Aterm[nr]<-1
        }else{
          Aterm[nr]<-0
        }

        nr<-nr+1
      }

##########################
# heterogenous variances #
##########################

      if(length(c(grep("id\\(", rmodel.terms[r]), grep("idh\\(", rmodel.terms[r])))>0){

	if(length(grep("id\\(", rmodel.terms[r]))>0){fixed.variance=TRUE}else{fixed.variance=FALSE}
        components<-strsplit(rmodel.terms[r], "\\(|\\)\\:")[[1]][2:3]                               # find component terms
        if(is.factor(data[,components[1]])==FALSE){stop(paste(components[1], " is not a factor"))}
        if(is.factor(data[,components[2]])==FALSE){stop(paste(components[2], " is not a factor"))}
        if(any(levels(data[,components[1]])%in%unique(data[,components[1]])==FALSE)){stop(paste("some levels not observed in ", components[1]))}
	if(any(levels(data[,components[2]])%in%unique(data[,components[2]])==FALSE)){stop(paste("some levels not observed in ", components[2]))}

        nfl[nr+1:nlevels(data[,components[1]])-1]<-1
        nrl[nr+1:nlevels(data[,components[1]])-1]<-colSums(table(data[,components[2]],  data[,components[1]])>0)
        rmodel.terms[r]<-paste(rev(components), collapse=":")
        variance.names<-c(variance.names, paste(components[2],components[1], levels(data[,components[1]]), sep="."))
        if(NOpriorG==FALSE){
          if(length(prior$G)<r){stop("priorG/priorR have the wrong number of structures")}
          if(is.null(prior$G[[r]]$V)){stop(paste("V does not exist for G/R structure", r))}
	  if(is.matrix(prior$G[[r]]$V)==FALSE){prior$G[[r]]$V<-as.matrix(prior$G[[r]]$V)}	
          if(dim(prior$G[[r]]$V)[1]!=nlevels(data[,components[1]]) | dim(prior$G[[r]]$V)[2]!=nlevels(data[,components[1]])){stop(paste("V is the wrong dimension for G/R structure", r))}
          if(is.null(prior$G[[r]]$n)){stop(paste("n does not exist for G/R structure", r))}
        }
        if(NOstartG==FALSE){
          if(length(start$G)<r){stop("starting G/R has the wrong number of structures")}
	  if(is.matrix(start$G[[r]]$V)==FALSE){start$G[[r]]$V<-as.matrix(start$G[[r]]$V)}	
          if(dim(start$G[[r]])[1]!=nlevels(data[,components[1]])  | dim(start$G[[r]])[2]!=nlevels(data[,components[1]])){stop(paste("starting G/R structure", r, "is the wrong dimension"))}
          if(is.positive.definite(start$G[[r]])==FALSE){stop(paste("starting G/R structure", r, " is not positive definite"))}
        }

        for(l in 1:nlevels(data[,components[1]])){

          if(NOpriorG==TRUE){
            GRprior[[nr]]<-list(V=matrix(1), n=0)
          }else{
            GRprior[[nr]]<-list(V=as.matrix(diag(prior$G[[r]]$V)[l]), n=prior$G[[r]]$n-(nlevels(data[,components[1]])-1))
            if(is.null(prior$G[[r]]$fix)==FALSE){
              if(l>=prior$G[[r]]$fix){
                GRprior[[nr]]$fix<-1
              }
            }
          }
          if(NOstartG==TRUE){
            GR[[nr]]<-GRprior[[nr]]$V
          }else{
            GR[[nr]]<-as.matrix(diag(start$G[[r]])[l])
          }
 
          if(length(grep("animal", rmodel.terms[r]))>0){
            Aterm[nr]<-1
          }else{
            Aterm[nr]<-0
          }
          nr<-nr+1
        }
        if(r==length(rmodel.terms)){         
          data<-data[order(100000000*as.numeric(data[,components[1]])+as.numeric(data[,components[2]])),]   # rearrange data to match R-structure
        }
      }

#########################
# forming pol structure #
#########################

      if(length(grep("leg\\(", rmodel.terms[r]))>0){
        rr.covariate<-strsplit(rmodel.terms[r], "\\(|\\)|\\,|\\:")[[1]][c(2:3,5)]            # find component terms and order of polynomial
        legvar.name<-paste("leg", rr.covariate[1], abs(as.numeric(rr.covariate[2])), sep=".")
        if(legvar.name%in%names(data)==FALSE){
          lp<-legendre.polynomials(abs(as.numeric(rr.covariate[2])), TRUE)	
          if(as.numeric(rr.covariate[2])<0){
            lp<-lp[-1]
          }			
          legvar<-sapply(lp,function(lp){as.function(lp)(data[,rr.covariate[1]])})
          data$legvar<-legvar
          names(data)[which(names(data)=="legvar")]<-legvar.name
        }
        if(is.factor(data[,rr.covariate[3]])==FALSE){stop(paste(rr.covariate[3], " is not a factor"))}
	if(any(levels(data[,rr.covariate[3]])%in%unique(data[,rr.covariate[3]])==FALSE)){stop(paste("some levels not observed in ",rr.covariate[3]))}

        rmodel.terms[r]<-paste(rr.covariate[3], ":", legvar.name, sep="")
        nfl[nr]<-abs(as.numeric(rr.covariate[2]))+(as.numeric(rr.covariate[2])>0)
        nrl[nr]<-nlevels(data[,rr.covariate[3]])
        pol.order<-(0+(as.numeric(rr.covariate[2])<0)):abs(as.numeric(rr.covariate[2]))
        variance.names<-c(variance.names, paste(rr.covariate[3],rr.covariate[1], expand.grid(pol.order, pol.order)[,1], expand.grid(pol.order, pol.order)[,2], sep="."))

        if(NOpriorG==TRUE){
          if(rr.covariate[1]=="MCMC_mev"){
             GRprior[[nr]]<-list(V=diag(nfl[nr]), n=0, fix=1)
          }else{
             GRprior[[nr]]<-list(V=diag(nfl[nr]), n=0)
          }
        }else{
          if(length(prior$G)<r){stop("priorG/priorR have the wrong number of structures")}
          if(is.null(prior$G[[r]]$V)){stop(paste("V does not exist for G/R structure", r))}
          if(dim(prior$G[[r]]$V)[1]!=nfl[nr] | dim(prior$G[[r]]$V)[2]!=nfl[nr]){stop(paste("V is the wrong dimension for G/R structure", r))}
          if(is.null(prior$G[[r]]$n)){stop(paste("n does not exist for G/R structure", r))}
          GRprior[[nr]]<-prior$G[[r]]
        }
        if(NOstartG==TRUE){
          GR[[nr]]<-GRprior[[nr]]$V
        }else{
          if(length(start$G)<r){stop("starting G/R has the wrong number of structures")}
          if(dim(start$G[[r]])[1]!=nfl[nr]  | dim(start$G[[r]])[2]!=nfl[nr]){stop(paste("starting G/R structure", r, "is the wrong dimension"))}
          GR[[nr]]<-start$G[[r]]
        }

        if(length(grep("animal", rmodel.terms[r]))>0){
          Aterm[nr]<-1
        }else{
          Aterm[nr]<-0
        }
        nr<-nr+1
      }

###########################
# forming group structure #
###########################

      if(length(grep("group\\(", rmodel.terms[r]))>0){

       components<-strsplit(rmodel.terms[r],"group\\( *| *, *| *\\)")[[1]][2:3]

       if(is.factor(data[,components[1]])==FALSE){stop(paste(components[1], " is not a factor"))}
       if(is.factor(data[,components[2]])==FALSE){stop(paste(components[2], " is not a factor"))}
       if(any(levels(data[,components[1]])%in%unique(data[,components[1]])==FALSE)){stop(paste("some levels not observed in ", components[1]))}
       if(any(levels(data[,components[2]])%in%unique(data[,components[2]])==FALSE)){stop(paste("some levels not observed in ", components[2]))}
       rmodel.terms[r]<-paste(components, collapse=".group_MCMC.")
 		  
       nfl[nr]<-1
       nrl[nr]<-nlevels(data[,components[2]])

       variance.names<-c(variance.names, rmodel.terms[r])
		  
       if(NOpriorG==TRUE){
         GRprior[[nr]]<-list(V=matrix(1), n=0)
       }else{
         if(length(prior$G)<r){stop("priorG/priorR have the wrong number of structures")}
         if(is.null(prior$G[[r]]$V)){stop(paste("V does not exist for G/R structure", r))}
	 if(is.matrix(prior$G[[r]]$V)==FALSE){prior$G[[r]]$V<-as.matrix(prior$G[[r]]$V)}	
         if(dim(prior$G[[r]]$V)[1]!=1  | dim(prior$G[[r]]$V)[2]!=1){stop(paste("V is the wrong dimension for G/R structure", r))}
         if(is.null(prior$G[[r]]$n)){stop(paste("n does not exist for G/R structure", r))}
         GRprior[[nr]]<-prior$G[[r]]
       }

       if(NOstartG==TRUE){

          GR[[nr]]<-GRprior[[nr]]$V

        }else{

          if(length(start$G)<r){stop("starting G/R has the wrong number of structures")}
	  if(is.matrix(start$G[[r]]$V)==FALSE){start$G[[r]]$V<-as.matrix(start$G[[r]]$V)}	
          if(dim(start$G[[r]])[1]!=1  | dim(start$G[[r]])[2]!=1){stop(paste("starting G/R structure", r, "is the wrong dimension"))}
          if(is.positive.definite(start$G[[r]])==FALSE){stop(paste("starting G/R structure", r, " is not positive definite"))}

          GR[[nr]]<-start$G[[r]]
        }

        if(length(grep("animal", rmodel.terms[r]))>0){
          Aterm[nr]<-1
        }else{
          Aterm[nr]<-0
        }
        nr<-nr+1
      }  		 
    }

     if(length(rmodel.terms)>1){   # random effects exist
       for(i in 1:(length(rmodel.terms)-1)){
         if(length(grep("\\.group_MCMC\\.", rmodel.terms[i]))>0){
           components<-strsplit(rmodel.terms[i],"\\.group_MCMC\\.")[[1]]
           data_comb<-data[,components[1]]
           xfactor<-rep(1,length(na.omit(data_comb)))
           rlevels<-levels(data_comb)
           data_pos<-match(data_comb,rlevels)
           Ztmp<-Matrix(0,dim(data)[1],length(table(data_pos)), dimnames=list(as.character(1:dim(data)[1]), paste(rmodel.terms[i], rlevels[as.numeric(names(table(data_pos)))], sep=".")))
           Ztmp[,1][2]<-1              # force it out of vbeing upper triangle!!!!
           Ztmp@p<-as.integer(c(0,cumsum(table(data_pos))))    
           cnt<-0
           for(j in 1:length(rlevels)){
             hit<-which(data_pos==j)
             hit<-hit-dim(data)[1]*(ceiling(hit/dim(data)[1])-1)
             if(length(hit)>0){
               Ztmp@i[cnt+1:length(hit)]<-as.integer(hit-1)
               cnt<-cnt+length(hit)
             }
           }  
           Ztmp@x<-xfactor
 
           Ztmp2<-Matrix(0,nlevels(data[,components[2]]),nlevels(data[,components[1]]))
           Ztmp2[,1][2]<-1              # force it out of vbeing upper triangle!!!!
           data_comb2<-data[,components[2]]
           xfactor2<-rep(1,length(na.omit(data_comb2)))
           rlevels2<-levels(data_comb2)

           cnt<-0
           cnt2<-1
           for(j in rlevels){
             hit<-which(rlevels2%in%data_comb2[which(data_comb==j)])
             if(length(hit)>0){
               Ztmp2@i[cnt+1:length(hit)]<-as.integer(hit-1)
               cnt<-cnt+length(hit)
             }
             cnt2<-cnt2+1
             Ztmp2@p[cnt2]<-as.integer(cnt)
           }  
           Ztmp2@x<-xfactor
           Ztmp<-Ztmp%*%t(Ztmp2)
         }else{
           if(length(grep("\\:", rmodel.terms[i]))>0){
             if(is.factor(data[,strsplit(rmodel.terms[i], "\\:")[[1]][2]])){
               data_comb<-data[,strsplit(rmodel.terms[i], "\\:")[[1]][2]]:data[,strsplit(rmodel.terms[i], "\\:")[[1]][1]]
               xfactor<-rep(1,length(na.omit(data_comb)))
             }else{
               dpol<-dim(data[,strsplit(rmodel.terms[i], "\\:")[[1]][2]])[2]
               data_comb<-as.factor(rep(1:dpol, each=dim(data)[1])):as.factor(rep(data[,strsplit(rmodel.terms[i], "\\:")[[1]][1]], dpol))
               xfactor<-c(data[,strsplit(rmodel.terms[i], "\\:")[[1]][2]])
               xfactor[which(is.na(xfactor))]<-0
             }
           }else{
             data_comb<-data[,rmodel.terms[i]]
             xfactor<-rep(1,length(na.omit(data_comb)))
           }
           rlevels<-levels(data_comb)
           data_pos<-match(data_comb,rlevels)
           Ztmp<-Matrix(0,dim(data)[1],length(table(data_pos)), dimnames=list(as.character(1:dim(data)[1]), paste(rmodel.terms[i], rlevels[as.numeric(names(table(data_pos)))], sep=".")))
           Ztmp[,1][2]<-1              # force it out of vbeing upper triangle!!!!
           Ztmp@p<-as.integer(c(0,cumsum(table(data_pos))))    
           cnt<-0
           for(j in 1:length(rlevels)){
             hit<-which(data_pos==j)
             hit<-hit-dim(data)[1]*(ceiling(hit/dim(data)[1])-1)
             if(length(hit)>0){
               Ztmp@i[cnt+1:length(hit)]<-as.integer(hit-1)
               cnt<-cnt+length(hit)
             }
           }  
           Ztmp@x<-xfactor
         }
         if(i==1){
           Z<-Ztmp
         }else{
           Z<-cBind(Z, Ztmp)          
         }
       }
     }else{
       Z<-as(matrix(0,1,0), "sparseMatrix")
     }

     nR<-nr-nG-1  # number of R structures
     if(sum(nfl[nG+1:nR]*nrl[nG+1:nR])!=dim(data)[1]){stop("R-structure does not define unique residual for each data point")}
     if(is.null(tune)){
       AMtune=c(rep(FALSE, nG), rep(TRUE, nR))
       for(i in 1:nR){
         tune[[i]] = diag(nfl[nG+i])
       }
     }else{
       AMtune=rep(FALSE, nR+nG)
       if(nR>1){
         tune<-sapply(diag(tune), as.matrix, simplify=FALSE)
       }else{
         tune<-list(tune)
       }
       for(i in 1:nR){
         if(dim(tune[[i]])[1]!= dim(tune[[i]])[2] |  dim(tune[[i]])[2]!= nfl[nG+i]){stop(paste("proposal distribution ", i, " is the wrong dimension"))}
         if(is.positive.definite(tune[[i]])==FALSE){stop(paste("proposal distribution ", i, " is not positive definite"))}
       }
     }	

############################
# Build Fixed Effect Model #
############################

       fmodel.terms<-strsplit(as.character(fixed)[3], " *\\+ *")[[1]]
       for(r in 1:length(fmodel.terms)){          
         if(length(grep("leg\\(", fmodel.terms[r]))>0){
           rr.covariate<-strsplit(fmodel.terms[r], "\\(|\\)|\\,|\\):")[[1]][c(2:4)]     
           legvar.name<-paste("leg", rr.covariate[1],  abs(as.numeric(rr.covariate[2])),sep=".")
           if(legvar.name%in%names(data)==FALSE){
             lp<-legendre.polynomials(abs(as.numeric(rr.covariate[2])), TRUE)
             if(as.numeric(rr.covariate[2])<0){
               lp<-lp[-1]
             }			
             legvar<-sapply(lp,function(lp){as.function(lp)(data[,rr.covariate[1]])})
             data$legvar<-legvar
             names(data)[which(names(data)=="legvar")]<-legvar.name
           }
           if(is.na(rr.covariate[3])){
             fmodel.terms[r]<-legvar.name
           }else{
             fmodel.terms[r]<-paste(legvar.name,rr.covariate[3],sep=":")
           }
         }
       }

   X<-model.matrix(as.formula(paste("~", paste(fmodel.terms, collapse="+"), sep="")),data)
   X[which(is.na(X))]<-0 
   if(any(apply(X,2, function(x){all(x==0)}))){
     warning(paste("fixed effects ", paste(colnames(X)[which(apply(X,2, function(x){all(x==0)}))], collapse=" "),  " not represented in data and have been deleted"))
     X<-X[,-which(apply(X,2, function(x){all(x==0)}))] 
   }		 
   X<-as(X, "sparseMatrix")

   if(is.null(prior$B)){
      prior$B=list(V=diag(dim(X)[2])*1e+10, mu=matrix(0,dim(X)[2],1))
   }else{
	   if(any(prior$B$mu!=0)){stop("sorry; non-zero mean vector for fixed effect prior not yet implemented")}
	   if(is.matrix(prior$B$mu)==FALSE){
		   prior$B$mu<-matrix(prior$B$mu,  length(prior$B$mu), 1)
	   }
	   if(is.matrix(prior$B$V)==FALSE){
		   prior$B$V<-as.matrix(prior$B$V)
	   }	   
	   if(dim(X)[2]!=dim(prior$B$mu)[1]){stop("fixed effect mu prior is the wrong dimension")}
	   if(dim(X)[2]!=dim(prior$B$V)[1] | dim(X)[2]!=dim(prior$B$V)[2]){stop("fixed effect V prior is the wrong dimension")}
   }	   

################
# Missing data #
################

    cnt<-1
    mvtype<-c()
    proposal<-c()
    for(i in 1:nR){
       mp<-matrix(data$MCMC_y[cnt+1:(nfl[i+nG]*nrl[i+nG])-1], nrl[i+nG], nfl[i+nG])
       fp<-matrix(match(data$MCMC_family.names, family.types)[cnt+1:(nfl[i+nG]*nrl[i+nG])-1], nrl[i+nG], nfl[i+nG])
       proposal<-c(proposal, AMtune[i+nG]*((fp==11)*(mp==0 & is.na(mp)==FALSE))[,1])
       missing.pattern<-fp*(is.na(mp)==FALSE) # 0 for missing data 1 for gaussian >1 for other
       # missing data codes: 2 complete gaussian (ignore);
       #                     1 completely missing (unconditional Gibbs)
       #                     0 partial missing with observed being Guassian (conditional Gibbs)
       #                    -1 partial or fully observed with non-gaussian (MH)  currently 0 & -1 are both MHed. 

       mvtype_tmp<-apply(missing.pattern, 1,function(x){any(x!=1 & x!=0)})*-1                 
       mvtype_tmp[which(apply(missing.pattern, 1,function(x){any(x!=1)})==FALSE)]<-2          
       mvtype_tmp[which(apply(missing.pattern, 1,function(x){any(x!=0)})==FALSE)]<-1       
       if(all(mvtype_tmp>0)){AMtune[i+nG]=FALSE}  # If everything can be gibbsed do not tune
       mvtype<-c(mvtype,mvtype_tmp) 
       cnt<-cnt+nfl[i+nG]*nrl[i+nG]
    }
	
    if(is.null(start$liab)){
       data$MCMC_liab<-rep(NA, length(data$MCMC_y)) 
       if(QUASI==TRUE){ 
        for(i in as.character(unique(data$trait:as.factor(data$MCMC_family.names)))){
          trait_set<-which((data$trait:as.factor(data$MCMC_family.names))==i)
          data_tmp<-data[trait_set,]
          missing_set<-which(is.na(data_tmp$MCMC_y))
          if(length(missing_set)==length(trait_set)){
            mu<-0
            v<-1
            data$MCMC_liab[trait_set]<-rnorm(length(trait_set), mu,v)
          }else{
            if(data_tmp$MCMC_family.names[1]=="poisson" | data_tmp$MCMC_family.names[1]=="cenpoisson" ){
              mu<-mean((data_tmp$MCMC_y+1), na.rm=TRUE)
              v<-log(((var(data_tmp$MCMC_y+1, na.rm=TRUE)-mu)/(mu^2))+1)
              mu<-log(mu)-0.5*v
            }
            if(data_tmp$MCMC_family.names[1]=="multinomial"){
              if(length(table(data_tmp$MCMC_y))>2){
                 m1<-summary(glm(cbind(MCMC_y, MCMC_y.additional)~1, family="quasibinomial", data=data_tmp,subset=is.na(data_tmp$MCMC_y)==FALSE))
                 v<-abs(((as.numeric(m1$dispersion[1])-0.5)/2)^2)
                 mu<-as.numeric(m1$coef[1])
              }else{
                 v<-1
                 mu<-0
              }
            }
            if(data_tmp$MCMC_family.names[1]=="exponential" | data_tmp$MCMC_family.names[1]=="cenexponential"){
              m1<-summary(glm(MCMC_y~1, family="Gamma", data=data_tmp,subset=is.na(data_tmp$MCMC_y)==FALSE))
              v<-abs((as.numeric(m1$dispersion[1])-1)/2)
              mu<-as.numeric(m1$coef[1])
            }
            if(data_tmp$MCMC_family.names[1]=="cengaussian"){ 
              v<-var(apply(cbind(data_tmp$MCMC_y,data_tmp$MCMC_y.additional),1,function(x){min(abs(x))}) , na.rm=T)
              mu<-mean(apply(cbind(data_tmp$MCMC_y,data_tmp$MCMC_y.additional),1,function(x){min(abs(x))}) , na.rm=T)
            }
            if(data_tmp$MCMC_family.names[1]=="gaussian"){ 
              v<-var(data_tmp$MCMC_y, na.rm=T)
              mu<-mean(data_tmp$MCMC_y, na.rm=T)
            }
            if(data_tmp$MCMC_family.names[1]=="zipoisson"){
              if(max(data_tmp$MCMC_y)==1){
                mu<-logit(mean(data_tmp$MCMC_y==1))
                v<-diag(GRprior[[nR]]$V)[length(diag(GRprior[[nR]]$V))]
              }else{
                data_tmp<-data_tmp[-which(data_tmp$MCMC_y==0),]  
                mu<-mean(data_tmp$MCMC_y, na.rm=TRUE)
                v<-log(((var(data_tmp$MCMC_y, na.rm=TRUE)-mu)/(mu^2))+1)
                mu<-log(mu)-0.5*v
              }
            }
            if(length(missing_set)>0){
              data$MCMC_liab[trait_set][missing_set]<-sample(data$MCMC_liab[trait_set][-missing_set], length(missing_set), TRUE)
            }
            l_tmp<-sort(rnorm(length(trait_set), mu, sqrt(v)))
            data$MCMC_liab[trait_set][order(data$MCMC_liab[trait_set])]<-l_tmp
          }
        }
      }	
      if(any(is.na(data$MCMC_liab))){
         warning("good starting values not obtained: using Norm(0,1)")
         data$MCMC_liab<-rnorm(length(data$MCMC_liab),0,1)
      }
    }else{
      if(length(c(start$liab))!=length(data$MCMC_y)){stop("liabilities must have the same dimensions as the response")}
      if(any(is.na(start$liab))){stop("starting liabilities must not contain missing vlaues")}
      data$MCMC_liab<-c(start$liab)
    }

    observed<-as.numeric(is.na(data$MCMC_y)==FALSE)
    data$MCMC_y[is.na(data$MCMC_y)]<-0
    data$MCMC_y.additional[is.na(data$MCMC_y.additional)]<-0
    data$MCMC_y[which(data$MCMC_y==-Inf | data$MCMC_y==Inf)]<-sign(data$MCMC_y[which(data$MCMC_y==-Inf | data$MCMC_y==Inf)])*1e+32
    data$MCMC_y.additional[which(data$MCMC_y.additional==-Inf | data$MCMC_y.additional==Inf)]<-sign(data$MCMC_y.additional[which(data$MCMC_y.additional==-Inf | data$MCMC_y.additional==Inf)])*1e+32
	
    if(any(substr(data$MCMC_family.names, 1,3)=="cen" & data$MCMC_y==data$MCMC_y.additional)){      # replace liabilities of family of censored variables with y=y.additional          
      cen_areknown<-which(substr(data$MCMC_family.names, 1,3)=="cen" & data$MCMC_y==data$MCMC_y.additional)                          
      data$MCMC_family.names[cen_areknown]<-substr(data$MCMC_family.names[cen_areknown], 4, nchar(as.character(data$MCMC_family.names[cen_areknown])))
    }
    if(any(data$MCMC_family.names=="gaussian")){                                       # replace liabilities of ovserved gaussian data with data                                    
      data$MCMC_liab[which(data$MCMC_family.names=="gaussian" & observed)]<-data$MCMC_y[which(data$MCMC_family.names=="gaussian" & observed)]
    }

    split<-unlist(lapply(GRprior, function(x){if(is.null(x$fix)){return(-998)}else{return(x$fix)}}))
    if(any(split>nfl | (split<1 & split!=-998))){stop("fix term in priorG/priorR must be at least one less than the dimension of V")}
    update<-as.numeric(split!=1)

    GRinv<-unlist(lapply(GR, function(x){c(solve(x))}))
    GRvpP<-lapply(GRprior, function(x){(x$V)*(x$n)})
    for(i in 1:length(GRprior)){
      if(is.null(GRprior[[i]]$fix)==FALSE){
        Crows<-GRprior[[i]]$fix:dim(GRprior[[i]]$V)[1]
        GRvpP[[i]][Crows,Crows]<-GRprior[[i]]$V[Crows,Crows,drop=FALSE]*(nrl[i]+GRprior[[i]]$n)
      }
    }

    GRvpP<-unlist(GRvpP)
    GRnpP<-unlist(lapply(GRprior, function(x){c(x$n)}))
    BvpP<-c(solve(prior$B$V))
    BmupP<-c(prior$B$mu)
	
    diagP<-rep(0, length(nfl))
    if(diagR){
      diagP[length(nfl)]<-1
    }
    data$MCMC_family.names<-match(data$MCMC_family.names, family.types)     # add measurement error variances and y.additional 

    if(nitt%%1!=0){stop("nitt must be integer")}
    if(thin%%1!=0){stop("thin must be integer")}
    if(burnin%%1!=0){stop("burnin must be integer")}

    nkeep<-ceiling((nitt-burnin)/thin)

    if(nkeep<1){stop("burnin is equal to, or greater than number of iterations")}

    Loc<-1:((sum((nfl*nrl)[1:nG])*pr+dim(X)[2])*nkeep)
    dbar<-1:(2+nkeep)

    if(pl==TRUE){
      PLiab<-1:(length(data$MCMC_y)*nkeep)
    }else{
      PLiab<-1
    }
    Var<-1:(length(GRinv)*nkeep)

    if(sum(Aterm)!=0){
      id<-match(Ai$node.names,Ai$pedigree[,1])         
      dam<-match(Ai$pedigree[,2], Ai$pedigree[,1])
      dam[which(is.na(dam))]<--998
      sire<-match(Ai$pedigree[,3], Ai$pedigree[,1])
      sire[which(is.na(sire))]<--998
      MSsd<-sqrt(Ai$dii)
      PedDim<-c(length(dam), 2-Ai$phylogeny)
      A<-Ai$Ainv
    }else{
      A<- as(diag(1), "sparseMatrix")
      MSsd<-sqrt(0.5)
      id<--998
      dam<--998
      sire<--998
      PedDim<-c(1,1)
    }

    if(length(mfac)==0){mfac=-999}
    	
	output<-.C("MCMCglmm",
        as.double(data$MCMC_y),   
        as.double(data$MCMC_y.additional),
        as.double(data$MCMC_liab), 
        as.integer(mvtype),   
        as.integer(length(data$MCMC_y)),
        as.integer(X@i),         
        as.integer(X@p),        
        as.double(X@x),         
        as.integer(X@Dim),         
    	as.integer(length(X@x)),	  
        as.integer(Z@i),         
        as.integer(Z@p),          
        as.double(Z@x),         
        as.integer(Z@Dim),         
    	as.integer(length(Z@x)),
        as.integer(A@i),         
        as.integer(A@p),       
        as.double(A@x),         
        as.integer(A@Dim),         
    	as.integer(length(A@x)),
    	as.double(MSsd),
      	as.integer(id-1),
    	as.integer(dam-1),
    	as.integer(sire-1),
        as.integer(PedDim),
    	as.integer(Aterm),	  	  	  	  
        as.integer(nfl),
        as.integer(nrl),
        as.integer(update),
        as.integer(split-1),
        as.integer(c(nG,nR)), 
        as.double(GRinv), 
        as.double(GRvpP),            
        as.double(GRnpP),            
        as.integer(nitt),
        as.integer(thin),
        as.integer(burnin),
        as.integer(c(pr, pl)),
        as.double(Loc),
        as.double(Var),
        as.double(PLiab),
        as.integer(data$MCMC_family.names),
        as.double(c(unlist(tune))),
        as.integer(verbose),
        as.double(BvpP),
        as.double(BmupP),
        as.integer(mfac), 
	as.integer(observed),
        as.integer(diagP),
        as.integer(AMtune),
	as.integer(DIC),
        as.double(dbar),	  
        as.integer(proposal) 
        )

        Sol<-t(matrix(output[[39]], sum((nfl*nrl)[1:nG])*pr+dim(X)[2], nkeep))
        if(pr){      
          colnames(Sol)<-c(colnames(X),colnames(Z))
        }else{
          colnames(Sol)<-c(colnames(X))
        }
        colnames(Sol)<-gsub("MCMC_", "", colnames(Sol))

        VCV<-t(matrix(output[[40]], length(GRinv), nkeep))
        colnames(VCV)<-variance.names
        colnames(VCV)<-gsub("MCMC_", "", colnames(VCV))

        if(DIC==TRUE){
         deviance<--2*output[[52]][1:nkeep]
         DIC<--4*output[[52]][nkeep+1]+2*output[[52]][nkeep+2]
        }else{
         deviance<-NULL
         DIC<-NULL
        }
        if(pl==TRUE){
          Liab<-mcmc(t(matrix(output[[41]], length(data$MCMC_y), nkeep)))
        }else{
          Liab<-NULL
        }
    	options("na.action"="na.omit")

        list(Sol=mcmc(Sol), VCV=mcmc(VCV), Liab=Liab, Fixed=original.fixed, Random=original.random, Residual=original.rcov, Deviance=mcmc(deviance),DIC=DIC)
	
}

