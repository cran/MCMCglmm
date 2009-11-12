"MCMCglmm"<-function(fixed, random=NULL, rcov=~units, family="gaussian", mev=NULL, data, start=NULL, prior=NULL, tune=NULL, pedigree=NULL, nodes="ALL",scale=TRUE,  nitt=13000, thin=10, burnin=3000, pr=FALSE, pl=FALSE, verbose=TRUE, DIC=TRUE, singular.ok=FALSE, saveX=FALSE, saveZ=FALSE){
   
    orig.na.action<-options("na.action")[[1]]
    options("na.action"="na.pass")	

    if(is.null(data)){stop("data argument is NULL")}
    if(is.null(fixed)){stop("fixed is NULL")}

    if(any(names(data)%in%c("units", "MCMC_y", "MCMC_y.additional","MCMC_liab","MCMC_meta", "MCMC_mev", "MCMC_family.names"))){
      stop(paste(names(data)[which(names(data)%in%c("units", "MCMC_y", "MCMC_y.additional","MCMC_liab","MCMC_meta", "MCMC_mev", "MCMC_family.names"))], " is a reserved variable please rename it"))
    }

    family.types<-c("gaussian", "poisson", "multinomial", "notyet_weibull", "exponential", "cengaussian", "cenpoisson", "notyet_cenweibull", "cenexponential",  "notyet_zigaussian", "zipoisson", "notyet_ziweibull", "notyet_ziexponential", "ordinal")

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


######################################################################################
# for phyloegnetic/pedigree analyses form A and augment with missing nodes if needed #
###################################################################################### 	
	
	if(is.null(pedigree)==FALSE){
	  if(is.null(data$animal) |length(grep("animal", random))==0){stop("pedigree/phylogeny has been passed, but no 'animal' variable")}
	  Ai<-inverseA(pedigree, nodes=nodes, scale=scale)
          if(any(data$animal%in%Ai$node.names==FALSE)){stop("individuals in data not appearing in pedigree/phylogeny")}
	  data$animal<-factor(data$animal, levels=Ai$node.names)                           # add factor levels to animal in data for additional nodes
	}else{
          if(any(colnames(data)=="animal") & length(grep("animal", as.character(random))>0)){
            stop("animal is a reserved variable for pedigrees/phylogenies only, please rename")
          }
        }

##############################################################################################
# if R-structure is of form idh(!=trait) assign new records with arbitrary levels of !-trait #
##############################################################################################

      nadded<-0  # number of records added

      rterms<-c(split.direct.sum(as.character(random)[2]), split.direct.sum(as.character(rcov)[2]))
      data$MCMC_dummy<-rep(1,dim(data)[1])

      for(i in 1:length(rterms)){

	components<-find.components(rterms[i], data)
      
	if(is.null(components[[1]])==FALSE | length(components[[2]])>0){

          for(j in 1:length(components[[2]])){                                                                                

            MCMC_components1<-interaction(data[,components[[1]]], sep=".MCMC.", drop=TRUE*(length(components[[1]])>1))  	# random effect factors 
            MCMC_components2<-interaction(data[,components[[2]][[j]]], sep=".MCMC.", drop=FALSE)  	                        # fixed effect factors
            allc<-expand.grid(levels(MCMC_components1),levels(MCMC_components2))  						# all pairwise combinations
            missing.combinations<-allc[which(paste(allc[,1], allc[,2])%in%paste(MCMC_components1, MCMC_components2)==FALSE),]   # missing pairwise combinations

            if(dim(missing.combinations)[1]>0){
			  
     	      warning(paste("some combinations in", rterms[i], "do not exist and", dim(missing.combinations)[1], "missing records have been generated"))

              ###################################
              # find if there are already holes #
              #    in the missing records       #
              ###################################
			  			  
	      missing.comb1<-which(is.na(MCMC_components1) & is.na(MCMC_components2)==FALSE)  # na's for random but have fixed
	      missing.comb2<-which(is.na(MCMC_components1)==FALSE & is.na(MCMC_components2))  # na's for fixed bt have random
	      missing.comb12<-which(is.na(MCMC_components1) & is.na(MCMC_components2))        # na's for both

	      matching<-match(unique(missing.combinations[,2]), MCMC_components2[missing.comb1])  # some places that can be filled
			
	      while(any(is.na(matching)==FALSE)){
                a.match<-match(MCMC_components2[missing.comb1[na.omit(matching)]],missing.combinations[,2])  # missing combinations to be placed
	        data[missing.comb1[na.omit(matching)],unlist(components[[1]])]<-as.matrix(apply(missing.combinations[a.match,1,drop=FALSE], 1, function(x){unlist(strsplit(x, "\\.MCMC\\."))}))
                MCMC_components1[missing.comb1[na.omit(matching)]]<-missing.combinations[a.match,1]   # fill in data columns and MCMC_componenets1
       	        missing.comb1<-missing.comb1[-na.omit(matching)]	
	        missing.combinations<-missing.combinations[-a.match,, drop=FALSE]
	        matching<-match(unique(missing.combinations[,2]), MCMC_components2[missing.comb1])
	      }

	      matching<-match(unique(missing.combinations[,1]), MCMC_components1[missing.comb2])  # some places that can be filled
			
	      while(any(is.na(matching)==FALSE)){
                a.match<-match(MCMC_components1[missing.comb2[na.omit(matching)]],missing.combinations[,1])  # missing combinations to be placed
	        data[missing.comb2[na.omit(matching)],unlist(components[[2]][[j]])]<-as.matrix(apply(missing.combinations[a.match,2,drop=FALSE], 1, function(x){unlist(strsplit(x, "\\.MCMC\\."))}))
                MCMC_components2[missing.comb2[na.omit(matching)]]<-missing.combinations[a.match,2]   # fill in data columns and MCMC_componenets1
       	        missing.comb2<-missing.comb2[-na.omit(matching)]	
   		missing.combinations<-missing.combinations[-a.match,, drop=FALSE]
	        matching<-match(unique(missing.combinations[,1]), MCMC_components1[missing.comb2])
	      }

              n.rm<-min(length(missing.comb12), dim(missing.combinations)[1])

	      if(n.rm>0){
	        data[missing.comb12,c(components[[1]], components[[2]][[j]])]<-t(apply(missing.combinations[1:n.rm,], 1, function(x){unlist(strsplit(x, "\\.MCMC\\."))}))
                MCMC_components1[missing.comb12][1:n.rm]<-missing.combinations[,1][1:n.rm]
                MCMC_components2[missing.comb12][1:n.rm]<-missing.combinations[,2][1:n.rm]
	      }	

	      if(length(missing.combinations)>0){       
                nadded<-nadded+dim(missing.combinations)[1]                                               # add dummy records if still needed
		data[dim(data)[1]+1:dim(missing.combinations)[1],]<-NA                  
		data[,c(components[[1]], components[[2]][[j]])][dim(data)[1]-(dim(missing.combinations)[1]-1):0,]<-t(apply(missing.combinations, 1, function(x){unlist(strsplit(x, "\\.MCMC\\."))}))
                if(is.null(mev)==FALSE){
                  mev<-c(mev,rep(1, dim(missing.combinations)[1]))
                }
	      }
	    }
          }		  
	}
      }
      data$MCMC_dummy<-rep(0,dim(data)[1])
      if(nadded>0){
        data$MCMC_dummy[(dim(data)[1]-nadded+1):dim(data)[1]]<-1   
      }
  
    MVasUV=FALSE  # Multivaraite as Univariate
    if(is.null(family)){
       if(is.null(data$family)){
         stop("no family specified")
       }else{
         if(is.null(data$trait)){
           stop("data.frame needs to have a column indexing traits if fitting multivaraite models as univariate")
         }else{
           if(is.factor(data$trait)==FALSE){
             stop("trait must be a factor")
           }
           if(any(tapply(data$family, data$trait, function(x){length(unique(x))})!=1)){
             stop("all data from the same trait must come from the same distribution")
           }
         }           
         if(length(unique(data$family))==1){
           family.names<-as.character(data$family[1])
         }else{     
           family.names<-as.character(data$family)   
           MVasUV=TRUE
         }
         if(length(grep("cen|multinomial|zi", family.names))>0){ 
           stop("For setting up multi-trait models as univariate the responses cannot come from distributions that require more than one data column or have more than one liability (i.e. censored, multinomial, zero-inflated, categorical with k>2): set it up as multivariate using cbind(...)")
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

    diagR=FALSE  # should the residual structure be diagonal even if us is sued?
	
    if(length(grep("zi|multinomial[3:19]", family.names))>0){ 
      if(length(grep("idh\\(|us\\(", rcov))==0){
        stop("please use idh() or us() error structure")
      }else{
        if(length(grep("idh\\(", rcov))>0){
          diagR=TRUE
        }
        rcov=as.formula(paste(gsub("idh\\(", "us\\(", rcov), collapse=""))
      }	  
    }
	
    nS<-dim(data)[1]                               # number of subjects
    y.additional<-matrix(NA, nS,0)                 # matrix equal in dimension to y holding the additional parameters of the distribution (n, upper interval etc.)
    nt<-1                                          # number of traits (to be iterated because y may change dimension with multinomial/categorical/censored distributions)

    if(MVasUV){
      mfac<-rep(0, nlevels(data$trait))
      if(any(family.names=="categorical")){
        ncat<-tapply(data[,response.names][which(family.names=="categorical")], data$trait[which(family.names=="categorical")], function(x){length(unique(x))})
        if(any(ncat)>2){
           stop("For setting up multi-trait models as univariate the responses cannot come from distributions that require more than one data column or have more than one liability (i.e. censored, multinomial, zero-inflated, categorical with k>2): set it up as multivariate using cbind(...)")
         }
        y.additional<-matrix(NA, nS,1)
        y.additional[which(family.names=="categorical")]<-1
        family.names[which(family.names=="categorical")]<-"multinomial"
        family.names[which(is.na(family.names))]<-"gaussian"
      }
      if(any(family.names=="ordinal")){
        ordinal.traits<-unique(as.numeric(data$trait)[which(family.names=="ordinal")])
        mfac[ordinal.traits]<-0:(length(ordinal.traits)-1)
        ncutpoints<-tapply(data[,response.names][which(family.names=="ordinal")], data$trait[which(family.names=="ordinal")], function(x){length(unique(x))})+2
      }else{
        ncutpoints<-c()
      }
    }else{
      if(any(names(data)=="trait")){
        stop("trait is a reserved variable please remove or rename this column in the data.frame")
      }
      mfac<-c()                                    # stores additional information for each R-structure term. For multinomial models it stores the number of k-2 categories
                                                   # for zero inflated models it indicates whether the R-structure term is for the Poisson part (0) or zero infaltion part (1)
                                                   # for ordinal responses it indicates whether its the i^th ordinal response. For example R-structure terms for data A:E with 
                                                   # multinomialA, multinomialA, multinomialB, ordinalC, zero-inflationD, zero-infaltionD, ordinalE: mfac=c(1,1,0,0,0,1,1) other 
                                                   # traits all have 0
      ncutpoints<-c()                              # number of cutpoints for ordinal variables = k+1
      ordinal.names<-c()

      for(i in 1:length(family)){

        dist.preffix<-substr(family[i],1,2)                  
        if(nt>length(response.names)){stop("family is the wrong length")}
        if(any(dist.preffix%in%c("ce", "mu", "ca", "tr", "zi"))){

######################
# categorical traits #
######################

          if(dist.preffix=="ca"){
            cont<-as.matrix(model.matrix(~ as.factor(data[[response.names[nt]]]))[,-1])                            # form new J-1 variable   
            nJ<-dim(cont)[2]                                                                                       # number of J-1 categories    
            if(length(grep("idh\\(|us\\(", rcov))==0 & nJ>1){
               stop("please use idh() or us() error structure for categorical traits with more than 2 categories")
            }else{
               if(length(grep("idh\\(", rcov))>0){
                 rcov=as.formula(paste(gsub("idh\\(", "us\\(", rcov), collapse=""))
                 diagR=TRUE
               }
            }
            mfac<-c(mfac, rep(nJ-1,nJ))             
            new.names<-paste(response.names[nt], ".", levels(as.factor(data[[response.names[nt]]]))[-1], sep="")   # give new variables names
            colnames(cont)<-new.names  
            data<-data[,-which(names(data)==response.names[nt]), drop=FALSE]    # remove original variable
	    row.names(data)<-row.names(cont)
            data<-cbind(data, cont)     # add new variables to data.frame
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
           mfac<-c(mfac, rep(nJ-1,nJ))  
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
           mfac<-c(mfac, 0)  
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
           mfac<-c(mfac, 0)  
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
	   mfac<-c(mfac, c(0,1))
	   ones<-rep(1, length(response.names))
	   ones[which(response.names==response.names[nt])]<-2
	   response.names<-rep(response.names, ones)
	   family.names<-rep(family.names, ones)
	   family.names[which(response.names==response.names[nt])]<-family.names[nt]
	   response.names[which(response.names==response.names[nt])]<-c(response.names[nt], paste("zi", response.names[nt], sep="."))
	   nt<-nt+2			  
	 }

###############################################
# gaussian/poisson/exponential/ordinal traits #
###############################################

        }else{
          mfac<-c(mfac, 0)  
          if(family.names[nt]=="poisson"){ 
	    if(all(data[,response.names[nt]]%%1==0, na.rm=T)==FALSE | all(data[,response.names[nt]]>=0, na.rm=T)==FALSE){stop("Poisson data must be positive integers")}
          }
	  if(family.names[nt]=="exponential"){ 
	    if(any(data[,response.names[nt]]<0, na.rm=T)){stop("Exponential data must be positive")}
	  }	
	  if(family.names[nt]=="ordinal"){
            mfac[length(mfac)]<-length(ncutpoints)  
            data[,response.names[nt]]<-as.numeric(as.factor(data[,response.names[nt]]))
            ncutpoints<-c(ncutpoints, max(data[,response.names[nt]])+1)         
            ordinal.names<-c(ordinal.names, response.names[nt])     
	  }	
          y.additional<-cbind(y.additional,matrix(NA,nS,1))     
          nt<-nt+1
        }
      }	
      nt<-nt-1
    }
    if(sum((family.names%in%family.types)==FALSE)!=0){stop(paste(unique(family[which((family.names%in%family.types)==FALSE)]), "not a supported distribution"))}

###**************************************########################

    if(MVasUV){
      data<-reshape(data, varying=response.names, v.names="MCMC_y", direction="long", idvar="units")       # reshape the data into long format 
      data$MCMC_family.names<-family.names
    }else{
      data<-reshape(data, varying=response.names, v.names="MCMC_y", direction="long", timevar="trait", idvar="units")       # reshape the data into long format 
      data$trait<-factor(response.names[data$trait], response.names)
      if(length(response.names)!=length(family.names)){stop("family must have the same length as the number of responses")}
      data$MCMC_family.names<-rep(family.names, each=nS)       
    }
    data$units<-as.factor(data$units)
    data$MCMC_y.additional<-c(y.additional) 

######################################################
# for (random) meta-analysis add weights/model terms #
###################################################### 	

    if(is.null(mev)==FALSE){
      if(any(dim(mev)!=dim(y.additional))){stop("mev has to be the same dimension as y")}
      data$MCMC_mev<-sqrt(mev)	
      data$MCMC_meta<-factor(1:dim(data)[1], levels=1:dim(data)[1])
      if(is.null(random)){
        random = ~us(leg(MCMC_mev, -1, FALSE)):MCMC_meta
        if(is.null(prior$R)==FALSE){
          prior$G<-list(G1=list(V=as.matrix(1), nu=1, fix=1))
        }
      }else{
        random<-update(random,~.+us(leg(MCMC_mev, -1, FALSE)):MCMC_meta)
        if(is.null(start$G)==FALSE){
          start$G[[length(start$G)+1]]<-as.matrix(1)
        }
        if(is.null(prior$G)==FALSE){
          prior$G[[length(prior$G)+1]]<-list(V=as.matrix(1), nu=1, fix=1)
        }
      } 
    }

    rmodel.terms<-split.direct.sum(as.character(random)[2])
    ngstructures<-length(rmodel.terms)
    rmodel.terms<-c(rmodel.terms, strsplit(as.character(rcov)[2], " *\\+ *")[[1]]) 
 
    nfl<-c()                                                       # number of fixed levels the random term is structured by
    nrl<-c()                                                       # number of random levels
    variance.names<-c()
    GR<-list()
    GRprior<-list()
    Aterm<-c()
    ordering<-c()
    trait.ordering<-c()
    missing<-c()

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

##########################
# Build Z ################
##########################

    nr<-1

    for(r in 1:length(rmodel.terms)){

       if(r==(ngstructures+1)){nG<-nr-1}  # number of (new) G structures

       if(r<length(rmodel.terms)){
         Zlist<-buildZ(rmodel.terms[r], data=data)
       }else{
         Zlist<-buildZ(rmodel.terms[r], data=data, formZ=FALSE, add.missing=FALSE)

         if(Zlist$missing!=0){                 # haven't been using this any more - seems to slow the cholesky factorisation up
           missing.rows<-data[dim(data)[1],]
           missing.rows$MCMC_y<-NA
           data[dim(data)[1]+1:Zlist$missing,]<-missing.rows
           Z<-rBind(Z, t(Matrix(missing, ncol(Z), Zlist$missing)))
           nadded<-TRUE
         }
         Zlist$missing<-NULL
       }

       nfl<-c(nfl,Zlist$nfl)
       nrl<-c(nrl,Zlist$nrl)
       Aterm<-c(Aterm, Zlist$Aterm)
       variance.names<-c(variance.names, Zlist$vnames)
       ordering<-c(ordering ,Zlist$ordering)
       trait.ordering<-c(trait.ordering ,Zlist$trait.ordering)
       missing<-c(missing ,Zlist$missing)

       if(NOpriorG==TRUE){
         GRprior[[nr]]<-list(V=diag(sum(Zlist$nfl)), nu=0)
         if(length(grep("MCMC_meta", rmodel.terms[r]))>0){
            GRprior[[nr]]<-list(V=as.matrix(1), fix=1)
         }
       }else{
         if(r<=ngstructures){
           if(length(prior$G)<r){stop("priorG/priorR have the wrong number of structures")}
           GRprior[[nr]]<-prior$G[[r]]
         }else{
           GRprior[[nr]]<-prior$R
         }
         if(is.null(GRprior[[nr]]$V)){stop("V not specified for some priorG/priorR elements")}
         if(is.matrix(GRprior[[nr]]$V)==FALSE){GRprior[[nr]]$V<-as.matrix(GRprior[[nr]]$V)}
         if(dim(GRprior[[nr]]$V)[1]!=sum(Zlist$nfl)  | dim(GRprior[[nr]]$V)[2]!=sum(Zlist$nfl)){stop("V is the wrong dimension for some priorG/priorR elements")}
         if(is.positive.definite(GRprior[[nr]]$V)==FALSE){stop("V is not positive definite for some priorG/priorR elements")}
         if(is.null(GRprior[[nr]]$fix)==FALSE){
           CM<-GRprior[[nr]]$V[GRprior[[nr]]$fix:dim(GRprior[[nr]]$V)[1],GRprior[[nr]]$fix:dim(GRprior[[nr]]$V)[1]]
           if(sum(CM!=0)>dim(GRprior[[nr]]$V)[1] & GRprior[[nr]]$fix>1){stop("sorry - matrices to be conditioned on must be diagonal")}           
           if(GRprior[[nr]]$fix!=1){
             if(is.null(GRprior[[nr]]$n)){stop("nu not specified for some priorG/priorR elements")}
           }else{
             GRprior[[nr]]$nu=1
           }
         }else{
           if(is.null(GRprior[[nr]]$n)){stop("nu not specified for some priorG/priorR elements")}
         }
       }
       if(NOstartG==TRUE){
         if(det(GRprior[[nr]]$V)<1e-8 & is.null(GRprior[[nr]]$fix)){
           GR[[nr]]<-GRprior[[nr]]$V+diag(dim(GRprior[[nr]]$V)[1])
         }else{
           GR[[nr]]<-GRprior[[nr]]$V
         }
       }else{
         GR[[nr]]<-start$G[[r]]
         if(r<=ngstructures){
           if(length(start$G)<r){stop("starting G/R has the wrong number of structures")}
           GR[[nr]]<-start$G[[r]]
         }else{
           GR[[nr]]<-start$R
         }
	 if(is.matrix(GR[[r]]$V)==FALSE){GR[[r]]$V<-as.matrix(GR[[r]]$V)}	
         if(dim(GR[[r]])[1]!=sum(Zlist$nfl)  | dim(GR[[r]])[2]!=sum(Zlist$nfl)){stop("V is the wrong dimension for some startG/startR elements")}
         if(is.positive.definite(GR[[r]])==FALSE){stop(paste("starting G/R structure", r, " is not positive definite"))}
       }
       if(Zlist$vtype[1]=="idh"){

         for(l in 1:length(Zlist$nfl)){

           if(NOpriorG==TRUE){
             GRprior[[nr]]<-list(V=matrix(1), nu=0)
           }else{
             if(r<=ngstructures){
               if(length(prior$G)<r){stop("priorG/priorR have the wrong number of structures")}
               GRprior[[nr]]<-list(V=as.matrix(prior$G[[r]]$V)[l,l,drop=FALSE], nu=prior$G[[r]]$n-(sum(Zlist$nfl)-1), fix=prior$G[[r]]$fix)
             }else{
               GRprior[[nr]]<-list(V=as.matrix(prior$R$V)[l,l,drop=FALSE], nu=prior$R$n-(sum(Zlist$nfl)-1), fix=prior$R$fix)
             }
             if(is.null(GRprior[[nr]]$fix)==FALSE){
               if(l>=GRprior[[nr]]$fix){
                 GRprior[[nr]]$fix<-1
                 GRprior[[nr]]$nu<-1
               }else{
                 GRprior[[nr]]$fix<-NULL
               }
             }
           }

           if(NOstartG==TRUE){
             if(det(GRprior[[nr]]$V)<1e-8 & is.null(GRprior[[nr]]$fix)){
               GR[[nr]]<-GRprior[[nr]]$V+diag(dim(GRprior[[nr]]$V)[1])
             }else{
               GR[[nr]]<-GRprior[[nr]]$V
             }
           }else{
             if(r<=ngstructures){
               GR[[nr]]<-as.matrix(diag(start$G[[r]])[l])
             }else{
               GR[[nr]]<-as.matrix(diag(start$R[[r]])[l])
             }
           } 
           nr<-nr+1
         }
       }else{
         nr<-nr+1
       }
       if(r<=ngstructures){
         if(r==1){
           Z<-Zlist$Z
         }else{
           Z<-cBind(Z, Zlist$Z)     
         }     
       }
     }

    if(any(duplicated(ordering))){stop("R-structure miss-specified: each residual must be unique to a data point")}
    if(any(range(ordering)!=c(1, dim(data)[1]))){stop("R-structure miss-specified: each data point must have a residual")}

    data<-data[ordering,]         

    if(length(rmodel.terms)==1){
       Z<-as(matrix(0,1,0), "sparseMatrix")
    }else{                                                     # rearrange data to match R-structure
       Z<-Z[ordering,]                                                     
    }

    mfac<-mfac[trait.ordering]

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

   fixed<-as.formula(paste("~",as.character(fixed)[3]))

   X<-model.matrix(fixed,data)

   if(nadded>0){
     X[which(data$MCMC_dummy==1),]<-0
#      X[ncol(X),]<-0
   }
   if(any(is.na(X))){stop("missing values in the fixed predictors")}
   if(singular.ok==FALSE){
     if(all(is.na(data$MCMC_y))){stop("all data are missing. Use singular.ok=TRUE to sample these effects, but use an informative prior!")}

     sing.rm<-lm(data$MCMC_y~X-1, subset=is.na(data$MCMC_y)==FALSE)

     sing.rm<-which(is.na(sing.rm$coef))
     if(length(sing.rm)){
       warning("some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!")
       X<-X[,-sing.rm]
     }
   }	
   X<-as(X, "sparseMatrix")
   if(any(rowSums(X)==0)){
     X[,1][which(rowSums(X)==0)]<-1e-18
   }

#####################################
# Check prior for the fixed effects #
#####################################

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
       #             o        1 completely missing (unconditional Gibbs)
       #                     0 partial missing with observed being Guassian (conditional Gibbs)
       #                    -1 partial or fully observed with non-gaussian (MH)  currently 0 & -1 are both MHed. 

       mvtype_tmp<-apply(missing.pattern, 1,function(x){any(x!=1 & x!=0)})*-1                 
       mvtype_tmp[which(apply(missing.pattern, 1,function(x){any(x!=1)})==FALSE)]<-2          
       mvtype_tmp[which(apply(missing.pattern, 1,function(x){any(x!=0)})==FALSE)]<-1       
       if(all(mvtype_tmp>0)){AMtune[i+nG]=FALSE}  # If everything can be gibbsed do not tune
       mvtype<-c(mvtype,mvtype_tmp) 
       cnt<-cnt+nfl[i+nG]*nrl[i+nG]
    }

    stcutpoints<-c()	
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
            if(data_tmp$MCMC_family.names[1]=="poisson" | data_tmp$MCMC_family.names[1]=="cenpoisson"){
              mu<-mean((data_tmp$MCMC_y+1), na.rm=TRUE)
              v<-abs(log(((var(data_tmp$MCMC_y+1, na.rm=TRUE)-mu)/(mu^2))+1))
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
            if(data_tmp$MCMC_family.names[1]=="ordinal"){
              v<-1
              cps<-qnorm(cumsum(c(0,table(data_tmp$MCMC_y)/length(data_tmp$MCMC_y))), 0, sqrt(2))
              mu<-cps[2]
              cps<-cps-cps[2]
              cps[1]<--1e+64
              cps[length(cps)]<-1e+64
              stcutpoints<-c(stcutpoints, cps)              
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
              if(max(data_tmp$MCMC_y, na.rm=T)==1){
                mu<-logit(mean(data_tmp$MCMC_y==1))
                v<-diag(GRprior[[nR]]$V)[length(diag(GRprior[[nR]]$V))]
              }else{
                data_tmp<-data_tmp[-which(data_tmp$MCMC_y==0),]  
                mu<-mean(data_tmp$MCMC_y, na.rm=TRUE)
                v<-abs(log(((var(data_tmp$MCMC_y, na.rm=TRUE)-mu)/(mu^2))+1))
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
    GRvpP<-unlist(GRvpP)
    GRnpP<-unlist(lapply(GRprior, function(x){c(x$n)}))
    BvpP<-c(solve(prior$B$V))
    BmupP<-c(prior$B$mu)
	
    diagP<-rep(0, length(nfl))
    if(diagR){
      diagP[length(nfl)]<-1
    }

    nordinal<-length(ncutpoints)
    if(nordinal==0){
      ncutpoints<-1
      stcutpoints<-1
    }  # no cutpoints need to be estimated if cutpoints=3 (standard binary)
    ncutpoints_store<-sum((ncutpoints-3)*(ncutpoints>3))

    data$MCMC_family.names<-match(data$MCMC_family.names, family.types)     # add measurement error variances and y.additional 

    if(nitt%%1!=0){stop("nitt must be integer")}
    if(thin%%1!=0){stop("thin must be integer")}
    if(burnin%%1!=0){stop("burnin must be integer")}

    nkeep<-floor((nitt-burnin)/thin)

    if(nkeep<1){stop("burnin is equal to, or greater than number of iterations")}

    Loc<-1:((sum((nfl*nrl)[1:nG])*pr+dim(X)[2])*nkeep)

    if(ncutpoints_store>0){
      CP<-1:(ncutpoints_store*nkeep)
    }else{
      CP<-1
    }
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
        as.integer(proposal),
        as.integer(ncutpoints),
        as.integer(nordinal),
        as.double(stcutpoints),
        as.double(CP)
        )

        if(ncutpoints_store==0){
          Sol<-t(matrix(output[[39]], sum((nfl*nrl)[1:nG])*pr+dim(X)[2], nkeep))
          if(pr){      
            colnames(Sol)<-c(colnames(X),colnames(Z))
          }else{
            colnames(Sol)<-c(colnames(X))
          }
          colnames(Sol)<-gsub("MCMC_", "", colnames(Sol))
       }else{
          Sol<-matrix(0,nkeep, sum((nfl*nrl)[1:nG])*pr+dim(X)[2]+ncutpoints_store)
          Sol[,1:ncutpoints_store]<-t(matrix(output[[57]],ncutpoints_store, nkeep))
          Sol[,ncutpoints_store+1:(sum((nfl*nrl)[1:nG])*pr+dim(X)[2])]<-t(matrix(output[[39]],sum((nfl*nrl)[1:nG])*pr+dim(X)[2], nkeep))
          if(pr){      
            colnames(Sol)<-c(paste("cutpoint.trait", rep(ordinal.names, ncutpoints-3), ".", unlist(sapply(ncutpoints-3, function(x){1:x})), sep=""), colnames(X),colnames(Z))
          }else{
            colnames(Sol)<-c(paste("cutpoint.trait", rep(ordinal.names, ncutpoints-3), ".", unlist(sapply(ncutpoints-3, function(x){1:x})), sep=""), colnames(X))
          }
          colnames(Sol)<-gsub("MCMC_", "", colnames(Sol))
       }

        VCV<-t(matrix(output[[40]], length(GRinv), nkeep))
        colnames(VCV)<-variance.names
        colnames(VCV)<-gsub("MCMC_", "", colnames(VCV))

        if(DIC==TRUE){
         deviance<-mcmc(-2*output[[52]][1:nkeep], start=burnin+1, end=nitt-(nitt-burnin)%%thin, thin=thin)
         DIC<--4*output[[52]][nkeep+1]+2*output[[52]][nkeep+2]
        }else{
         deviance<-NULL
         DIC<-NULL
        }
        if(pl==TRUE){
          Liab<-mcmc(t(matrix(output[[41]], length(data$MCMC_y), nkeep)), start=burnin+1, end=nitt-(nitt-burnin)%%thin, thin=thin)
        }else{
          Liab<-NULL
        }
        if(saveX==FALSE){
          X<-NULL
        }
        if(saveZ==FALSE){
          Z<-NULL
        }

    	options("na.action"=orig.na.action)

        list(Sol=mcmc(Sol, start=burnin+1, end=nitt-(nitt-burnin)%%thin, thin=thin), VCV=mcmc(VCV, start=burnin+1, end=nitt-(nitt-burnin)%%thin, thin=thin), Liab=Liab, Fixed=original.fixed, Random=original.random, Residual=original.rcov, Deviance=deviance,DIC=DIC, X=X, Z=Z)
	
}

