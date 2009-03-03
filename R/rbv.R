rbv<-function(pedigree=NULL, G, nodes="ALL", scale=TRUE, ggroups=NULL, gmeans=NULL){

  ped=TRUE
  if(length(attr(pedigree, "class"))!=0){
    if(attr(pedigree, "class")=="phylo"){ped=FALSE}
  }

  if(is.matrix(G)==FALSE){G<-as.matrix(G)}
    
  if(ped){

    if(dim(pedigree)[2]!=3){stop("pedigree must have three columns: id, dam and sire")}
    if(sum((na.omit(pedigree[,2])%in%pedigree[,1])==FALSE)>0){stop("individuals appearing as dams but not in pedigree")}
    if(sum((na.omit(pedigree[,3])%in%pedigree[,1])==FALSE)>0){stop("individuals appearing as sire but not in pedigree")}

    node.names<-pedigree[,1]

    id<-match(pedigree[,1], pedigree[,1], nomatch=-998)
    dam<-match(pedigree[,2], pedigree[,1], nomatch=-998)
    sire<-match(pedigree[,3], pedigree[,1], nomatch=-998)

    if(is.null(ggroups)==FALSE){
       if(is.factor(ggroups)==FALSE){stop("ggroups must be a factor")}
       if(dim(pedigree)[1]!=length(ggroups)){stop("ggroups must be the same length as number of individuals in the pedigree")}
       if(any(is.na(ggroups[which(dam==-998 | sire==-998)]))){stop("all individuals without both parents should have a genetic group")}
       if(is.null(gmeans)){stop("genetic groups specified but no gmeans")}else{gmeans<-as.matrix(gmeans)}
       ngroups<-nlevels(ggroups)
       if(ngroups!=dim(gmeans)[1]){stop("number of rows in gmeans should have length nlevels(ggroups)")}
       if(dim(G)[2]!=dim(gmeans)[2]){stop("number of columns in gmeans should have length dim(G)[1]")}
       ggroups<-match(ggroups, levels(ggroups), nomatch=-998)
    }else{
       ngroups<-1
       ggroups<-rep(1, dim(pedigree)[1])
       gmeans<-rep(0, dim(G)[1])
    }

    dnmiss<-which(dam!=-998)   		        # dams not missing
    snmiss<-which(sire!=-998)                   # sires not missing
    bnmiss<-which(dam!=-999 & sire!=-998)       # neither missing

    if(length(intersect(dam[dnmiss], sire[snmiss]))>0 & (length(dnmiss)>0) & (length(snmiss)>0)){stop("dams appearing as sires")}
    if(sum(dam[dnmiss]-id[dnmiss])>=0 & (length(dnmiss)>0)){stop("dams appearing before their offspring: try orderPed from MasterBayes")}
    if(sum(sire[snmiss]-id[snmiss])>=0 & (length(snmiss)>0)){stop("sires appearing before their offspring: try orderPed from MasterBayes")}

    d<-rep(1,length(id))
    
  }else{

    roots<-setdiff(pedigree$edge[,1], pedigree$edge[,2])  
    reorder<-order((pedigree$edge[,1]%in%roots==FALSE)*pedigree$edge[,2]+10000000*(pedigree$edge[,2]<=length(pedigree$tip.label)))
        
    id<-pedigree$edge[,2][reorder]
    tips<-which(id<=length(pedigree$tip.label))

    node.names<-id
    node.names[tips]<-pedigree$tip.label
    node.names[-tips]<-paste("a", 1:(length(id)-length(tips)), sep="")

    dam<-match(pedigree$edge[,1][reorder], id)
    dam[which(is.na(dam))]<--998
    d<-pedigree$edge.length[reorder]
    id<-match(id, id)
    sire<--998

    if(scale){
      root2tip<-0
      ind<-id[length(id)]
      while(ind!=-998){
        root2tip<-root2tip+d[ind]
        ind<-dam[ind]
      }
      d<-d/root2tip
    }
    ngroups<-1
    ggroups<-rep(1, length(id))
    gmeans<-rep(0, dim(G)[1])
  }

 
  
   rbv<-rep(0, dim(G)[1]*length(id))

   output<-.C("rbv",
        as.integer(id-1),      
        as.integer(dam-1),       
        as.integer(sire-1),         
        as.double(d),
        as.double(c(rbv)),     
        as.integer(length(id)),
        as.integer(dim(G)[1]),        
        as.double(c(solve(G))),
        as.integer(ped),
        as.integer(ggroups-1),
        as.double(c(gmeans)),
        as.integer(ngroups)
   )       
 
   rbv<-matrix(output[[5]], length(id),  dim(G)[1])
   rownames(rbv)<-node.names

   if(nodes=="ALL"){
      return(rbv)      
   }
   if(nodes=="TIPS"){
      return(as.matrix(rbv[tips,]))
   }
}
