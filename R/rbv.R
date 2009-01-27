rbv<-function(pedigree=NULL, G, nodes="ALL", scale=TRUE){

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

    dnmiss<-which(dam!=-998)   		        # dams not missing
    snmiss<-which(sire!=-998)                   # sires not missing
    bnmiss<-which(dam!=-999 & sire!=-998)       # neither missing

    if(length(intersect(dam[dnmiss], sire[snmiss]))>0){stop("dams appearing as sires")}
    if(sum(dam[dnmiss]-id[dnmiss])>=0){stop("dams appearing before their offspring: try order.ped")}
    if(sum(sire[snmiss]-id[snmiss])>=0){stop("sires appearing before their offspring: try order.ped")}

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
        as.logical(ped))      
 
   rbv<-matrix(output[[5]], length(id),  dim(G)[1])
   rownames(rbv)<-node.names

   if(nodes=="ALL"){
      return(rbv)      
   }
   if(nodes=="TIPS"){
      return(as.matrix(rbv[tips,]))
   }
}
