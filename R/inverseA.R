"inverseA"<-function(pedigree=NULL, nodes="ALL", scale=TRUE){

  ped=TRUE
  if(length(attr(pedigree, "class"))!=0){
    if(attr(pedigree, "class")=="phylo"){ped=FALSE}
  }

  if(ped){

    if(all(is.na(pedigree[,2])) & all(is.na(pedigree[,3]))){stop("All dams and sires are missing")}
    if(dim(pedigree)[2]!=3){stop("pedigree must have three columns: id, dam and sire")}
    if(sum((na.omit(pedigree[,2])%in%pedigree[,1])==FALSE)>0 & any(is.na(pedigree[,2])==FALSE)){stop("individuals appearing as dams but not in pedigree")}
    if(sum((na.omit(pedigree[,3])%in%pedigree[,1])==FALSE)>0 & any(is.na(pedigree[,3])==FALSE)){stop("individuals appearing as sire but not in pedigree")}
    if(sum(duplicated(pedigree[,1]))>0){stop("some individuals appear more than once in the pedigree")}

    numeric.pedigree<-matrix(-998, dim(pedigree)[1], dim(pedigree)[2])

    numeric.pedigree[,1]<-match(pedigree[,1], pedigree[,1], nomatch=-998)
    numeric.pedigree[,2]<-match(pedigree[,2], pedigree[,1], nomatch=-998)
    numeric.pedigree[,3]<-match(pedigree[,3], pedigree[,1], nomatch=-998)
 
    dnmiss<-which(numeric.pedigree[,2]!=-998)   		                 # dams not missing
    snmiss<-which(numeric.pedigree[,3]!=-998)                                    # sires not missing
    bnmiss<-which(numeric.pedigree[,2]!=-998 &  numeric.pedigree[,3]!=-998)      # one present 

    if(length(intersect(numeric.pedigree[,2][dnmiss], numeric.pedigree[,3][snmiss]))>0 & (length(dnmiss)>0) & (length(snmiss)>0)){warning("dams appearing as sires")}
    if(any(numeric.pedigree[,2][dnmiss]>numeric.pedigree[,1][dnmiss]) & (length(dnmiss)>0)){stop("dams appearing before their offspring: try orderPed form MasterBayes")}
    if(any(numeric.pedigree[,3][snmiss]>numeric.pedigree[,1][snmiss]) & (length(snmiss)>0) ){stop("sires appearing before their offspring: try orderPed from MasterBayes")}

    n<-dim(numeric.pedigree)[1]
    nA<-n+2*length(dnmiss)+2*length(snmiss)                                                               # parent-offspring elements
    nA<-nA+2*sum(duplicated(paste(numeric.pedigree[,2], numeric.pedigree[,3])[bnmiss])==FALSE)            # add spouse elements

    inbreeding<-rep(0, n)
    dii<-rep(0, n)
    iA<-rep(0, nA)
    pA<-rep(0, n+1)
    xA<-rep(0, nA)

    Tinv.row<-c(numeric.pedigree[,1][dnmiss], numeric.pedigree[,1][snmiss], 1:n)
    Tinv.col<-c(numeric.pedigree[,2][dnmiss], numeric.pedigree[,3][snmiss], 1:n)
    Tinv.x<-c(rep(-0.5, length(dnmiss)+length(snmiss)),rep(1, n))
    el.order<-order(Tinv.col+Tinv.row/(n+1), decreasing=FALSE)
    Tinv<-Matrix(0,n,n)
    Tinv[1,2]<-1
    Tinv@i<-as.integer(Tinv.row[el.order]-1)
    Tinv@p<-as.integer(c(match(1:n, Tinv.col[el.order]),length(el.order)+1)-1)
    Tinv@x<-as.double(Tinv.x[el.order])

    output<-.C("inverseA",
        as.integer(numeric.pedigree[,1]-1),      
        as.integer(numeric.pedigree[,2]-1),      
        as.integer(numeric.pedigree[,3]-1),             
        as.double(inbreeding),  
        as.double(dii),       
        as.integer(iA),                  
	as.integer(pA), 	         
	as.double(xA), 
        as.integer(n), 	
        as.integer(nA), 
        as.integer(Tinv@i),                  
	as.integer(c(Tinv@p, length(Tinv@x))), 	         
	as.double(Tinv@x), 	
        as.integer(n),                     
	as.integer(length(Tinv@x))
    )

    inbreeding<-output[[4]]
    dii<-output[[5]]
    Ainv<-Matrix(0,n,n)
    Ainv[1,2]<-1
    Ainv@i<-output[[6]][1:output[[10]]]
    Ainv@p<-output[[7]]
    Ainv@x<-output[[8]][1:output[[10]]]

    if(nodes[1]!="ALL"){
       tips<-match(nodes, pedigree[,1])
       if(is.na(tips[1])){stop("some nodes do not appear in pedigree")}
       Ainv<-Ainv[,tips][tips,]-Ainv[,-tips][tips,]%*%as(solve(Ainv[,-tips][-tips,]), "sparseMatrix")%*%Ainv[,tips][-tips,]
       node.names<-pedigree[,1][tips]
    }else{
       node.names<-pedigree[,1]
    }

  }else{
    if(is.rooted(pedigree)==FALSE){
	   stop("phyloegny needs to be rooted")
	}	
    if(is.null(pedigree$edge.length)){
      warning("no branch lengths: compute.brlen from ape has been used")
      pedigree$edge.length<-compute.brlen(pedigree)$edge.length
      scale=TRUE
    }	  
    if(is.ultrametric(pedigree)==FALSE & scale==TRUE){stop("can't scale non-ultrametric trees")}

    roots<-setdiff(pedigree$edge[,1], pedigree$edge[,2])  
    reorder<-order((pedigree$edge[,1]%in%roots==FALSE)*pedigree$edge[,2]+10000000*(pedigree$edge[,2]<=length(pedigree$tip.label)))
        
    id<-pedigree$edge[,2][reorder]
    tips<-which(id<=length(pedigree$tip.label))

    node.names<-id
    node.names[tips]<-pedigree$tip.label
	  
    if(is.null(pedigree$node.label)){
       pedigree<-makeNodeLabel(pedigree) 
    }
	  
    node.names[-tips]<-pedigree$node.label[-1]
		
    inbreeding<-pedigree$edge.length[reorder]
    if(any(inbreeding<1e-16)){stop("some phylogeny edge.lengths are zero")}

    dam<-match(pedigree$edge[,1][reorder], id)
    id<-match(id, id)

    pedigree<-cbind(node.names, node.names[dam], rep(NA, length(node.names)))

    dam[which(is.na(dam))]<--998

    if(scale){
      root2tip<-0
      ind<-id[length(id)]
      while(ind!=-998){
        root2tip<-root2tip+inbreeding[ind]
        ind<-dam[ind]
      }
      inbreeding<-inbreeding/root2tip
    }

    Ainv<-Matrix(0, length(dam), length(dam))
    off<-tapply(id, dam,function(x){x})
    off<-off[-which(names(off)=="-998")]
    off<-lapply(off, function(x){sum(1/inbreeding[x], na.rm=T)})
    diag(Ainv)<-as.numeric(1/inbreeding)
    diag(Ainv)[as.numeric(names(off))]<-diag(Ainv)[as.numeric(names(off))]+unlist(off)
    Ainv[(dam[which(dam!="-998")]-1)*dim(Ainv)[1]+id[which(dam!="-998")]]<-(-1/inbreeding[id[which(dam!="-998")]])
    Ainv[(id[which(dam!="-998")]-1)*dim(Ainv)[1]+dam[which(dam!="-998")]]<-(-1/inbreeding[id[which(dam!="-998")]])
    if(nodes[1]!="ALL"){
      if(nodes[1]=="TIPS"){
        Ainv<-Ainv[,tips][tips,]-Ainv[,-tips][tips,]%*%as(solve(as.matrix(Ainv[,-tips][-tips,])), "sparseMatrix")%*%Ainv[,tips][-tips,]
        node.names<-node.names[tips]
      }
    }

    dii<-inbreeding  # diagonal elements of D equal branch lengths in the case of phylogenies
  }

  return(list(Ainv=Ainv, inbreeding=inbreeding, dii=dii, node.names=node.names, pedigree=pedigree, phylogeny=ped!=TRUE))
}

