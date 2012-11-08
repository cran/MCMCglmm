"prunePed"<-function(pedigree, keep, make.base=FALSE){

   ind.keep<-keep
   nind<-length(ind.keep)+1
   while(length(ind.keep)!=nind){
     nind<-length(ind.keep)
     ind.keep<-union(na.omit(c(unlist(pedigree[,2:3][match(ind.keep,pedigree[,1]),]))), ind.keep)
   }
   pedigree<-pedigree[sort(match(ind.keep, pedigree[,1])),]

   if(make.base){

     if(any(match(pedigree[,2], pedigree[,1])>match(pedigree[,1], pedigree[,1]), na.rm=T)){stop("dams appearing before their offspring: try orderPed form MasterBayes")}
     if(any(match(pedigree[,3], pedigree[,1])>match(pedigree[,1], pedigree[,1]), na.rm=T)){stop("sires appearing before their offspring: try orderPed from MasterBayes")}

     phenotyped<-pedigree[,1]%in%keep
     delete<-rep(FALSE, dim(pedigree)[1])

     for(i in 1:dim(pedigree)[1]){
       nlinks<-phenotyped[i]+sum(pedigree[,2]%in%pedigree[,1][i])+sum(pedigree[,3]%in%pedigree[,1][i])+sum(is.na(pedigree[i,][2:3])==FALSE)    
       if(nlinks<2 & phenotyped[i]==FALSE){                  
         pedigree[,2][which(as.character(pedigree[,2])==as.character(pedigree[,1][i]))]<-NA
         pedigree[,3][which(as.character(pedigree[,3])==as.character(pedigree[,1][i]))]<-NA
         delete[i]<-TRUE                                                            
       }
     } 
     if(any(delete)){
       pedigree<-pedigree[-which(delete),]
     }
   }
   pedigree
   
}
