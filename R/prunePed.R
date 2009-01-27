"prunePed"<-function(pedigree, keep){

   ind.keep<-keep
   nind<-length(ind.keep)+1
   while(length(ind.keep)!=nind){
     nind<-length(ind.keep)
     ind.keep<-union(na.omit(unlist(pedigree[,2:3][match(ind.keep,pedigree[,1]),])), ind.keep)
   }
   pedigree[sort(match(ind.keep, pedigree[,1])),]
}
