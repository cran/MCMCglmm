test<-function(G, n=1, nG=1){

   rbv<-1:(n*dim(G)[1]*nG)
   
output<-.C("test",
        as.double(rbv),      
        as.double(c(solve(G))),
        as.integer(dim(G)[1]),     
        as.integer(length(c(G))),
        as.integer(n),
        as.integer(nG)

   )       
 
   return(matrix(output[[1]], nG*dim(G)[1],  n))
 }
