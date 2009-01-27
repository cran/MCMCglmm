\name{inverseA}
\alias{inverseA}
\title{Inverse Relatedness Matrix and Phylogenetic Covariance Matrix}
\description{
 Henderson (1976) and Meuwissen and Luo (1992) algorithm for inverting 
relatedness matrices, and Hadfield and Nakagawa (2009) algorithm for inverting 
phylogenetic covariance matrices. 
}
\usage{
  inverseA(pedigree=NULL, nodes="ALL", scale=TRUE)
}
\arguments{
  \item{pedigree}{ordered pedigree with 3 columns: id, dam and sire, or a 
\code{phylo} object.}
  \item{nodes}{\code{"ALL"} calculates the inverse for all individuals/nodes. 
For phylogenies \code{"TIPS"} calculates the inverse for the species tips only, 
and for pedigrees a vector of id's can be passed which inverts the relatedness 
matrix for that subset.}
  \item{scale}{logical: should a phylogeny (needs to be ultrametric) be scaled 
to unit length (distance from root to tip)?}}
\details{
}
\value{
  \item{Ainv}{inverse as \code{sparseMatrix}}
  \item{inbreeding}{inbreeding coefficients/branch lengths}
  \item{pedigree}{pedigree/pedigree representation of phylogeny}
}
\references{
Henderson (1976) Biometrics 32 (1) 69:83
Meuwissen and Luo (1992) 24 (4) 305:313
Hadfield and Nakagawa (2009) submitted
}
\author{Jarrod Hadfield j.hadfield@ed.ac.uk}
\examples{
data(bird.families)
Ainv<-inverseA(bird.families)
}
\keyword{array}



