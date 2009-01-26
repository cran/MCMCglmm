\name{prunePed}
\alias{prunePed}
\title{Pedigree pruning}
\description{
 Creates a subset of a pedigree by retaining the ancestors of a specified subset 
of individuals
}
\usage{
prunePed(pedigree, keep)
}
\arguments{
  \item{pedigree}{pedigree with id in column 1 dam in column 2 and sire in 
column 3}
  \item{keep}{individuals in pedigree for which the ancestors should be 
retained}
}
\value{
  subsetted pedigree
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\examples{
}
\keyword{manip}


