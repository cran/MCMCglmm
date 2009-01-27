\name{list2bdiag}
\alias{list2bdiag}
\title{Forms the direct product from a list of matrices}
\description{
 Forms a block-diagonal matrix through the direct product from a list of 
matrices
 }
\usage{
list2bdiag(x)
}
\arguments{
  \item{x}{list of square matrices}
}
\details{
}
\value{
    matrix
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{}
\examples{
M<-list(rIW(10, diag(3)), rIW(10, diag(2)))
list2bdiag(M)
}
\keyword{manip}

