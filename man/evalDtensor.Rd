\name{evalDtensor}
\alias{evalDtensor}
\title{Evaluates a list of (mixed) partial derivatives}
\description{
 Evaluates a list of (mixed) partial derivatives
}
\usage{
evalDtensor(x, mu)
}
\arguments{
  \item{x}{unevaluated (list) of expression(s)}
  \item{mu}{values at which the derivatives are evaluated: names need to match 
terms in x}
}
\details{
}
\value{
   tensor 
}
\author{Jarrod Hadfield j.hadfield@ed.ac.uk}
\seealso{Dtensor, D}
\examples{
f<-expression(beta_1 + time*beta_2+u)
Df<-Dtensor(f, eval=FALSE, m=2)
evalDtensor(Df, mu=data.frame(beta_1=0.5, beta_2=1, time=3, u=2.3))
Dtensor(f, mu=c(1,3,1,2.3), m=2)
}
\keyword{math}
\keyword{array}
