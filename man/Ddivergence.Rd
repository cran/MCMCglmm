\name{Ddivergence}
\alias{Ddivergence}
\title{d-divergence}
\description{Calculates Ovaskainen's (2008) d-divergence between 2 zero-mean multivariate normal distributions.}
\usage{
Ddivergence(CA=NULL, CB=NULL, n=10000)
}
\arguments{
  \item{CA}{Matrix A}
  \item{CB}{Matrix B}
  \item{n}{number of Monte Carlo samples for approximating the integral}
}
\details{
 }
\value{
  \item{}{d-divergence}
}
\references{Ovaskainen et. al. (2008) Proc. Roy. Soc - B (275) 1635 593-750}
\author{Jarrod Hadfield j.hadfield@ed.ac.uk}
\note{}
\seealso{}
\examples{
CA<-rIW(10, diag(2), n=1)
CB<-rIW(10, diag(2), n=1)
Ddivergence(CA, CB)
}
\keyword{multivariate}
