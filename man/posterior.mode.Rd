\name{posterior.mode}
\alias{posterior.mode}
\title{Estimates the mode of a vector using kernel density estimation}
\description{Estimates the mode of a vector using kernel density estimation}
\usage{
posterior.mode(x)
}
\arguments{
 \item{x}{numeric vector}
}
\value{
  mode of the kernel density estimate
 }
\author{Jarrod Hadfield j.hadfield@ed.ac.uk}
\seealso{density}
\examples{
v<-rIW(10, as.matrix(1), n=1000)
hist(v)
abline(v=posterior.mode(v), col="red")
}
\keyword{distribution}

