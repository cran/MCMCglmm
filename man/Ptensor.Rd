\name{Ptensor}
\alias{Ptensor}
\title{Tensor of Sample (Mixed) Central Moments}
\description{Forms a tensor of sample (mixed) central moments}
\usage{
Ptensor(x, k)
}
\arguments{
  \item{x}{matrix; traits in columns samples in rows}
  \item{k}{kth central moment}
}
\details{
}
\value{
   \item{}{tensor}
 }
\author{Jarrod Hadfield j.hadfield@ed.ac.uk}
\seealso{}
\examples{
n<-1000
y<-mvrnorm(n, c(0,0), diag(2))
Ptensor(y,2)
cov(y)*((n-1)/n)
}
\keyword{array}

