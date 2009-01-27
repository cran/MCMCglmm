\name{rIW}
\alias{rIW}
\title{Random Generation from the Conditional Inverse Wishart Distribution}
\description{Samples from the inverse Wishart distribution, with the possibility 
of conditioning on a diagonal submatrix}
\usage{
rIW(nu,V,split=NULL, n=1)
}
\arguments{
 \item{nu}{degrees of freedom}
 \item{V}{inverse scale matrix}
 \item{split}{optional integer indexing the partition to be conditioned on}
 \item{n}{integer: number of samples to be drawn}
}
\details{
 \code{split} indexes the diagonal element of \code{V} which partitions \code{V} 
into 4 submatrices. \code{split} indexes the upper left corner of the lower 
diagonal matrix and it is this matrix that is conditioned on.
 
 For example partioning \code{V} such that 

 \deqn{
 {\bf V} = \left[
             \begin{array}{cc}
		       {\bf V}_{11}&{\bf V}_{12}\\
			   {\bf V}_{21}&{\bf V}_{12}\\
		     \end{array}
		   \right]
  }{V = V11 V12}
  \deqn{}{        V21 V22}

the split indexes the upper left corner of \eqn{{\bf V}_{22}}{V22}. If 
\code{dim(V)}=4 and \code{split=2} then \eqn{{\bf V}_{11}}{V11} is a 1X1 matrix 
and V22 is a 3X3 matrix.
}
\value{
  if \code{n} = 1 a matrix equal in dimension to \code{V}, if \code{n}>1 a 
matrix of dimension \code{n} x \code{length(V)}
}
\author{Jarrod Hadfield j.hadfield@ed.ac.uk}
\references{Korsgaard IR et. al. 1999 Genetics Selection Evolution 31 (2) 
177:181}
\seealso{rwishart}
\examples{
n<-10
V<-diag(4)
rIW(n, solve(n*V), split=2)
}
\keyword{distribution}

