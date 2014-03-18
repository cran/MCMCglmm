\name{leg}
\alias{leg}
\title{Evaluates Orthogonal Legendre Polynomials}
\description{Evaluates orthogonal Legendre polynomials of degree 1 to \code{degree} over the specified set of points \code{x}}
\usage{
leg(x, degree, normalized=TRUE)
}
\arguments{
  \item{x}{a numeric vector at which to evaluate the polynomial}
  \item{degree}{the degree of the polynomial. If negative, degree 0 is not evaluated}
  \item{normalized}{logical: shoudlthe polynomial be normalised}
}
\details{
  Please refer to \code{\link[orthopolynom]{legendre.polynomials}} in the package orthopolynom for more detils.
}
\value{
	Matrix of evaluated polynomials
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link[orthopolynom]{legendre.polynomials}}}
\examples{
x<-rnorm(30)
model.matrix(~leg(x,2)-1)
}
\keyword{manip}

