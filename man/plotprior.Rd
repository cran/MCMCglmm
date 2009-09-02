\name{plotprior}
\alias{plotprior}
\title{Histogram of draws from the prior distribution for the variance structure of an MCMCglmm model}
\description{Samples from the prior distribution of the R-structure and G-structure used in MCMCglmm models, and plots the histograms}
\usage{
plotprior(prior, cor=FALSE, n=1000, ask = dev.interactive(), ...)
}
\arguments{
 \item{prior}{list of prior variance-structure specifications having 2 possible elements: 
the \code{R} (R-structure) and \code{G} (G-structure) specifying the expected (co)variances (\code{V}) and degree 
of belief parameter}
 \item{cor}{logical; if \code{TRUE}, the prior correlation matrices are generated}
 \item{n}{integer: number of samples to be drawn}
 \item{ask}{logical: prompt user before each page of plots}
 \item{...}{Further arguments}
}
\author{Jarrod Hadfield j.hadfield@ed.ac.uk}
\seealso{\code{\link{rIW}}}
\examples{
prior=list(R=list(V=diag(2), n=2),
           G=list(G1=list(V=diag(2), n=3),
                  G2=list(V=diag(2), n=4)))
plotprior(prior, cor=TRUE)
plotprior(prior, cor=FALSE, xlim=c(-1,5))
}
\keyword{hplot}

