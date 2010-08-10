\name{predict.MCMCglmm}
\alias{predict.MCMCglmm}
\title{Predict method for GLMMs fitted with MCMCglmm}
\description{Predicted values for GLMMs fitted with MCMCglmm}
\usage{
\method{predict}{MCMCglmm}(object, newdata=NULL, marginal=object$Random$formula, type="response", interval="none", level=0.95, \dots)
}
\arguments{
 \item{object}{an object of class \code{"MCMCglmm"}}
 \item{newdata}{An optional data frame in which to look for variables with which to predict}
 \item{marginal}{formula defining random effects to be maginalised}
 \item{type}{character; either "terms" (link scale) or "response" (data scale)}
 \item{interval}{character; either "none", "confidence" or "prediction"}
 \item{level}{A numeric scalar in the interval (0,1) giving the target probability content of the intervals.}
 \item{\dots}{Further arguments to be passed}
}
\value{Expectation and credible interval}  
\author{Jarrod Hadfield j.hadfield@ed.ac.uk}
\seealso{\code{\link{MCMCglmm}}}
\keyword{models}

