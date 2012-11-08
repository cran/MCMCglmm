\name{gelman.prior}
\alias{gelman.prior}
\title{Prior Covariance Matrix for Fixed Effects.}
\description{Prior Covariance Matrix for Fixed Effects.}
\usage{
gelman.prior(formula, data, scale=1, intercept=scale)
}
\arguments{
 \item{formula}{\code{\link{formula}} for the fixed effects.}
 \item{data}{\code{\link{data.frame}}.}
 \item{intercept}{prior variance for intercept}
 \item{scale}{prior variance for regression parameters}
 }
\value{
  prior covariance matrix
 }
\references{
        Gelman, A. et al. (2008) The Annals of Appled Statistics 2 4  1360-1383
}
\author{Jarrod Hadfield \email{jarrod.hadfield@zoo.ox.ac.uk}}
\details{
 Gelman et al. (2008) suggest that the input variables of a categorical regression are standardised and that the associated regression parameters are assumed independent in the prior.  Gelman et al. (2008) recommend a scaled t-distribution with a single degree of freedom (scaled Cauchy) and a scale of 10 for the intercept and 2.5 for the regression parameters. If the degree of freedom is infinity (i.e. a normal distribution) then a prior covariance matrix \code{B$V} can be defined for the regression parameters without input standardisation that corresponds to a diagonal prior \eqn{{\bf D}}{D} for the regression parameters had the inputs been standardised.  The diagonal elements of \eqn{{\bf D}}{D} are set to \code{scale} except the first which is set to \code{intercept}.  With logistic regression \eqn{\pi^{2}/3+\sigma^{2}}{pi^2/3+v} gives a prior that is approximately flat on the probability scale, where \eqn{\sigma^{2}}{v} is the total variance due to the random effects. For probit regression it is \eqn{1+\sigma^{2}}{1+v}.
}
\examples{
dat<-data.frame(y=c(0,0,1,1), x=gl(2,2))
# data with complete separation

prior<-list(B=list(mu=c(0,0),V=gelman.prior(~x, data=dat, scale=1+pi^2/3)), R=list(V=1,fix=1))

m1<-MCMCglmm(y~x, prior=prior, data=dat, family="ordinal", verbose=FALSE)

c2 <- (16 * sqrt(3)/(15 * pi))^2

plot(plogis(m1$Sol[,1]/sqrt(1+c2))) # marginal probability when x=1

}
\keyword{distribution}

