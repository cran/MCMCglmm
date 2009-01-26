\name{MCMCglmm}
\alias{MCMCglmm}
\title{Multivariate Generalised Linear Mixed Models}
\description{
 Markov chain Monte Carlo Sampler for Multivariate Generalised Linear Mixed 
Models with special emphasis on correlated random effects arising from pedigrees 
and phylogenies. Please read the tutorial: \code{vignette("Tutorial", 
"MCMCglmm")}
}
\usage{
MCMCglmm(fixed, random=NULL, rcov=~units, family="gaussian", mev=NULL, 
data=NULL,
      start=NULL, prior=NULL, tune=NULL, pedigree=NULL, nodes="ALL", scale=TRUE,
      nitt=13000, thin=10, burnin=3000, pr=FALSE, pl=FALSE, verbose=TRUE)
}
\arguments{
  \item{fixed}{\code{\link{formula}} for the fixed effects, multiple responses 
are passed as a matrix using cbind}
  \item{random}{\code{\link{formula}} for the random effects. There are three 
reserved variables: \code{units} which indexes rows of the response variable, 
\code{trait} which indexes columns of the response variable and \code{animal} 
which is associated with a pedigree or phylogeny. \code{~rfactor} fits a single 
variance across all levels of \code{rfactor}, \code{idh(factor):rfactor} fits 
different \code{rfactor} variances for each level of \code{factor}, and 
\code{us(factor):rfactor} fits different variances and covariances terms. Random 
regression using Legendre polynomials is possible using \code{leg(covariate, 
n):rfactor} where \code{n} is the order of the polynomial. \code{-n} fits an n-
order polynomial but without the intercept.}
  \item{rcov}{\code{\link{formula}} for residual covariance structure. This has 
to be set up so that each data point is associated with a unique residual. For 
example a multi-trait model might have the R-structure defined by 
\code{~us(trait):units}}
  \item{family}{optional character vector of trait distributions. Currently, 
\code{"gaussian"}, \code{"poisson"},  \code{"categorical"},  
\code{"multinomial"},  \code{"exponential"},  \code{"cengaussian"},  
\code{"cenpoisson"},  \code{"cenexponential"},  \code{"zipoisson"} are 
supported, where the prefix \code{"cen"} means censored and the prefix 
\code{"zi"} means zero inflated. If \code{NULL} \code{data} needs to contain a 
\code{family} column.}
  \item{mev}{optional vector of measurement error variances for each data point 
for random effect meta-analysis.}
  \item{data}{\code{data.frame}}
  \item{start}{optional list of starting parameters having 3 possible elements: 
\code{R} (R-structure) \code{G} (G-structure) and \code{liab} (latent variables 
or liabilities). \code{G} itself is also a list with as many elements as random 
effect components}
  \item{prior}{optional list of prior specifications having 3 possible elements: 
\code{R} (R-structure) \code{G} (G-structure) and \code{B} (fixed effects). Each 
element is a list containing the expected (co)variances (\code{V}) and a degree 
of belief parameter (\code{n}) for the variance structures (\code{R} and 
\code{G}). \code{B} is a list containing the expected value (\code{mu}) and a 
(co)variance matrix (\code{V}) representing the strength of belief. Defaults to 
\code{n=0} and \code{B$V<-I*1e+8}, where I is an identity matrix of appropriate 
dimension.}
  \item{tune}{optional (co)variance matrix defining the proposal distribution 
for the latent variables. If NULL an adaptive algorithm is used which ceases to 
adapt once the burnin phase has finished.}
  \item{pedigree}{ordered pedigree with 3 columns id, dam and sire or a 
\code{phylo} object.}
  \item{nodes}{pedigree/phylogeny nodes to be estimated. The default, 
\code{"ALL"} estimates effects for all individuals in a pedigree or nodes in a 
phylogeny (including ancestral nodes). For phylogenies \code{"TIPS"} estimates 
effects for the tips only, and for pedigrees a vector of ids can be passed to 
\code{nodes} specifying the subset of individuals for which animal effects are 
estimated.  Note that all analyses are equivalent if omitted nodes have missing 
data but by absorbing these nodes the chain max mix better. However, the 
algorithm may be less numerically stable and may iterate slower, especially for 
large phylogenies.}
  \item{scale}{logical: should the phylogeny (needs to be ultrametric) be scaled 
to unit length (distance from root to tip)?}
  \item{nitt}{number of MCMC iterations}
  \item{thin}{thinning interval}
  \item{burnin}{burnin}
  \item{pr}{logical: should the posterior distribution of random effects be 
saved?}
  \item{pl}{logical: should the posterior distribution of latent variables be 
saved?}
  \item{verbose}{logical: if \code{TRUE} MH diagnostics are printed to screen}
}
\details{
}
\value{
  \item{Sol}{Posterior Distribution of MME solutions, including fixed effects}
  \item{VCV}{Posterior Distribution of (co)variance matrices}
  \item{Liab}{Posterior Distribution of latent variables}
  \item{Fixed}{formula: fixed terms}
  \item{Random}{formula: random terms}
  \item{Residual}{formula: residual terms}
}
\references{
For phylogenetic analyses: Hadfield & Nakagawa (2009) \emph{submitted}
For other analyses: Hadfield \emph{in prep}
Sorensen \& Gianola (2004) Springer
}
\author{Jarrod Hadfield j.hadfield@ed.ac.uk}
\seealso{mcmc}
\examples{

# Example 1: univariate gaussian model with standard random effect
 
data(PlodiaPO)  
model1<-MCMCglmm(PO~1, random=~FSfamily, data=PlodiaPO, verbose=FALSE)
plot(model1$VCV)

# Example 2: univariate gaussian model with phylogenetically correlated random effect

data(bird.families) 

phylo.effect<-rbv(bird.families, 1, nodes="TIPS") # simulate phylogenetic effects with unit variance
phenotype<-phylo.effect+rnorm(dim(phylo.effect)[1], 0, 1)  # add residual with unit variance

test.data<-data.frame(phenotype=phenotype, animal=as.factor(row.names(phenotype)))

model2<-MCMCglmm(phenotype~1, random=~animal, data=test.data, pedigree=bird.families, verbose=FALSE)
plot(model2$VCV)

}
\keyword{models}


