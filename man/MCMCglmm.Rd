\name{MCMCglmm}
\alias{MCMCglmm}
\title{Multivariate Generalised Linear Mixed Models}
\description{
 Markov chain Monte Carlo Sampler for Multivariate Generalised Linear Mixed 
Models with special emphasis on correlated random effects arising from pedigrees 
and phylogenies (Hadfield 2009). Please read the tutorial: \code{vignette("Tutorial", 
"MCMCglmm")}
}
\usage{
MCMCglmm(fixed, random=NULL, rcov=~units, family="gaussian", mev=NULL, 
    data,start=NULL, prior=NULL, tune=NULL, pedigree=NULL,
    nodes="ALL", scale=TRUE, nitt=13000, thin=10, burnin=3000, pr=FALSE,
    pl=FALSE, verbose=TRUE, DIC=TRUE, singular.ok=FALSE, saveX=FALSE, saveZ=FALSE)
}
\arguments{
  \item{fixed}{\code{\link{formula}} for the fixed effects, multiple responses 
are passed as a matrix using cbind}
  \item{random}{\code{\link{formula}} for the random effects. Multiple random terms
can be passed using the \code{+} operator, and in the most general case each random term has the form \code{variance.function(formula):random.term}.
Currently, the only \code{variance.functions} available are \code{idh} and \code{us}. Both fit different variances across each component in \code{formula}, but \code{us} will also fit the covariances. The \code{formula} can contain both factors and numeric terms (i.e. random regression) although it should be noted that the intercept term is suppressed.  The (co)variances are the (co)variances of the \code{random.term} effects. For simple random effects the \code{variance.function(formula)} can be omitted and the model syntax has the simpler form \code{~random1+random2+...}.  There are three reserved variables: \code{units} which index rows of the response variable, \code{trait} which index columns of the response variable and \code{animal} which is associated with a pedigree or phylogeny.} 
  \item{rcov}{\code{\link{formula}} for residual covariance structure. This has 
to be set up so that each data point is associated with a unique residual. For 
example a multi-response model might have the R-structure defined by 
\code{~us(trait):units}}
  \item{family}{optional character vector of trait distributions. Currently, 
\code{"gaussian"}, \code{"poisson"},  \code{"categorical"},  
\code{"multinomial"},  \code{"ordinal"}, \code{"ordinal"}, \code{"exponential"},  \code{"cengaussian"},  
\code{"cenpoisson"},  \code{"cenexponential"},  \code{"zipoisson"} are 
supported, where the prefix \code{"cen"} means censored and the prefix 
\code{"zi"} means zero inflated. If \code{NULL} \code{data} needs to contain a 
\code{family} column.}
  \item{mev}{optional vector of measurement error variances for each data point 
for random effect meta-analysis.}
  \item{data}{\code{data.frame}}
  \item{start}{optional list having 4 possible elements: 
\code{R} (R-structure) \code{G} (G-structure) and \code{liab} (latent variables 
or liabilities) should contain the starting values where \code{G} itself is also a list with as many elements as random 
effect components. The fourth element \code{QUASI} should be logical: if \code{TRUE} starting latent variables are obtained heuristically, if \code{FALSE} then they are sampled from a Z-distribution}
  \item{prior}{optional list of prior specifications having 3 possible elements: 
\code{R} (R-structure) \code{G} (G-structure) and \code{B} (fixed effects). Each 
element is a list containing the expected (co)variances (\code{V}) and a degree 
of belief parameter (\code{nu}) for the variance structures (\code{R} and 
\code{G}). \code{B} is a list containing the expected value (\code{mu}) and a 
(co)variance matrix (\code{V}) representing the strength of belief. Defaults to 
\code{n=0} and \code{B$V<-I*1e+8}, where I is an identity matrix of appropriate 
dimension.}
  \item{tune}{optional (co)variance matrix defining the proposal distribution 
for the latent variables. If NULL an adaptive algorithm is used which ceases to 
adapt once the burn-in phase has finished.}
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
  \item{DIC}{logical: if \code{TRUE} deviance and deviance information criterion are calculated}
  \item{singular.ok}{logical: if \code{FALSE} linear dependencies in the fixed effects are removed. if \code{TRUE} they are left in an estimated, although all information comes form the prior}
  \item{saveX}{logical: save fixed effect design matrix}
  \item{saveZ}{logical: save random effect design matrix}

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
  \item{Deviance}{deviance -2*log(p(y|...))}
  \item{DIC}{deviance information criterion}
  \item{X}{sparse fixed effect design matrix}
  \item{Z}{sparse random effect design matrix}
}
\references{
	General analyses: Hadfield, J.D. (2009) \emph{submitted}
	
	Phylogenetic analyses: Hadfield, J.D.  and Nakagawa, S. (2009) \emph{submitted}
	
	Background Sorensen, D. \& Gianola, D. (2002) Springer
}
\author{Jarrod Hadfield j.hadfield@ed.ac.uk}
\seealso{\code{\link[coda]{mcmc}}}
\examples{

# Example 1: univariate Gaussian model with standard random effect
 
data(PlodiaPO)  
model1<-MCMCglmm(PO~1, random=~FSfamily, data=PlodiaPO, verbose=FALSE)
plot(model1$VCV)

# Example 2: univariate Gaussian model with phylogenetically correlated
# random effect

data(bird.families) 

phylo.effect<-rbv(bird.families, 1, nodes="TIPS") 
phenotype<-phylo.effect+rnorm(dim(phylo.effect)[1], 0, 1)  

# simulate phylogenetic and residual effects with unit variance

test.data<-data.frame(phenotype=phenotype, animal=row.names(phenotype))

prior<-list(R=list(V=1, nu=1), G=list(G1=list(V=1, nu=1)))

model2<-MCMCglmm(phenotype~1, random=~animal, data=test.data, 
pedigree=bird.families, prior=prior, verbose=FALSE)

plot(model2$VCV)

}
\keyword{models}


