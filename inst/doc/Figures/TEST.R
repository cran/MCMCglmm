#source("~/Desktop/MCMCglmmTEST.R")
#source("~/Desktop/recovery//MCMCglmm_2.05/inst/doc/Figures/TEST.R")
library(VGAM)
library(MASS)
library(MCMCglmm)
verbose=FALSE
plotit=TRUE
leg=TRUE
nsim<-1
nitt<-130
thin<-1
burnin<-30

# poisson test 1.2 seconds OK
print("res1")
R<-diag(1)
res1<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1))
tpar<-c(1, 1)
for(i in 1:nsim){
l<-exp(rnorm(100,1,R))
y<-rpois(100,l)
data=data.frame(y1=y, y2=y)
m1<-MCMCglmm(cbind(y1)~1, family="poisson", data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res1 different from expected"))
}
res1[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# multinomial test J=1 test 0.9 seconds  OK
print("res2")
res2<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1))
tpar<-c(1,1)
for(i in 1:nsim){
y<-rbinom(100,10,plogis(rnorm(100,1,1)))
data=data.frame(y1=y, y2=10-y)
m1<-MCMCglmm(cbind(y1,y2)~1, family="multinomial2", data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res2 different from expected"))
}
res2[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# categorical test J=1 test 0.9 seconds OK
print("res3")
res3<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
tapr<-c(1,1)
for(i in 1:nsim){
y<-rbinom(100,1,plogis(rnorm(100,1,1)))
data=data.frame(y1=y, y2=1-y)
prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
m1<-MCMCglmm(y1~1, family="categorical", data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res3 different from expected"))
}
res3[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# categorical test J=1 test with slice sampling 0.9 seconds OK
print("res3b")
res3b<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
tpar<-c(1,1)
for(i in 1:nsim){
y<-rbinom(100,1,plogis(rnorm(100,1,1)))
data=data.frame(y1=y, y2=1-y)
prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
m1<-MCMCglmm(y1~1, family="categorical", data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin, slice=TRUE)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res3b different from expected"))
}
res3b[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


# gauss with blocked random  1.5 seconds OK
print("res4")
res4<-matrix(NA, nsim,3)
R<-as.matrix(2)
G<-as.matrix(1)
prior<-list(R=list(V=as.matrix(1), n=1), G=list(G1=list(V=G, n=1)))
tpar<-c(-1,1,2)
for(i in 1:nsim){
fac<-as.factor(sample(1:50,300,replace=TRUE))
ffac<-gl(2,150)
y<-mvrnorm(300, c(-1), R)+mvrnorm(50, c(0), G)[fac]
data=data.frame(y1=y, fac=fac, ffac=ffac)
m1<-MCMCglmm(y1~1,random=~fac,  data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res4 different from expected"))
}
res4[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# binary with blocked random  1.5 seconds OK

print("res4c")
res4c<-matrix(NA, nsim,3)
R<-as.matrix(2)
G<-as.matrix(1)
tpar<-c(-1,1,2)
prior=list(R=list(V=R, n=1, fix=1),G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-mvrnorm(300, c(-1), R)+mvrnorm(75, c(0), G)[fac]
data=data.frame(y1=rbinom(300, 1,plogis(y)), fac=fac, y2=y)
m1<-MCMCglmm(y1~1,random=~fac,  data=data, prior=prior, family="categorical",verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res4c different from expected"))
}
res4c[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# gauss with correlated random 9.6 seconds OK

print("res5")
res5<-matrix(NA, nsim,3)
R<-as.matrix(2)
G<-as.matrix(1)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=1)))
tpar<-c(-1,1,2)
for(i in 1:nsim){
Ped<-cbind(1:400, c(rep(NA,100), sample(1:50,300,TRUE)),c(rep(NA,100), sample(51:100,300,TRUE)))
y<-mvrnorm(300, c(-1), R)+rbv(Ped,G)[101:400]
data=data.frame(y1=y, animal=as.factor(Ped[,1][101:400]), fac=gl(2,150))
system.time(m1<-MCMCglmm(y1~1, random=~animal, pedigree=Ped, data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin))
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res5 different from expected"))
}
res5[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


if(leg){
# gauss with correlated random regression OK
print("res5b")
res5b<-matrix(NA, nsim,6)
R<-as.matrix(2)
G<-diag(2)
prior=list(R=list(V=R, nu=1),G=list(G1=list(V=G, nu=2)))
tpar<-c(-1, 1, 0,0,1,2)
for(i in 1:nsim){
Ped<-cbind(1:400, c(rep(NA,100), sample(1:50,300,TRUE)),c(rep(NA,100), sample(51:100,300,TRUE)))
x<-runif(300, -1,1)
bv<-rbv(Ped,G)[101:400,]
y<-mvrnorm(300, c(-1), R)+bv[,1]+bv[,2]*x
data=data.frame(y1=y, animal=as.factor(Ped[,1][101:400]), x=x)
m1<-MCMCglmm(y1~1, random=~us(leg(x,1,FALSE)):animal, pedigree=Ped, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res5b different from expected"))
}
res5b[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
}

# gauss with us random effect 3 seconds
print("res6")
res6<-matrix(NA, nsim,6)
G=matrix(c(1,0.5,0.5,2),2,2)
R<-as.matrix(2)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=1)))
tpar<-c(-1, 1, 0.5, 0.5, 2, 2)
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
facf<-as.factor(sample(1:2,300,replace=TRUE))
rear<-sample(1:300)
fac<-gl(75,4)[rear]
facf<-gl(2,2,300)[rear]
y<-mvrnorm(300, c(-1), R)+rowSums(mvrnorm(75, c(0,0), G)[fac,]*cbind(facf==1,facf==2))
data=data.frame(y1=y, fac=fac, facf=facf)
rmodel.terms<-"fac:facf"
m1<-MCMCglmm(y1~1,random=~us(facf):fac,  data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res6 different from expected"))
}
res6[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


# gauss with idh random effect 1.7 seconds
print("res7")
res7<-matrix(NA,nsim,4)
res7b<-matrix(NA,nsim,6)
R<-as.matrix(2)
G=matrix(c(1,0,0,2),2,2)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=2)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
facf<-as.factor(sample(1:2,300,replace=TRUE))
y<-mvrnorm(300, c(-1), R)+rowSums(mvrnorm(75, c(0,0), G)[fac,]*cbind(facf==1,facf==2))
data=data.frame(y1=y, fac=fac, facf=facf)
print(table(table(fac,facf)[,1]>0, table(fac,facf)[,2]>0))
m1<-MCMCglmm(y1~1,random=~idh(facf):fac,  data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
tpar<-c(-1, 1, 2, 2)
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res7 different from expected"))
}
m2<-MCMCglmm(y1~1,random=~us(facf):fac,  data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m2$Sol, m2$VCV)), ask=FALSE)
}
tpar<-c(-1, 1,0,0, 2, 2)
if(any(HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,2]<tpar)){
print("res7b different from expected")
}
res7[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
res7b[i,]<-posterior.mode(mcmc(cbind(m2$Sol, m2$VCV)))
print(i)
}


# bivariate gauss 1.9 seconds
print("res8")
res8<-matrix(NA, nsim,6)
R=matrix(c(1,0.5,0.5,2),2,2)
prior=list(R=list(V=R, n=1))
tpar<-c(-1,1,1,0.5,0.5,2)
for(i in 1:nsim){
y<-mvrnorm(300, c(-1,1), R)
data=data.frame(y1=y[,1], y2=y[,2])
data$y1[sample(1:300, 30)]<-NA
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","gaussian"), rcov=~us(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res8 different from expected"))
}
res8[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


# bivariate gauss idh residual 1.7 seconds
print("res9")
res9<-matrix(NA, nsim,4)
R=matrix(c(1,0,0,2),2,2)
prior=list(R=list(V=R, n=2))
tpar<-c(-1,1,1,2)
for(i in 1:nsim){
y<-mvrnorm(300, c(-1,1), R)
data=data.frame(y1=y[,1], y2=y[,2])
data$y1[sample(1:300, 50)]<-NA
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","gaussian"), rcov=~idh(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res9 different from expected"))
}
res9[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# bivariate gauss + random 2.8 seconds
print("res10")
res10<-matrix(NA, nsim,6)
R=matrix(c(1,0,0,2),2,2)
G=matrix(c(2,0,0,1),2,2)
tpar<-c(-1,1,2,1,1,2)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-mvrnorm(300, c(-1,1), R)+mvrnorm(75, c(0,0), G)[fac,]
data=data.frame(y1=y[,1], y2=y[,2], fac=fac)
m1<-MCMCglmm(cbind(y1,y2)~trait-1, random=~idh(trait):fac, family=c("gaussian","gaussian"), rcov=~idh(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res10 different from expected"))
}
res10[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# bivariate binoimal + random 5.9 seconds
res11<-matrix(NA, nsim,8)
R=matrix(c(1,0,0,2),2,2)
G=matrix(c(2,0.5,0.5,1),2,2)
tpar<-c(-1,1,2,0.5,0.5,1,1,2)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-mvrnorm(300, c(-1,1), R)+mvrnorm(75, c(0,0), G)[fac,]
y1<-rbinom(300,10,plogis(y[,1]))
y2<-rbinom(300,10,plogis(y[,2]))
data=data.frame(y1s=y1,y1f=10-y1,y2s=y2,y2f=10-y2, fac=fac)
m1<-MCMCglmm(cbind(y1s,y1f,y2s,y2f)~trait-1, random=~us(trait):fac, family=c("multinomial2","multinomial2"), rcov=~idh(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res11 different from expected"))
}
res11[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


# gauss with blocked random - all missing 1.9 seconds
print("res12")
res12<-matrix(NA, nsim,3)
R<-as.matrix(2)
G<-as.matrix(1)
prior=list(R=list(V=R, n=100),G=list(G1=list(V=G, n=100)), B=list(mu=0, V=0.000000001))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-rep(NA,300)
tpar<-c(0,1,2)
data=data.frame(y1=y, fac=fac)
m1<-MCMCglmm(y1~1,random=~fac,  data=data, prior=prior, singular.ok=TRUE,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)

if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res12 different from expected"))
}

res12[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# Jennys data - once as gaussian once as binomial 32 seconds

if(file.exists("~/Work/Jenny/Data/Intermediate/ThirdC.R")){

  m1R<-dget("~/Work/Jenny/Data/Intermediate/ThirdC.R")
  res13<-matrix(NA, nsim,39)
  res14<-matrix(NA, nsim,39)
  print("res13")
  print("res14")

  firstP<-read.table("~/Work/Jenny/Data/Raw/Third_paternal.txt", header=T)
  firstP$day<-as.factor(firstP$day)

  coef<-apply(cbind(m1R$Sol, m1R$VCV), 2, median)
  G<-list(G1=matrix(coef[9:33],5,5), G2=as.matrix(coef[34]))


  R<-diag(coef[35:39])
  prior=list(R=list(V=R, n=5),G=list(G1=list(V=G[[1]], n=1), G2=list(V=G[[2]], n=1)))


  for(i in 1:nsim){
    location<-coef[firstP$virus]+c(0,coef[1])[((firstP$virus!=levels(firstP$virus)[1]))+1]+c(0,coef[6:8])[firstP$day]
    resid<-rnorm(dim(firstP)[1], 0, sqrt(diag(R)[firstP$virus]))
    reffects<-rowSums(mvrnorm(nlevels(firstP$line),rep(0,5),G[[1]])[as.numeric(firstP$line),]*outer(firstP$virus, levels(firstP$virus), "=="))
    reffects<-reffects+rnorm(nlevels(firstP$f2rep),0,sqrt(G[[2]]))[as.numeric(firstP$f2rep)]
    l<-location+resid+reffects
    firstP$l<-l
    y<-rbinom(dim(firstP)[1], firstP$total, plogis(l))
    firstP$noalive<-y
    firstP$nodead<-firstP$total-y
    tpar<-coef
    m1test<-MCMCglmm(l~virus+day, random=~us(virus):line+f2rep, rcov=~idh(virus):units, family=c("gaussian"), data=firstP, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1test$Sol, m1test$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1test$Sol, m1test$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1test$Sol, m1test$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1test$Sol, m1test$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1test$Sol, m1test$VCV)))[,2]<tpar)/length(tpar), "res13 different from expected"))
}


    m1test2<-MCMCglmm(cbind(noalive, nodead)~virus+day, random=~us(virus):line+f2rep, rcov=~idh(virus):units, family=c("multinomial2"), data=firstP, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1test2$Sol, m1test2$VCV)), ask=FALSE)
}

if(any(HPDinterval(mcmc(cbind(m1test2$Sol, m1test2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1test2$Sol, m1test2$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1test2$Sol, m1test2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1test2$Sol, m1test2$VCV)))[,2]<tpar)/length(tpar), "res14 different from expected"))
}

    res13[i,]<-posterior.mode(mcmc(cbind(m1test$Sol, m1test$VCV)))
    res14[i,]<-posterior.mode(mcmc(cbind(m1test2$Sol, m1test2$VCV)))
    print(i)
  }
if(plotit){
  plot(apply(res14, 2, median)~coef)
  points(coef, apply(res13, 2, median), col="red")
}
}else{
  print("file ~/Work/Jenny/Data/Intermediate/ThirdC.R does not exist: skipping res 13 and 14")
}

res15<-matrix(NA, nsim,3)
print("res15")
R<-as.matrix(2)
prior=list(R=list(V=R, n=1))
for(i in 1:nsim){
mev<-rchisq(300,10)
y<-mvrnorm(300, c(-1), R)+rnorm(300, c(0), sqrt(mev))
data=data.frame(y1=y)
tpar<-c(-1,1,2)
m1<-MCMCglmm(y1~1,data=data, prior=prior, mev=mev,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res15 different from expected"))
}

res15[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, median)
print(i)
}



# censored gaussian data

res17a<-matrix(NA, nsim,2)
res17b<-matrix(NA, nsim,2)
print("res17")
for(i in 1:nsim){
l<-rnorm(500,0,1)
y<-cut(l,c(-Inf,seq(-4,4,1), Inf), , include.lowest=TRUE)
y1<-c(c(-Inf,seq(-4,4,1), Inf),Inf)[y]
y2<-c(c(-Inf,seq(-4,4,1), Inf),Inf)[as.numeric(y)+1]
ym<-(y1+y2)/2
data=data.frame(y1=y1, y2=y2, ym=ym)
if(any(y1==-Inf) | any(y2==Inf)){
  m1$VCV<-rep(-1,1000)
  m1$Sol<-rep(-1,1000)
}else{
  m1<-MCMCglmm(ym~1, data=data, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
}
  m2<-MCMCglmm(cbind(y1, y2)~1, data=data, family="cengaussian",verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
  if(plotit){
    plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
  }
  if(plotit){
    plot(mcmc(cbind(m2$Sol, m2$VCV)), ask=FALSE)
  }
tpar<-c(0,1)
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res17a different from expected"))
}
if(any(HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,2]<tpar)/length(tpar), "res17b different from expected"))
}
res17a[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
res17b[i,]<-posterior.mode(mcmc(cbind(m2$Sol, m1$VCV)))
print(i)
}

# censored Possion data
res18<-matrix(NA, nsim,2)
print("res18")
tpar<-c(0,1)
prior=list(R=list(V=diag(1), n=1))
for(i in 1:nsim){
l<-rpois(100, exp(rnorm(100,0,1)))
y<-cut(l,c(seq(0,10,2), Inf), include.lowest=TRUE, right=FALSE)
y1<-c(c(seq(0,10,2), Inf),Inf)[y]
y2<-c(c(seq(0,10,2), Inf),Inf)[as.numeric(y)+1]-1
y1[100]<-l[100]
y2[100]<-l[100]
data=data.frame(y1=y1, y2=y2)
m1<-MCMCglmm(cbind(y1, y2)~1, data=data, family="cenpoisson", prior=prior, pl=TRUE,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res18 different from expected"))
}
res18[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

print("res19")
res19<-matrix(0, nsim, 4)
R<-diag(2)
prior=list(R=list(V=R, n=1, fix=2), B=list(mu=c(0,0), V=matrix(c(1000,0,0,pi^2/3),2,2)))
tune=list(diag(2))
tpar<-c(1,-1,1,0,0,1)
for(i in 1:nsim){
l<-mvrnorm(300,c(0,-1), R)
y<-rzipois(300, exp(1+l[,1]), plogis(l[,2]))
data=data.frame(y1=y)
m1<-MCMCglmm(y1~trait-1, rcov=~idh(trait):units, data=data, family="zipoisson",prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res19 different from expected"))
}
res19[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

print("res19b")
res19b<-matrix(0, nsim, 6)
R<-diag(2)
prior=list(R=list(V=R, n=1, fix=2),G=list(G1=list(V=1, n=1)))
tune=list(diag(2))
tpar<-c(1, -0.5, 0.2, 1,1,0,0,1)
for(i in 1:nsim){
	x<-rnorm(300)
	l<-mvrnorm(300,c(0,-0.5), R)
        fac<-gl(50,6)
        r<-rnorm(50)
	y<-rzipois(300, exp(1+l[,1]+r[fac]+0.2*x), plogis(l[,2]))
	data=data.frame(y1=y, x=x, fac=fac)
	m1<-MCMCglmm(y1~trait+at.level(trait, 1):x-1, random=~us(at.level(trait,1)):fac, rcov=~idh(trait):units, data=data, family="zipoisson",prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res19b different from expected"))
}

	res19b[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
	print(i)
}

# bivariate gauss + categorical  residual 1.7 seconds

res20<-matrix(NA, nsim,6)
print("res20")
tune<-diag(2)
R=matrix(c(2,0.25,0.25,1),2,2)
prior=list(R=list(V=R, n=2, fix=2))
tpar<-c(-1,1,2,0.25, 0.25,1)
for(i in 1:nsim){
y<-mvrnorm(300, c(-1,1), R)
data=data.frame(y1=y[,1], y2=rbinom(300, 1, plogis(y[,2])), y3=y[,2])
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","categorical"), rcov=~us(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res20 different from expected"))
}
res20[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# gauss random regression

if(leg){
res21<-matrix(NA, nsim,5)
res21b<-matrix(NA, nsim,5)
res21c<-matrix(NA, nsim,7)
print("res21")
G=matrix(c(2,0,0,1),2,2)
R=matrix(1,1,1)
prior=list(R=list(V=R, n=1), G=list(G1=list(V=G[1,1,drop=FALSE], n=1), G2=list(V=G[2,2,drop=FALSE], n=0.1)))
prior2=list(R=list(V=R, n=1), G=list(G1=list(V=G, n=1)))

int.slope<-mvrnorm(300, c(0,0), G)
ind<-gl(300,3)
time<-rnorm(900)
y<-int.slope[,1][ind]+time*int.slope[,2][ind]
y<-y+rnorm(900,-1,R)
data=data.frame(y1=y, time=time, ind=ind)

for(i in 1:nsim){
int.slope<-mvrnorm(300, c(0,0), G)
ind<-gl(300,3)
time<-rnorm(900)
y<-int.slope[,1][ind]+time*int.slope[,2][ind]
y<-y+rnorm(900,-1,R)
data=data.frame(y1=y, time=time, ind=ind)
	m1<-MCMCglmm(y1~time, random=~us(leg(time,0,FALSE)):ind+us(leg(time,-1,FALSE)):ind, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
	m2<-MCMCglmm(y1~time, random=~idh(leg(time,1,FALSE)):ind, data=data, prior=prior2,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
        m3<-MCMCglmm(y1~time, random=~us(1+poly(time,1, raw=TRUE)):ind, data=data, prior=prior2,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
tpar<-c(-1, 0, 2,1,1)
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res21-m1 different from expected"))
}
if(any(HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,2]<tpar)/length(tpar), "res21-m2 different from expected"))
}
tpar<-c(-1, 0, 2,0,0,1,1)
if(any(HPDinterval(mcmc(cbind(m3$Sol, m3$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m3$Sol, m3$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m3$Sol, m3$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m3$Sol, m3$VCV)))[,2]<tpar)/length(tpar), "res21-m3 different from expected"))
}
res21[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
res21b[i,]<-posterior.mode(mcmc(cbind(m2$Sol, m2$VCV)))
res21c[i,]<-posterior.mode(mcmc(cbind(m3$Sol, m3$VCV)))
print(i)
}
}

if(file.exists("~/Work/Boots/Data/Raw/PO.csv")){
print("MVasUV")
  POdata<-read.csv("~/Work/Boots/Data/Raw/PO.csv") 
  bc<-boxcox(POdata[,"PO"]~1)  
  bc.par<-bc$x[which.max(bc$y)]     
  POdata$PO<-c(scale(POdata$PO^bc.par))

  Rdata<-read.csv("~/Work/Boots/Data/Raw/Resistant.csv")

  Rdata2<-Rdata[rep(1:length(Rdata[,1]), Rdata[,4]),]
  Rdata2[,4]<-0
  Rdata3<-Rdata[rep(1:length(Rdata[,1]), Rdata[,5]),]
  Rdata3[,4]<-1
  Rdata4<-rbind(Rdata2, Rdata3)
  Rdata4$PO<-rep(NA, dim(Rdata4)[1])
  Rdata4<-Rdata4[,-c(2,3,5)]

  data<-rbind(Rdata4, cbind(Family=POdata[,1], Infected=rep(NA, dim(POdata)[1]), PO=POdata$PO))
  data$Family<-as.factor(data$Family)

  prior=list(R=list(V=diag(2), n=2, fix=2), G=list(G1=list(V=diag(2), n=1)))

  system.time(model1<-MCMCglmm(cbind(PO, Infected)~trait, random=~us(trait):Family, rcov=~idh(trait):units, family=c("gaussian", "categorical"), data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin))
if(plotit){
plot(mcmc(cbind(model1$Sol, model1$VCV)), ask=FALSE)
}

  names(POdata)[2]<-"variable"
  names(Rdata4)[2]<-"variable"

  ndata<-cbind(rbind(POdata[,1:2], Rdata4[,c(1:2)]), family=factor(c(rep("gaussian", dim(POdata)[1]), rep("categorical", dim(Rdata4)[1])), c("gaussian","categorical" )))
  ndata$Family<-as.factor(ndata$Family)
  ndata$trait<-as.factor(ndata$family)
  prior=list(R=list(V=diag(2), n=2, fix=2), G=list(G1=list(V=diag(2), n=1)))

  system.time(model2<-MCMCglmm(variable~family, random=~us(family):Family, rcov=~idh(family):units, data=ndata, prior=prior, family=NULL,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin))
}else{
print("file ~/Work/Boots/Data/Raw/PO.csv does not exist")
}

# caetgorical k=3 + random 5.9 seconds
print("res22")
res22<-matrix(NA, nsim,10)
R=diag(2)
G=matrix(c(2,0.5,0.5,1),2,2)
prior=list(R=list(V=R, n=1, fix=1),G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:100,900,replace=TRUE))
l<-mvrnorm(900, c(-1,0), R)+mvrnorm(100, c(0,0), G)[fac,]

y<-1:900
for(j in 1:900){
y[j]<-which(rmultinom(1,1,prob=c(1,exp(l[j,][1]), exp(l[j,][2])))==1)
}

data=data.frame(y=y, fac=fac)
m1<-MCMCglmm(y~trait-1, random=~us(trait):fac, family=c("categorical"), rcov=~us(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
tpar<-c(-1,0,2,0.5,0.5,1,1,0,0,1)
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res22 different from expected"))
}
res22[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# binary probit

print("res23")
res23<-matrix(NA, nsim,3)
R<-as.matrix(1)
G<-as.matrix(2)
prior=list(R=list(V=R, n=1, fix=1))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
x<-runif(300)
g<-rnorm(300, -1+0.5*x, R)
cp<-0.5
pr<-cbind(pnorm(-g),1-pnorm(-g))

y1<-1:300
for(j in 1:300){
  y1[j]<-which(rmultinom(1, 1, pr[j,])==1)
}

table(y1)
cp<-qnorm(cumsum(c(0,table(y1)/length(y1))), 0, sqrt(1+R))
cp-cp[2]

data=data.frame(y1=y1, fac=fac, xcov=x)
m1<-MCMCglmm(y1~xcov,data=data, prior=prior, family="ordinal", verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
tpar<-c(-1,0.5,1)
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res23 different from expected"))
}
res23[i,]<-c(posterior.mode(mcmc(cbind(m1$Sol, m1$VCV))))
print(i)
}

# 4 category ordinal
print("res24")
res24<-matrix(NA, nsim,5)
R<-as.matrix(1)
G<-as.matrix(2)
prior=list(R=list(V=R, n=1, fix=1)) #,G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
x<-runif(300)
g<-rnorm(300, -1+1*x, R)
cp<-c(0.5,1)
pr<-cbind(pnorm(-g),pnorm(cp[1]-g)-pnorm(-g),  pnorm(cp[2]-g)-pnorm(cp[1]-g), 1-pnorm(cp[2]-g))

y1<-1:300
for(j in 1:300){
  y1[j]<-which(rmultinom(1, 1, pr[j,])==1)
}

data=data.frame(y1=y1, fac=fac, xcov=x)
m1<-MCMCglmm(y1~xcov,data=data, prior=prior, family="ordinal", verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)), ask=FALSE)
}
tpar<-c(cp,-1,1,1)
if(any(HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res24 different from expected"))
}
res24[i,]<-c(posterior.mode(mcmc(cbind(m1$CP, m1$Sol, m1$VCV))))
print(i)
}


# bivariate ordinal

print("res25")
res25<-matrix(NA, nsim,13)
R<-diag(2)
G<-matrix(c(1,0.5,0.5,1),2,2)
prior=list(R=list(V=R, n=1, fix=1),G=list(G1=list(V=G, n=2)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
x<-runif(300)
g<-mvrnorm(300, c(0,0), R)+mvrnorm(75, c(0,0), G)[fac,]
g[,1]<-g[,1]-1+1*x
g[,2]<-g[,2]-0.5-1*x

cp1<-c(0.5,1)
cp2<-c(0.75)
pr1<-cbind(pnorm(-g[,1]),pnorm(cp1[1]-g[,1])-pnorm(-g[,1]),  pnorm(cp1[2]-g[,1])-pnorm(cp1[1]-g[,1]), 1-pnorm(cp1[2]-g[,1]))
pr2<-cbind(pnorm(-g[,2]),pnorm(cp2[1]-g[,2])-pnorm(-g[,2]),  1-pnorm(cp2[1]-g[,2]))


y1<-1:300
y2<-1:300
for(j in 1:300){
  y1[j]<-which(rmultinom(1, 1, pr1[j,])==1)
  y2[j]<-which(rmultinom(1, 1, pr2[j,])==1)
}
table(y1,y2)

data=data.frame(y1=y1, y2=y2, fac=fac, xcov=x)
m1<-MCMCglmm(cbind(y1,y2)~trait+trait:xcov-1,random=~us(trait):fac, rcov=~idh(trait):units, data=data, prior=prior, family=cbind("ordinal", "ordinal"),verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)), ask=FALSE)
}
tpar<-c(cp1, cp2,-1,-0.5,1,-1,1,0.5,0.5,1,1,1)
if(any(HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res25 different from expected"))
}
res25[i,]<-posterior.mode(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))
print(i)
}

print("res26")
res26<-matrix(NA, nsim,3)
R<-diag(1)
G<-diag(1)
prior=list(R=list(V=R, nu=1),G=list(G1=list(V=G, nu=1)))
for(i in 1:nsim){
fac1<-as.factor(sample(1:75,300,replace=TRUE))
fac2<-as.factor(sample(1:75,300,replace=TRUE))
y<-mvrnorm(300, 0, R)+mvrnorm(75, 0, G)[fac1]+mvrnorm(75, 0, G)[fac2]
data=data.frame(y=y, fac1=fac1,fac2=fac2)
m1<-MCMCglmm(y~1,random=~idv(as.numeric(fac1)+fac2), data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}

res26[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
tpar<-c(0,1,1)
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res26 different from expected"))
}


print("res27")
res27<-matrix(NA, nsim,6)
R<-diag(1)
G<-diag(1)
G2<-diag(1)
prior=list(R=list(V=R, nu=1),G=list(G1=list(V=G, nu=1, alpha.mu=1, alpha.V=100), G2=list(V=G2, nu=1, alpha.mu=1, alpha.V=100)))

for(i in 1:nsim){
fac1<-as.factor(sample(1:75,300,replace=TRUE))
fac2<-as.factor(sample(1:75,300,replace=TRUE))
id<-as.factor(sample(1:75,300,replace=TRUE))
fac3<-as.factor(sample(1:3, 300, T))
y<-mvrnorm(300, 0, R)+mvrnorm(75, 0, G)[fac1]+mvrnorm(75, 0, G)[fac2]+rnorm(75,0,sqrt(G2))[id]
data=data.frame(y=y, fac1=fac1,fac2=fac2, fac3=fac3, id=id)
m1<-MCMCglmm(y~fac3,random=~idv(fac1+fac2)+id, data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)

if(plotit){
plot(mcmc(cbind(mcmc(m1$Sol[,1:3]), m1$VCV)), ask=FALSE)
}
tpar<-c(0,0,0,1,1,1)
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res27 different from expected"))
}

res27[i,]<-posterior.mode(mcmc(cbind(m1$Sol[,1:3], m1$VCV)))
print(i)
}

# bivariate binary
print("res28")
T<-2
res28<-matrix(NA, nsim,T*(1+T))
R<-cbind(c(1,-0.25),c(-0.25,1))
mu<-matrix(c(0.3,0.9),T,1)
prior=list(R=list(V=diag(T), nu=T+1))

for(i in 1:nsim){
	y<-mvrnorm(300, mu, R)
	data=data.frame(apply(y,2, function(x){rbinom(length(x), 1, plogis(x))}))
	m1<-MCMCglmm(cbind(X1, X2)~trait-1,rcov=~cor(trait):units, family=rep("categorical", T), data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
	
	if(plotit){
		plot(mcmc(cbind(mcmc(m1$Sol), m1$VCV)), ask=FALSE)
	}
	tpar<-c(mu,c(R))
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res28 different from expected"))
}

	res28[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
	print(i)
}


# bivaraite 4 category ordinal
print("res29")
res29<-matrix(NA, nsim,12)
R<-cbind(c(1,-0.35), c(-0.35,1))

prior=list(R=list(V=diag(2), n=3)) #,G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
	fac<-as.factor(sample(1:75,300,replace=TRUE))
	x<-runif(300)
	g<-mvrnorm(300, c(-1,0), R)+cbind(x*-1, x*0.5) #+mvrnorm(75, c(0), G)[fac]
	cp<-c(0.5,1)
	pr<-cbind(pnorm(-g[,1]),pnorm(cp[1]-g[,1])-pnorm(-g[,1]),  pnorm(cp[2]-g[,1])-pnorm(cp[1]-g[,1]), 1-pnorm(cp[2]-g[,1]))
	cp2<-c(1,2)
	pr2<-cbind(pnorm(-g[,2]),pnorm(cp2[1]-g[,2])-pnorm(-g[,2]),  pnorm(cp2[2]-g[,2])-pnorm(cp2[1]-g[,2]), 1-pnorm(cp2[2]-g[,2]))
	
	y1<-1:300
	for(j in 1:300){
		y1[j]<-which(rmultinom(1, 1, pr[j,])==1)
	}
	y2<-1:300
	for(j in 1:300){
		y2[j]<-which(rmultinom(1, 1, pr2[j,])==1)
	}
	data=data.frame(y1=y1, y2=y2, fac=fac, xcov=x)
	m1<-MCMCglmm(cbind(y1, y2)~trait+trait:xcov-1,rcov=~cor(trait):units, data=data, prior=prior, family=c("ordinal", "ordinal"), verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
	if(plotit){
		plot(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)), ask=FALSE)
	}
	tpar<-c(cp,cp2, -1,0,-1,0.5, c(R))
        if(any(HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,2]<tpar)){
        print(paste(sum(HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res29 different from expected"))
        }
	res29[i,]<-c(posterior.mode(m1$CP),posterior.mode(m1$Sol),posterior.mode(m1$VCV))
	print(i)
}

# hurdle Poisson
print("res30")
res30<-matrix(NA, nsim,12)
data("bioChemists", package = "pscl")
prior = list(R = list(V = diag(2), nu = 1.002, fix = 2)) 
c2<-(16*sqrt(3)/(15*pi))^2

for(i in 1:nsim){
m5d.1 <- MCMCglmm(art ~ trait-1+trait:(fem+mar+kid5 +phd+ment), rcov = ~idh(trait):units,  data = bioChemists, prior = prior, family = "hupoisson", nitt=13000, thin=10, burnin=3000) 
res30[i,][seq(2,12,2)]<-colMeans(m5d.1$Sol[,seq(2,12,2)]/sqrt(1+c2))
res30[i,][seq(1,11,2)]<-colMeans(m5d.1$Sol[,seq(1,11,2)]+0.5*m5d.1$VCV[,1])
}

n<-100
l<-rnorm(n, -1, sqrt(1))
t<-(-log(1-runif(n)*(1-exp(-exp(l)))))
y<-rpois(n,exp(l)-t)+1
y<-c(rep(0, length(y)/2), y)
data=data.frame(y=y)
prior=list(R=list(V=diag(2), fix=2, nu=1))
m1<-MCMCglmm(y~trait-1, rcov=~idh(trait):units, data=data, family="hupoisson", prior=prior)
# truncated Poisson

print("res31")
res31<-matrix(NA, nsim,2)

n<-200

for(i in 1:nsim){
l<-rnorm(n, -1, sqrt(1))
t<-(-log(1-runif(n)*(1-exp(-exp(l)))))
y<-rpois(n,exp(l)-t)+1

dat<-data.frame(y=y)
m1<-MCMCglmm(y~1, family="ztpoisson", data=dat)

res31[i,]<-c(posterior.mode(m1$Sol),posterior.mode(m1$VCV))
}


# geometric

n<-200
print("res32")
res32<-matrix(NA, nsim,2)

for(i in 1:nsim){

y<-rgeom(n,plogis(rnorm(n,-1,sqrt(0.5))))

dat<-data.frame(y=y)
m1<-MCMCglmm(y~1, family="geometric", data=dat)

res32[i,]<-c(posterior.mode(m1$Sol),posterior.mode(m1$VCV))

}



print("res33")
res33<-matrix(0, nsim, 4)
R<-diag(2)
prior=list(R=list(V=R, n=1, fix=2), B=list(mu=c(0,0), V=matrix(c(1000,0,0,pi^2/3),2,2)))
tune=list(diag(2))
tpar<-c(0, -1, 1, 1)
for(i in 1:nsim){
l<-mvrnorm(300,c(0,-1), R)
y<-rzibinom(300, 20, exp(l[,1])/(1+exp(l[,1])), plogis(l[,2]))
data=data.frame(success=y,failure=20-y)
m1<-MCMCglmm(cbind(success,failure)~trait-1, rcov=~idh(trait):units, data=data, family="zibinomial",prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res19 different from expected"))
}
res33[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}



tpar<-c(1, 1,1,1,1,1,1,1,-1,1,2,-1,1,2,-1,1,2,-1, 1, 0,0,1,2,-1, 1, 0.5, 0.5, 2, 2,-1, 1, 2, 2,-1, 1,0,0, 2, 2,-1,1,1,0.5,0.5,2,-1,1,1,2,-1,1,2,1,1,2,-1,1,2,0.5,0.5,1,1,2,0,1,2, coef, coef,-1,1,2,0,1,0,1,0,1,1,-1,1,0,0,1,1, -0.5, 0.2, 1,1,0,0,1,-1,1,2,0.25, 0.25,1,-1, 0, 2,1,1,-1, 0, 2,1,1,-1, 0, 2,0,0,1,1,-1,0,2,0.5,0.5,1,1,0,0,1,-1,0.5,1,0.5,1,-1,1,1,0.5,1,0.75,-1,-0.5,1,-1,1,0.5,0.5,1,1,1,0,1,1,0,0,0,1,1,1,0.3,0.9,1,-0.25,-0.25,1,0.5,1,1,2,-1,0,-1,0.5,1,-0.35, -0.35,1)

est<-colMeans(cbind(res1, res2, res3, res3b, res4, res4c, res5,res5b,res6, res7,res7b, res8,res9, res10, res11, res12, res13, res14, res15, res17a,res17b, res18, res19, res19b, res20, res21, res21b, res21c, res22, res23, res24, res25, res26, res27, res28, res29))

plot(est~tpar)
abline(0,1)


