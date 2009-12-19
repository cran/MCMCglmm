#source("~/Desktop/MCMCglmmTEST.R")
#source("~/Work/AManal/MCMCglmm_1.11/inst/doc/Figures/TEST.R")
library(MCMCglmm)
library(VGAM)
verbose=FALSE
plotit=TRUE
leg=TRUE
nsim<-1
nitt<-5000
thin<-4
burnin<-1000

# poisson test 1.2 seconds OK
print("res1")
R<-diag(1)
res1<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1))
for(i in 1:nsim){
l<-exp(rnorm(100,1,R))
y<-rpois(100,l)
data=data.frame(y1=y, y2=y)
m1<-MCMCglmm(cbind(y1)~1, family="poisson", data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
res1[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


# multinomial test J=1 test 0.9 seconds  OK
print("res2")
res2<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1))
for(i in 1:nsim){
y<-rbinom(100,10,inv.logit(rnorm(100,1,1)))
data=data.frame(y1=y, y2=10-y)
m1<-MCMCglmm(cbind(y1,y2)~1, family="multinomial2", data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
res2[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# categorical test J=1 test 0.9 seconds OK
print("res3")
res3<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
for(i in 1:nsim){
y<-rbinom(100,1,inv.logit(rnorm(100,1,1)))
data=data.frame(y1=y, y2=1-y)
prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
m1<-MCMCglmm(y1~1, family="categorical", data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
res3[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# categorical test J=1 test with slice sampling 0.9 seconds OK
print("res3b")
res3<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
for(i in 1:nsim){
y<-rbinom(100,1,inv.logit(rnorm(100,1,1)))
data=data.frame(y1=y, y2=1-y)
prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
m1<-MCMCglmm(y1~1, family="categorical", data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin, slice=TRUE)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
res3[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


# gauss with blocked random  1.5 seconds OK
print("res4")
res4<-matrix(NA, nsim,3)
R<-as.matrix(2)
G<-as.matrix(1)
prior<-list(R=list(V=as.matrix(1), n=1), G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:50,300,replace=TRUE))
ffac<-gl(2,150)
y<-mvrnorm(300, c(-1), R)+mvrnorm(50, c(0), G)[fac]
data=data.frame(y1=y, fac=fac, ffac=ffac)
m1<-MCMCglmm(y1~1,random=~fac,  data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
res4[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# binary with blocked random  1.5 seconds OK

print("res4")
res4c<-matrix(NA, nsim,3)
R<-as.matrix(2)
G<-as.matrix(1)
prior=list(R=list(V=R, n=1, fix=1),G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-mvrnorm(300, c(-1), R)+mvrnorm(75, c(0), G)[fac]
data=data.frame(y1=rbinom(300, 1,inv.logit(y)), fac=fac, y2=y)
m1<-MCMCglmm(y1~1,random=~fac,  data=data, prior=prior, family="categorical",verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
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
for(i in 1:nsim){
Ped<-cbind(1:400, c(rep(NA,100), sample(1:50,300,TRUE)),c(rep(NA,100), sample(51:100,300,TRUE)))
y<-mvrnorm(300, c(-1), R)+rbv(Ped,G)[101:400]
data=data.frame(y1=y, animal=as.factor(Ped[,1][101:400]))
system.time(m1<-MCMCglmm(y1~1, random=~animal, pedigree=Ped, data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin))
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
m2<-MCMCglmm(y1~1, random=~animal, pedigree=Ped, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m2$Sol, m2$VCV)), ask=FALSE)
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
m2<-MCMCglmm(y1~1,random=~us(facf):fac,  data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m2$Sol, m2$VCV)), ask=FALSE)
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
m2<-MCMCglmm(y1~1,random=~us(facf):fac,  data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m2$Sol, m2$VCV)), ask=FALSE)
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
for(i in 1:nsim){
y<-mvrnorm(300, c(-1,1), R)
data=data.frame(y1=y[,1], y2=y[,2])
data$y1[sample(1:300, 30)]<-NA
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","gaussian"), rcov=~us(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
res8[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# bivariate gauss idh residual 1.7 seconds
print("res9")
res9<-matrix(NA, nsim,4)
R=matrix(c(1,0,0,2),2,2)
prior=list(R=list(V=R, n=2))
for(i in 1:nsim){
y<-mvrnorm(300, c(-1,1), R)
data=data.frame(y1=y[,1], y2=y[,2])
data$y1[sample(1:300, 50)]<-NA
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","gaussian"), rcov=~idh(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
res9[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# bivariate gauss + random 2.8 seconds
print("res10")
res10<-matrix(NA, nsim,6)
R=matrix(c(1,0,0,2),2,2)
G=matrix(c(2,0.5,0.5,1),2,2)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-mvrnorm(300, c(-1,1), R)+mvrnorm(75, c(0,0), G)[fac,]
data=data.frame(y1=y[,1], y2=y[,2], fac=fac)
m1<-MCMCglmm(cbind(y1,y2)~trait-1, random=~idh(trait):fac, family=c("gaussian","gaussian"), rcov=~idh(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
res10[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# bivariate binoimal + random 5.9 seconds
res11<-matrix(NA, nsim,8)
R=matrix(c(1,0,0,2),2,2)
G=matrix(c(2,0.5,0.5,1),2,2)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-mvrnorm(300, c(-1,1), R)+mvrnorm(75, c(0,0), G)[fac,]
y1<-rbinom(300,10,inv.logit(y[,1]))
y2<-rbinom(300,10,inv.logit(y[,2]))
data=data.frame(y1s=y1,y1f=10-y1,y2s=y2,y2f=10-y2, fac=fac)
m1<-MCMCglmm(cbind(y1s,y1f,y2s,y2f)~trait-1, random=~us(trait):fac, family=c("multinomial2","multinomial2"), rcov=~idh(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
res11[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


# gauss with blocked random - all missing 1.9 seconds
print("res12")
R<-as.matrix(2)
G<-as.matrix(1)
prior=list(R=list(V=R, n=100),G=list(G1=list(V=G, n=100)), B=list(mu=0, V=0.000000001))
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-rep(NA,300)
data=data.frame(y1=y, fac=fac)
res12<-MCMCglmm(y1~1,random=~fac,  data=data, prior=prior, singular.ok=TRUE,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)


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
    y<-rbinom(dim(firstP)[1], firstP$total, inv.logit(l))
    firstP$noalive<-y
    firstP$nodead<-firstP$total-y
    m1test<-MCMCglmm(l~virus+day, random=~us(virus):line+f2rep, rcov=~idh(virus):units, family=c("gaussian"), data=firstP, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1test$Sol, m1test$VCV)), ask=FALSE)
}
    m1test2<-MCMCglmm(cbind(noalive, nodead)~virus+day, random=~us(virus):line+f2rep, rcov=~idh(virus):units, family=c("multinomial2"), data=firstP, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1test2$Sol, m1test2$VCV)), ask=FALSE)
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
m1<-MCMCglmm(y1~1,data=data, prior=prior, mev=mev,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
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

res17a[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
res17b[i,]<-posterior.mode(mcmc(cbind(m2$Sol, m1$VCV)))
print(i)
}

# censored Possion data
res18<-matrix(NA, nsim,2)
print("res18")
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
res18[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


print("res19")
res19<-matrix(0, nsim, 6)
R<-diag(2)
prior=list(R=list(V=R, n=1, fix=2), B=list(mu=c(0,0), V=matrix(c(1000,0,0,pi^2/3),2,2)))
tune=list(diag(2))
for(i in 1:nsim){
l<-mvrnorm(300,c(0,-1), R)
y<-rzipois(300, exp(1+l[,1]), inv.logit(l[,2]))
data=data.frame(y1=y)
m1<-MCMCglmm(y1~trait-1, rcov=~idh(trait):units, data=data, family="zipoisson",prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
res19[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


print("res19b")
res19b<-matrix(0, nsim, 8)
R<-diag(2)
prior=list(R=list(V=R, n=1, fix=2),G=list(G1=list(V=1, n=1)))
tune=list(diag(2))
for(i in 1:nsim){
	x<-rnorm(300)
	l<-mvrnorm(300,c(0,-0.5), R)
        fac<-gl(50,6)
        r<-rnorm(50)
	y<-rzipois(300, exp(1+l[,1]+r[fac]+0.2*x), inv.logit(l[,2]))
	data=data.frame(y1=y, x=x, fac=fac)
	m1<-MCMCglmm(y1~trait+at.level(trait, 1):x-1, random=~us(at.level(trait,1)):fac, rcov=~idh(trait):units, data=data, family="zipoisson",prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
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
for(i in 1:nsim){
y<-mvrnorm(300, c(-1,1), R)
data=data.frame(y1=y[,1], y2=rbinom(300, 1, inv.logit(y[,2])), y3=y[,2])
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","categorical"), rcov=~us(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
res20[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# gauss random regression

if(leg){
res21<-matrix(NA, nsim,5)
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


m2<-MCMCglmm(y1~time, random=~us(1+poly(time,1, raw=TRUE)):ind, data=data, prior=prior2,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)

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
res21[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
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

res22[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


# binary probit

print("res23")
res23<-matrix(NA, nsim,2)
R<-as.matrix(1)
G<-as.matrix(2)
prior=list(R=list(V=R, n=1, fix=1)) #,G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
x<-runif(300)
g<-rnorm(300, -1+0.5*x, R) #+mvrnorm(75, c(0), G)[fac]
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
plot(mcmc(cbind(m1$Sol)), ask=FALSE)
}
res23[i,]<-posterior.mode(m1$Sol)
print(i)
}


# 4 category ordinal
print("res24")
res24<-matrix(NA, nsim,4)
R<-as.matrix(1)
G<-as.matrix(2)
prior=list(R=list(V=R, n=1, fix=1)) #,G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
x<-runif(300)
g<-rnorm(300, -1+1*x, R) #+mvrnorm(75, c(0), G)[fac]
cp<-c(0.5,1)
pr<-cbind(pnorm(-g),pnorm(cp[1]-g)-pnorm(-g),  pnorm(cp[2]-g)-pnorm(cp[1]-g), 1-pnorm(cp[2]-g))

y1<-1:300
for(j in 1:300){
  y1[j]<-which(rmultinom(1, 1, pr[j,])==1)
}

data=data.frame(y1=y1, fac=fac, xcov=x)
m1<-MCMCglmm(y1~xcov,data=data, prior=prior, family="ordinal", verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol)), ask=FALSE)
}
res24[i,]<-posterior.mode(m1$Sol)
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
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}

res25[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
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
m1<-MCMCglmm(y~1,random=~idv(fac1+fac2), data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}

res26[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

#source("~/Desktop/MCMCglmmTEST.R")
#source("~/Work/AManal/MCMCglmm_1.11/inst/doc/Figures/TEST.R")

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
m1<-MCMCglmm(y~fac3,random=~idv(fac1+fac2)+id, data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin, pr=TRUE)

if(plotit){
plot(mcmc(cbind(mcmc(m1$Sol[,1:3]), m1$VCV)), ask=FALSE)
}

res27[i,]<-posterior.mode(mcmc(cbind(m1$Sol[,1:3], m1$VCV)))
print(i)
}


