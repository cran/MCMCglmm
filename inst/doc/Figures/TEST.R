#source("~/Desktop/MCMCglmmTEST.R")
#source("~/Work/AManal/MCMCglmm_1.02/R/MCMCglmm.R")
library(MCMCglmm)
library(VGAM)

nsim<-1

# poisson test 1.2 seconds OK
print("res1")
R<-diag(1)
res1<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1), G=list(G1=list(V=as.matrix(1), n=1)))
for(i in 1:nsim){
l<-exp(rnorm(100,1,R))
y<-rpois(100,l)
data=data.frame(y1=y, y2=y)
m1<-MCMCglmm(cbind(y1)~1, family="poisson", data=data, prior=prior)
res1[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, posterior.mode)
print(i)
}


# multinomial test J=1 test 0.9 seconds  OK
print("res2")
res2<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1))
for(i in 1:nsim){
y<-rbinom(100,10,inv.logit(rnorm(100,1,1)))
data=data.frame(y1=y, y2=10-y)
m1<-MCMCglmm(cbind(y1,y2)~1, family="multinomial2", data=data, prior=prior)
res2[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, posterior.mode)
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
m1<-MCMCglmm(y1~1, family="categorical", data=data, prior=prior)
res3[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, posterior.mode)
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
m1<-MCMCglmm(y1~1,random=~fac,  data=data, prior=prior)
res4[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, posterior.mode)
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
m1<-MCMCglmm(y1~1,random=~fac,  data=data, prior=prior, family="categorical")
res4c[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, posterior.mode)
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
m1<-MCMCglmm(y1~1, random=~animal, pedigree=Ped, data=data, prior=prior, nitt=25000, thin=20, burnin=5000)
m2<-MCMCglmm(y1~1, random=~animal, pedigree=Ped, data=data, prior=prior, nitt=25000, thin=20, burnin=5000)
res5[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, posterior.mode)
print(i)
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
m1<-MCMCglmm(y1~1,random=~us(facf):fac,  data=data, prior=prior, nitt=100000, thin=8, burnin=20000)
m2<-MCMCglmm(y1~1,random=~us(facf):fac,  data=data, prior=prior, nitt=100000, thin=8, burnin=20000)

res6[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, posterior.mode)
print(i)
}

# gauss with idh random effect 1.7 seconds
print("res7")
res7<-matrix(NA,nsim,4)
R<-as.matrix(2)
G=matrix(c(1,0,0,2),2,2)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=2)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
facf<-as.factor(sample(1:2,300,replace=TRUE))
y<-mvrnorm(300, c(-1), R)+rowSums(mvrnorm(75, c(0,0), G)[fac,]*cbind(facf==1,facf==2))
data=data.frame(y1=y, fac=fac, facf=facf)
m1<-MCMCglmm(y1~1,random=~idh(facf):fac,  data=data, prior=prior)
res7[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, posterior.mode)
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
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","gaussian"), rcov=~us(trait):units, data=data, prior=prior)
res8[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, posterior.mode)
print(i)
}

# bivariate gauss idh residual 1.7 seconds
print("res9")
res9<-matrix(NA, nsim,4)
R=matrix(c(1,0,0,2),2,2)
prior=list(R=list(V=R, n=1))
for(i in 1:nsim){
y<-mvrnorm(300, c(-1,1), R)
data=data.frame(y1=y[,1], y2=y[,2])
data$y1[sample(1:300, 50)]<-NA
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","gaussian"), rcov=~idh(trait):units, data=data, prior=prior)
res9[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, posterior.mode)
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
m1<-MCMCglmm(cbind(y1,y2)~trait-1, random=~idh(trait):fac, family=c("gaussian","gaussian"), rcov=~idh(trait):units, data=data, prior=prior)
res10[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, posterior.mode)
print(i)
}

# bivariate binoimal + random 5.9 seconds
print("res11")
res11<-matrix(NA, nsim,8)
R=matrix(c(1,0,0,2),2,2)
G=matrix(c(2,0.5,0.5,1),2,2)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-mvrnorm(300, c(-1,1), R)+mvrnorm(300, c(0,0), G)[fac,]
y1<-rbinom(300,10,inv.logit(y[,1]))
y2<-rbinom(300,10,inv.logit(y[,2]))
data=data.frame(y1s=y1,y1f=10-y1,y2s=y2,y2f=10-y2, fac=fac)
m1<-MCMCglmm(cbind(y1s,y1f,y2s,y2f)~trait-1, random=~us(trait):fac, family=c("multinomial2","multinomial2"), rcov=~idh(trait):units, data=data, prior=prior)
res11[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, posterior.mode)
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
res12<-MCMCglmm(y1~1,random=~fac,  data=data, prior=prior)


# Jennys data - once as gaussian once as binomial 32 seconds

if(file.exists("~/Work/Jenny/Data/Intermediate/ThirdC.R")){
  nsim<-1
  m1R<-dget("~/Work/Jenny/Data/Intermediate/ThirdC.R")
  res13<-matrix(NA, nsim,39)
  res14<-matrix(NA, nsim,39)
  print("res13")
  print("res14")

  firstP<-read.table("~/Work/Jenny/Data/Raw/Third_paternal.txt", header=T)
  firstP$day<-as.factor(firstP$day)

  coef<-apply(cbind(m1R$Sol, m1R$VCV), 2, median)[c(1:8,34:59, 9,15,21,27,33)]
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
    m1test<-MCMCglmm(l~virus+day, random=~us(virus):line+f2rep, rcov=~idh(virus):units, family=c("gaussian"), data=firstP, prior=prior)
    m1test2<-MCMCglmm(cbind(nodead, noalive)~virus+day, random=~us(virus):line+f2rep, rcov=~idh(virus):units, family=c("multinomial2"), data=firstP, prior=prior)
    res13[i,]<-apply(cbind(m1test$Sol, m1test$VCV),2, posterior.mode)
    res14[i,]<-apply(cbind(m1test2$Sol, m1test2$VCV),2, posterior.mode)
    print(i)
  }

  plot(apply(res13, 2, median)~coef)
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
m1<-MCMCglmm(y1~1,random=NULL,  data=data, prior=prior, mev=mev)
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
  m1<-MCMCglmm(ym~1, data=data)
}
m2<-MCMCglmm(cbind(y1, y2)~1, data=data, family="cengaussian")
res17a[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, median, na.rm=T)
res17b[i,]<-apply(mcmc(cbind(m2$Sol, m2$VCV)), 2, median, na.rm=T)
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
m1<-MCMCglmm(cbind(y1, y2)~1, data=data, family="cenpoisson", prior=prior, pl=TRUE)
res18[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, median, na.rm=T)
print(i)
}


print("res19")
res19<-matrix(0, nsim, 6)
R<-diag(2)
prior=list(R=list(V=R, n=1, fix=2), B=list(mu=c(0,0), V=matrix(c(10000000000,0,0,3.06),2,2)))
tune=list(diag(2))
for(i in 1:nsim){
l<-mvrnorm(300,c(0,-0.5), R)
y<-rzipois(300, exp(1+l[,1]), inv.logit(l[,2]))
data=data.frame(y1=y)
m1<-MCMCglmm(y1~trait-1, rcov=~idh(trait):units, data=data, family="zipoisson",prior=prior)
res19[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, median)
print(i)
}

# bivariate gauss + categorical  residual 1.7 seconds

res20<-matrix(NA, nsim,6)
print("res16")
tune<-diag(2)
R=matrix(c(2,0.25,0.25,1),2,2)
prior=list(R=list(V=R, n=2, fix=2))
for(i in 1:nsim){
y<-mvrnorm(300, c(-1,1), R)
data=data.frame(y1=y[,1], y2=rbinom(300, 1, inv.logit(y[,2])), y3=y[,2])
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","categorical"), rcov=~us(trait):units, data=data, prior=prior)
res20[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, median)
print(i)
}



# gauss random regression

res21<-matrix(NA, nsim,7)
print("res21")
G=matrix(c(2,0.25,0.25,1),2,2)
R=matrix(1,1,1)
prior=list(R=list(V=R, n=1), G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
int.slope<-mvrnorm(300, c(0,0), G)
ind<-gl(300,3)
time<-rnorm(900)
y<-(-1+int.slope[,1][ind]*0.7071068)+(1+1.224745*time*int.slope[,2][ind])
y<-y+rnorm(900,0,R)
data=data.frame(y1=y, time=time, ind=ind)
m1<-MCMCglmm(y1~time, random=~leg(time,1):ind, data=data, prior=prior)
res21[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, median)
print(i)
}
if(file.exists("~/Work/Boots/Data/Raw/PO.csv")){

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

  system.time(model1<-MCMCglmm(cbind(PO, Infected)~trait, random=~us(trait):Family, rcov=~idh(trait):units, family=c("gaussian", "categorical"), data=data, prior=prior))

  names(POdata)[2]<-"variable"
  names(Rdata4)[2]<-"variable"

  ndata<-cbind(rbind(POdata[,1:2], Rdata4[,c(1:2)]), family=factor(c(rep("gaussian", dim(POdata)[1]), rep("categorical", dim(Rdata4)[1])), c("gaussian","categorical" )))
  ndata$Family<-as.factor(ndata$Family)
  prior=list(R=list(V=diag(2), n=2, fix=2), G=list(G1=list(V=diag(2), n=1)))

  system.time(model2<-MCMCglmm(variable~family, random=~us(family):Family, rcov=~idh(family):units, data=ndata, prior=prior, family=NULL))
}else{
print("file ~/Work/Boots/Data/Raw/PO.csv does not exist")
}





