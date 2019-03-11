#question1 (b)
pmean<- seq(-20,0,0.002)
n<-16
mu<- -7.9875
ssquare<-30.2665
sd<-30.2665^(1/2)
z<-sqrt(n)*(pmean-mu)/sd
n1<- n-1
den1<-dt(z,n1)
plot(pmean,den1,xlab = "mu", ylab = "density", main = "marginal posterior distribution of mu")

#question1(d)
ssq<-seq(1,100,0.001)
library(LaplacesDemon)
den2<-dinvchisq(ssq,df=15,scale=30.2665)
plot(ssq,den2,xlab = "sigma2", ylab = "density", main = "marginal posterior distribution of sgima2")

#question2(b)
k0<-0.25
pmean1<- seq(-40,20,0.001)
mun<-(n*mu+k0*0)/(k0+n)
mun
sun<-(k0*mu*mun+(n-1)*30.2665+2*0.5)/(n+2*0.5)
sun
z2<-(pmean1-mun)/sun*sqrt(n+k0)
k0<-0.25
df<-n+2*0.5
df
den3<- dt(z2,17)
den3
plot(pmean1,den3,type ="l",xlab = "mu", ylab = "density", main = "marginal posterior distribution of mu")
length(pmean1)
length(den3)

#question2(d)
ssq2<- seq(1,60,0.001)
den4<- dinvchisq(ssq2,df=17,scale=27.68836)
plot(ssq2,den4,xlab = "sigma2", ylab = "density", main = "marginal posterior distribution of sgima2")

#question3(a)
y=c(3.7,-6.7,-10.5,-6.1,-17.6,2.3,-7.9,-8.9,-4.5,-7.7,-9.4,-10.4,-10.9,-9.3,-16.7,-7.2)
tau=25
l=100
m=1000
a11=0.5
a111=n/2+a11
b11=0.5
mugb=rep(NA,l+m)
sigb=rep(NA,1+m)
sigb[1]=1
for(i in 1:(l+m)){
  musigma=1/(1/tau+n/sigb[i])
  mumu=n*mu*musigma/sigb[i]
  mugb[i]=rnorm(1,mumu,sqrt(musigma))
  b111=(n*mugb[i]^2-2*sum(y)*mugb[i]+sum(y^2)+2*b11)/2
  invsigb=rgamma(1,a111,beta1)
  sigb[i+1]=1/invsigb
}
musigma=1/(1/tau+n/sigb[l+m])
mumu=n*mu*musigma/sigb[l+m]
mugb[l+m]=rnorm(1,mumu,sqrt(musigma))
lastmu=mugb[(l+1):(l+m)]
lastmu

#problem3(c)
pmean3<-seq(0,1000,0.01)
plot(lastmu,type = "l", xlab = "mu", ylab = "density", main = "marginal posterior density of mu")
hist(lastmu,xlab = "mu", ylab = "density", main = "marginal posterior density of mu")

#problem3(e)
y=c(3.7,-6.7,-10.5,-6.1,-17.6,2.3,-7.9,-8.9,-4.5,-7.7,-9.4,-10.4,-10.9,-9.3,-16.7,-7.2)
tau=25
l=100
m=1000
a11=0.5
a111=n/2+a11
b11=0.5
mugb=rep(NA,l+m)
sigb=rep(NA,1+m)
sigb[1]=100
for(i in 1:(l+m)){
  musigma=1/(1/tau+n/sigb[i])
  mumu=n*mu*musigma/sigb[i]
  mugb[i]=rnorm(1,mumu,sqrt(musigma))
  b111=(n*mugb[i]^2-2*sum(y)*mugb[i]+sum(y^2)+2*b11)/2
  invsigb=rgamma(1,a111,beta1)
  sigb[i+1]=1/invsigb
}
musigma=1/(1/tau+n/sigb[l+m])
mumu=n*mu*musigma/sigb[l+m]
mugb[l+m]=rnorm(1,mumu,sqrt(musigma))
lastmu=mugb[(l+1):(l+m)]
lastmu

#problem3(iv)
ciindependent<-qt(c(0.025,0.975),15)
ciindependent
q1<- ciindependent*sd/sqrt(16)-7.9875
q1
ciconjugate<- qt(c(0.025,0.975),17)
ciconjugate
q2<- ciconjugate*sqrt(sun)/sqrt(n+k0)+mun
q2
cigb<- quantile(lastmu,probs = c(0.025,0.975))
cigb




#problem4(a)
install.packages("OOmisc")
library(OOmisc)
x<- seq(0,1,0.01)
denbeta<-function(x){
  denbetaf=0.4*dbeta(x,4,2)+0.6*dbeta(x,2,6)
  denbetaf
}
den<-denbeta(x)
den
plot(x, den, type="l",ylim=c(0,3))


candidate<-function(x){
  candid<- dunif(x,0,1)
  candid
}
can<- candidate(x)
can
c = max(den/can)
c
lines(x,c*can,col=3)
dev.off()

prob1=can/sum(can*0.001)
prob1
samplecan<-runif(1000,0,1)
samplecan
install.packages("distr")
library(distr)
mix<- UnivarMixingDistribution(Beta(4,2),Beta(2,6),mixCoeff = c(0.4,0.6))
rmyMix <- r(mix)
sampledenf<-rmyMix(1000)
sampledenf
k=1
n1<-sampledenf/(c*samplecan)
n1
length(n1)
T=1000
re<-rep(NA,T)
xxx<-seq(0,1,0.001)
re
u<- runif(1)
u
while(k <= T){
  nre<-sample(xxx,1,replace=TRUE)
  cnew<-can(nre)
  dnew<-den(nre)
  u=runif(1)
  if(u<=dnew/c*cnew){
    re[k]=nre;
    k=k+1;
  }
}
length(u)
xxx
re
hist(re,prob=T,ylim = c(0,5))



l=500
m=1000
theta=rep(NA,l+m)
theta[1]=1
N<-1
for(k in 2:(l+m)){
  theanew=runif(1);
  theta[k] = theta[k-1]
  u=runif(1)
  ratio=f(theanew)/f(theta[k-1])
  r=min(ratio,1)
  if(u<r){
    theta[k]=theanew;
    N<-N+1
  }
}
N/(l+m)
plot(theta,type='l')
thetas<- theta[(l+1):(l+m)]
hist(thetas,probability = T)
points(density(thetas),lty=1)



    
    
