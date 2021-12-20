## -----------------------------------------------------------------------------
library(StatComp21027)
library(compiler)
library(MASS)
n<-c(20,20)
nu<-c(0.1,0.1)
omega<-c(0.1,0.1)
a<-c(1,0.6)
b<-c(0.6,1)
set.seed(3)
dat<-NULL
label<-NULL
for(i in 1:length(n)){
  ni<-rmultinom(1,size=n[i],prob=c(nu[i],omega[i],1-nu[i]-omega[i]))
  x<-c(rep(0,ni[1]),rep(1,ni[2]),rbeta(ni[3],shape1=a[i],shape2=b[i]))
  label1<-rep(letters[i],n[i])
  for(j in 1:4)
  { label2<-rep(letters[j+2],n[i]/4)
  label<-c(label,label2)}
  dat<-rbind(dat,data.frame(label1,label,x))
}
selr(dat)
boot.selr(dat)


