## -----------------------------------------------------------------------------
set.seed(12345) 
library(stats)
r<-runif(1000) 
x<-2/((r)^(1/2))#  Derive the inverse function
hist(x,prob=TRUE,main=expression(f(x)==8*x^(-3)))
y<-seq(2,20,.01)
lines(y,8*y^(-3),col="blue")

## -----------------------------------------------------------------------------
u<-runif(1000) 
x<-(-2*log(u))^(1/2)
hist(x,prob=TRUE,main=expression(f(x)==x*e^(-(x^2)/2)))
y<-seq(0,100,.01)
y1<-y*exp(-(y^2)/2)
lines(y,y1,col="green")

## -----------------------------------------------------------------------------
p<-c(0.1,0.2,0.2,0.2,0.3)
library(knitr)
res<-data.frame(x=c(0,1,2,3,4),p=c(0.1,0.2,0.2,0.2,0.3))
kable(res,digits = getOption("digits"),align = "lccrr")
x<-sample(c(0,1,2,3,4),size=1000,replace=TRUE,prob=c(0.1,0.2,0.2,0.2,0.3))
p0<-sum(x==0)/1000
p1<-sum(x==1)/1000
p2<-sum(x==2)/1000
p3<-sum(x==3)/1000
p4<-sum(x==4)/1000
data.frame(p0,p1,p2,p3,p4)
res1<-data.frame(x=c(0,1,2,3,4),t_p=c(0.1,0.2,0.2,0.2,0.3),e_p=c(p0,p1,p2,p3,p4))
kable(res1,digits = getOption("digits"),align = "lccrr",label = "Empirical with the theoretical probabilities ")

## -----------------------------------------------------------------------------
library(ggplot2)
set.seed(12345)
u<-runif(1000) #Generate 1000 independent and uniformly distributed random variables
# for sigma=1
x<-(-2*log(u))^(1/2)
#graph the histogram
hist(x,prob=TRUE,main=expression(f(x)==x*e^(-(x^2)/2) ))
#Graph the density curve
y<-seq(0,100,.01)
y1<-y*exp(-(y^2)/2)
lines(y,y1,col="green")
# for sigma=5
x1<-(-50*log(u))^(1/2)
#graph the histogram
hist(x1,prob=TRUE,main=expression(f(x)==x/25*e^(-(x^2)/50)))
#Graph the density curve
y<-seq(0,100,.01)
y1<-y/25*exp(-(y^2)/50)
lines(y,y1,col="blue")
# for sigma=10
x2<-(-200*log(u))^(1/2)
#graph the histogram
hist(x2,prob=TRUE,main=expression(f(x)==x/100*e^(-(x^2)/200)))
#Graph the density curve
y<-seq(0,100,.01)
y1<-y/100*exp(-(y^2)/200)
lines(y,y1,col="red")

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 1000
z<-c(0,length(n))
## for p1=0.75
for(i in 1:n)
{
 u<-runif(1)
  if(u<0.75)
   z[i]=rnorm(1)
  else 
    z[i]=rnorm(1,3,1)
}
hist(z,prob=TRUE,main=expression(mixture_distribution_p1=0.75))
y<-seq(-100,100,.01)
y1<-(0.75/((2*pi)^(1/2)))*exp(-(y^2)/2)+(0.25/((2*pi)^(1/2)))*exp(-((y-3)^2)/2)
lines(y,y1,col="purple")
## we will consider the other value of p1
# p1=0.9,then we have
for(i in 1:n)
{
 u<-runif(1)
  if(u<0.9)
   z[i]=rnorm(1)
  else 
    z[i]=rnorm(1,3,1)
}
hist(z,prob=TRUE,main=expression(mixture_distribution_p1=0.9))
y<-seq(-100,100,.01)
y1<-(0.9/((2*pi)^(1/2)))*exp(-(y^2)/2)+(0.1/((2*pi)^(1/2)))*exp(-((y-3)^2)/2)
lines(y,y1,col="blue")

# p1=0.5,then we have
for(i in 1:n)
{
 u<-runif(1)
  if(u<0.5)
   z[i]=rnorm(1)
  else 
    z[i]=rnorm(1,3,1)
}
hist(z,prob=TRUE,main=expression(mixture_distribution_p1=0.5))
y<-seq(-100,100,.01)
y1<-(0.5/((2*pi)^(1/2)))*exp(-(y^2)/2)+(0.5/((2*pi)^(1/2)))*exp(-((y-3)^2)/2)
lines(y,y1,col="red")

## -----------------------------------------------------------------------------
set.seed(12345)
lambda <-2
t0<-10
supper<-100
N_t<-numeric(10000)
for (i in 1:10000) {
N<-rpois(1,lambda*supper)
Un<-runif(N,0,supper) #unordered arrival times
Sn<-sort(Un) #arrival times
n <- min(which(Sn>t0)) #arrivals+1 in [0, t0]
N_t[i] <- n - 1 #arrivals in [0, t0]
}
N_10=ceiling(mean(N_t))
N_10
x<-c(0,length(1000))
for(i in 1:1000)
{
  y<- rgamma(N_10,4,5)
  x[i]<-sum(y)
}
Y<- rgamma(N_10,4,5)
res<-data.frame(mean(x),lambda*t0*0.8,var(x),lambda*t0*4/25)
kable(res,digits = getOption("digits"),align = "lccrr")

## -----------------------------------------------------------------------------
beta_cdf<-function(x){
  n<-10000
  t<-runif(n,0,x)
  theta.hat<-mean(t^2*(1-t)^2/beta(3,3)*x)
  return(theta.hat)
}
Monte_Carlo_estimate<-numeric(9)
Beta_cdf<-numeric(9)
for (i in 1:9){
  Monte_Carlo_estimate[i]=beta_cdf(i*0.1)
  Beta_cdf[i]=pbeta(i*0.1,3,3)
}
res<-data.frame(Monte_Carlo_estimate,Beta_cdf)
kable(res,digits = getOption("digits"),align = "lccrr")

## -----------------------------------------------------------------------------
###using antithetic variables
library(knitr)
Mc<-function(x,m,R=10000,antithetic=TRUE)##,and the m is parameter.
  {
  u<-runif(R/2)
  if(!antithetic)v<-runif(R/2) 
  else v<-1-u
  t<-c(u,v)
  theta_hat<-mean(x^2*t/m^2*exp(-x^2*t^2/2/m^2))
}
k<-1000
Mc1<-numeric(k)
Mc2<-numeric(k)
Mc3<-numeric(k)
Mc4<-numeric(k)
Mc5<-numeric(k)
x<-c(1,2,3,4,5)
m<-c(1,2,3,4,5)
sigma<-c(1,2,3,4,5)
for(j in 1:k)
 {
      Mc1[j]<-Mc(x[1],m[1])
      Mc2[j]<-Mc(x[2],m[2])
      Mc3[j]<-Mc(x[3],m[3])
      Mc4[j]<-Mc(x[4],m[4])
      Mc5[j]<-Mc(x[5],m[5])
}

antithetic_variables<-c(var(Mc1),var(Mc2),var(Mc3),var(Mc4),var(Mc5))
####for the independent variables
MC_independent<-function(x,m,R=10000)##,and the m is parameter.
  {
  u<-runif(R/2)
  v<-runif(R/2)
  theta_hat<-mean(x^2*u/m^2*exp(-x^2*u^2/2/m^2)+x^2*v/m^2*exp(-x^2*v^2/2/m^2))
}
Mc6<-numeric(k)
Mc7<-numeric(k)
Mc8<-numeric(k)
Mc9<-numeric(k)
Mc10<-numeric(k)
for(j in 1:k)
 {
      Mc6[j]<-MC_independent(x[1],m[1])
      Mc7[j]<-MC_independent(x[2],m[2])
      Mc8[j]<-MC_independent(x[3],m[3])
      Mc9[j]<-MC_independent(x[4],m[4])
      Mc10[j]<-MC_independent(x[5],m[5])
}
independent_variables<-c(var(Mc6),var(Mc7),var(Mc8),var(Mc9),var(Mc10))
variance_reduction_percent<-c((var(Mc6)-var(Mc1))/var(Mc6),(var(Mc7)-var(Mc2))/var(Mc7),(var(Mc8)-var(Mc3))/var(Mc8),(var(Mc9)-var(Mc4))/var(Mc9),(var(Mc10)-var(Mc5))/var(Mc10))
res<-data.frame(x,sigma,antithetic_variables,independent_variables,variance_reduction_percent)
kable(res,digits = getOption("digits"),align = "lccrr")

## -----------------------------------------------------------------------------
 x <- seq(1, 100, .01)
g<-x^2/sqrt(2*pi)*exp(-x^2/2)
f1<-1/sqrt(2*pi)*exp(-x^2/2)
f2<-x*exp(-x^2/2)
gs <- c(
            expression(g(x)==x^2/sqrt(2*pi)*exp(-x^2/2)),
            expression(f[1](x)==1/sqrt(2*pi)*exp(-x^2/2)),
            expression(f[2](x)==x*exp(-x^2/2))
            )
 plot(x,g,type = "l", ylab = "",
         ylim = c(0,2), lwd=2,col=1,main='(A)')
    lines(x, f1, lty = 2, lwd = 2,col=2)
    lines(x, f2, lty = 3, lwd = 2,col=3)
legend("topright", legend = gs,
           lty = 1:3,lwd=2, inset=0.02,col=1:3)


## -----------------------------------------------------------------------------
### now we will compute the intergation 
set.seed(12345)
m<-10000
  est<-sd<-numeric(2)
  g<-function(x) {
 x^2/sqrt(2*pi)*exp(-x^2/2)*(x>1)
  }
### for the f1(x)
x<-rnorm(m)
  fg<-g(x)/(1/sqrt(2*pi)*exp(-x^2/2))
  est[1]<-mean(fg)
  sd[1]<-sd(fg)
### for the f2(x)
u<-runif(m)
x<-sqrt(-2*log(u))
fg<-g(x)/(x*exp(-x^2/2))
  est[2]<-mean(fg)
  sd[2]<-sd(fg)
  res<-rbind(est=round(est,3), sd=round(sd,3))
  colnames(res)<-paste0('f',1:2)
 knitr::kable(res,digits = getOption("digits"),align = "lccrr") 


## -----------------------------------------------------------------------------
set.seed(12345)
####Now we will compute the probability that the confidence interval covers in example 6.4
## we use sample which are normal
n<-20
alpha <-.05
UCL<-replicate(1000, expr = {
x<-rnorm(n, mean = 0, sd = 2)
(n-1)*var(x)/qchisq(alpha,df = n-1)
} )
#or compute the mean to get the confidence level
print("the coverage probability of the t-interval for variance which sample are normal")
mean(UCL>4)

###we use sample which are non-normal
n<-20
alpha <-.05
UCL <- replicate(1000, expr = {
x<-rchisq(n,2)
(n-1)*var(x)/qchisq(alpha,df = n-1)
} )
#or compute the mean to get the confidence level
print("the coverage probability of the t-interval for variance which sample are non-normal")
mean(UCL>4)

####Now we will compute the probability that the confidence interval covers in this case
n<-20
CL<-replicate(1000, expr = {
   x<-rchisq(n,2) 
   sqrt(n)*(mean(x)-2)/2
} )
In<-mean((CL>(-qt(0.025,n-1,lower.tail=F)))&(CL<qt(0.025,n-1,lower.tail = F)))
print("the coverage probability of the t-interval for mean which sample are non-normal:")
print(In)

## -----------------------------------------------------------------------------
set.seed(1234)
library(knitr)
##(i)chi-square
n <-100
alpha <-.05
mu0<-1
sigma <-2
x_hat<-x_se <- numeric()
p.hat<-se.hat<-p_value<-Mc_Simulation<-numeric(3)
m<-10000 #number of replicates
p<-numeric(m) #storage for p-values
for(j in 1:m){
x<-rchisq(n,1)
ttest <- t.test(x, mu=mu0)
p[j]<-ttest$p.value
x_hat[j] <-mean(x)
x_se[j] <- sd(x)
}
p2<-2*(1-pt(abs(sqrt(n)*(x_hat-mu0)/x_se),df=n-1))
Mc_Simulation[1]<-mean(p2<0.05)#we obtain the tle by Monte Carlo simulation
p_value[1]<-mean(p)##p-value, help judge whether we should reject null hypothesis
p.hat[1]<-mean(p<alpha)#we obtain the tle by using t.test 
se.hat[1]<-sqrt(p.hat[1]*(1-p.hat[1])/m)#standard error

##(ii)Uniform(0,2)
n<-100
alpha <-.05
mu0<-1
sigma<-2
m<-10000 #number of replicates
p<-numeric(m) #storage for p-values
for(j in 1:m){
x<-runif(n,0,2)
ttest<-t.test(x,alternative = "greater", mu=mu0)
p[j]<-ttest$p.value
x_hat[j]<-mean(x)
x_se[j]<- sd(x)

}
p_value[2]<-mean(p)##p-value, help judge whether we should reject null hypothesis
p.hat[2]<-mean(p<alpha)#we obtain the tle by using t.test
se.hat[2]<-sqrt(p.hat[2]*(1-p.hat[2])/m)#standard error
p2<-2*(1-pt(abs(sqrt(n)*(x_hat-mu0)/x_se),df=n-1))
Mc_Simulation[2]<-mean(p2<0.05)#we obtain the tle by Monte Carlo simulation

##(iii)exp(1)
n <-100
alpha <-.05
mu0<-1
sigma <-2
m<-10000 #number of replicates
p<-numeric(m) #storage for p-values
for(j in 1:m) {
x<-rexp(n,1)
ttest <- t.test(x, mu=mu0)
p[j]<-ttest$p.value
x_hat[j] <-mean(x)
x_se[j] <- sd(x)
}
p_value[3]<-mean(p)##p-value, help judge whether we should reject null hypothesis
p.hat[3]<-mean(p<alpha)#we obtain the tle by using t.test
se.hat[3]<-sqrt(p.hat[3]*(1-p.hat[3])/m)#standard error
p2<-2*(1-pt(abs(sqrt(n)*(x_hat-mu0)/x_se),df=n-1))
Mc_Simulation[3]<-mean(p2<0.05)#we obtain the tle by Monte Carlo simulation

res<-rbind(Mc_Simulation=round(Mc_Simulation,4),p.hat=round(p.hat,4),se.hat=round(se.hat,4),p_value=round(p_value,4))
colnames(res)<-paste0(c("chi-square","Uniform","Exponential"))
kable(res,digits = getOption("digits"),align = "lccrr") 
#"chi-square","Uniform","Exponential"

## -----------------------------------------------------------------------------
n<-c(10,20,30,50,100,500) #sample sizes
d<-3
v1<-qchisq(0.95,d*(d+1)*(d+2)/6) ##Chi squared quantiles

## -----------------------------------------------------------------------------
### sk funcion  return 0 when we accept H0, return 1 when we reject H0
#computes the sample skewness coeff.
sk<-function(x) {
r<-nrow(x) ## we use r to note the sample size
x_center<-x # we will center x
for(i in 1:d)  
{
  x_center[,i]<-x[,i]-mean(x[,i])
}
## sigma_hat is the maximum likelihood estimator of covariance
sigma_hat<-cov(x)*(r-1)/r
b0<-x_center%*%solve(sigma_hat)%*%t(x_center)
b1<-sum(colSums(b0^3))/(r^2)
b<-r*b1/6 # Approximately follows the Chi-square distribution
as.integer(b>v1)##return 0(accept H0) or 1(reject H0)
}

## -----------------------------------------------------------------------------
#n is a vector of sample sizes
#we are doing length(n) different simulations
# set.seed(123)
# library(MASS)
# p.reject<-numeric(length(n))    #to store sim.results
# mean1<-c(0,0,0)
# sigma<-matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
# m<-1000 
#  for(i in 1:length(n)) {
#    p.reject[i]<-mean(replicate(m,expr ={
#      x<-mvrnorm(n[i],mean1,sigma)
#      sk(x)
#    })) #proportion rejected
#     }

## -----------------------------------------------------------------------------
# res<-rbind(p.reject=p.reject)
# colnames(res)<-paste0(n)
# kable(res,digits=getOption("digits"),align ="lccrr") 

## -----------------------------------------------------------------------------
alpha <-.1
n<-30
m<-2500
x<-matrix(0,nrow=n,ncol=3)
epsilon<-c(seq(0, 0.15, 0.01),seq(0.15, 1, 0.05))
N<-length(epsilon)
power<-numeric(N)
sigma1<-matrix(c(100,0,0,0,100,0,0,0,100),nrow=3,ncol=3)

## -----------------------------------------------------------------------------
#critical value for the skewness test
# set.seed(1234)
# for(j in 1:N) { #for each epsilon
# e<-epsilon[j] ###
# sktest<-numeric(m)
# for(i in 1:m)
# {
#   for(k in 1:n)
#   {
#     u<-runif(1)
#     if(u>e)
#       x[k,]<-mvrnorm(1,mean1,sigma)
#     else x[k,]<-mvrnorm(1,mean1,sigma1)
#   }
# sktest[i]<-sk(x)
# }
# power[j]<-mean(sktest)
# }

## -----------------------------------------------------------------------------
#plot power vs epsilon
# plot(epsilon,power,type ="b",col="blue",
# xlab =bquote(epsilon),ylim =c(0,1))
# abline(h =.1,lty =3)
# se<-sqrt(power*(1-power)/m) #add standard errors
# lines(epsilon,power+se,lty=3,col="red")
# lines(epsilon,power-se,lty=3,col="green")

## -----------------------------------------------------------------------------
##Bootstrap
library(knitr)
library(bootstrap)
library(boot)
set.seed(1234)
r<-nrow(scor)
sigma_hat1<-cov(scor)*(r-1)/r
theta_hat<-(eigen(sigma_hat1)$values[1])/sum(eigen(sigma_hat1)$values)
B<-2000

theta<-function(x,i){
data<-x[i,]
y<-cov(data)*(r-1)/r
theta<-eigen(y)$value[1]/sum(eigen(y)$value)
return(theta)
}
results<-boot(data=cbind(scor$mec,scor$vec,scor$alg,scor$ana,scor$sta),statistic = theta,R=B)
theta_b<-results$t
res<-rbind(c(theta_hat=theta_hat,theta_hat_star=mean(theta_b),boot.bias=mean(theta_b)-theta_hat,boot.se=sd(theta_b)))
kable(res,digits=getOption("digits"),align ="lccrr") 

## -----------------------------------------------------------------------------
#jackknife method
theta.jack<-numeric(r)
for(i in 1:r){
   data<-cov(scor[-i,])*(r-2)/(r-1)
   theta.jack[i]<-(eigen(data)$values[1])/sum(eigen(data)$values)
}
res<-rbind(c(theta_hat=theta_hat,theta_jack=mean(theta.jack),jack.bias=mean(theta.jack)-theta_hat,jack.se=sd(theta.jack)))
kable(res,digits=getOption("digits"),align ="lccrr") 

## -----------------------------------------------------------------------------
#  m<-1e2
# ci.perc<-ci.bca<-matrix(NA,m,2)
# for(i in 1:m){
# ci<-boot.ci(results,type=c("perc","bca"))
# ci.perc[i,]<-ci$percent[4:5]
# ci.bca[i,]<-ci$bca[4:5] }
# print("the 95% percentile confidence intervals for theta_hat")
# colMeans(ci.perc)
# print("the 95% BCa confidence intervals for theta_hat")
# colMeans(ci.bca)

## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)
library(knitr)
n<-1e1
m<-1e2
sk<-function(x,i) mean((x[i]-mean(x[i]))^3)/mean((x[i]-mean(x[i]))^2)^1.5 
#computes the sample skewness coeff.
ci.norm<-ci.basic<-ci.perc<-matrix(NA,m,2)
norm<-basic<-perc<-numeric(3)

## -----------------------------------------------------------------------------
####  normal populations  (skewness 0)
set.seed(100000)
sk1<-0
for(i in 1:m){
   x<-rnorm(n)
 result1<-boot(data=x,statistic=sk,R=999)
 ci<-boot.ci(result1,type=c("norm","basic","perc"))
 ci.norm[i,]<-ci$norm[2:3]
 ci.basic[i,]<-ci$basic[4:5]
 ci.perc[i,]<-ci$percent[4:5]
 #ci.bca[i,]<-ci$bca[4:5],"bca"
}
#empirical coverage rates for the sample skewness
norm[1]<-mean(ci.norm[,1]<=sk1& ci.norm[,2]>=sk1)
basic[1]<-mean(ci.basic[,1]<=sk1&ci.basic[,2]>=sk1)
perc[1]<-mean(ci.perc[,1]<=sk1& ci.perc[,2]>=sk1)
#BCa[1]<-mean(ci.bca[,1]<=sk1& ci.bca[,2]>=sk1)

##the proportion of times that the confidence intervals miss on the left
norm[2]<-sum(ci.norm[,1]>=sk1)/m
basic[2]<-sum(ci.basic[,1]>=sk1)/m
perc[2]<-sum(ci.perc[,1]>=sk1)/m
#BCa[2]<-sum(ci.bca[,2]>=sk1)/m

##the porportion of times that the confidence intervals miss on the right
norm[3]<-sum(ci.norm[,2]<=sk1)/m
basic[3]<-sum(ci.basic[,2]<=sk1)/m
perc[3]<-sum(ci.perc[,2]<=sk1)/m
#BCa[3]<-sum(ci.bca[,2]<=sk1)/m,BCa=round(BCa,3)

res<-rbind(norm=round(norm,3),basic=round(basic,3),perc=round(perc,3))
colnames(res)<-paste0(c("coverage probabilities","proportion_miss on the left","porportion_miss on the right"))
kable(res,digits = getOption("digits"),align = "lccrr")


## -----------------------------------------------------------------------------
####chi-square(5) skewness is sk=gamma(11/2)/gamma(5/2)
set.seed(12345)
sk2<-sqrt(8/5)
for(i in 1:m){
y<-rchisq(n,5)
result2<-boot(data=y,statistic=sk,R=999)
ci<-boot.ci(result2,type=c("norm","basic","perc"))
ci.norm[i,]<-ci$norm[2:3];ci.basic[i,]<-ci$basic[4:5]
ci.perc[i,]<-ci$percent[4:5] 
}

#empirical coverage rates for the sample skewness
norm[1]<-mean(ci.norm[,1]<=sk2& ci.norm[,2]>=sk2)
basic[1]<-mean(ci.basic[,1]<=sk2&ci.basic[,2]>=sk2)
perc[1]<-mean(ci.perc[,1]<=sk2& ci.perc[,2]>=sk2)

##the proportion of times that the confidence intervals miss on the left
norm[2]<-sum(ci.norm[,1]>=sk2)/m
basic[2]<-sum(ci.basic[,1]>=sk2)/m
perc[2]<-sum(ci.perc[,1]>=sk2)/m

##the porportion of times that the confidence intervals miss on the right
norm[3]<-sum(ci.norm[,2]<=sk2)/m
basic[3]<-sum(ci.basic[,2]<=sk2)/m
perc[3]<-sum(ci.perc[,2]<=sk2)/m

res<-rbind(norm=round(norm,3),basic=round(basic,3),perc=round(perc,3))
colnames(res)<-paste0(c("coverage probabilities","proportion_miss on the left","porportion_miss on the right"))
kable(res,digits = getOption("digits"),align = "lccrr")

## -----------------------------------------------------------------------------
### Spearman rank correlation test for  independence
set.seed(12345)
x<-rnorm(20,1,2)
y<-rnorm(20,3,4)
R<-1e3 #number of replicates
z<-c(x,y) #pooled sample
K<-1:40
reps<-numeric(R) #storage for replicates
c<-cor(x,y,method="spearman") # the statistic we need
for (i in 1:R) {
#generate indices k for the first sample
k<-sample(K,size=20,replace=FALSE)
x1<-z[k]
y1<-z[-k] #complement of x1
reps[i]<-cor(x1,y1,method="spearman")
}
##significance level as follow
p0<-mean(abs(c(c,reps))>=abs(c))
p0
####  p-value reported by cor.test
p<-cor.test(x,y,method = "spearman")
p$p.value
## The histogram of replicates  is shown as follow
hist(reps,main="",freq=FALSE,xlab="T(p = 0.164)",
breaks="scott")
points(c,0,cex=1,col=2,pch =16)

## -----------------------------------------------------------------------------
### NN
library(RANN) # implementing a fast algorithm
library(boot)
library(energy)
library(Ball)
library(knitr)

## -----------------------------------------------------------------------------
### power compare
m<-1e3; k<-3; p<-2; 
R<-999;
### Tn is the kth nearest neighbor statistic
Tn<-function(z,ix,sizes,k) 
  {
   n1<-sizes[1]  
   n2<-sizes[2]
   n<-n1+n2
   if(is.vector(z))z<-data.frame(z)
   z<-z[ix,]
   NN<-nn2(data=z, k=k+1) 
   block1<-NN$nn.idx[1:n1,-1]
   block2<-NN$nn.idx[(n1+1):n,-1]
   i1<-sum(block1<=n1)
   i2<-sum(block2>n1)
   (i1+i2)/(k*n)
  }

### we can obtain the p-value of NN test from the eqdist.nn
eqdist.nn<-function(z,sizes,k)
  {
   boot.obj<-boot(data=z,statistic=Tn,R=R,sim="permutation",sizes=sizes,k=k)
   ts<-c(boot.obj$t0,boot.obj$t)
   p.value<-mean(ts>=ts[1])
   list(statistic=ts[1],p.value=p.value)
 }

### we compute the NN, energy, and ball's power when the random varible are normal
power_compare<-function(n1,n2,mu1,mu2,sigma1,sigma2){
   p_values<-matrix(NA,m,3)
   N<-c(n1,n2)
   for(i in 1:m)
  {
    x<-matrix(rnorm(n1*p,mu1,sigma1),ncol=p)
    y<-matrix(rnorm(n2*p,mu2,sigma2),ncol=p)
    z<-rbind(x,y)
    p_values[i,1]<-eqdist.nn(z,N,k)$p.value
    p_values[i,2]<-eqdist.etest(z,sizes=N,R=R)$p.value
    p_values[i,3]<-bd.test(x=x,y=y,num.permutations=R,seed=i*12345)$p.value
  }
  alpha<-0.05
  pow<-colMeans(p_values<alpha)
  pow ##the power vector
}
#we compute the NN, energy, and ball's power when the random varible are normal and student.
power_compare1<-function(n1,n2){
   p_values<-matrix(NA,m,3)
   N<-c(n1,n2)
  for(i in 1:m)
  {
    x<-matrix(rt(n1*p,1),ncol=p);
    y<-cbind(rnorm(n2),rnorm(n2,mean=1))
    z<-rbind(x,y)
    p_values[i,1]<-eqdist.nn(z,N,k)$p.value
    p_values[i,2]<-eqdist.etest(z,sizes=N,R=R)$p.value
    p_values[i,3]<-bd.test(x=x,y=y,num.permutations=R,seed=i*12345)$p.value
  }
  alpha<-0.05
  pow<-colMeans(p_values<alpha)
  pow ##the power vector
}

## -----------------------------------------------------------------------------
## we use vector pow1,2,3,4 denote the power of three test under four cases
## In the following four cases, n1 and n2 were uniformly used to record the sample size of two samples, mu1 and mu2 were used to record the mean values of two normal samples, and sigma1 and sigma2 were used to record the standard deviations of two normal samples.

# pow1<-pow2<-pow3<-pow4<-numeric(3)  
# ###(1)Unequal variances and equal expectations
# set.seed(12345)
# n1<-n2<-50; mu1<-mu2<-0; sigma1<-0.1; sigma2<-0.17;
# pow1<-power_compare(n1,n2,mu1,mu2,sigma1,sigma2)
# 
# ###(2)Unequal variances and unequal expectations
# n1<-n2<-50; mu1<-0.15; mu2<-0.1; sigma1<-0.11; sigma2<-0.12
# pow2<-power_compare(n1,n2,mu1,mu2,sigma1,sigma2)
# 
# ###(3)Non-normal distributions: t distribution with 1 df (heavy-tailed distribution), bimodel distribution (mixture of two normal distributions)
# n1<-n2<-30 ;
# pow3<-power_compare1(n1,n2)
# 
# ###(4)Unbalanced samples (say, 1 case versus 10 controls)
# n1<-10; n2<-100 ; mu1<-mu2<-0; sigma1<-1; sigma2<-3
# pow4<-power_compare(n1,n2,mu1,mu2,sigma1,sigma2)
# 
# res<-rbind(pow1=round(pow1,4),pow2=round(pow2,4),pow3=round(pow3,4),pow4=round(pow4,4))
# colnames(res)<-paste0(c("NN","energy","Ball"))
# kable(res,digits = getOption("digits"),align = "lccrr")

## -----------------------------------------------------------------------------
library(knitr)
set.seed(23466)
#####
y<-rnorm(1)
#####
m<-1e2
x<-numeric(m)   ### Markov chains
x[1]<-rnorm(1)
k<-0  #### The chain length
u<-runif(m)
 for (i in 2:m) {
            y<-rnorm(1,x[i-1],5)
                if(u[i] <= (dt(y,1)/dt(x[i-1],1)))
                    x[i] <- y  
                else {
                    x[i] <- x[i-1]
                    k <-k+1
                } }
index<-11:100
    y1<-x[index]
Q_sample<-quantile(y1,0.1)
a<-ppoints(100)
Q_Cauchy<-qt(0.1,1)
Quantile<-data.frame(sample_quantile=Q_sample,Cauchy_quantile=Q_Cauchy)
  kable(Quantile)


## -----------------------------------------------------------------------------
a<-ppoints(100)
Q_Cauchy1<-qt(a,1)
Q_sample1<-quantile(y1,a)
par(mfrow=c(1,2))
    qqplot(Q_Cauchy1,Q_sample1, main="",
        xlab="Cauchy Quantiles", ylab="Sample Quantiles")
    abline(0,1,col='blue',lwd=2)
    hist(y1, breaks="scott", main="", xlab="",ylim=c(0,0.4),xlim=c(-20,20),freq=FALSE)
    lines(Q_Cauchy1,dt(Q_Cauchy1,1))

## -----------------------------------------------------------------------------
#initialize constants and parameters
set.seed(11093)
N<-5000 #length of chain
burn<-1000 #burn-in length
X<-matrix(0, N, 2) #the chain, a bivariate sample
a<-2
b<-3
n<-100
###### generate the chain #####
X[1,]<-c(1, 0.1) #initialize
for(i in 2:N) {
y<-X[i-1,2]
X[i,1]<-rbinom(1,n,y)
x<-X[i,1]
X[i,2]<-rbeta(1,x+a,n-x+b)
}
b<-burn+1
M<-X[b:N,]
# compare sample statistics to parameters
colMeans(M)
plot(M, main="", cex=.5, xlab=bquote(X),ylab=bquote(Y), ylim=range(M[,2]))

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
                                # psi[i,j] is the statistic psi(X[i,1:j])
                               # for chain in i-th row of X
    psi<-as.matrix(psi)
    n<-ncol(psi)
    k<-nrow(psi)
    psi.means<-rowMeans(psi)    #row means
    B<-n*var(psi.means)         #between variance est.
    psi.w<-apply(psi,1,"var")   #within variances
    W<-mean(psi.w)              #within est.
    v.hat<-W*(n-1)/n+(B/n)      #upper variance est.
    r.hat<-v.hat/W              #G-R statistic
    return(r.hat)
}

## -----------------------------------------------------------------------------
Cauchy.chain <- function(N,X1) 
  {
 #with Normal(X[t], 5) proposal distribution
 #and starting value X1
    x<-rep(0,N)
    x[1]<-X1
    u<-runif(N)
 for (i in 2:1000) {
            y <- rnorm(1,x[i-1],5)
                if (u[i]<=(dt(y,1)/dt(x[i-1],1)))
                    x[i] <-y  
                else {
                    x[i] <- x[i-1]
                    k<-k+1 } }
    return(x)
}

## -----------------------------------------------------------------------------
k<-4 #number of chains to generate
n<-2000 #length of chains
bu<-500 #burn-in length
#choose overdispersed initial values
x0<-c(-2, -1, 0.1, 0.5)
#generate the chains
set.seed(123455)
X<-matrix(0,nrow=k, ncol=n)
for (i in 1:k)
X[i,]<-Cauchy.chain(n, x0[i])
#compute diagnostic statistics
psi<-t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,]<-psi[i,]/(1:ncol(psi))
print(Gelman.Rubin(psi))
#plot psi for the four chains
# par(mfrow=c(2,2))
# for (i in 1:k)
# plot(psi[i,(bu+1):n], type="l",
# xlab=i, ylab=bquote(psi))

rhat <- rep(0, n)
for(j in (bu+1):n)
rhat[j]<-Gelman.Rubin(psi[,1:j])
plot(rhat[(bu+1):n],type="l",xlab="",ylab="R")
abline(h=1.2, lty=2,col=2)

## -----------------------------------------------------------------------------
k<-4 #number of chains to generate
N<-1500 #length of chains
bu<-500 #burn-in length
#generate the chains
MC.chain<-function(N,x0){
X<-matrix(0, N, 2) #the chain, a bivariate sample
a<-2
b<-3
n<-100
###### generate the chain #####
X[1,]<-x0 #initialize
for(i in 2:N) {
y<-X[i-1,2]
X[i,1]<-rbinom(1,n,y)
x<-X[i,1]
X[i,2]<-rbeta(1,x+a,n-x+b)
}
return(X)
}

## -----------------------------------------------------------------------------
set.seed(12345)
x0<-matrix(c(1,0.1,5,0.5),ncol=2)
X1<-MC.chain(N,x0[,1])
#X1
X2<-MC.chain(N,x0[,2])
Y<-cbind(X1,X2)
#compute diagnostic statistics
psi1<-apply(Y, 1, cumsum)
for (i in 1:nrow(psi1))
psi1[i,]<-psi1[i,]/(1:ncol(psi1))
print(Gelman.Rubin(psi1))

## -----------------------------------------------------------------------------
#plot psi for the four chains
# par(mfrow=c(2,2))
# for (i in 1:4)
# plot(psi1[i,(bu+1):N], type="l",
# xlab=i)
# par(mfrow=c(1,1)) #restore default
#plot the sequence of R-hat statistics
rhat <- rep(0, N)
for(j in (bu+1):N)
 rhat[j]<-Gelman.Rubin(psi1[,1:j])
plot(rhat[(bu+1):N],type="l",xlab="",ylab="R")
abline(h=1.2, lty=2,col=2)
abline(h=1.1, lty=2,col=3)

## -----------------------------------------------------------------------------
##(a)
### the kth term,k!=gamma(k+1)
f<-function(k,a,d)
  {(-1)^k*exp(((2*k+2)*log(norm(a,type="F"))+lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+1)-k*log(2)-log(2*k+1)-log(2*k+2)-lgamma(k+d/2+1)) )
}
k<-8
d<-3
a<-as.matrix(c(1,2,3))
f(k,a,d)


## -----------------------------------------------------------------------------
###(b)
# K<-0:1e5
# sum(f(k,a,d))
summ<-function(a,d)
{
  k=0
  x<-y<-f(k,a,d)
  while(x>1e-5){
   k=k+1
   x<-f(k,a,d)
   y<-y+x
  }
  y
}

## -----------------------------------------------------------------------------
### (c)
a<-matrix(c(1,2))
d<-2
summ(a,d)

## -----------------------------------------------------------------------------
k<-c(4:25,100,500,1000)
s<-function(k,a){
  return(pt(sqrt(a^2*k/(k+1-a^2)),df=k,lower.tail = F))
}
delta.s <- function(k,a){
  s(k,a)-s(k-1,a)
}
n<-length(k)
A<-numeric(n)
for(i in 1:n){
  A[i]<-uniroot(delta.s,interval=c(1,2),k=k[i])$root
}
result<-cbind(k,A)
knitr::kable(result)

## -----------------------------------------------------------------------------
#11.5
f1<-function(u)(1+u^2/(k-1))^(-k/2)
f2<-function(u)(1+u^2/(k))^(-(k+1)/2)
m<-function(k)exp(log(2)+lgamma(k/2)-0.5*(log(pi)+log(k-1))-lgamma((k-1)/2))
fu<-function(k){
  f<-function(a){
  m(k)*integrate(f1,lower=0,upper=sqrt(a^2*(k-1)/(k-a^2)),rel.tol=.Machine$double.eps^0.25)$value-m(k+1)*integrate(f2,lower=0,upper=sqrt(a^2*k/(k+1-a^2)),rel.tol=.Machine$double.eps^0.25)$value
}
uniroot(f,c(-1e-7,1.99))$root
}
res<-numeric(25)
for(k in 4:25)
  res[k-3]<-fu(k)

## -----------------------------------------------------------------------------
f<-function(a){
  m(k)*integrate(f1,lower=0,upper=sqrt(a^2*(k-1)/(k-a^2)),rel.tol=.Machine$double.eps^0.25)$value-m(k+1)*integrate(f2,lower=0,upper=sqrt(a^2*k/(k+1-a^2)),rel.tol=.Machine$double.eps^0.25)$value
}
k<-100
res[23]<-uniroot(f,c(-1e-7,1.99))$root
k<-500
res[24]<-uniroot(f,c(-1e-7,1.99))$root
k<-1000
res[25]<-uniroot(f,c(-1e-7,1.99))$root
res

## -----------------------------------------------------------------------------
trims<-c(0, 0.1, 0.2, 0.5)
x<-rcauchy(100)
lapply(trims,function(trim)mean(x,trim=trim))
lapply(trims, mean, x=x)


## -----------------------------------------------------------------------------
##(5)
rsq<-function(mod){summary(mod)$r.squared}
###3. Use both for loops and lapply() to fit linear models to the
#mtcars using the formulas stored in this list:
data1<-mtcars
formulas<-list(
 mpg~disp,
 mpg~I(1/disp),
 mpg~disp+wt,
 mpg~I(1/disp)+wt
)
l1<-lapply(formulas,lm,data=data1)
print("In Exercise 3, the R^2 for each model is below")
lapply(l1,rsq)

### 4. Fit the model mpg~disp to each of the bootstrap replicatesof mtcars in #the list below by using a for loop and lapply(). Can you do it without an #anonymous function?
bootstraps<-lapply(1:10,function(i){
    rows<-sample(1:nrow(mtcars),rep=TRUE)
    mtcars[rows,]
  })
l2<-lapply(bootstraps,lm,formula=mpg~disp)
print("In  Exercise 4, the R^2 for each model is below")
lapply(l2, rsq)


## -----------------------------------------------------------------------------
#(a)the standard deviation of every column in a numeric data frame.
vapply(data1,sd,FUN.VALUE=c(sd=0))

## -----------------------------------------------------------------------------
#(b)the standard deviation of every numeric column in a mixed data frame.
data2<-esoph
summary1<-function(x){
vapply(x,as.numeric,FUN.VALUE=c(value1=0))
}
round(vapply(data2,function(x)(sd(summary1(x))),FUN.VALUE=c(sd=0)),3)

## -----------------------------------------------------------------------------
library(parallel)
n<-10
cl<-makeCluster(mc <-getOption("cl.cores", 2))
parLapply(cl,1:n,sqrt)

#### we write the mcsapply based parLapply
mcsapply<-function(cl=NULL,X, fun, ...) {
  
     res<-parLapply(cl,X, fun, ...)
     simplify2array(res)
}
mcvapply<-function(cl,X, fun, fun.value,...){
    
     fun1<-function(x)x
     vapply(parLapply(cl,X,fun),fun1,FUN.VALUE=fun.value,...)
}
mcsapply(cl,1:n,sqrt)

## ----eval=FALSE---------------------------------------------------------------
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  //[[Rcpp::export]]
#  NumericMatrix gibbsC(int a,int b,int n) {
#    int N = 10000;
#    NumericMatrix X(N,2);
#    double x=0,y=0;
#    X(0,0)=0;
#    X(0,1)=0.3;
#    for(int i=1; i<N; ++i) {
#      y=X(i-1,1);
#      X(i,0)=rbinom(1,n,y)[0];
#      x=X(i,0);
#      X(i,1)=rbeta(1,x+a,n-x+b)[0];
#    }
#    return (X);
#  }

## -----------------------------------------------------------------------------
library(Rcpp)
library(StatComp21027)
 m<-gibbsC(2,3,100)
 plot(m[,1],m[,2])

## -----------------------------------------------------------------------------
    set.seed(12)
    x_C<-gibbsC(2,3,100)
    y_R<-gibbsR(2,3,100)
## we will draw the qqplot
c<-ppoints(100)
QC1<-quantile(x_C[,1],c)
QR1<-quantile(y_R[,1],c)
QC2<-quantile(x_C[,2],c)
QR2<-quantile(y_R[,2],c)
 qqplot(QC1,QR1,main="",xlab="QC1", ylab="QR1")
 abline(0,1,col = "red")
 qqplot(QC2,QR2,main="",xlab="QC2", ylab="QR2")
 abline(0,1,col = "blue")

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(gibbR=gibbsR(2,3,100),
gibbC=gibbsC(2,3,100))
summary(ts)[,c(1,3,5,6)]

