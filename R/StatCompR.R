#' @title the statement of R packages we use
#' @name statement
#' @description Use R package \code{microbenchmark} to compare the performance of C functions (\code{gibbsR} and \code{vaccR}) and Cpp functions (\code{gibbsC} and \code{vaccC}).
#' @examples
#' \dontrun{
#' tm1<-microbenchmark::microbenchmark(
#'   rnR = gibbsR(2,3,100),
#'   rnC = gibbsC(2,3,100)
#' )
#' print(summary(tm1)[,c(1,3,5,6)])
#' }
#' @import parallel
#' @import MASS
#' @import qpdf
#' @import nnet 
#' @import compiler
#' @import microbenchmark
#' @import stats
#' @import knitr
#' @import ggplot2
#' @import bootstrap
#' @import boot
#' @import RANN
#' @import energy
#' @import Ball
#' @import Rcpp
#' @useDynLib StatComp21027
NULL

selr.part0<-function(n){
  n0<-n[,1] 
  n1<-n[,2] 
  n2<-n[,3] 
  nu<-n0/(n0+n1+n2)
  omega<-n1/(n0+n1+n2)
  R0.alt<-sum(n0*log(nu+1e-50))+sum(n1*log(omega +1e-50))+sum(n2*log(1-nu-omega+1e-50))
  nu.nul<-sum(n0)/(sum(n0)+sum(n1)+sum(n2))
  omega.nul<-sum(n1)/(sum(n0)+sum(n1)+sum(n2))
  R0.nul<-sum(n0)*log(nu.nul+1e-50)+sum(n1)*log(omega.nul+1e-50)+ sum(n2)*log(1-nu.nul-omega.nul+1e-50)
  R0<-2*(R0.alt-R0.nul)
  R0
}

selr.part1<-function(x,n){
  x[,1]<-factor(x[,1])
  group1<-unique(x[,1])
  m<-length(group1)
  x[,2]<-factor(x[,2])
  group2<-unique(x[,2])
  n<-c() 
  for(i in 1:m){
    n<-c(n, length(x[x[,1]==group1[i],2]))
  }
  rho<-n/sum(n)
  x1=x[,3]
  x2=log(x1)
  x3=log(1-x1)
  x4=log(x1/(1-x1))
  x5=x4^2
  result1<-summary(multinom(x[,1]~x1,trace=F))
  result2<-summary(multinom(x[,1]~x2,trace=F)) 
  result3<-summary(multinom(x[,1]~x3,trace=F))
  result4<-summary(multinom(x[,1]~x2+x3,trace=F)) 
  result5<-summary(multinom(x[,1]~x4,trace=F))
  result6<-summary(multinom(x[,1]~x5,trace=F))
  result7<-summary(multinom(x[,1]~x4+x5,trace=F)) 
  results<-c(result1$value,result2$value,result3$value,result4$value,result5$value,result6$value,result7$value)
  loglik<-(-results-sum(n*log(rho)))
  R1<-2*loglik
  R1
}

#' @title Composite semiparametric empirical likelihood ratio 
#' @name CSELR
#' @description compute the composite semiparametric empirical likelihood value
#' @param dat the data
#' @return the value of Composite semiparametric empirical likelihood ratio 
#' @examples
#' \dontrun{
#' s<-selr(data)
#' }
#' @export
selr<-function(dat){
  dat[,1]<-factor(dat[,1])
  group1<-unique(dat[,1])
  dat[,2]<-factor(dat[,2])
  m<-length(group1)
  n<-numeric(m) 
  ncount=c()
  for(i in 1:m){
    x=dat[dat[,1]==group1[i],3]
    n0=sum(x==0)
    n1=sum(x==1)
    n01=sum(x>0&x<1)
    ncount=rbind(ncount,c(n0,n1,n01))
  }
  dat.continuous<-dat[dat[,3]>0&dat[,3]<1,]
  part0<-selr.part0(ncount)
  part1<-selr.part1(dat.continuous)
  list(selrt=part0+part1)
}

#' @title Use bootstrap method to compute the p_vaule of Composite semiparametric empirical likelihood ratio 
#' @name boot_CSELR
#' @description use bootstrap method to compute the p value
#' @param dat the data
#' @param B is the number of bootstrap replication, with default B=999
#' @param n.cores is the number of cores used for parelell computing, with default n.cores=1
#' @return the p_value of Composite semiparametric empirical likelihood ratio 
#' @examples
#' \dontrun{
#' bs<-boot.selr(data)
#' }
#' @export
boot.selr<-function(dat,B=999,n.cores=1)
{
  arguments<-as.list(match.call())
  dat[,1]<-factor(dat[,1])
  group<-unique(dat[,1])
  m<-length(group)
  dat[,2]<-factor(dat[,2])
  group2<-unique(dat[,2])
  N<-nrow(dat)
  output.list<-mclapply(seq(1:B),function(i)
  { 
    y<-sample(dat[,3], size =N, replace = TRUE)
    boot.sample<-data.frame(dat[,1:2],y)
    selr.boot<-selr(boot.sample)$selrt
    return(selr.boot)
  }, mc.cores = n.cores)
  output <- unlist(output.list)
  boot.Rn <- matrix(output, 7, B, byrow = FALSE)
  obs.teststat <- selr(dat)$selrt
  pvalues <- apply(cbind(obs.teststat, boot.Rn), 1, function(y){ mean(y[-1]>y[1]) })
  res <- cbind(obs.teststat, pvalues)
  rnames <- c("qx=x", "qx=logx", "qx=log(1-x)", "qx=(logx,log(1-x))", "qx=log(x/1-x)",
              "qx=log(x/1-x)^2", "qx=(log(x/1-x),log(x/1-x)^2")
  cnames <- c("obs.SELR", "boot.pvalue")
  dimnames(res)<-list(rnames,cnames)
  list(boot.SELR.test=res)
}

#' @title A Gibbs sampler using R
#' @name GibbsR
#' @description A Gibbs sampler using R
#' @param n the parameter of the binomial distribution
#' @param a the component elements of the first parameter of beta
#' @param b the component elements of the two parameter of beta
#' @return a random sample  
#' @examples
#' \dontrun{
#' rnR <- gibbsR(100,2,3)
#' par(mfrow=c(2,1))
#' plot(rnR[,1],type='l')
#' plot(rnR[,2],type='l')
#' }
#' @export
gibbsR<-function(a,b,n){
  N<-5000 #length of chain
  X<-matrix(0, N, 2) #the chain, a bivariate sample
  X[1,]<-c(0, 0.3) #initialize
  for(i in 2:N) {
    y<-X[i-1,2]
    X[i,1]<-rbinom(1,n,y)
    x<-X[i,1]
    X[i,2]<-rbeta(1,x+a,n-x+b)
  }
  return (X)
}
