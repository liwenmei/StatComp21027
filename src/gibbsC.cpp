#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param a the component elements of the first parameter of beta
//' @param b the component elements of the two parameter of beta
//' @param n the parameter of the binomial distribution
//' @return a random sample 
//' @examples
//' \dontrun{
//' rnC<-gibbsC(2,3,100)
//' plot(rnC[,1],type='l')
//' plot(rnC[,2],type='l')
//' }
//' @export
//[[Rcpp::export]]
NumericMatrix gibbsC(int a,int b,int n) {
  int N=10000;
  NumericMatrix X(N,2);
  double x=0,y=0;
  X(0,0)=0;
  X(0,1)=0.3;
  for(int i=1;i<N;++i) {
    y=X(i-1,1);
    X(i,0)=rbinom(1,n,y)[0];
    x=X(i,0);
    X(i,1)=rbeta(1,x+a,n-x+b)[0];
  }
  return (X);
} 
