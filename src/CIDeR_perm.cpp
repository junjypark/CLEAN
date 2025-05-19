#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

// [[Rcpp::export]]
// ymat is V*N
// xmat is N*p
Rcpp::List Ciderperm(arma::mat& ymat,  arma::mat& xmat, arma::sp_mat& NNmatrix, int nperm, int s){
  int q=NNmatrix.n_rows;
  int n=ymat.n_cols;
  int V=ymat.n_rows;
  arma::vec onevec(n); onevec.fill(1);
  arma::vec U(q);
  arma::mat T(V,1);
  arma::mat ymattemp(V,n);
  double sd;
  double mn;
  arma::mat permU(q, nperm); 
  
  arma::vec xsum=arma::sum(xmat, 1);
  T=ymat*xsum;
  U=NNmatrix*T;  
  
  //Rcout << "perm";
  
  arma::mat permT(V, nperm);
  
  set_seed(s);
  
  for (int i=0; i < nperm; ++i) {
    ymattemp=shuffle(ymat,1);
    permT.col(i)=ymattemp*xsum;
  }
  
  permU=NNmatrix*permT;
  
  //Rcout << "stand";
  for (int k=0; k<q; ++k){
    mn=mean(permU.row(k));
    sd=stddev(permU.row(k));
    permU.row(k)=(permU.row(k)-mn)/sd;
    U(k)=(U(k)-mn)/sd;
  }
  
  arma::vec permMax(nperm);
  arma::vec permMin(nperm);
  //Rcout << "max" <<endl;
  for (int i=0; i<nperm; ++i){
    permMax(i)=permU.col(i).max();
    permMin(i)=permU.col(i).min();
  }
  
  return Rcpp::List::create(Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("permMin")=permMin,
                            Rcpp::Named("nperm")=nperm);
}