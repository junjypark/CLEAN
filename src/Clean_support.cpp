
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void set_seed(unsigned int seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(seed);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List CleanMeanC(arma::mat& ymat, arma::sp_mat& NNmatrix, int nperm, int s){
  int q=NNmatrix.n_rows;
  int n=ymat.n_cols;
  arma::vec onevec(n); onevec.fill(1);
  arma::vec U(q);
  double sd;
  double mn;
  
  arma::mat permU(q, nperm); 
  
  U=NNmatrix*ymat*onevec;  

  arma::mat flipmat(n,nperm); 
  set_seed(s);
  flipmat.randn(); 
  flipmat=sign(flipmat);
  permU=NNmatrix*ymat*flipmat;

  for (int k=0; k<q; ++k){
    mn=mean(permU.row(k));
    sd=stddev(permU.row(k));
    permU.row(k)=(permU.row(k)-mn)/sd;
    U(k)=(U(k)-mn)/sd;
  }
  
  arma::vec permMax(nperm);
  arma::vec permMin(nperm);
    
  for (int i=0; i<nperm; ++i){
    permMin(i)=permU.col(i).min();
    permMax(i)=permU.col(i).max();
  }
  
  return Rcpp::List::create(Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("permMin")=permMin,
                            Rcpp::Named("nperm")=nperm);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List CleanDiffC(arma::mat& ymat, arma::sp_mat& NNmatrix, arma::vec group, int nperm, int s){
  int q=NNmatrix.n_rows;
  int n=group.size();
  arma::vec U(q);
  double sd;
  double mn;

  arma::mat permU(q, nperm);

  U=NNmatrix*ymat*group;  
  
  arma::mat permmat(n,nperm); permmat.fill(0);
  set_seed(s);
  for (int i=0; i<nperm; ++i){
    permmat.col(i)=shuffle(group);
  }

  permU=NNmatrix*ymat*permmat;
    
  for (int k=0; k<q; ++k){
    mn=mean(permU.row(k));
    sd=stddev(permU.row(k));
    permU.row(k)=(permU.row(k)-mn)/sd;
    U(k)=(U(k)-mn)/sd;
  }

  arma::vec permMax(nperm);
  arma::vec permMin(nperm);
    
  for (int i=0; i<nperm; ++i){
    permMin(i)=permU.col(i).min();
    permMax(i)=permU.col(i).max();
  }
  
  return Rcpp::List::create(Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("permMin")=permMin,
                            Rcpp::Named("nperm")=nperm);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List MassiveMeanC(arma::mat ymat, int nperm,  int s){
  int q=ymat.n_rows;
  int n=ymat.n_cols;
  arma::vec onevec(n); onevec.fill(1);
  arma::vec U(q);
  double sd;
  double mn;  

  arma::mat permU(q, nperm); 
  
  U=ymat*onevec;  

  arma::mat flipmat(n,nperm); 
  set_seed(s);
  flipmat.randn(); 
  flipmat=sign(flipmat);
  permU=ymat*flipmat;
  
  for (int k=0; k<q; ++k){
    mn=mean(permU.row(k));
    sd=stddev(permU.row(k));
    permU.row(k)=(permU.row(k)-mn)/sd;
    U(k)=(U(k)-mn)/sd;
  }
  
  arma::vec permMax(nperm);
  arma::vec permMin(nperm);
    
  for (int i=0; i<nperm; ++i){
    permMin(i)=permU.col(i).min();
    permMax(i)=permU.col(i).max();
  }
  
  return Rcpp::List::create(Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("permMin")=permMin,
                            Rcpp::Named("nperm")=nperm);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List MassiveDiffC(arma::mat ymat, arma::vec group, int nperm,  int s){
  int q=ymat.n_rows;
  int n=group.size();
  arma::vec U(q);
  double sd;
  double mn;
  arma::mat permU(q, nperm);

  U=ymat*group;  

  arma::mat permmat(n,nperm);
  set_seed(s);
  for (int i=0; i<nperm; ++i){
    permmat.col(i)=shuffle(group);
  }
  permU=ymat*permmat;

  for (int k=0; k<q; ++k){
    mn=mean(permU.row(k));
    sd=stddev(permU.row(k));
    permU.row(k)=(permU.row(k)-mn)/sd;
    U(k)=(U(k)-mn)/sd;
  }

  arma::vec permMax(nperm);
  arma::vec permMin(nperm);
    
  for (int i=0; i<nperm; ++i){
    permMin(i)=permU.col(i).min();
    permMax(i)=permU.col(i).max();
  }
  
  return Rcpp::List::create(Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("permMin")=permMin,
                            Rcpp::Named("nperm")=nperm);
}

// [[Rcpp::export]]
double computetraceABA(arma::mat A, arma::mat B){
  arma::mat BA=B*A;
  double out=accu(A%BA);
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double CovRegOptimC(double rho, arma::mat& epsilon, arma::mat corMat_base){
  int p=epsilon.n_rows;
  int n=epsilon.n_cols;
  double sigma2;
  double tau2;
  double ss;

  arma::mat corMat = pow(corMat_base, rho);
  double corMat_norm=pow(norm(corMat, "fro"),2);
  if (corMat_norm>p+pow(10,-10)){
    double y1=computetraceABA(epsilon, corMat)/n;
    double y2=pow(norm(epsilon,"fro"),2)/n;
    sigma2= (y1-y2)/(corMat_norm-p);
    tau2= (-p*y1+corMat_norm*y2)/(corMat_norm*p-p*p);
    ss= -2*(sigma2*y1*n+tau2*y2*n)+p*tau2*tau2+corMat_norm*sigma2*sigma2+2*sigma2*tau2*p;
  } else{
    sigma2=-1;
    tau2=-1;
    ss=pow(10, 10);
  }
  return ss;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List ObtainVarCompsC(double rho, arma::mat& epsilon, arma::mat corMat_base){
  int p=epsilon.n_rows;
  int n=epsilon.n_cols;
  double sigma2;
  double tau2;
  arma::mat corMat = pow(corMat_base, rho);
  double corMat_norm=pow(norm(corMat, "fro"),2);
  if (corMat_norm>p+pow(10,-10)){
    double y1=computetraceABA(epsilon, corMat)/n;
    double y2=pow(norm(epsilon,"fro"),2)/n;
    sigma2= (y1-y2)/(corMat_norm-p);
    tau2= (-p*y1+corMat_norm*y2)/(corMat_norm*p-p*p);
  } else{
    sigma2=-1;
    tau2=-1;
  }
  return Rcpp::List::create(Rcpp::Named("sigma2")=sigma2,
                            Rcpp::Named("tau2")=tau2);
}

