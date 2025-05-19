// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double CovRegOptim2C(NumericVector phi, const arma::mat& epsilon,
                     const arma::sp_mat& corMat_base1,
                     const arma::sp_mat& corMat_base2) {
  
  double phi1 = phi[0];
  double phi2 = phi[1];
  
  int p = epsilon.n_rows;
  int n = epsilon.n_cols;
  
  arma::sp_mat corMat1 = corMat_base1;
  arma::sp_mat corMat2 = corMat_base2;
  
  corMat1.for_each([=](double& val) { val = std::pow(val, phi1); });
  corMat2.for_each([=](double& val) { val = std::pow(val, phi2); });
  
  
  double corMat1_norm  = accu(corMat1 % corMat1);
  double corMat2_norm  = accu(corMat2 % corMat2);
  double corMat12_norm = accu(corMat1 % corMat2);
  
  arma::mat M(3, 3);
  M(0, 0) = corMat1_norm;
  M(0, 1) = corMat12_norm;
  M(0, 2) = p;
  M(1, 0) = corMat12_norm;
  M(1, 1) = corMat2_norm;
  M(1, 2) = p;
  M(2, 0) = p;
  M(2, 1) = p;
  M(2, 2) = p;
  
  double ss;
  
  if (rcond(M) > 1e-12) {
    arma::mat cor1_eps = corMat1 * epsilon;
    arma::mat cor2_eps = corMat2 * epsilon;
    
    double y1 = accu(epsilon % cor1_eps) / n;
    double y2 = accu(epsilon % cor2_eps) / n;
    double y3 = accu(epsilon % epsilon) / n;
    
    arma::vec y = {y1, y2, y3};
    arma::vec theta = solve(M, y);
    
    double sigma1 = theta(0);
    double sigma2 = theta(1);
    double tau2   = theta(2);
    
    ss = sigma1 * sigma1 * corMat1_norm +
      sigma2 * sigma2 * corMat2_norm +
      tau2 * tau2 * p +
      2 * (sigma1 * sigma2 * corMat12_norm +
      sigma1 * tau2 * p + sigma2 * tau2 * p) -
      2 * (sigma1 * y1 * n + sigma2 * y2 * n + tau2 * y3 * n);
  } else {
    ss = 1e10;
  }
  
  //Rcpp::Rcout << "phi = (" << phi1 << ", " << phi2 << "), loss = " << ss << std::endl;
  return ss;
}
