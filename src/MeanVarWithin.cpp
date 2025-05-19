#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat MeanVarWithin(const arma::mat& m, const arma::sp_mat& NNmatrix) {
  int N = m.n_rows;  // number of subjects
  int V = m.n_cols;  // number of vertices
  
  arma::mat m_adjusted(N, V, arma::fill::zeros);
  
  for (int v = 0; v < V; ++v) {
    arma::uvec neighbors = arma::find(NNmatrix.row(v));
    int K = neighbors.n_elem;
    
    if (K == 1) {
      m_adjusted.col(v) = m.col(neighbors[0]);
      continue;
    }
    
    arma::mat m_neighbors = m.cols(neighbors);
    arma::mat row_mean = arma::mean(m_neighbors, 1);  // N x 1
    arma::mat centered = m_neighbors.each_col() - row_mean;
    
    // Fixed: compute row-wise sample variance manually
    arma::vec row_var = arma::sum(centered % centered, 1) / (K - 1.0);
    row_var.transform([](double val) { return val == 0.0 ? 1.0 : val; });
    
    arma::mat m_std = (m.col(v) - row_mean) / arma::sqrt(row_var);
    m_adjusted.col(v) = m_std;
  }
  
  return m_adjusted;
}


