#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat adjust_local_variance(const arma::mat& m1,
                                const arma::mat& m2,
                                const arma::sp_mat& NNmatrix) {
  int N = m1.n_rows;  // number of individuals
  int V = m1.n_cols;  // number of vertices
  
  arma::mat rho_mat(N, V, fill::zeros);
  
  for (int v = 0; v < V; ++v) {
    // Extract neighbors of vertex v
    arma::uvec neighbors = arma::find(NNmatrix.row(v));
    int K = neighbors.n_elem;
    
    if (K == 1) {
     arma::mat prod = m1.cols(neighbors) % m2.cols(neighbors);  // N x 1
     rho_mat.col(v) = prod.col(0);
     continue;
    }
    
    // Subset data: N x K
    arma::mat X = m1.cols(neighbors);
    arma::mat Y = m2.cols(neighbors);
    
    // Apply weights row-wise (broadcasted across each row)
   
    //arma::vec x_ss = arma::sum(X % X, 1);
    //arma::vec y_ss = arma::sum(Y % Y, 1);  // sum of squares
    
    //arma::vec denom = arma::sqrt(x_ss % y_ss);  // element-wise product, then sqrt
    // Number of columns (observations per row)
    double n = static_cast<double>(X.n_cols);
    // Center X and Y by subtracting row means
    arma::mat X_centered = X.each_col() - arma::mean(X, 1);
    arma::mat Y_centered = Y.each_col() - arma::mean(Y, 1);
    
    // Compute row-wise sum of squares of centered data
    arma::vec x_ss = arma::sum(X_centered % X_centered, 1);
    arma::vec y_ss = arma::sum(Y_centered % Y_centered, 1);
    
    // Calculate sample variance explicitly by dividing by (n - 1)
    arma::vec x_var = x_ss / (n - 1);
    arma::vec y_var = y_ss / (n - 1);
    
    // Compute normalization denominator: sqrt of product of sample variances
    arma::vec denom = arma::sqrt(x_var % y_var);
    
    // // Compute normalization denominator: sqrt of product of variances
    // arma::vec denom = arma::sqrt(x_ss % y_ss);
    
    for (int i = 0; i < N; ++i) {
      rho_mat(i, v) = m1(i, v) * m2(i, v) / denom(i);
    }
  }
  
  return rho_mat;
}


// [[Rcpp::export]]
arma::mat adjust_data_local_sd(const arma::mat& m1,
                               const arma::sp_mat& NNmatrix) {
  int N = m1.n_rows;
  int V = m1.n_cols;
  arma::mat m1_norm(N, V, arma::fill::zeros);
  
  for (int v = 0; v < V; ++v) {
    arma::uvec neighbors = arma::find(NNmatrix.row(v));
    int K = neighbors.n_elem;
    
    if (K == 0) {
      m1_norm.col(v).fill(arma::datum::nan);
      continue;
    }
    
    // extract submatrix: N × K (values across all subjects for v's neighbors)
    arma::mat neighbor_vals = m1.cols(neighbors);  // faster than per-row extraction
    
    // standard deviation across each row (subject)
    arma::vec row_sd = arma::stddev(neighbor_vals, 0, 1);  // 0 = normalize by N-1, 1 = row-wise
    
    // avoid division by 0 (keep original values)
    for (int i = 0; i < N; ++i) {
      double sd = row_sd(i);
      m1_norm(i, v) = (sd > 0) ? m1(i, v) / sd : m1(i, v);
    }
    
    // Print progress every 100 vertices
   // if (v % 100 == 0) {
   //    
   //    
   //    Rcpp::Rcout << "Processed vertex " << v << "/" << V;
   //  
   //  }
  }
  
  return m1_norm;
}



// [[Rcpp::export]]
arma::mat adjust_data_local_sd_both(const arma::mat& m1,
                                         const arma::sp_mat& NNmatrix) {
  int N = m1.n_rows;  // number of subjects
  int V = m1.n_cols;  // number of vertices
  
  arma::mat m1_norm(N, V, arma::fill::zeros);
  
  for (int v = 0; v < V; ++v) {
    arma::uvec neighbors = arma::find(NNmatrix.row(v));
    int K = neighbors.n_elem;
    
    if (K == 0) {
      m1_norm.col(v).fill(arma::datum::nan);
      continue;
    }
    
    // Subset: N × K matrix of values for neighbors of v
    arma::mat m1_neigh = m1.cols(neighbors);
    
    // Flatten to vector of length N * K
    arma::vec m1_vec = arma::vectorise(m1_neigh);
    
    // Compute standard deviation
    double sd = arma::stddev(m1_vec);
    
    if (sd > 0) {
      m1_norm.col(v) = m1.col(v) / sd;
    } else {
      m1_norm.col(v).zeros();  // or NA
    }
  }
  
  return m1_norm;
}

