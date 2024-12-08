//' @useDynLib SA24204149
//' @import Rcpp
 
#include <Rcpp.h>
using namespace Rcpp;
//' @title Gibbs
//' @description Gibbs sampler
//' @param n_samples samples
//' @param n a int value
//' @param a a double value1
//' @param b a double value2
//' @return List 
//' @export
// [[Rcpp::export]]
List Gibbs(int n_samples, int n, double a, double b) {
  // Initialize the chain with arbitrary starting values
  int x = 0;
  double y = 0.5;
  
  // Vectors to store the samples
  IntegerVector x_samples(n_samples);
  NumericVector y_samples(n_samples);
  
  for (int i = 0; i < n_samples; ++i) {
    // Sample from the conditional distribution of x given y
    x = R::rbinom(n, y);
    
    // Sample from the conditional distribution of y given x
    y = R::rbeta(x + a, n - x + b);
    
    // Store the samples
    x_samples[i] = x;
    y_samples[i] = y;
  }
  
  return List::create(Named("x") = x_samples, Named("y") = y_samples);
}