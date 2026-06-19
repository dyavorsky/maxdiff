// llmnl_rcpp.cpp — MNL log-likelihood for plain (unconstrained) beta.
//
// Original code by Wayne Taylor (2014) from the bayesm package (Peter Rossi).
// Ragged choice-set extension (p_vec parameter) by Dan Yavorsky.
//
// Reference: Rossi, Allenby, and McCulloch (2024), Bayesian Statistics and
//   Marketing (2nd ed.), Wiley.

#include "maxdiff.h"

double llmnl(vec const& beta, vec const& y, mat const& X,
             IntegerVector p_vec) {

  // Evaluates log-likelihood for the multinomial logit model.
  // p_vec (optional): length-n integer vector giving the choice-set size at
  //   each occasion. Empty => legacy uniform-p path (j = X.n_rows / n).

  int  n     = y.size();
  mat  Xbeta = X * beta;
  vec  xby   = zeros<vec>(n);
  vec  denom = zeros<vec>(n);

  if (p_vec.size() == 0) {
    int j = X.n_rows / n;
    for (int i = 0; i < n; i++) {
      for (int p = 0; p < j; p++) denom[i] += exp(Xbeta[i * j + p]);
      xby[i] = Xbeta[i * j + y[i] - 1];
    }
  } else {
    int offset = 0;
    for (int i = 0; i < n; i++) {
      int j = p_vec[i];
      for (int p = 0; p < j; p++) denom[i] += exp(Xbeta[offset + p]);
      xby[i]  = Xbeta[offset + y[i] - 1];
      offset += j;
    }
  }
  return sum(xby - log(denom));
}
