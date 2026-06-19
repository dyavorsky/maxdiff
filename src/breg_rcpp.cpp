// breg_rcpp.cpp — Bayesian linear regression draw (sigma^2 = 1).
//
// Original code by Keunwoo Kim (2014) from the bayesm package (Peter Rossi).
// Internal C++ function (not exported to R).

#include "maxdiff.h"

vec breg(vec const& y, mat const& X, vec const& betabar, mat const& A) {

  int k  = betabar.size();
  mat RA = chol(A);
  mat W  = join_cols(X, RA);
  vec z  = join_cols(y, RA * betabar);
  mat IR = solve(trimatu(chol(trans(W) * W)), eye(k, k));

  return (IR * trans(IR)) * (trans(W) * z) + IR * vec(rnorm(k));
}
