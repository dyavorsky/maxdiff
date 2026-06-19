// rmultireg_rcpp.cpp — Bayesian multivariate regression draw.
//
// Original code by Keunwoo Kim (2014) from the bayesm package (Peter Rossi).
// Internal C++ function (not exported to R).

#include "maxdiff.h"

List rmultireg(mat const& Y, mat const& X, mat const& Bbar, mat const& A,
               double nu, mat const& V) {

  int n = Y.n_rows;
  int m = Y.n_cols;
  int k = X.n_cols;

  mat RA     = chol(A);
  mat W      = join_cols(X, RA);
  mat Z      = join_cols(Y, RA * Bbar);
  mat IR     = solve(trimatu(chol(trans(W) * W)), eye(k, k));
  mat Btilde = (IR * trans(IR)) * (trans(W) * Z);
  mat E      = Z - W * Btilde;
  mat S      = trans(E) * E;

  mat ucholinv = solve(trimatu(chol(V + S)), eye(m, m));
  mat VSinv    = ucholinv * trans(ucholinv);
  List rwout   = rwishart(nu + n, VSinv);

  mat CI   = rwout["CI"];
  mat draw = mat(rnorm(k * m));
  draw.reshape(k, m);
  mat B = Btilde + IR * draw * trans(CI);

  return List::create(Named("B")     = B,
                      Named("Sigma") = rwout["IW"]);
}
