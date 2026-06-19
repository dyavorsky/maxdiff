// rdirichlet_rcpp.cpp — Dirichlet draw.
//
// Original code by Wayne Taylor (2015) from the bayesm package (Peter Rossi).
// Internal C++ function (not exported to R).

#include "maxdiff.h"

vec rdirichlet(vec const& alpha) {

  int dim = alpha.size();
  vec y   = zeros<vec>(dim);

  for (int i = 0; i < dim; i++)
    y[i] = rgamma(1, alpha[i])[0];

  return y / sum(y);
}
