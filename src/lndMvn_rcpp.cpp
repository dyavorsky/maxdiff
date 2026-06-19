// lndMvn_rcpp.cpp — log multivariate normal density.
//
// Original code by Wayne Taylor (2014) from the bayesm package (Peter Rossi).
// Internal C++ function (not exported to R).

#include "maxdiff.h"

double lndMvn(vec const& x, vec const& mu, mat const& rooti) {

  // rooti is upper-triangular Cholesky root of Sigma^-1 (UL decomposition).
  vec z = vectorise(trans(rooti) * (x - mu));
  return (-(x.size() / 2.0) * log(2 * M_PI) - .5 * (trans(z) * z) +
          sum(log(diagvec(rooti))))[0];
}
