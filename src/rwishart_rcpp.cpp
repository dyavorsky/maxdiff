// rwishart_rcpp.cpp — Wishart and inverse-Wishart draws.
//
// Original code by Wayne Taylor (2015) from the bayesm package (Peter Rossi).
// Internal C++ function (not exported to R).

#include "maxdiff.h"

List rwishart(double nu, mat const& V) {

  int m = V.n_rows;
  mat T = zeros(m, m);

  for (int i = 0; i < m; i++)
    T(i, i) = sqrt(rchisq(1, nu - i)[0]);

  for (int j = 0; j < m; j++)
    for (int i = j + 1; i < m; i++)
      T(i, j) = rnorm(1)[0];

  mat C  = trans(T) * chol(V);
  mat CI = solve(trimatu(C), eye(m, m));

  return List::create(Named("W")  = trans(C) * C,
                      Named("IW") = CI * trans(CI),
                      Named("C")  = C,
                      Named("CI") = CI);
}
