// rmixGibbs_rcpp.cpp — Gibbs sampler for mixture of normals.
//
// Original implementation by Wayne Taylor (2014, 2015) and Peter Rossi from
// the bayesm package.
//
// Reference: Rossi, Allenby, and McCulloch (2024), Bayesian Statistics and
//   Marketing (2nd ed.), Wiley.

#include "maxdiff.h"

List drawCompsFromLabels(mat const& y, mat const& Bbar,
                         mat const& A, double nu,
                         mat const& V, int ncomp,
                         vec const& z) {

  vec b, r, mu;
  mat yk, Xk, Ck, sigma, rooti, S, IW, CI;
  List temp, rw, comps(ncomp);

  int n = z.n_rows;
  vec nobincomp = zeros<vec>(ncomp);

  for (int i = 0; i < n; i++) nobincomp[z[i] - 1]++;

  for (int k = 0; k < ncomp; k++) {
    if (nobincomp[k] > 0) {
      yk   = y.rows(find(z == (k + 1)));
      Xk   = ones(nobincomp[k], 1);
      temp = rmultireg(yk, Xk, Bbar, A, nu, V);
      sigma = as<mat>(temp["Sigma"]);
      rooti = solve(trimatu(chol(sigma)), eye(sigma.n_rows, sigma.n_cols));
      mu    = as<vec>(temp["B"]);
      comps(k) = List::create(
        Named("mu")    = NumericVector(mu.begin(), mu.end()),
        Named("rooti") = rooti);
    } else {
      S  = solve(trimatu(chol(V)), eye(V.n_rows, V.n_cols));
      S  = S * trans(S);
      rw = rwishart(nu, S);
      IW = as<mat>(rw["IW"]);
      CI = as<mat>(rw["CI"]);
      rooti = solve(trimatu(chol(IW)), eye(IW.n_rows, IW.n_cols));
      b = vectorise(Bbar);
      r = rnorm(b.n_rows, 0, 1);
      mu = b + (CI * r) / sqrt(A(0, 0));
      comps(k) = List::create(
        Named("mu")    = NumericVector(mu.begin(), mu.end()),
        Named("rooti") = rooti);
    }
  }
  return comps;
}

vec drawLabelsFromComps(mat const& y, vec const& p, List comps) {

  double logprod;
  vec mu, u;
  mat rooti;
  List compsk;

  int n     = y.n_rows;
  vec res   = zeros<vec>(n);
  int ncomp = comps.size();
  mat prob(n, ncomp);

  for (int k = 0; k < ncomp; k++) {
    compsk = comps[k];
    mu     = as<vec>(compsk["mu"]);
    rooti  = as<mat>(compsk["rooti"]);
    logprod = log(prod(diagvec(rooti)));
    mat z(y);
    z.each_row() -= trans(mu);
    z = trans(rooti) * trans(z);
    z = -(y.n_cols / 2.0) * log(2 * M_PI) + logprod - .5 * sum(z % z, 0);
    prob.col(k) = trans(z);
  }

  prob = exp(prob);
  prob.each_row() %= trans(p);
  prob = cumsum(prob, 1);
  u    = as<vec>(runif(n)) % prob.col(ncomp - 1);

  for (int i = 0; i < n; i++) {
    while (u[i] > prob(i, res[i]++));
  }
  return res;
}

vec drawPFromLabels(vec const& a, vec const& z) {

  vec a2 = a;
  int n  = z.n_rows;
  for (int i = 0; i < n; i++) a2[z[i] - 1]++;
  return rdirichlet(a2);
}

// rmixGibbs is C++-only (not exported to R) — called from rhierMnlRwMixture_rcpp_loop.
List rmixGibbs(mat const& y, mat const& Bbar,
               mat const& A, double nu,
               mat const& V, vec const& a,
               vec const& p, vec const& z) {

  List comps = drawCompsFromLabels(y, Bbar, A, nu, V, a.size(), z);
  vec  z2    = drawLabelsFromComps(y, p, comps);
  vec  p2    = drawPFromLabels(a, z2);

  return List::create(Named("p")     = p2,
                      Named("z")     = z2,
                      Named("comps") = comps);
}
