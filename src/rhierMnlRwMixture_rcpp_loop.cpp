// rhierMnlRwMixture_rcpp_loop.cpp — MNL log-likelihood (sign-constrained)
// and the main Gibbs loop for hierarchical MNL with mixture-of-normals
// heterogeneity.
//
// Original implementation by Wayne Taylor (2014, 2016) and Peter Rossi from
// the bayesm package. Ragged choice-set extensions (p_vec, p_it_list) by
// Dan Yavorsky.
//
// Reference: Rossi, Allenby, and McCulloch (2024), Bayesian Statistics and
//   Marketing (2nd ed.), Wiley.

#include "maxdiff.h"

// [[Rcpp::export]]
double llmnl_con(vec const& betastar, vec const& y, mat const& X,
                 vec const& SignRes = NumericVector::create(0),
                 IntegerVector p_vec = IntegerVector::create()) {

  // Evaluates log-likelihood for MNL with sign constraints.
  // p_vec (optional): length-n integer vector giving the choice-set size at
  //   each occasion. Empty => legacy uniform-p path (j = X.n_rows / n).

  vec beta = betastar;
  if (any(SignRes)) {
    uvec signInd = find(SignRes != 0);
    beta.elem(signInd) = SignRes.elem(signInd) % exp(beta.elem(signInd));
  }

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

// Internal: RW Metropolis step for sign-constrained MNL.
// p_vec forwarded to llmnl_con.
mnlMetropOnceOut mnlMetropOnce_con(vec const& y, mat const& X,
                                   vec const& oldbeta, double oldll,
                                   double s, mat const& incroot,
                                   vec const& betabar, mat const& rootpi,
                                   vec const& SignRes = NumericVector::create(2),
                                   IntegerVector p_vec = IntegerVector::create()) {

  mnlMetropOnceOut out_struct;
  double unif;
  vec betadraw, alphaminv(2);
  int stay = 0;

  vec betac     = oldbeta + s * trans(incroot) * as<vec>(rnorm(X.n_cols));
  double cll    = llmnl_con(betac, y, X, SignRes, p_vec);
  double clpost = cll + lndMvn(betac, betabar, rootpi);
  double ldiff  = clpost - oldll - lndMvn(oldbeta, betabar, rootpi);
  alphaminv = { 1, exp(ldiff) };
  double alpha  = min(alphaminv);

  if (alpha < 1) { unif = as_scalar(vec(runif(1))); } else { unif = 0; }
  if (unif <= alpha) { betadraw = betac; oldll = cll; }
  else               { betadraw = oldbeta; stay = 1; }

  out_struct.betadraw = betadraw;
  out_struct.stay     = stay;
  out_struct.oldll    = oldll;
  return out_struct;
}

// [[Rcpp::export]]
List rhierMnlRwMixture_rcpp_loop(List const& lgtdata, mat const& Z,
                                  vec const& deltabar, mat const& Ad,
                                  mat const& mubar,    mat const& Amu,
                                  double nu,           mat const& V,
                                  double s,
                                  int R, int keep, int nprint, bool drawdelta,
                                  mat olddelta, vec const& a, vec oldprob,
                                  mat oldbetas, vec ind, vec const& SignRes,
                                  List const& p_it_list) {

  // Main Gibbs loop for hierarchical MNL with mixture-of-normals heterogeneity.
  // p_it_list: list of length nlgt; element lgt is an IntegerVector of length
  //   T_lgt giving the choice-set size at each occasion.

  int nlgt = lgtdata.size();
  int nvar = V.n_cols;
  int nz   = Z.n_cols;

  mat rootpi, betabar_loop, ucholinv, incroot;
  int mkeep;
  mnlMetropOnceOut metropout_struct;
  List lgtdatai, nmix;

  std::vector<moments>       lgtdata_vector;
  std::vector<IntegerVector> p_vec_vector;
  moments lgtdatai_struct;
  for (int lgt = 0; lgt < nlgt; lgt++) {
    lgtdatai = lgtdata[lgt];
    lgtdatai_struct.y    = as<vec>(lgtdatai["y"]);
    lgtdatai_struct.X    = as<mat>(lgtdatai["X"]);
    lgtdatai_struct.hess = as<mat>(lgtdatai["hess"]);
    lgtdata_vector.push_back(lgtdatai_struct);
    p_vec_vector.push_back(as<IntegerVector>(p_it_list[lgt]));
  }

  vec  oldll    = zeros<vec>(nlgt);
  cube betadraw(nlgt, nvar, R / keep);
  mat  probdraw(R / keep, oldprob.size());
  vec  loglike(R / keep);
  mat  Deltadraw(1, 1);
  if (drawdelta) Deltadraw.zeros(R / keep, nz * nvar);
  List compdraw(R / keep);

  if (nprint > 0) startMcmcTimer();

  for (int rep = 0; rep < R; rep++) {
    List mgout;
    if (drawdelta) {
      olddelta.reshape(nvar, nz);
      mgout = rmixGibbs(oldbetas - Z * trans(olddelta), mubar, Amu, nu, V, a, oldprob, ind);
    } else {
      mgout = rmixGibbs(oldbetas, mubar, Amu, nu, V, a, oldprob, ind);
    }
    List oldcomp = mgout["comps"];
    oldprob = as<vec>(mgout["p"]);
    ind     = as<vec>(mgout["z"]);

    if (drawdelta) olddelta = drawDelta(Z, oldbetas, ind, oldcomp, deltabar, Ad);

    for (int lgt = 0; lgt < nlgt; lgt++) {
      List oldcomplgt = oldcomp[ind[lgt] - 1];
      rootpi = as<mat>(oldcomplgt[1]);
      if (drawdelta) {
        olddelta.reshape(nvar, nz);
        betabar_loop = as<vec>(oldcomplgt[0]) + olddelta * vectorise(Z(lgt, span::all));
      } else {
        betabar_loop = as<vec>(oldcomplgt[0]);
      }
      if (rep == 0)
        oldll[lgt] = llmnl_con(vectorise(oldbetas(lgt, span::all)),
                               lgtdata_vector[lgt].y, lgtdata_vector[lgt].X,
                               SignRes, p_vec_vector[lgt]);
      ucholinv = solve(trimatu(chol(lgtdata_vector[lgt].hess + rootpi * trans(rootpi))),
                       eye(nvar, nvar));
      incroot  = chol(ucholinv * trans(ucholinv));
      metropout_struct = mnlMetropOnce_con(
        lgtdata_vector[lgt].y, lgtdata_vector[lgt].X,
        vectorise(oldbetas(lgt, span::all)), oldll[lgt],
        s, incroot, betabar_loop, rootpi, SignRes, p_vec_vector[lgt]);
      oldbetas(lgt, span::all) = trans(metropout_struct.betadraw);
      oldll[lgt]               = metropout_struct.oldll;
    }

    if (nprint > 0) if ((rep + 1) % nprint == 0) infoMcmcTimer(rep, R);
    if ((rep + 1) % keep == 0) {
      mkeep = (rep + 1) / keep;
      betadraw.slice(mkeep - 1) = oldbetas;
      probdraw(mkeep - 1, span::all) = trans(oldprob);
      loglike[mkeep - 1] = sum(oldll);
      if (drawdelta) Deltadraw(mkeep - 1, span::all) = trans(vectorise(olddelta));
      compdraw[mkeep - 1] = oldcomp;
    }
  }

  if (nprint > 0) endMcmcTimer();

  nmix = List::create(Named("probdraw") = probdraw,
                      Named("zdraw")    = R_NilValue,
                      Named("compdraw") = compdraw);

  bool conStatus = any(SignRes);
  if (conStatus) {
    int SignResSize = SignRes.size();
    for (int i = 0; i < SignResSize; i++) {
      if (SignRes[i] != 0) {
        for (int s = 0; s < R / keep; s++)
          betadraw(span(), span(i), span(s)) =
            SignRes[i] * exp(betadraw(span(), span(i), span(s)));
      }
    }
  }

  if (drawdelta) {
    return List::create(Named("Deltadraw") = Deltadraw,
                        Named("betadraw")  = betadraw,
                        Named("nmix")      = nmix,
                        Named("loglike")   = loglike,
                        Named("SignRes")   = SignRes);
  } else {
    return List::create(Named("betadraw") = betadraw,
                        Named("nmix")     = nmix,
                        Named("loglike")  = loglike,
                        Named("SignRes")  = SignRes);
  }
}
