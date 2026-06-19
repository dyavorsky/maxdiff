// utilityFunctions.cpp — shared utility functions for the maxdiff C++ layer.
//
// Original implementations by Wayne Taylor, Keunwoo Kim, and Peter Rossi from
// the bayesm package. All statistical algorithms are their work.
// p_vec extension to mnlMetropOnce by Dan Yavorsky.
//
// Reference: Rossi, Allenby, and McCulloch (2024), Bayesian Statistics and
//   Marketing (2nd ed.), Wiley.

#include "maxdiff.h"

// Used in rmvpGibbs and rmnpGibbs -----------------------------------------
vec condmom(vec const& x, vec const& mu, mat const& sigmai, int p, int j) {

  vec out(2);
  int jm1 = j - 1;
  int ind  = p * jm1;
  double csigsq = 1 / sigmai(ind + jm1);
  double m = 0.0;
  for (int i = 0; i < p; i++)
    if (i != jm1) m += -csigsq * sigmai(ind + i) * (x[i] - mu[i]);
  out[0] = mu[jm1] + m;
  out[1] = sqrt(csigsq);
  return out;
}

double rtrun1(double mu, double sigma, double trunpt, int above) {

  double FA, FB, rnd, result, arg;
  if (above) {
    FA = 0.0;
    FB = R::pnorm((trunpt - mu) / sigma, 0.0, 1.0, 1, 0);
  } else {
    FB = 1.0;
    FA = R::pnorm((trunpt - mu) / sigma, 0.0, 1.0, 1, 0);
  }
  rnd = runif(1)[0];
  arg = rnd * (FB - FA) + FA;
  if (arg > .999999999)  arg = .999999999;
  if (arg < .0000000001) arg = .0000000001;
  result = mu + sigma * R::qnorm(arg, 0.0, 1.0, 1, 0);
  return result;
}

// Used in rhierLinearModel, rhierLinearMixture, rhierMnlRwMixture -----------
mat drawDelta(mat const& x, mat const& y, vec const& z, List const& comps,
              vec const& deltabar, mat const& Ad) {

  int p     = y.n_cols;
  int k     = x.n_cols;
  int ncomp = comps.length();
  mat xtx   = zeros<mat>(k * p, k * p);
  mat xty   = zeros<mat>(p, k);

  uvec colAlly(p), colAllx(k);
  for (int i = 0; i < p; i++) colAlly(i) = i;
  for (int i = 0; i < k; i++) colAllx(i) = i;

  for (int compi = 0; compi < ncomp; compi++) {
    uvec ind = find(z == (compi + 1));
    if (ind.size() > 0) {
      mat yi   = y.submat(ind, colAlly);
      mat xi   = x.submat(ind, colAllx);
      List compsi  = comps[compi];
      rowvec mui   = as<rowvec>(compsi[0]);
      mat rootii   = trimatu(as<mat>(compsi[1]));
      yi.each_row() -= mui;
      mat sigi = rootii * trans(rootii);
      xtx      = xtx + kron(trans(xi) * xi, sigi);
      xty      = xty + (sigi * (trans(yi) * xi));
    }
  }
  xty.reshape(xty.n_rows * xty.n_cols, 1);
  mat ucholinv = solve(trimatu(chol(xtx + Ad)), eye(k * p, k * p));
  mat Vinv     = ucholinv * trans(ucholinv);
  return Vinv * (xty + Ad * deltabar) +
         trans(chol(Vinv)) * as<vec>(rnorm(deltabar.size()));
}

unireg runiregG(vec const& y, mat const& X, mat const& XpX, vec const& Xpy,
                double sigmasq, mat const& A, vec const& Abetabar,
                double nu, double ssq) {

  unireg out_struct;
  int n = y.size();
  int k = XpX.n_cols;

  mat IR     = solve(trimatu(chol(XpX / sigmasq + A)), eye(k, k));
  vec btilde = (IR * trans(IR)) * (Xpy / sigmasq + Abetabar);
  vec beta   = btilde + IR * vec(rnorm(k));
  double s   = sum(square(y - X * beta));
  sigmasq    = (s + nu * ssq) / rchisq(1, nu + n)[0];

  out_struct.beta    = beta;
  out_struct.sigmasq = sigmasq;
  return out_struct;
}

// Used in rnegbinRW and rhierNegbinRw --------------------------------------
double llnegbin(vec const& y, vec const& lambda, double alpha, bool constant) {

  int   i;
  int   nobs = y.size();
  vec   prob = alpha / (alpha + lambda);
  vec   logp(nobs);
  if (constant) {
    for (i = 0; i < nobs; i++)
      logp[i] = R::dnbinom(y[i], alpha, prob[i], 1);
  } else {
    logp = sum(alpha * log(prob) + y % log(1 - prob));
  }
  return sum(logp);
}

double lpostbeta(double alpha, vec const& beta, mat const& X, vec const& y,
                 vec const& betabar, mat const& rootA) {

  vec lambda = exp(X * beta);
  double ll  = llnegbin(y, lambda, alpha, FALSE);
  vec z      = rootA * (beta - betabar);
  return ll - 0.5 * sum(z % z);
}

double lpostalpha(double alpha, vec const& beta, mat const& X, vec const& y,
                  double a, double b) {

  vec lambda = exp(X * beta);
  double ll  = llnegbin(y, lambda, alpha, TRUE);
  return ll + (a - 1) * log(alpha) - b * alpha;
}

// Used in rbprobitGibbs and rordprobitGibbs --------------------------------
vec breg1(mat const& root, mat const& X, vec const& y, vec const& Abetabar) {

  mat cov = trans(root) * root;
  return cov * (trans(X) * y + Abetabar) + trans(root) * vec(rnorm(root.n_cols));
}

vec rtrunVec(vec const& mu, vec const& sigma, vec const& a, vec const& b) {

  int n = mu.size();
  vec FA(n), FB(n), out(n);
  for (int i = 0; i < n; i++) {
    FA[i]  = R::pnorm((a[i] - mu[i]) / sigma[i], 0, 1, 1, 0);
    FB[i]  = R::pnorm((b[i] - mu[i]) / sigma[i], 0, 1, 1, 0);
    out[i] = mu[i] + sigma[i] *
             R::qnorm(R::runif(0, 1) * (FB[i] - FA[i]) + FA[i], 0, 1, 1, 0);
  }
  return out;
}

// Used in rhierMnlDP and rhierMnlRwMixture ---------------------------------
mnlMetropOnceOut mnlMetropOnce(vec const& y, mat const& X, vec const& oldbeta,
                               double oldll, double s, mat const& incroot,
                               vec const& betabar, mat const& rootpi,
                               IntegerVector p_vec) {

  mnlMetropOnceOut metropout_struct;
  double unif;
  vec betadraw, alphaminv(2);
  int stay = 0;

  vec betac     = oldbeta + s * trans(incroot) * as<vec>(rnorm(X.n_cols));
  double cll    = llmnl(betac, y, X, p_vec);
  double clpost = cll + lndMvn(betac, betabar, rootpi);
  double ldiff  = clpost - oldll - lndMvn(oldbeta, betabar, rootpi);
  alphaminv = { 1, exp(ldiff) };
  double alpha  = min(alphaminv);

  if (alpha < 1) { unif = runif(1)[0]; } else { unif = 0; }
  if (unif <= alpha) { betadraw = betac; oldll = cll; }
  else               { betadraw = oldbeta; stay = 1; }

  metropout_struct.betadraw = betadraw;
  metropout_struct.stay     = stay;
  metropout_struct.oldll    = oldll;
  return metropout_struct;
}

// Used in rDPGibbs, rhierMnlDP, rivDP -------------------------------------
int rmultinomF(vec const& p) {

  vec csp = cumsum(p);
  double rnd = runif(1)[0];
  int res    = 0;
  int psize  = p.size();
  for (int i = 0; i < psize; i++)
    if (rnd > csp[i]) res++;
  return res + 1;
}

mat yden(std::vector<murooti> const& thetaStar_vector, mat const& y) {

  int nunique = thetaStar_vector.size();
  int n = y.n_rows;
  int k = y.n_cols;
  mat ydenmat = zeros<mat>(nunique, n);
  vec mu;
  mat rooti, transy, quads;

  for (int i = 0; i < nunique; i++) {
    mu     = thetaStar_vector[i].mu;
    rooti  = thetaStar_vector[i].rooti;
    transy = trans(y);
    transy.each_col() -= mu;
    quads = sum(square(trans(rooti) * transy), 0);
    ydenmat(i, span::all) =
      exp(-(k / 2.0) * log(2 * M_PI) + sum(log(rooti.diag())) - .5 * quads);
  }
  return ydenmat;
}

ivec numcomp(ivec const& indic, int k) {

  ivec ncomp(k);
  for (int comp = 0; comp < k; comp++)
    ncomp[comp] = sum(indic == (comp + 1));
  return ncomp;
}

murooti thetaD(mat const& y, lambda const& lambda_struct) {

  mat  X = ones<mat>(y.n_rows, 1);
  mat  A(1, 1); A.fill(lambda_struct.Amu);
  List rout = rmultireg(y, X, trans(lambda_struct.mubar), A,
                        lambda_struct.nu, lambda_struct.V);

  murooti out_struct;
  out_struct.mu    = as<vec>(rout["B"]);
  out_struct.rooti = solve(chol(trimatu(as<mat>(rout["Sigma"]))),
                           eye(y.n_cols, y.n_cols));
  return out_struct;
}

thetaStarIndex thetaStarDraw(ivec indic,
                              std::vector<murooti> thetaStar_vector,
                              mat const& y, mat ydenmat,
                              vec const& q0v, double alpha,
                              lambda const& lambda_struct, int maxuniq) {

  int n = indic.size();
  ivec ncomp, indicC;
  int k, inc, cntNonzero;
  std::vector<murooti> listofone_vector(1);
  std::vector<murooti> thetaStarC_vector;

  for (int i = 0; i < n; i++) {
    k = thetaStar_vector.size();
    vec probs(k + 1);
    probs[k]  = q0v[i] * (alpha / (alpha + (n - 1)));

    ivec indicmi = zeros<ivec>(n - 1);
    inc = 0;
    for (int j = 0; j < (n - 1); j++) {
      if (j == i) inc++;
      indicmi[j] = indic[inc];
      inc++;
    }
    ncomp = numcomp(indicmi, k);
    for (int comp = 0; comp < k; comp++)
      probs[comp] = ydenmat(comp, i) * ncomp[comp] / (alpha + (n - 1));
    probs = probs / sum(probs);
    indic[i] = rmultinomF(probs);

    if (indic[i] == (k + 1)) {
      if ((k + 1) > maxuniq)
        stop("max number of comps exceeded");
      else {
        listofone_vector[0] = thetaD(y(i, span::all), lambda_struct);
        thetaStar_vector.push_back(listofone_vector[0]);
        ydenmat(k, span::all) = yden(listofone_vector, y);
      }
    }
  }

  k = thetaStar_vector.size();
  indicC = zeros<ivec>(n);
  ncomp  = numcomp(indic, k);
  cntNonzero = 0;
  for (int comp = 0; comp < k; comp++) {
    if (ncomp[comp] != 0) {
      thetaStarC_vector.push_back(thetaStar_vector[comp]);
      cntNonzero++;
      for (int i = 0; i < n; i++)
        if (indic[i] == (comp + 1)) indicC[i] = cntNonzero;
    }
  }

  thetaStarIndex out_struct;
  out_struct.indic             = indicC;
  out_struct.thetaStar_vector  = thetaStarC_vector;
  return out_struct;
}

vec q0(mat const& y, lambda const& lambda_struct) {

  int k = y.n_cols;
  mat R = chol(lambda_struct.V);
  double logdetR = sum(log(R.diag()));
  double lnk1k2, constant;
  mat transy, m, vivi, lnq0v;

  if (k > 1) {
    vec km1(k - 1);
    for (int i = 0; i < (k - 1); i++) km1[i] = i + 1;
    lnk1k2 = (k / 2.0) * log(2.0) + log((lambda_struct.nu - k) / 2) +
              lgamma((lambda_struct.nu - k) / 2) - lgamma(lambda_struct.nu / 2) +
              sum(log(lambda_struct.nu / 2 - km1 / 2));
  } else {
    lnk1k2 = (k / 2.0) * log(2.0) + log((lambda_struct.nu - k) / 2) +
              lgamma((lambda_struct.nu - k) / 2) - lgamma(lambda_struct.nu / 2);
  }
  constant = -(k / 2.0) * log(2 * M_PI) +
             (k / 2.0) * log(lambda_struct.Amu / (1 + lambda_struct.Amu)) +
             lnk1k2 + lambda_struct.nu * logdetR;

  transy = trans(y);
  transy.each_col() -= lambda_struct.mubar;
  m    = sqrt(lambda_struct.Amu / (1 + lambda_struct.Amu)) *
         trans(solve(trimatu(R), eye(y.n_cols, y.n_cols))) * transy;
  vivi = sum(square(m), 0);
  lnq0v = constant - ((lambda_struct.nu + 1) / 2) * (2 * logdetR + log(1 + vivi));
  return trans(exp(lnq0v));
}

vec seq_rcpp(double from, double to, int len) {

  vec res(len);
  res[len - 1] = to; res[0] = from;
  double increment = (res[len - 1] - res[0]) / (len - 1);
  for (int i = 1; i < (len - 1); i++) res[i] = res[i - 1] + increment;
  return res;
}

double alphaD(priorAlpha const& priorAlpha_struct, int Istar, int gridsize) {

  vec alpha = seq_rcpp(priorAlpha_struct.alphamin,
                       priorAlpha_struct.alphamax - .000001, gridsize);
  vec lnprob(gridsize);
  for (int i = 0; i < gridsize; i++)
    lnprob[i] = Istar * log(alpha[i]) + lgamma(alpha[i]) -
                lgamma(priorAlpha_struct.n + alpha[i]) +
                priorAlpha_struct.power *
                log(1 - (alpha[i] - priorAlpha_struct.alphamin) /
                        (priorAlpha_struct.alphamax - priorAlpha_struct.alphamin));
  lnprob = lnprob - median(lnprob);
  vec probs = exp(lnprob);
  probs = probs / sum(probs);
  return alpha(rmultinomF(probs) - 1);
}

murooti GD(lambda const& lambda_struct) {

  int k = lambda_struct.mubar.size();
  List Rout = rwishart(lambda_struct.nu,
                       solve(trimatu(lambda_struct.V), eye(k, k)));
  mat Sigma = as<mat>(Rout["IW"]);
  mat root  = chol(Sigma);
  mat draws = rnorm(k);
  mat mu    = lambda_struct.mubar +
              (1 / sqrt(lambda_struct.Amu)) * trans(root) * draws;

  murooti out_struct;
  out_struct.mu    = mu;
  out_struct.rooti = solve(trimatu(root), eye(k, k));
  return out_struct;
}

lambda lambdaD(lambda const& lambda_struct,
               std::vector<murooti> const& thetaStar_vector,
               vec const& alim, vec const& nulim, vec const& vlim,
               int gridsize) {

  vec  lnprob, probs, rowSumslgammaarg;
  int  ind;
  murooti thetaStari_struct;
  mat  rootii, mout, rimu, arg, lgammaarg;
  vec  mui;
  double sumdiagriri, sumlogdiag, sumquads, adraw, nudraw, vdraw;

  murooti thetaStar0_struct = thetaStar_vector[0];
  int d     = thetaStar0_struct.mu.size();
  int Istar = thetaStar_vector.size();

  vec aseq  = seq_rcpp(alim[0],  alim[1],  gridsize);
  vec nuseq = d - 1 + exp(seq_rcpp(nulim[0], nulim[1], gridsize));
  vec vseq  = seq_rcpp(vlim[0],  vlim[1],  gridsize);

  mout = zeros<mat>(d, Istar * d);
  ind  = 0;
  for (int i = 0; i < Istar; i++) {
    thetaStari_struct = thetaStar_vector[i];
    rootii = thetaStari_struct.rooti;
    ind    = i * d;
    mout.submat(0, ind, d - 1, ind + d - 1) = trans(rootii);
  }
  sumdiagriri = sum(sum(square(mout), 0));

  sumlogdiag = 0.0;
  for (int i = 0; i < Istar; i++) {
    ind = i * d;
    for (int j = 0; j < d; j++) sumlogdiag += log(mout(j, ind + j));
  }

  rimu = zeros<mat>(d, Istar);
  for (int i = 0; i < Istar; i++) {
    thetaStari_struct = thetaStar_vector[i];
    mui    = thetaStari_struct.mu;
    rootii = thetaStari_struct.rooti;
    rimu(span::all, i) = trans(rootii) * mui;
  }
  sumquads = sum(sum(square(rimu), 0));

  // draw a
  lnprob = Istar * (-(d / 2.0) * log(2 * M_PI)) -
           .5 * aseq * sumquads + Istar * d * log(sqrt(aseq)) + sumlogdiag;
  lnprob = lnprob - max(lnprob) + 200;
  probs  = exp(lnprob);
  probs  = probs / sum(probs);
  adraw  = aseq[rmultinomF(probs) - 1];

  // draw nu
  lnprob = zeros<vec>(nuseq.size());
  arg    = zeros<mat>(gridsize, d);
  for (int i = 0; i < d; i++) {
    vec indvec(gridsize);
    indvec.fill(-(i + 1) + 1);
    arg(span::all, i) = indvec;
  }
  arg.each_col() += nuseq;
  arg = arg / 2.0;
  lgammaarg = zeros<mat>(gridsize, d);
  for (int i = 0; i < gridsize; i++)
    for (int j = 0; j < d; j++)
      lgammaarg(i, j) = lgamma(arg(i, j));
  rowSumslgammaarg = sum(lgammaarg, 1);
  lnprob = zeros<vec>(gridsize);
  for (int i = 0; i < gridsize; i++)
    lnprob[i] = -Istar * log(2.0) * d / 2.0 * nuseq[i] -
                Istar * rowSumslgammaarg[i] +
                Istar * d * log(sqrt(lambda_struct.V(0, 0))) * nuseq[i] +
                sumlogdiag * nuseq[i];
  lnprob = lnprob - max(lnprob) + 200;
  probs  = exp(lnprob);
  probs  = probs / sum(probs);
  nudraw = nuseq[rmultinomF(probs) - 1];

  // draw v
  lnprob = Istar * nudraw * d * log(sqrt(vseq * nudraw)) -
           .5 * sumdiagriri * vseq * nudraw;
  lnprob = lnprob - max(lnprob) + 200;
  probs  = exp(lnprob);
  probs  = probs / sum(probs);
  vdraw  = vseq[rmultinomF(probs) - 1];

  lambda out_struct;
  out_struct.mubar = zeros<vec>(d);
  out_struct.Amu   = adraw;
  out_struct.nu    = nudraw;
  out_struct.V     = nudraw * vdraw * eye(d, d);
  return out_struct;
}

// Root-finder (internal only) ---------------------------------------------
static double root_find(double c1, double c2, double tol, int iterlim) {

  int iter = 0;
  double uold = .1, unew = .00001;
  while (iter <= iterlim && fabs(uold - unew) > tol) {
    uold = unew;
    unew = uold + (uold * (c1 - c2 * uold - log(uold))) / (1. + c2 * uold);
    if (unew < 1.0e-50) unew = 1.0e-50;
    iter++;
  }
  return unew;
}

// callroot: no longer exported to R in maxdiff ----------------------------
vec callroot(vec const& c1, vec const& c2, double tol, int iterlim) {

  int n = c1.size();
  vec u = zeros<vec>(n);
  for (int i = 0; i < n; i++)
    u[i] = root_find(c1[i], c2[i], tol, iterlim);
  return u;
}
