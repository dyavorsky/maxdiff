// maxdiff.h — shared types and declarations for the maxdiff C++ layer.
//
// Adapted from bayesm/inst/include/bayesm.h by Peter Rossi, Wayne Taylor,
// and Robert McCulloch. All core statistical algorithms and C++ implementations
// are their work.
//
// This header retains only the structs and function declarations needed by
// rhierMnlRwMixture and its dependencies.
//
// Reference: Rossi, Allenby, and McCulloch (2024), Bayesian Statistics and
//   Marketing (2nd ed.), Wiley. ISBN 978-1-118-97218-3.

#ifndef __MAXDIFF_H__
#define __MAXDIFF_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <time.h>

using namespace arma;
using namespace Rcpp;

// Structs ----------------------------------------------------------------

struct moments {
  vec y;
  mat X;
  mat XpX;
  vec Xpy;
  mat hess;
};

struct unireg {
  vec    beta;
  double sigmasq;
};

struct mnlMetropOnceOut {
  vec    betadraw;
  int    stay;
  double oldll;
};

struct lambda {
  vec    mubar;
  double Amu;
  double nu;
  mat    V;
};

struct priorAlpha {
  double power;
  double alphamin;
  double alphamax;
  int    n;
};

struct murooti {
  vec mu;
  mat rooti;
};

struct thetaStarIndex {
  ivec indic;
  std::vector<murooti> thetaStar_vector;
};

struct DPOut {
  ivec indic;
  std::vector<murooti> thetaStar_vector;
  std::vector<murooti> thetaNp1_vector;
  double alpha;
  int    Istar;
  lambda lambda_struct;
};

// Exposed C++ functions (not wrapped to R) --------------------------------

List rwishart(double nu, mat const& V);
List rmultireg(mat const& Y, mat const& X, mat const& Bbar, mat const& A, double nu, mat const& V);
vec  rdirichlet(vec const& alpha);

double llmnl(vec const& beta, vec const& y, mat const& X,
             IntegerVector p_vec = IntegerVector::create());

mat    lndIChisq(double nu, double ssq, mat const& X);
double lndMvn(vec const& x, vec const& mu, mat const& rooti);
double lndIWishart(double nu, mat const& V, mat const& IW);
vec    rmvst(double nu, vec const& mu, mat const& root);
vec    breg(vec const& y, mat const& X, vec const& betabar, mat const& A);

List rmixGibbs(mat const& y, mat const& Bbar, mat const& A, double nu,
               mat const& V, vec const& a, vec const& p, vec const& z);

vec condmom(vec const& x, vec const& mu, mat const& sigmai, int p, int j);
double trunNorm(double mu, double sig, double trunpt, int above);

mat drawDelta(mat const& x, mat const& y, vec const& z, List const& comps,
              vec const& deltabar, mat const& Ad);

unireg runiregG(vec const& y, mat const& X, mat const& XpX, vec const& Xpy,
                double sigmasq, mat const& A, vec const& Abetabar,
                double nu, double ssq);

double llnegbin(vec const& y, vec const& lambda, double alpha, bool constant);
double lpostbeta(double alpha, vec const& beta, mat const& X, vec const& y,
                 vec const& betabar, mat const& rootA);
double lpostalpha(double alpha, vec const& beta, mat const& X, vec const& y,
                  double a, double b);

vec breg1(mat const& root, mat const& X, vec const& y, vec const& Abetabar);
vec rtrunVec(vec const& mu, vec const& sigma, vec const& a, vec const& b);
vec trunNorm_vec(vec const& mu, vec const& sig, vec const& trunpt, vec const& above);

mnlMetropOnceOut mnlMetropOnce(vec const& y, mat const& X, vec const& oldbeta,
                               double oldll, double s, mat const& incroot,
                               vec const& betabar, mat const& rootpi,
                               IntegerVector p_vec = IntegerVector::create());

int  rmultinomF(vec const& p);
mat  yden(std::vector<murooti> const& thetaStar, mat const& y);
ivec numcomp(ivec const& indic, int k);
murooti thetaD(mat const& y, lambda const& lambda_struct);
thetaStarIndex thetaStarDraw(ivec indic, std::vector<murooti> thetaStar_vector,
                             mat const& y, mat ydenmat, vec const& q0v,
                             double alpha, lambda const& lambda_struct, int maxuniq);
vec    q0(mat const& y, lambda const& lambda_struct);
vec    seq_rcpp(double from, double to, int len);
double alphaD(priorAlpha const& priorAlpha_struct, int Istar, int gridsize);
murooti GD(lambda const& lambda_struct);
lambda lambdaD(lambda const& lambda_struct,
               std::vector<murooti> const& thetaStar_vector,
               vec const& alim, vec const& nulim, vec const& vlim, int gridsize);

// MCMC timing (functionTiming.cpp) ----------------------------------------
void startMcmcTimer();
void infoMcmcTimer(int rep, int R);
void endMcmcTimer();

#endif
