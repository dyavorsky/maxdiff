# Internal HB MNL Gibbs sampler.
#
# Adapted from Peter Rossi's rhierMnlRwMixture function in the bayesm package
# (version 3.1-7, https://cran.r-project.org/package=bayesm), with an
# extension by Dan Yavorsky to support ragged choice sets (list-of-X per
# occasion). All statistical development, the C++ kernels, and the core
# algorithm are Peter Rossi's work.
#
# Reference: Rossi, Allenby, and McCulloch (2024), Bayesian Statistics and
#   Marketing (2nd ed.), Wiley. ISBN 978-1-118-97218-3.
#
# The Data API differs from upstream bayesm:
#   lgtdata[[i]] = list(y = <length-T_i integer vector>,
#                       X = list(X_t1, X_t2, ...))   # list of p_it x nvar matrices
#   Data = list(lgtdata = lgtdata, Z = Z)             # no Data$p argument
#
# Choice-set sizes are derived from the data via sapply(X_list, nrow).

.rhierMnlRwMixture <- function(Data, Prior, Mcmc) {

  if (missing(Data)) pandterm("Requires Data argument -- list of lgtdata and (possibly) Z")
  if (is.null(Data$lgtdata)) pandterm("Requires Data element lgtdata (list of data for each unit)")
  lgtdata <- Data$lgtdata
  nlgt    <- length(lgtdata)
  drawdelta <- TRUE
  if (is.null(Data$Z)) {
    cat("Z not specified", fill = TRUE); fsh(); drawdelta <- FALSE
  } else {
    if (!is.matrix(Data$Z)) pandterm("Z must be a matrix")
    if (nrow(Data$Z) != nlgt) pandterm(paste("Nrow(Z)", nrow(Data$Z), "ne number logits", nlgt))
    Z <- Data$Z
  }
  if (drawdelta) {
    nz       <- ncol(Z)
    colmeans <- apply(Z, 2, mean)
    if (sum(colmeans) > 0.00001)
      pandterm(paste("Z does not appear to be de-meaned: colmeans=", colmeans))
  }

  ypooled <- NULL
  Xpooled <- NULL
  p_it    <- vector("list", nlgt)
  oldncol <- NULL
  for (i in 1:nlgt) {
    if (is.null(lgtdata[[i]]$y)) pandterm(paste0("Requires element y of lgtdata[[", i, "]]"))
    if (is.null(lgtdata[[i]]$X)) pandterm(paste0("Requires element X of lgtdata[[", i, "]]"))
    if (!is.list(lgtdata[[i]]$X) || is.data.frame(lgtdata[[i]]$X))
      pandterm(paste0("lgtdata[[", i, "]]$X must be a list of per-occasion design matrices"))
    yi    <- as.vector(lgtdata[[i]]$y)
    Xlist <- lgtdata[[i]]$X
    T_i   <- length(Xlist)
    if (T_i < 1) pandterm(paste0("lgtdata[[", i, "]]$X must be non-empty"))
    if (length(yi) != T_i)
      pandterm(paste0("length(lgtdata[[", i, "]]$y) (=", length(yi), ") must equal length(lgtdata[[", i, "]]$X) (=", T_i, ")"))
    p_i <- integer(T_i)
    for (t in seq_len(T_i)) {
      if (!is.matrix(Xlist[[t]]))
        pandterm(paste0("lgtdata[[", i, "]]$X[[", t, "]] must be a matrix"))
      p_i[t] <- nrow(Xlist[[t]])
      if (is.null(oldncol)) oldncol <- ncol(Xlist[[t]])
      if (ncol(Xlist[[t]]) != oldncol)
        pandterm(paste0("All X matrices must have same # of cols; exception at unit ", i, " occasion ", t))
    }
    if (any(p_i < 2)) pandterm(paste0("Every choice set must have >= 2 alternatives; exception at unit ", i))
    if (any(yi < 1) || any(yi > p_i))
      pandterm(paste0("lgtdata[[", i, "]]$y has an element outside [1, nrow(X[[t]])] at unit ", i))
    p_it[[i]]      <- p_i
    Xstacked       <- do.call(rbind, Xlist)
    lgtdata[[i]]$X <- Xstacked
    ypooled        <- c(ypooled, yi)
    Xpooled        <- rbind(Xpooled, Xstacked)
  }
  nvar     <- ncol(Xpooled)
  p_pooled <- unlist(p_it, use.names = FALSE)
  p_max    <- max(p_pooled)
  cat("Table of Y values pooled over all units", fill = TRUE)
  print(table(ypooled))

  if (missing(Prior)) pandterm("Requires Prior list argument (at least ncomp)")
  if (is.null(Prior$ncomp)) pandterm("Requires Prior element ncomp") else ncomp <- Prior$ncomp
  if (is.null(Prior$SignRes)) SignRes <- rep(0, nvar) else SignRes <- Prior$SignRes
  if (length(SignRes) != nvar) pandterm("Length of SignRes must equal ncol(X)")
  if (sum(!(SignRes %in% c(-1, 0, 1)) > 0)) pandterm("All elements of SignRes must be -1, 0, or 1")

  if (is.null(Prior$mubar) & sum(abs(SignRes)) == 0) {
    mubar <- matrix(rep(0, nvar), nrow = 1)
  } else if (is.null(Prior$mubar) & sum(abs(SignRes)) > 0) {
    mubar <- matrix(rep(0, nvar) + 2 * abs(SignRes), nrow = 1)
  } else {
    mubar <- matrix(Prior$mubar, nrow = 1)
  }
  if (ncol(mubar) != nvar) pandterm(paste("mubar must have nvar cols, ncol(mubar)=", ncol(mubar)))

  if (is.null(Prior$Amu) & sum(abs(SignRes)) == 0) {
    Amu <- matrix(.mdConst$A, ncol = 1)
  } else if (is.null(Prior$Amu) & sum(abs(SignRes)) > 0) {
    Amu <- matrix(.mdConst$A * 10, ncol = 1)
  } else {
    Amu <- matrix(Prior$Amu, ncol = 1)
  }
  if (ncol(Amu) != 1 | nrow(Amu) != 1) pandterm("Amu must be a 1 x 1 array")

  if (is.null(Prior$nu) & sum(abs(SignRes)) == 0) {
    nu <- nvar + .mdConst$nuInc
  } else if (is.null(Prior$nu) & sum(abs(SignRes)) > 0) {
    nu <- nvar + .mdConst$nuInc + 12
  } else {
    nu <- Prior$nu
  }
  if (nu < 1) pandterm("invalid nu value")

  if (is.null(Prior$V) & sum(abs(SignRes)) == 0) {
    V <- nu * diag(nvar)
  } else if (is.null(Prior$V) & sum(abs(SignRes)) > 0) {
    V <- nu * diag(abs(SignRes) * 0.1 + (!abs(SignRes)) * 4)
  } else {
    V <- Prior$V
  }
  if (sum(dim(V) == c(nvar, nvar)) != 2) pandterm("Invalid V in prior")

  if (is.null(Prior$Ad) & drawdelta) Ad <- .mdConst$A * diag(nvar * nz) else Ad <- Prior$Ad
  if (drawdelta && (ncol(Ad) != nvar * nz | nrow(Ad) != nvar * nz))
    pandterm("Ad must be nvar*nz x nvar*nz")
  if (is.null(Prior$deltabar) & drawdelta) deltabar <- rep(0, nz * nvar) else deltabar <- Prior$deltabar
  if (drawdelta && length(deltabar) != nz * nvar) pandterm("deltabar must be of length nvar*nz")
  if (is.null(Prior$a)) a <- rep(.mdConst$a, ncomp) else a <- Prior$a
  if (length(a) != ncomp) pandterm("Requires dim(a) = ncomp")
  if (any(a < 0)) pandterm("invalid values in a vector")

  if (is.null(Prior$nu) & sum(abs(SignRes)) > 0) nu <- nvar + 15
  if (is.null(Prior$Amu) & sum(abs(SignRes)) > 0) Amu <- matrix(0.1)
  if (is.null(Prior$V) & sum(abs(SignRes)) > 0) V <- nu * (diag(nvar) - diag(abs(SignRes) > 0) * 0.8)

  if (missing(Mcmc)) pandterm("Requires Mcmc list argument")
  if (is.null(Mcmc$s)) s <- .mdConst$RRScaling / sqrt(nvar) else s <- Mcmc$s
  if (is.null(Mcmc$w)) w <- .mdConst$w else w <- Mcmc$w
  if (is.null(Mcmc$keep)) keep <- .mdConst$keep else keep <- Mcmc$keep
  if (is.null(Mcmc$R)) pandterm("Requires R argument in Mcmc list") else R <- Mcmc$R
  if (is.null(Mcmc$nprint)) nprint <- .mdConst$nprint else nprint <- Mcmc$nprint
  if (nprint < 0) pandterm("nprint must be >= 0")

  cat(" ", fill = TRUE)
  cat("Starting MCMC Inference for Hierarchical Logit:", fill = TRUE)
  cat("   Normal Mixture with", ncomp, "components for first stage prior", fill = TRUE)
  cat(paste("   up to", p_max, "alternatives per occasion;", nvar, "variables in X"), fill = TRUE)
  cat(paste("   for", nlgt, "cross-sectional units"), fill = TRUE)
  cat(" ", fill = TRUE)
  cat("Prior Parms: ", fill = TRUE)
  cat("nu =", nu, fill = TRUE)
  cat("V ", fill = TRUE); print(V)
  cat("mubar ", fill = TRUE); print(mubar)
  cat("Amu ", fill = TRUE); print(Amu)
  cat("a ", fill = TRUE); print(a)
  if (drawdelta) {
    cat("deltabar", fill = TRUE); print(deltabar)
    cat("Ad", fill = TRUE); print(Ad)
  }
  if (sum(abs(SignRes)) != 0) {
    cat("Sign Restrictions Vector (0: unconstrained, 1: positive, -1: negative)", fill = TRUE)
    print(matrix(SignRes, ncol = 1))
  }
  cat(" ", fill = TRUE)
  cat("MCMC Parms: ", fill = TRUE)
  cat(paste("s=", round(s, 3), " w=", w, " R=", R, " keep=", keep, " nprint=", nprint), fill = TRUE)
  cat("", fill = TRUE)

  oldbetas <- matrix(double(nlgt * nvar), ncol = nvar)

  llmnlFract <- function(beta, y, X, betapooled, rootH, w, wgt,
                         SignRes = rep(0, ncol(X)), p_vec = integer(0)) {
    if (is.null(p_vec)) p_vec <- integer(0)
    z <- as.vector(rootH %*% (beta - betapooled))
    (1 - w) * llmnl_con(beta, y, X, SignRes, p_vec) + w * wgt * (-.5 * (z %*% z))
  }

  .mnlHess_con <- function(betastar, y, X, SignRes = rep(0, ncol(X)), p_vec = NULL) {
    beta <- betastar
    beta[SignRes != 0] <- SignRes[SignRes != 0] * exp(betastar[SignRes != 0])
    n    <- length(y)
    k    <- ncol(X)
    Xbeta <- X %*% beta
    Hess  <- matrix(double(k * k), ncol = k)
    if (is.null(p_vec) || length(p_vec) == 0) {
      j    <- nrow(X) / n
      Xbm  <- matrix(Xbeta, byrow = TRUE, ncol = j)
      Xbm  <- exp(Xbm)
      iota <- c(rep(1, j))
      denom <- Xbm %*% iota
      Prob  <- Xbm / as.vector(denom)
      for (i in 1:n) {
        p    <- as.vector(Prob[i, ])
        A    <- diag(p) - outer(p, p)
        Xt   <- X[(j * (i - 1) + 1):(j * i), , drop = FALSE]
        Hess <- Hess + crossprod(Xt, A) %*% Xt
      }
    } else {
      offset <- 0L
      for (i in 1:n) {
        j   <- p_vec[i]
        idx <- (offset + 1L):(offset + j)
        eb  <- exp(Xbeta[idx])
        p   <- as.vector(eb / sum(eb))
        A   <- diag(p) - outer(p, p)
        Xt  <- X[idx, , drop = FALSE]
        Hess <- Hess + crossprod(Xt, A) %*% Xt
        offset <- offset + j
      }
    }
    lambda <- c(rep(1, length(SignRes)))
    lambda[SignRes == 1]  <- beta[SignRes == 1]
    lambda[SignRes == -1] <- -beta[SignRes == -1]
    Hess <- Hess * crossprod(t(lambda))
    return(Hess)
  }

  cat("initializing Metropolis candidate densities for", nlgt, "units ...", fill = TRUE)
  fsh()

  betainit <- c(rep(0, nvar))
  noRes    <- c(rep(0, nvar))
  out <- optim(betainit, llmnl_con, method = "BFGS",
               control = list(fnscale = -1, trace = 0, reltol = 1e-6),
               X = Xpooled, y = ypooled, SignRes = noRes, p_vec = p_pooled)
  betainit <- out$par
  betainit[SignRes != 0] <- 0
  out <- optim(betainit, llmnl_con,
               control = list(fnscale = -1, trace = 0, reltol = 1e-6),
               X = Xpooled, y = ypooled, SignRes = SignRes, p_vec = p_pooled)
  betapooled <- out$par
  if (sum(abs(betapooled[as.logical(SignRes)]) > 10)) {
    cat("In tuning Metropolis, constrained pooled estimates contain very large values", fill = TRUE)
    print(cbind(betapooled, SignRes))
  }
  H     <- .mnlHess_con(betapooled, ypooled, Xpooled, SignRes, p_vec = p_pooled)
  rootH <- chol(H)
  for (i in 1:nlgt) {
    wgt <- length(lgtdata[[i]]$y) / length(ypooled)
    out <- optim(betapooled, llmnlFract, method = "BFGS",
                 control = list(fnscale = -1, trace = 0, reltol = 1e-4),
                 X = lgtdata[[i]]$X, y = lgtdata[[i]]$y,
                 betapooled = betapooled, rootH = rootH, w = w, wgt = wgt,
                 SignRes = SignRes, p_vec = p_it[[i]])
    if (out$convergence == 0) {
      hess <- .mnlHess_con(out$par, lgtdata[[i]]$y, lgtdata[[i]]$X, SignRes, p_vec = p_it[[i]])
      lgtdata[[i]] <- c(lgtdata[[i]], list(converge = 1, betafmle = out$par, hess = hess))
    } else {
      lgtdata[[i]] <- c(lgtdata[[i]], list(converge = 0, betafmle = c(rep(0, nvar)), hess = diag(nvar)))
    }
    oldbetas[i, ] <- lgtdata[[i]]$betafmle
    if (i %% 50 == 0) cat("  completed unit #", i, fill = TRUE)
    fsh()
  }

  ind    <- NULL
  ninc   <- floor(nlgt / ncomp)
  for (i in 1:(ncomp - 1)) ind <- c(ind, rep(i, ninc))
  if (ncomp != 1) ind <- c(ind, rep(ncomp, nlgt - length(ind))) else ind <- rep(1, nlgt)
  oldprob <- rep(1 / ncomp, ncomp)
  if (drawdelta) olddelta <- rep(0, nz * nvar) else { olddelta <- 0; Z <- matrix(0); deltabar <- 0; Ad <- matrix(0) }

  draws <- rhierMnlRwMixture_rcpp_loop(
    lgtdata, Z, deltabar, Ad, mubar, Amu,
    nu, V, s, R, keep, nprint, drawdelta,
    as.matrix(olddelta), a, oldprob, oldbetas, ind, SignRes, p_it
  )

  if (drawdelta) {
    attributes(draws$Deltadraw)$class <- c("bayesm.mat", "mcmc")
    attributes(draws$Deltadraw)$mcpar <- c(1, R, keep)
  }
  attributes(draws$betadraw)$class <- c("bayesm.hcoef")
  attributes(draws$nmix)$class     <- "bayesm.nmix"
  draws
}
