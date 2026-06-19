# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package overview

`maxdiff` is a self-contained R package for MaxDiff (Best-Worst Scaling) analysis using a hierarchical Bayesian MNL (HB MNL) model. It is a personal package by Dan Yavorsky, hosted at github.com/dyavorsky/maxdiff.

Key design principles:
- **Self-contained**: no dependency on `bayesm` or `gbkmaxdiff`; the C++ sampler is copied directly into `src/`
- **Exact best-then-worst likelihood**: each MaxDiff task encodes as two MNL observations (K-alt best, K-1-alt worst) via a ragged choice-set structure
- **HB MNL algorithm**: Rossi/Allenby/McCulloch Gibbs sampler with mixture-of-normals heterogeneity; all statistical credit to Peter Rossi

Install from GitHub via `pak::pak("dyavorsky/maxdiff")`.

## Commands

```r
# Standard R package workflow
devtools::load_all()        # Load interactively during development
devtools::document()        # Regenerate NAMESPACE and man/ from roxygen2
devtools::check()           # Full R CMD check
devtools::build()           # Build .tar.gz

# After editing any src/*.cpp signature that has //[[Rcpp::export]]:
Rcpp::compileAttributes()   # Regenerates src/RcppExports.cpp and R/RcppExports.R
# IMPORTANT: after compileAttributes(), manually add this line back to src/RcppExports.cpp
# (right before #include <RcppArmadillo.h>):
#   #include "maxdiff.h"
# compileAttributes does not include it; without it, vec/mat types won't resolve.
# NOTE: currently only llmnl_con and rhierMnlRwMixture_rcpp_loop are exported

# Install with compilation
devtools::install()

# Run the vignette locally
quarto::quarto_render("vignettes/getting-started.qmd")
```

## Architecture

### R layer

```
fit_hbmnl()          ← main user entry point for HB estimation
  └── assemble_data()   ← converts survey data to lgtdata (ragged X)
  └── .rhierMnlRwMixture()   ← internal R wrapper (R/mnl_sampler.R)
        └── rhierMnlRwMixture_rcpp_loop()  ← C++ Gibbs loop

fit_simple()         ← fast B-W tally scoring (no MCMC)
make_design()        ← experimental design via choiceDes::tradeoff.des()
subgroup_summary()   ← per-subgroup partworth averages
compare_groups()     ← Welch or Bayesian inter-group tests
to_probability_scale()  ← Sawtooth probability transform
summary_table()      ← gt publication table
plot.*()             ← ggplot2 visuals for all three S3 classes
```

### C++ layer (src/)

Only two functions are exported to R (via RcppExports):
- `llmnl_con` — MNL log-likelihood with sign constraints and ragged p_vec
- `rhierMnlRwMixture_rcpp_loop` — main Gibbs loop

Everything else is C++-internal:
- `llmnl` — basic MNL log-likelihood (called by `mnlMetropOnce` in utilityFunctions.cpp)
- `rmixGibbs`, `rmultireg`, `rwishart`, `rdirichlet`, `lndMvn`, `breg` — sampler building blocks
- `drawDelta`, `runiregG`, `mnlMetropOnce`, DP utility functions — shared helpers
- Timer functions in `functionTiming.cpp`

The shared header `maxdiff.h` declares all C++ functions used across translation units.

### Data structure (ragged choice sets)

```
lgtdata[[i]] = list(
  y    = integer vector length 2*T,     # choice index within each occasion
  X    = stacked design matrix,         # (sum_t p_it) x M
  hess = M x M Hessian matrix           # pre-computed at MAP for proposal scaling
)
p_it_list[[i]] = integer vector length 2*T  # per-occasion choice-set sizes
```

Each MaxDiff task t for respondent i contributes two occasions:
- Best pick: K alternatives, dummy X with +1 in chosen-item column
- Worst pick: K-1 alternatives (best item removed), dummy X with -1

### S3 classes

| Class | Constructor | Methods |
|-------|-------------|---------|
| `maxdiff_hbmnl` | `fit_hbmnl()` | print, summary, coef, plot |
| `maxdiff_simple` | `fit_simple()` | print, plot |
| `maxdiff_subgroup` | `subgroup_summary()` | print, plot |

### Sum-to-zero centering

For unanchored studies, the sampler estimates M-1 free parameters (reference
coded). After burnin discard, `fit_hbmnl()` augments to M and applies
per-respondent, per-draw sum-to-zero centering:
```r
for (d in seq_len(n_kept)) {
  slab <- betadraw[, , d]
  betadraw[, , d] <- slab - rowMeans(slab)
}
```

## Key files

| File | Purpose |
|------|---------|
| `R/mnl_sampler.R` | Internal `.rhierMnlRwMixture()` — R wrapper around the C++ loop |
| `R/assemble_data.R` | Builds `lgtdata` with ragged best+worst encoding |
| `R/fit_hbmnl.R` | User-facing `fit_hbmnl()`, S3 methods |
| `R/constants.R` | `.mdConst` — MCMC/prior defaults |
| `src/rhierMnlRwMixture_rcpp_loop.cpp` | Main Gibbs loop + `llmnl_con` |
| `src/maxdiff.h` | Shared C++ header with structs and declarations |

## Credits

The HB MNL sampler is adapted from `bayesm::rhierMnlRwMixture` (Peter Rossi,
Wayne Taylor, Robert McCulloch). The ragged choice-set extension was added by
Dan Yavorsky for exact best-then-worst likelihood encoding.
