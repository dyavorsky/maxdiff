# maxdiff

An R package for MaxDiff (Best-Worst Scaling) analysis using a hierarchical
Bayesian MNL model.

## Installation

```r
pak::pak("dyavorsky/maxdiff")
```

Or from source:

```r
devtools::install_github("dyavorsky/maxdiff")
```

## Overview

MaxDiff surveys ask respondents to repeatedly pick the **best** and **worst**
item from a small shown set. The `maxdiff` package estimates individual-level
partworth utilities via a Gibbs sampler (Rossi/Allenby/McCulloch HB MNL) with
an exact joint best-then-worst likelihood. 

## Quick start

```r
library(maxdiff)

# 1. Generate a design (10 items, 12 tasks, 4 items per task)
design <- make_design(n_items = 10, n_tasks = 12, n_items_per_task = 4)

# 2. Fast B-W scoring (no MCMC)
simple_fit <- fit_simple(
  responses   = responses,
  design      = design,
  version_col = "version",
  best_cols   = setNames(paste0("best",  1:12), as.character(1:12)),
  worst_cols  = setNames(paste0("worst", 1:12), as.character(1:12))
)
plot(simple_fit)

# 3. Full HB MNL estimation
hb_fit <- fit_hbmnl(
  responses   = responses,
  design      = design,
  version_col = "version",
  best_cols   = setNames(paste0("best",  1:12), as.character(1:12)),
  worst_cols  = setNames(paste0("worst", 1:12), as.character(1:12)),
  item_labels = c("Price", "Quality", "Speed", "Support", "Design",
                  "Battery", "Weight", "Screen", "Camera", "Storage"),
  R = 20000, burnin = 10000, keep = 10
)
summary(hb_fit)
plot(hb_fit)
plot(hb_fit, scale = "probability")

# 4. Subgroup analysis
sg <- subgroup_summary(hb_fit, list(
  Young = responses$age < 35,
  Old   = responses$age >= 35
))
plot(sg)
compare_groups(hb_fit, responses$age < 35, responses$age >= 35)

# 5. Publication table (requires gt)
summary_table(hb_fit)
```

## Functions

| Function | Description |
|----------|-------------|
| `make_design()` | Generate a balanced MaxDiff experimental design |
| `fit_simple()` | Best-minus-worst tally scoring (fast, no MCMC) |
| `fit_hbmnl()` | Full HB MNL Gibbs sampler estimation |
| `assemble_data()` | Low-level: build `lgtdata` from survey data |
| `subgroup_summary()` | Average partworths by subgroup |
| `compare_groups()` | Welch or Bayesian inter-group significance tests |
| `to_probability_scale()` | Sawtooth probability transform |
| `summary_table()` | `gt` publication table |
| `plot.*()` | Visualizations for all result types |

## Methodology

- **HB MNL**: random-walk Metropolis-within-Gibbs with mixture-of-normals
  heterogeneity (Rossi, Allenby, and McCulloch 2024).
- **Ragged choice sets**: enables the exact joint best-then-worst likelihood.
  Each task → 2 MNL observations: K-alternative best pick, then K-1-alternative
  worst pick (best item removed).
- **Sum-to-zero centering**: unanchored fits use M-1 reference coding; results
  are post-hoc augmented to M items and centered per respondent per draw.

## Credits

The HB MNL Gibbs sampler is adapted from Peter Rossi's
[bayesm](https://cran.r-project.org/package=bayesm) package. 

> Rossi, P. E., Allenby, G. M., & McCulloch, R. (2024).
> *Bayesian Statistics and Marketing* (2nd ed.). Wiley.

## License

MIT © Dan Yavorsky
