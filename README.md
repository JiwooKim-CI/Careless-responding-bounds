# Partial-Identification Bounds for Careless Responding

This repository contains an exact-population analysis of partial-identification methods for recovering the population mean of a latent careful response when observed survey responses may be contaminated by careless responding.

## Overview

The analysis treats the latent careful response, `V_star`, as a five-point Likert variable and evaluates what can be learned about

\[
E(V^*)
\]

when true careful-responding status is not fully observed. Rather than drawing finite samples, the script enumerates the exact population distribution under each data-generating condition. The resulting intervals are population identified sets, not confidence intervals.

The analysis compares three careless-responding mechanisms:

- **CRCR:** Careless responding completely at random.
- **CRAR:** Careless responding at random conditional on an observed covariate, `C`.
- **CRNAR:** Careless responding not at random because careful-responding status depends partly on the latent response, `V_star`.

It evaluates three partial-identification strategies:

- **Manski support bounds:** Use only the bounded support of the Likert response.
- **Monotonicity bounds:** Add an ordering assumption based on the auxiliary variable `Z`.
- **Proxy bounds:** Use an additional proxy, `W`, that is related to latent careful-responding status.

The script also reports the performance of retaining all observed responses, deleting respondents not certified as careful, and standardizing the certified-response mean over observed covariates.

## Coding conventions and maintained assumptions

The variables are coded as follows:

| Variable | Meaning |
|---|---|
| `R_CR = 1` | Truly careful response |
| `R_CR = 0` | Truly careless response |
| `D = 1` | Certified or classified as careful |
| `D = 0` | Not certified as careful |
| `V_star` | Latent careful response |
| `V_obs` | Observed response |
| `U` | Carelessly generated replacement response |
| `C` | Observed predictor of `V_star` and careful responding |
| `Z` | Auxiliary predictor used for the monotonicity restriction |
| `W` | Proxy for latent careful-responding status |

The analysis maintains two central assumptions:

1. `D = 1` implies `R_CR = 1`, corresponding to one-sided detection error.
2. `R_CR = 1` implies `V_obs = V_star`, corresponding to consistency for truly careful responses.

## Repository files

| File | Description |
|---|---|
| [`careless_responding_bounds_exact_population.R`](careless_responding_bounds_exact_population.R) | Defines the exact population law, computes all bounds and point estimates, conducts the factorial and differential-certification analyses, and produces paper-oriented summaries. |
| [`record_session_info.R`](record_session_info.R) | Records the R version, platform, locale, and loaded package versions in `session-info.txt`. |
| [`session-info.txt`](session-info.txt) | Environment information recorded for the current local run. |

## Software requirements

The analysis uses functions from base and recommended R packages only; no additional package installation is currently required. It was checked locally using:

- R 4.5.3
- macOS on Apple Silicon (`aarch64-apple-darwin20`)

Because the script calls `View()`, running it interactively in RStudio is the most convenient option. The proxy analysis uses a grid of 20,001 values in the full analysis and may require some computation time.

## Running the analysis

Clone or download the repository, set the repository directory as the working directory, and run:

```r
source("careless_responding_bounds_exact_population.R")
```

The script creates the following principal objects:

| Object | Contents |
|---|---|
| `population_results` | Full factorial results across mechanisms, prevalence levels, certification rates, proxy strengths, and bounding methods |
| `paper_condition_results` | Condition-level bounds and point-estimate results |
| `paper_width_summary` | Mean, minimum, and maximum bound widths by mechanism, method, and proxy setting |
| `paper_point_summary` | Mean absolute bias of the retain-all, deletion, and standardized-deletion estimators |
| `differential_results` | Results when certification rates differ across levels of `Z` |
| `differential_table` | Paper-oriented summary of the differential-certification analysis |
| `grid_check` | Proxy-bound results across alternative grid sizes |

## Recording the R and package versions

After running the main analysis, run the following command in the **same R session**:

```r
source("record_session_info.R")
```

This creates `session-info.txt`, which records the R version and packages loaded during the analysis. The same information can be displayed without creating a file by running:

```r
sessionInfo()
```

`sessionInfo()` reports packages that are attached or loaded in the current session; it does not list every package installed on the computer. If the project later acquires external package dependencies, an `renv` lockfile can be used to record and restore their exact versions:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

After cloning the repository, another user can restore that environment with:

```r
renv::restore()
```

## Interpretation

`contains_population_truth` indicates whether a population identified set contains the known population value used to generate the exact distribution. It should not be interpreted as frequentist confidence-interval coverage. Finite-sample estimation and sampling uncertainty are outside the scope of the current script.
