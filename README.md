# toda.yamamoto

A simple R implementation of the modified Wald test for Granger non-causality in possibly integrated or cointegrated VAR systems following Toda & Yamamoto (1995).

## Description

This function implements the Toda-Yamamoto procedure for testing Granger causality in levels VAR models. The approach estimates a (*p* + *d_max*)-th order lag-augmented VAR (LA-VAR) and computes the Wald statistic over the first *p* coefficient matrices only, ignoring the last *d_max* augmentation lags. This guarantees an asymptotic chi-square distribution with *p* degrees of freedom regardless of the integration or cointegration properties of the data.

## Reference

Toda, H. Y., & Yamamoto, T. (1995). Statistical inference in vector autoregressions with possibly integrated processes. *Journal of Econometrics*, 66(1-2), 225-250.

## Installation

```r
# Install dependencies
install.packages(c("vars", "aod"))

# Load the function
source("toda_yamamoto.R")
```

## Usage

```r
toda.yamamoto(var.model, d.max = 1L, verbose = FALSE)
```

## Arguments

| Argument | Description |
|---|---|
| `var.model` | A `varest` object estimated with `vars::VAR()`. See details below. |
| `d.max` | Maximum order of integration suspected in the process (default 1). Must be set by the user based on prior unit root testing. |
| `verbose` | Logical. If `TRUE`, prints the exact restrictions applied in each test (default `FALSE`). |

## Details

### How to prepare the input model

The function receives a VAR model estimated with `p - 1` lags — that is, one lag fewer than the true VAR order — and internally estimates the LA-VAR with `(p - 1) + d_max` total lags. This design allows the user to select the optimal lag length via information criteria and then pass that model directly to the function.

**Step 1.** Select the optimal lag length using an information criterion:

```r
library(vars)

lag_sel <- VARselect(Y, lag.max = 12, type = "const")
p_opt <- lag_sel$selection["AIC(n)"]
```

**Step 2.** Estimate the VAR with `p_opt` lags and verify that the residuals show no serial correlation:

```r
var_mod <- VAR(Y, p = p_opt, type = "const")
serial.test(var_mod, lags.pt = 12, type = "BG")
```

**Step 3.** Pass the model to `toda.yamamoto()` and specify `d.max` based on your unit root tests. The function will automatically estimate the LA-VAR with `p_opt + d.max` total lags:

```r
res <- toda.yamamoto(var_mod, d.max = 1)
```

### What d.max to use

`d.max` is the maximum order of integration suspected in any variable in the system. It must be determined externally via unit root tests (ADF, KPSS, PP) before calling the function:

- If all variables are I(1): `d.max = 1`
- If any variable may be I(2): `d.max = 2`

### Restriction logic

The Wald statistic is computed over `seq_len(p_opt)` — exactly the first `p_opt` coefficient lags of the causal variable. The augmentation lags `p_opt + 1` through `p_opt + d_max` are included in the LA-VAR estimation but are not restricted, consistent with Toda & Yamamoto (1995, p. 230).

## Output

The function prints to console:

- Wald test results for all directional pairs with significance stars
- Adjusted R-squared per equation
- Multivariate Breusch-Godfrey LM test with 12 lags

And returns invisibly a named list:

| Element | Description |
|---|---|
| `ty.results` | `data.frame` with `cause`, `effect`, `chisq`, `pvalue` |
| `ty.augmented_var` | `varest` object of the estimated LA-VAR(*p* + *d_max*) |
| `ty.wald` | List of `wald.test` objects, one per directional pair |
| `ty.bg` | `serial.test` object (BG test, 12 lags) |
| `ty.regressors` | Character vector of LA-VAR regressor names |
| `ty.r2` | Named vector of adjusted R-squared per equation |

## Example

```r
library(vars)
source("toda_yamamoto.R")

# Simulate two I(1) series
set.seed(42)
Y <- apply(matrix(rnorm(200 * 2), 200, 2), 2, cumsum)
colnames(Y) <- c("x1", "x2")

# Step 1: select lag length
lag_sel <- VARselect(Y, lag.max = 12, type = "const")
p_opt   <- lag_sel$selection["AIC(n)"]

# Step 2: estimate VAR and check serial correlation
var_mod <- VAR(Y, p = p_opt, type = "const")
serial.test(var_mod, lags.pt = 12, type = "BG")

# Step 3: run Toda-Yamamoto with d_max = 1
res <- toda.yamamoto(var_mod, d.max = 1)

# Access results
res$ty.results
res$ty.augmented_var
res$ty.r2
```

## Dependencies

- [`vars`](https://CRAN.R-project.org/package=vars) — VAR estimation and serial correlation test
- [`aod`](https://CRAN.R-project.org/package=aod) — Wald test

---

> **Note:** This function is provided as-is for research and educational purposes. The user is solely responsible for the correct specification of the model, including the selection of the lag length *p* and the maximum order of integration *d_max*. Any errors arising from incorrect model specification, inappropriate use of the function, or misinterpretation of results are the sole responsibility of the user.
