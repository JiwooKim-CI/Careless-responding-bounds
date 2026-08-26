# ============================================================
# Exact-population analysis for bounds on a latent careful
# Likert response
#
# Coding:
#   R_CR = 1: truly careful
#   R_CR = 0: truly careless
#   D    = 1: certified/classified as careful
#   D    = 0: not certified as careful
#
# Maintained assumptions:
#   D = 1 implies R_CR = 1 (one-sided detection error)
#   R_CR = 1 implies V_obs = V_star (consistency)
#
# Target:
#   E(V_star)
#
# This script evaluates exact population identified sets. It does
# not draw finite samples and does not estimate confidence intervals.
# ============================================================


# ------------------------------------------------------------
# 1. Basic helpers
# ------------------------------------------------------------

calibrate_weighted_logit_intercept <- function(
    linear_part,
    weights,
    target_probability) {

  stopifnot(target_probability > 0, target_probability < 1)
  stopifnot(length(linear_part) == length(weights))

  weights <- weights / sum(weights)

  objective <- function(intercept) {
    sum(weights * plogis(intercept + linear_part)) -
      target_probability
  }

  uniroot(objective, interval = c(-30, 30))$root
}


# C and Z both predict V_star. The normal density is evaluated at
# the Likert categories and normalized to obtain categorical
# probabilities.
latent_response_probabilities <- function(values, c_value, z_value) {
  # Keep E(V*) above the center of U so that contamination does not appear
  # harmless merely because the two distributions have the same mean.
  location <- 3.5 +
    0.35 * (z_value - 0.5) +
    0.30 * (c_value - 0.5)

  probabilities <- dnorm(values, mean = location, sd = 1)
  probabilities / sum(probabilities)
}


replacement_response_probabilities <- function(values) {
  probabilities <- dnorm(values, mean = mean(values), sd = 0.90)
  probabilities / sum(probabilities)
}


# ------------------------------------------------------------
# 2. Exact population law
# ------------------------------------------------------------

exact_likert_cr_law <- function(
    lower = 1,
    upper = 5,
    careful_prevalence = 0.70,
    careful_certification_rate = 0.70,
    certification_rate_z0 = NULL,
    certification_rate_z1 = NULL,
    p_w1_careful = 0.30,
    p_w1_careless = 0.70,
    mechanism = c("CRCR", "CRAR", "CRNAR")) {

  mechanism <- match.arg(mechanism)
  values <- lower:upper

  stopifnot(length(values) >= 2L)
  stopifnot(careful_prevalence > 0, careful_prevalence < 1)
  stopifnot(p_w1_careful > 0, p_w1_careful < 1)
  stopifnot(p_w1_careless > 0, p_w1_careless < 1)
  stopifnot(p_w1_careful != p_w1_careless)

  if (is.null(certification_rate_z0)) {
    certification_rate_z0 <- careful_certification_rate
  }
  if (is.null(certification_rate_z1)) {
    certification_rate_z1 <- careful_certification_rate
  }

  stopifnot(
    certification_rate_z0 >= 0,
    certification_rate_z0 <= 1,
    certification_rate_z1 >= 0,
    certification_rate_z1 <= 1
  )

  p_c <- c(`0` = 0.50, `1` = 0.50)
  p_z <- c(`0` = 0.50, `1` = 0.50)

  # Joint population law of C, Z, and V_star.
  state_cv <- expand.grid(
    C = 0:1,
    Z = 0:1,
    V_star = values
  )

  state_cv$weight <- mapply(
    function(c_i, z_i, v_i) {
      p_c[as.character(c_i)] *
        p_z[as.character(z_i)] *
        latent_response_probabilities(
          values = values,
          c_value = c_i,
          z_value = z_i
        )[match(v_i, values)]
    },
    state_cv$C,
    state_cv$Z,
    state_cv$V_star
  )

  state_cv$weight <- state_cv$weight / sum(state_cv$weight)

  mean_v <- sum(state_cv$weight * state_cv$V_star)
  sd_v <- sqrt(
    sum(state_cv$weight * (state_cv$V_star - mean_v)^2)
  )

  # Careful-response mechanism.
  linear_part <- switch(
    mechanism,
    CRCR = rep(0, nrow(state_cv)),
    CRAR = 0.90 * (state_cv$C - 0.50),
    CRNAR =
      0.60 * (state_cv$C - 0.50) +
      0.90 * (state_cv$V_star - mean_v) / sd_v
  )

  intercept <- calibrate_weighted_logit_intercept(
    linear_part = linear_part,
    weights = state_cv$weight,
    target_probability = careful_prevalence
  )

  state_cv$p_careful <- plogis(intercept + linear_part)

  replacement_prob <- replacement_response_probabilities(values)
  names(replacement_prob) <- as.character(values)

  # Enumerate the full law. U is integrated over even when R_CR=1;
  # it then has no effect on V_obs.
  states <- merge(
    state_cv,
    expand.grid(
      R_CR = 0:1,
      U = values,
      D = 0:1,
      W = 0:1
    ),
    by = NULL
  )

  states$p_r <- ifelse(
    states$R_CR == 1L,
    states$p_careful,
    1 - states$p_careful
  )

  states$p_u <- replacement_prob[as.character(states$U)]

  # One-sided detector: D=1 is impossible for truly careless
  # responses. Among truly careful responses, certification may
  # depend on Z.
  certification_probability <- ifelse(
    states$Z == 0L,
    certification_rate_z0,
    certification_rate_z1
  )

  states$p_d <- ifelse(
    states$R_CR == 0L,
    as.numeric(states$D == 0L),
    ifelse(
      states$D == 1L,
      certification_probability,
      1 - certification_probability
    )
  )

  # W depends only on R_CR by construction.
  p_w1 <- ifelse(
    states$R_CR == 1L,
    p_w1_careful,
    p_w1_careless
  )

  states$p_w <- ifelse(
    states$W == 1L,
    p_w1,
    1 - p_w1
  )

  states$V_obs <- ifelse(
    states$R_CR == 1L,
    states$V_star,
    states$U
  )

  states$probability <- with(
    states,
    weight * p_r * p_u * p_d * p_w
  )

  states <- states[states$probability > 0, ]
  states$probability <- states$probability /
    sum(states$probability)

  stopifnot(abs(sum(states$probability) - 1) < 1e-12)
  stopifnot(all(states$D == 0L | states$R_CR == 1L))
  stopifnot(all(states$R_CR == 0L | states$V_obs == states$V_star))

  true_mean_by_z <- sapply(0:1, function(z) {
    in_z <- states$Z == z
    sum(states$probability[in_z] * states$V_star[in_z]) /
      sum(states$probability[in_z])
  })

  true_mean_by_c <- sapply(0:1, function(c_value) {
    in_c <- states$C == c_value
    sum(states$probability[in_c] * states$V_star[in_c]) /
      sum(states$probability[in_c])
  })

  list(
    full_law = states,
    true_mean = sum(states$probability * states$V_star),
    true_mean_by_z = true_mean_by_z,
    true_mean_by_c = true_mean_by_c,
    true_q0 = 1 - p_w1_careless,
    true_q1 = 1 - p_w1_careful,
    p_c = unname(p_c),
    values = values,
    mechanism = mechanism,
    careful_prevalence = careful_prevalence,
    certification_rate_z0 = certification_rate_z0,
    certification_rate_z1 = certification_rate_z1,
    proxy_parameters = c(
      p_w1_careful = p_w1_careful,
      p_w1_careless = p_w1_careless
    )
  )
}


# ------------------------------------------------------------
# 3. Support and monotonicity bounds
# ------------------------------------------------------------

population_support_bounds <- function(
    exact_law,
    lower = 1,
    upper = 5) {

  states <- exact_law$full_law

  known <- states$D == 1L
  known_contribution <- sum(
    states$probability[known] * states$V_obs[known]
  )
  p_d0 <- sum(states$probability[states$D == 0L])

  c(
    LB = known_contribution + p_d0 * lower,
    UB = known_contribution + p_d0 * upper
  )
}


population_monotone_bounds <- function(
    exact_law,
    lower = 1,
    upper = 5) {

  states <- exact_law$full_law

  stratum_bounds <- function(z) {
    in_z <- states$Z == z
    known <- in_z & states$D == 1L

    p_z <- sum(states$probability[in_z])
    p_d1_z <- sum(states$probability[known]) / p_z

    if (sum(states$probability[known]) <= 0) {
      return(c(a = lower, b = upper, p_z = p_z))
    }

    anchor_mean <- sum(
      states$probability[known] * states$V_obs[known]
    ) / sum(states$probability[known])

    c(
      a = p_d1_z * anchor_mean + (1 - p_d1_z) * lower,
      b = p_d1_z * anchor_mean + (1 - p_d1_z) * upper,
      p_z = p_z
    )
  }

  z0 <- stratum_bounds(0L)
  z1 <- stratum_bounds(1L)

  # No values satisfy theta_1 >= theta_0.
  if (z0["a"] > z1["b"]) {
    return(c(LB = NA_real_, UB = NA_real_, feasible = 0))
  }

  c(
    LB = unname(
      z0["p_z"] * z0["a"] +
        z1["p_z"] * max(z1["a"], z0["a"])
    ),
    UB = unname(
      z0["p_z"] * min(z0["b"], z1["b"]) +
        z1["p_z"] * z1["b"]
    ),
    feasible = 1
  )
}


# ------------------------------------------------------------
# 4. Observed point estimates
# ------------------------------------------------------------

population_point_estimates <- function(exact_law) {
  states <- exact_law$full_law

  true_mean <- sum(states$probability * states$V_star)
  retain_all <- sum(states$probability * states$V_obs)

  in_d1 <- states$D == 1L
  p_d1 <- sum(states$probability[in_d1])

  deletion <- if (p_d1 > 0) {
    sum(
      states$probability[in_d1] * states$V_obs[in_d1]
    ) / p_d1
  } else {
    NA_real_
  }

  # Standardize the mean among respondents certified as careful over the
  # population distribution of the fully observed covariates C and Z:
  #
  #   sum_{c,z} E(V_obs | D = 1, C = c, Z = z) P(C = c, Z = z).
  #
  # Including Z makes this estimator valid not only when certification is
  # constant, but also in the differential-certification analysis where the
  # probability of D = 1 varies across Z. Under CRAR, this recovers E(V_star)
  # when conditional exchangeability, positivity, consistency, and one-sided
  # detection error hold. It need not recover the target under CRNAR.
  population_cz <- aggregate(
    probability ~ C + Z,
    data = states,
    FUN = sum
  )

  certified <- states[in_d1, , drop = FALSE]
  certified$weighted_response <-
    certified$probability * certified$V_obs

  certified_cz <- aggregate(
    cbind(probability, weighted_response) ~ C + Z,
    data = certified,
    FUN = sum
  )

  certified_cz$conditional_mean <- with(
    certified_cz,
    weighted_response / probability
  )

  standardized_cells <- merge(
    population_cz,
    certified_cz[, c("C", "Z", "conditional_mean")],
    by = c("C", "Z"),
    all.x = TRUE
  )

  standardized_deletion <- if (
    all(is.finite(standardized_cells$conditional_mean))
  ) {
    sum(
      standardized_cells$probability *
        standardized_cells$conditional_mean
    )
  } else {
    NA_real_
  }

  data.frame(
    true_mean = true_mean,
    retain_all_estimate = retain_all,
    deletion_estimate = deletion,
    standardized_deletion_estimate = standardized_deletion,
    retain_all_bias = retain_all - true_mean,
    deletion_bias = deletion - true_mean,
    standardized_deletion_bias = standardized_deletion - true_mean
  )
}


# ------------------------------------------------------------
# 5. Proxy restoration
# ------------------------------------------------------------

exact_joint_w_c_d_v <- function(exact_law) {
  states <- exact_law$full_law
  values <- exact_law$values

  output <- array(
    0,
    dim = c(2, 2, 2, length(values)),
    dimnames = list(
      W = c("0", "1"),
      C = c("0", "1"),
      D = c("0", "1"),
      V_obs = as.character(values)
    )
  )

  for (i in seq_len(nrow(states))) {
    output[
      as.character(states$W[i]),
      as.character(states$C[i]),
      as.character(states$D[i]),
      as.character(states$V_obs[i])
    ] <- output[
      as.character(states$W[i]),
      as.character(states$C[i]),
      as.character(states$D[i]),
      as.character(states$V_obs[i])
    ] + states$probability[i]
  }

  output
}


proxy_relationship_bounds_from_law <- function(p_wcdv) {
  p_w_d <- apply(p_wcdv, c(1, 3), sum)

  p_d1 <- sum(p_w_d[, "1"])
  p_d0 <- sum(p_w_d[, "0"])

  if (p_d1 <= 0 || p_d0 <= 0) {
    stop("Both D levels must have positive probability.")
  }

  a <- p_w_d["0", "1"] / p_d1
  b <- p_w_d["0", "0"] / p_d0

  q0_bound <- if (b > a) {
    c(LB = b, UB = 1)
  } else if (b < a) {
    c(LB = 0, UB = b)
  } else {
    c(LB = 0, UB = 1)
  }

  list(
    q1 = a,
    q0_bound = q0_bound,
    observed_p_w0_d1 = a,
    observed_p_w0_d0 = b
  )
}


restore_r_c_d_v <- function(p_wcdv, q0, q1) {
  proxy_matrix <- matrix(
    c(q0, q1, 1 - q0, 1 - q1),
    nrow = 2,
    byrow = TRUE
  )

  if (abs(det(proxy_matrix)) < 1e-12) {
    return(NULL)
  }

  p_w_columns <- matrix(p_wcdv, nrow = 2)
  p_r_columns <- solve(proxy_matrix, p_w_columns)

  array(
    p_r_columns,
    dim = c(
      2,
      dim(p_wcdv)[2],
      dim(p_wcdv)[3],
      dim(p_wcdv)[4]
    ),
    dimnames = list(
      R_CR = c("0", "1"),
      C = c("0", "1"),
      D = c("0", "1"),
      V_obs = dimnames(p_wcdv)[[4]]
    )
  )
}


candidate_target <- function(
    p_rcdv,
    values,
    mechanism_assumption = c("general", "CRCR", "CRAR"),
    p_c,
    lower,
    upper,
    tolerance = 1e-10) {

  mechanism_assumption <- match.arg(mechanism_assumption)

  if (min(p_rcdv) < -tolerance) return(NULL)
  if (sum(p_rcdv["0", , "1", ]) > tolerance) return(NULL)

  p_rcdv[p_rcdv < 0 & p_rcdv >= -tolerance] <- 0

  if (abs(sum(p_rcdv) - 1) > 1e-8) return(NULL)

  p_r0 <- sum(p_rcdv["0", , , ])
  p_r1 <- sum(p_rcdv["1", , , ])
  if (p_r1 <= 0) return(NULL)

  p_r1_v <- apply(
    p_rcdv["1", , , , drop = FALSE],
    4,
    sum
  )

  careful_joint_contribution <- sum(values * p_r1_v)

  if (mechanism_assumption == "general") {
    return(c(
      LB = careful_joint_contribution + p_r0 * lower,
      UB = careful_joint_contribution + p_r0 * upper,
      p_r0 = p_r0
    ))
  }

  if (mechanism_assumption == "CRCR") {
    mean_value <- careful_joint_contribution / p_r1
    return(c(
      LB = mean_value,
      UB = mean_value,
      p_r0 = p_r0
    ))
  }

  # CRAR: recover the careful-response mean within C and
  # standardize over the population distribution of C.
  conditional_means <- numeric(2)

  for (c_value in 0:1) {
    block <- p_rcdv[
      "1",
      as.character(c_value),
      ,
      ,
      drop = FALSE
    ]

    p_r1_c <- sum(block)
    if (p_r1_c <= 0) return(NULL)

    p_v_r1_c <- apply(block, 4, sum)
    conditional_means[c_value + 1L] <-
      sum(values * p_v_r1_c) / p_r1_c
  }

  mean_value <- sum(p_c * conditional_means)

  c(
    LB = mean_value,
    UB = mean_value,
    p_r0 = p_r0
  )
}


proxy_mean_bounds_from_population_law <- function(
    exact_law,
    mechanism_assumption = c("general", "CRCR", "CRAR"),
    lower = 1,
    upper = 5,
    grid_n = 20001) {

  mechanism_assumption <- match.arg(mechanism_assumption)

  p_wcdv <- exact_joint_w_c_d_v(exact_law)
  relationship <- proxy_relationship_bounds_from_law(p_wcdv)
  q1 <- relationship$q1

  q0_grid <- seq(
    relationship$q0_bound["LB"],
    relationship$q0_bound["UB"],
    length.out = grid_n
  )

  q0_grid <- unique(c(q0_grid, relationship$q0_bound))
  q0_grid <- q0_grid[abs(q0_grid - q1) > 1e-12]

  candidates <- vector("list", length(q0_grid))
  keep <- 0L

  for (q0 in q0_grid) {
    restored <- restore_r_c_d_v(p_wcdv, q0, q1)
    if (is.null(restored)) next

    target <- candidate_target(
      p_rcdv = restored,
      values = exact_law$values,
      mechanism_assumption = mechanism_assumption,
      p_c = exact_law$p_c,
      lower = lower,
      upper = upper
    )

    if (is.null(target)) next

    keep <- keep + 1L
    candidates[[keep]] <- data.frame(
      q0 = q0,
      lower_bound = unname(target["LB"]),
      upper_bound = unname(target["UB"]),
      p_r0 = unname(target["p_r0"])
    )
  }

  if (keep == 0L) {
    stop("No population-compatible proxy relationship was found.")
  }

  candidates <- do.call(rbind, candidates[seq_len(keep)])

  # Verify that the true proxy relationship is compatible with the
  # observed population law. It is not inserted into the grid.
  restored_at_truth <- restore_r_c_d_v(
    p_wcdv,
    q0 = exact_law$true_q0,
    q1 = q1
  )

  target_at_truth <- candidate_target(
    p_rcdv = restored_at_truth,
    values = exact_law$values,
    mechanism_assumption = mechanism_assumption,
    p_c = exact_law$p_c,
    lower = lower,
    upper = upper
  )

  if (is.null(target_at_truth)) {
    stop("The true proxy relationship was incompatible.")
  }

  list(
    LB = min(candidates$lower_bound),
    UB = max(candidates$upper_bound),
    width = max(candidates$upper_bound) -
      min(candidates$lower_bound),
    compatible_q0_range = range(candidates$q0),
    candidates = candidates,
    proxy_relationship = relationship,
    true_q0_is_compatible = TRUE,
    target_at_true_q0 = target_at_truth[c("LB", "UB")]
  )
}


# ------------------------------------------------------------
# 6. Evaluate all methods for one exact population law
# ------------------------------------------------------------

evaluate_population_law <- function(
    exact_law,
    lower = 1,
    upper = 5,
    proxy_grid_n = 20001) {

  mechanism <- exact_law$mechanism

  proxy_assumption <- switch(
    mechanism,
    CRCR = "CRCR",
    CRAR = "CRAR",
    CRNAR = "general"
  )

  support <- population_support_bounds(
    exact_law,
    lower = lower,
    upper = upper
  )

  monotone <- population_monotone_bounds(
    exact_law,
    lower = lower,
    upper = upper
  )

  proxy <- proxy_mean_bounds_from_population_law(
    exact_law = exact_law,
    mechanism_assumption = proxy_assumption,
    lower = lower,
    upper = upper,
    grid_n = proxy_grid_n
  )

  point_estimates <- population_point_estimates(exact_law)
  truth <- exact_law$true_mean

  intervals <- data.frame(
    method = c("Manski support", "Monotonicity", "Proxy"),
    lower_bound = c(
      support["LB"],
      monotone["LB"],
      proxy$LB
    ),
    upper_bound = c(
      support["UB"],
      monotone["UB"],
      proxy$UB
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  intervals$width <-
    intervals$upper_bound - intervals$lower_bound

  intervals$contains_population_truth <- with(
    intervals,
    lower_bound - 1e-8 <= truth &
      truth <= upper_bound + 1e-8
  )

  intervals$true_mean <- truth
  intervals$retain_all_estimate <-
    point_estimates$retain_all_estimate
  intervals$deletion_estimate <-
    point_estimates$deletion_estimate
  intervals$standardized_deletion_estimate <-
    point_estimates$standardized_deletion_estimate
  intervals$retain_all_bias <-
    point_estimates$retain_all_bias
  intervals$deletion_bias <-
    point_estimates$deletion_bias
  intervals$standardized_deletion_bias <-
    point_estimates$standardized_deletion_bias
  intervals$proxy_assumption <- proxy_assumption
  intervals$true_q0_is_compatible <-
    proxy$true_q0_is_compatible

  intervals
}


# ------------------------------------------------------------
# 7. Baseline factorial design
# ------------------------------------------------------------

run_population_analysis <- function(
    mechanisms = c("CRCR", "CRAR", "CRNAR"),
    careful_prevalences = c(0.50, 0.70, 0.90),
    certification_rates = c(0.40, 0.70, 0.90),
    proxy_settings = list(
      weak = c(
        p_w1_careful = 0.45,
        p_w1_careless = 0.55
      ),
      moderate = c(
        p_w1_careful = 0.30,
        p_w1_careless = 0.70
      ),
      strong = c(
        p_w1_careful = 0.10,
        p_w1_careless = 0.90
      )
    ),
    lower = 1,
    upper = 5,
    proxy_grid_n = 20001) {

  design <- expand.grid(
    mechanism = mechanisms,
    careful_prevalence = careful_prevalences,
    certification_rate = certification_rates,
    proxy_setting = names(proxy_settings),
    stringsAsFactors = FALSE
  )

  output <- vector("list", nrow(design))

  for (row in seq_len(nrow(design))) {
    condition <- design[row, ]
    proxy_parameters <- proxy_settings[[condition$proxy_setting]]

    law <- exact_likert_cr_law(
      lower = lower,
      upper = upper,
      careful_prevalence = condition$careful_prevalence,
      careful_certification_rate = condition$certification_rate,
      p_w1_careful = proxy_parameters["p_w1_careful"],
      p_w1_careless = proxy_parameters["p_w1_careless"],
      mechanism = condition$mechanism
    )

    evaluation <- evaluate_population_law(
      exact_law = law,
      lower = lower,
      upper = upper,
      proxy_grid_n = proxy_grid_n
    )

    output[[row]] <- cbind(
      condition[rep(1, nrow(evaluation)), , drop = FALSE],
      evaluation
    )
  }

  results <- do.call(rbind, output)

  if (!all(results$contains_population_truth)) {
    stop("At least one identified set excluded the population truth.")
  }

  results
}


# ------------------------------------------------------------
# 8. Differential-certification analysis
# ------------------------------------------------------------

run_differential_certification_analysis <- function(
    mechanisms = c("CRCR", "CRAR", "CRNAR"),
    careful_prevalence = 0.70,
    proxy_setting = c(
      p_w1_careful = 0.30,
      p_w1_careless = 0.70
    ),
    lower = 1,
    upper = 5,
    proxy_grid_n = 20001) {

  scenarios <- data.frame(
    scenario = c(
      "Lower-bound binding",
      "Upper-bound binding"
    ),
    certification_rate_z0 = c(0.90, 0.30),
    certification_rate_z1 = c(0.30, 0.90),
    stringsAsFactors = FALSE
  )

  design <- merge(
    scenarios,
    data.frame(
      mechanism = mechanisms,
      stringsAsFactors = FALSE
    ),
    by = NULL
  )

  output <- vector("list", nrow(design))

  for (row in seq_len(nrow(design))) {
    condition <- design[row, ]

    law <- exact_likert_cr_law(
      lower = lower,
      upper = upper,
      careful_prevalence = careful_prevalence,
      certification_rate_z0 = condition$certification_rate_z0,
      certification_rate_z1 = condition$certification_rate_z1,
      p_w1_careful = proxy_setting["p_w1_careful"],
      p_w1_careless = proxy_setting["p_w1_careless"],
      mechanism = condition$mechanism
    )

    evaluation <- evaluate_population_law(
      exact_law = law,
      lower = lower,
      upper = upper,
      proxy_grid_n = proxy_grid_n
    )

    output[[row]] <- cbind(
      condition[rep(1, nrow(evaluation)), , drop = FALSE],
      careful_prevalence = careful_prevalence,
      proxy_setting = "moderate",
      evaluation
    )
  }

  do.call(rbind, output)
}


# ------------------------------------------------------------
# 9. Summaries for the paper
# ------------------------------------------------------------

condition_level_results <- function(results) {
  requested <- c(
    "mechanism",
    "careful_prevalence",
    "certification_rate",
    "proxy_setting",
    "method",
    "true_mean",
    "lower_bound",
    "upper_bound",
    "width",
    "contains_population_truth",
    "retain_all_estimate",
    "deletion_estimate",
    "standardized_deletion_estimate",
    "retain_all_bias",
    "deletion_bias",
    "standardized_deletion_bias"
  )

  output <- results[, requested]
  output[order(
    output$mechanism,
    output$careful_prevalence,
    output$certification_rate,
    output$proxy_setting,
    output$method
  ), ]
}


summarize_bound_widths <- function(results) {
  groups <- split(
    results,
    interaction(
      results$mechanism,
      results$method,
      results$proxy_setting,
      drop = TRUE
    )
  )

  output <- lapply(groups, function(group) {
    data.frame(
      mechanism = group$mechanism[1],
      method = group$method[1],
      proxy_setting = group$proxy_setting[1],
      mean_width = mean(group$width),
      min_width = min(group$width),
      max_width = max(group$width),
      truth_inclusion = mean(group$contains_population_truth)
    )
  })

  output <- do.call(rbind, output)
  rownames(output) <- NULL
  output[order(
    output$mechanism,
    output$method,
    output$proxy_setting
  ), ]
}


condition_level_point_estimates <- function(results) {
  # Point estimates do not depend on the bounding method or proxy strength.
  # Select one row for each mechanism/prevalence/certification condition rather
  # than relying on exact floating-point equality in unique().
  keys <- results[, c(
    "mechanism",
    "careful_prevalence",
    "certification_rate"
  )]

  output <- results[!duplicated(keys), c(
    "mechanism",
    "careful_prevalence",
    "certification_rate",
    "true_mean",
    "retain_all_estimate",
    "deletion_estimate",
    "standardized_deletion_estimate",
    "retain_all_bias",
    "deletion_bias",
    "standardized_deletion_bias"
  )]

  rownames(output) <- NULL
  output[order(
    output$mechanism,
    output$careful_prevalence,
    output$certification_rate
  ), ]
}


summarize_point_estimates <- function(results) {
  unique_rows <- condition_level_point_estimates(results)

  output <- aggregate(
    cbind(
      retain_all_bias,
      deletion_bias,
      standardized_deletion_bias
    ) ~ mechanism,
    data = unique_rows,
    FUN = function(x) mean(abs(x))
  )

  names(output)[names(output) == "retain_all_bias"] <-
    "mean_absolute_retain_all_bias"
  names(output)[names(output) == "deletion_bias"] <-
    "mean_absolute_deletion_bias"
  names(output)[names(output) == "standardized_deletion_bias"] <-
    "mean_absolute_standardized_deletion_bias"

  output
}


# ------------------------------------------------------------
# 10. Proxy-grid convergence check
# ------------------------------------------------------------

check_proxy_grid_convergence <- function(
    grid_sizes = c(1001, 5001, 20001),
    careful_prevalence = 0.70,
    certification_rate = 0.70,
    proxy_setting = c(
      p_w1_careful = 0.30,
      p_w1_careless = 0.70
    )) {

  output <- list()
  index <- 0L

  for (grid_n in grid_sizes) {
    for (mechanism in c("CRCR", "CRAR", "CRNAR")) {
      law <- exact_likert_cr_law(
        careful_prevalence = careful_prevalence,
        careful_certification_rate = certification_rate,
        p_w1_careful = proxy_setting["p_w1_careful"],
        p_w1_careless = proxy_setting["p_w1_careless"],
        mechanism = mechanism
      )

      evaluation <- evaluate_population_law(
        law,
        proxy_grid_n = grid_n
      )

      proxy_row <- evaluation[evaluation$method == "Proxy", ]

      index <- index + 1L
      output[[index]] <- data.frame(
        mechanism = mechanism,
        grid_n = grid_n,
        lower_bound = proxy_row$lower_bound,
        upper_bound = proxy_row$upper_bound,
        width = proxy_row$width
      )
    }
  }

  do.call(rbind, output)
}


# ------------------------------------------------------------
# 11. Example commands
# ------------------------------------------------------------

# Quick validation run:
# quick_results <- run_population_analysis(
#   careful_prevalences = 0.70,
#   certification_rates = 0.70,
#   proxy_settings = list(
#     moderate = c(
#       p_w1_careful = 0.30,
#       p_w1_careless = 0.70
#     )
#   ),
#   proxy_grid_n = 1001
# )
# condition_level_results(quick_results)

# Full baseline analysis for the paper:
# population_results <- run_population_analysis(
#   proxy_grid_n = 20001
# )
# paper_condition_results <- condition_level_results(population_results)
# paper_width_summary <- summarize_bound_widths(population_results)
# paper_point_summary <- summarize_point_estimates(population_results)

# Differential-certification analysis:
# differential_results <- run_differential_certification_analysis(
#   proxy_grid_n = 20001
# )

# Run the full analysis
# ------------------------------------------------------------

# Grid-convergence check
grid_check <- check_proxy_grid_convergence()

population_results <- run_population_analysis(
  proxy_grid_n = 20001
)

paper_condition_results <- condition_level_results(
  population_results
)

paper_width_summary <- summarize_bound_widths(
  population_results
)

paper_point_summary <- summarize_point_estimates(
  population_results
)

differential_results <- run_differential_certification_analysis(
  proxy_grid_n = 20001
)

# Inspect condition-level bounds
# ------------------------------------------------------------

View(paper_condition_results)
print(paper_condition_results, row.names = FALSE)

# Inspect mean bound widths by method
View(paper_width_summary)
print(paper_width_summary, row.names = FALSE)

# Inspect bias in the retain-all and deletion estimators
View(paper_point_summary)
print(paper_point_summary, row.names = FALSE)

# Summarize width and truth inclusion
aggregate(
  cbind(width, contains_population_truth) ~ mechanism + method,
  data = population_results,
  FUN = mean
)

aggregate(
  cbind(width, contains_population_truth) ~
    mechanism + proxy_setting + method,
  data = population_results,
  FUN = mean
)

# Construct condition-level point-estimate results
# ------------------------------------------------------------

point_results <- unique(
  population_results[
    ,
    c(
      "mechanism",
      "careful_prevalence",
      "certification_rate",
      "retain_all_estimate",
      "deletion_estimate",
      "true_mean",
      "retain_all_bias",
      "deletion_bias"
    )
  ]
)

View(point_results)
print(point_results, row.names = FALSE)

# Construct the differential-certification table
# ------------------------------------------------------------

differential_table <- differential_results[
  ,
  c(
    "scenario",
    "mechanism",
    "method",
    "lower_bound",
    "upper_bound",
    "width",
    "true_mean",
    "contains_population_truth"
  )
]
point_results$abs_retain_all_bias <-
  abs(point_results$retain_all_bias)

point_results$abs_deletion_bias <-
  abs(point_results$deletion_bias)

aggregate(
  cbind(abs_retain_all_bias, abs_deletion_bias) ~ mechanism,
  data = point_results,
  FUN = mean
)

differential_table <- differential_table[
  order(
    differential_table$scenario,
    differential_table$mechanism,
    differential_table$method
  ),
]

print(differential_table, row.names = FALSE)
View(differential_table)

# Recheck proxy-grid convergence at selected grid sizes
# ------------------------------------------------------------

grid_check <- check_proxy_grid_convergence(
  grid_sizes = c(1001, 5001, 20001)
)

print(grid_check, row.names = FALSE)
View(grid_check)
