# Work14 boundary diagnostics. This file does not implement a bootstrap or p-value.
# Official fwildclusterboot p-values are retained unchanged. Arithmetic on returned
# statistics is diagnostic only. The production data and models are out of scope.

work14_alpha <- 1 / 60

expected_rademacher_draws <- function(G, requested_B) {
  stopifnot(length(G) == 1L, G >= 1, G == floor(G),
            length(requested_B) == 1L, requested_B >= 1,
            requested_B == floor(requested_B))
  full <- 2^G <= requested_B
  list(full_enumeration = full,
       expected_draws = if (full) 2^G else requested_B,
       expected_columns_with_observed = if (full) 2^G + 1 else requested_B + 1)
}

cp99_tail_interval <- function(k, M) {
  stopifnot(length(k) == 1L, length(M) == 1L, is.finite(k), is.finite(M),
            M > 0, k >= 0, k <= M, k == floor(k), M == floor(M))
  unname(stats::binom.test(k, M, alternative = 'two.sided',
                         conf.level = 0.99)$conf.int)
}

interval_touches_alpha <- function(interval, alpha = work14_alpha) {
  stopifnot(length(interval) == 2L, all(is.finite(interval)))
  interval[1] <= alpha && alpha <= interval[2]
}

diagnose_boundary <- function(ans, id, fixture_kind, G, requested_B, direction,
                              delta_pp = NA_real_, observed_se = NULL,
                              official_division_contract = FALSE,
                              mc_independent_signs_confirmed = FALSE,
                              explicitly_degenerate = FALSE) {
  stopifnot(direction %in% c('>', '<'))
  expected <- expected_rademacher_draws(G, requested_B)
  t <- as.vector(ans$t_boot)
  observed <- as.numeric(ans$t_stat)
  p <- as.numeric(ans$p_val)
  M <- length(t)
  boot_iter <- if (is.null(ans$boot_iter)) NA_real_ else as.numeric(ans$boot_iter)
  observed_finite <- length(observed) == 1L && is.finite(observed)
  p_valid <- length(p) == 1L && is.finite(p) && p >= 0 && p <= 1
  raw_pass <- p_valid && p < work14_alpha
  returned_NA <- sum(is.na(t) & !is.nan(t))
  returned_NaN <- sum(is.nan(t))
  returned_posInf <- sum(t == Inf, na.rm = TRUE)
  returned_negInf <- sum(t == -Inf, na.rm = TRUE)
  returned_nonfinite <- sum(!is.finite(t))
  all_t_finite <- M > 0L && all(is.finite(t))
  exact_constant <- all_t_finite && length(unique(t)) == 1L
  valid_comparisons <- observed_finite && returned_NA == 0L && returned_NaN == 0L
  k <- if (valid_comparisons) {
    if (direction == '>') sum(observed < t) else sum(observed > t)
  } else NA_real_
  ties <- if (valid_comparisons) sum(observed == t) else NA_real_
  q <- if (M > 0L && is.finite(ties)) ties / M else NA_real_
  tie_hold <- raw_pass && is.finite(q) && p + q >= work14_alpha
  M_matches <- M == expected$expected_draws
  boot_iter_matches <- is.finite(boot_iter) && boot_iter == expected$expected_draws

  # A numerical SE is never reconstructed from an estimate/t ratio. The public
  # object does not expose SE. For the pinned official implementation only,
  # finite t = numer/sqrt(JJ_sum) excludes a zero denominator under R arithmetic.
  se_value <- if (is.null(observed_se)) NA_real_ else as.numeric(observed_se)
  se_status <- if (!is.null(observed_se) && length(se_value) == 1L) {
    if (!is.finite(se_value)) 'NONFINITE_MEASURED' else if (se_value == 0) 'ZERO_MEASURED' else if (se_value < 0) 'NEGATIVE_MEASURED' else 'POSITIVE_MEASURED'
  } else if (official_division_contract && observed_finite) {
    'ZERO_EXCLUDED_BY_OFFICIAL_DIVISION_AND_FINITE_T'
  } else 'UNAVAILABLE'

  reasons <- character()
  add_reason <- function(x) reasons <<- c(reasons, x)
  if (!observed_finite) add_reason('NONFINITE_OBSERVED_T')
  if (!p_valid) add_reason('INVALID_OFFICIAL_P')
  if (se_status %in% c('ZERO_MEASURED', 'NONFINITE_MEASURED', 'NEGATIVE_MEASURED')) add_reason(paste0('SE_', se_status))
  if (M == 0L) add_reason('EMPTY_BOOTSTRAP_STATISTICS')
  if (!M_matches) add_reason('UNEXPECTED_M')
  if (!boot_iter_matches) add_reason('UNEXPECTED_RETURNED_BOOT_ITER')
  if (returned_nonfinite > 0L) add_reason('NONFINITE_BOOTSTRAP_T')
  if (exact_constant || explicitly_degenerate) add_reason('DEGENERATE_STATISTICS')
  if (tie_hold) add_reason('TIE_SENSITIVITY_AT_ALPHA')

  cp <- c(NA_real_, NA_real_)
  cp_touches <- NA
  cp_status <- 'NOT_APPLICABLE_FULL_ENUMERATION'
  if (!expected$full_enumeration) {
    cp_status <- 'NOT_COMPUTED_INVALID_OR_INCOMPLETE_COUNTS'
    if (!mc_independent_signs_confirmed) {
      cp_status <- 'NOT_COMPUTED_INDEPENDENT_SIGN_DESIGN_UNVERIFIED'
      add_reason('MC_DRAW_DESIGN_UNVERIFIED')
    } else if (observed_finite && all_t_finite && M_matches && boot_iter_matches &&
               is.finite(k) && p_valid && !exact_constant && !explicitly_degenerate) {
      cp <- cp99_tail_interval(k, M)
      cp_touches <- interval_touches_alpha(cp)
      cp_status <- if (cp_touches) 'ALPHA_CONTACT_OR_CROSSING_HOLD' else 'CLEAR_OF_ALPHA'
      if (cp_touches) add_reason('MC99_INTERVAL_TOUCHES_ALPHA')
    }
  }

  # These fields refer to statistics BEFORE the package's NA/NaN deletion.
  # They are not publicly measured. Full expected M logically implies none were
  # deleted, but a deficit is retained as a deficit, not assigned a guessed cause.
  removed_combined_source_inference <- if (official_division_contract && M_matches) 0L else NA_integer_
  predelete_status <- if (!is.na(removed_combined_source_inference))
    'NO_NA_NAN_DELETION_DEDUCED_FROM_EXPECTED_M' else 'UNAVAILABLE_CAUSE_BREAKDOWN'
  hold <- length(reasons) > 0L
  status <- if (hold) 'HOLD' else if (raw_pass) 'LOCAL_DIAGNOSTICS_CLEAR_RAW_PASS' else 'LOCAL_DIAGNOSTICS_CLEAR_RAW_NONPASS'

  data.frame(
    id = id, fixture_kind = fixture_kind, G = G, requested_B = requested_B,
    full_enumeration = expected$full_enumeration,
    expected_draws = expected$expected_draws,
    expected_columns_with_observed = expected$expected_columns_with_observed,
    returned_boot_iter = boot_iter, actual_M = M,
    missing_vs_expected = max(0, expected$expected_draws - M),
    direction = direction, delta_pp = delta_pp,
    r = if (is.finite(delta_pp)) if (direction == '>') -delta_pp / 100 else delta_pp / 100 else NA_real_,
    alpha = work14_alpha, official_p_unchanged = p,
    official_p_strictly_below_alpha = raw_pass,
    official_p_equals_alpha = p_valid && p == work14_alpha,
    observed_t = if (length(observed) == 1L) observed else NA_real_,
    observed_se_value = se_value, observed_se_status = se_status,
    tail_count = k, tie_count = ties, tie_rate = q,
    diagnostic_p_plus_tie_rate = if (p_valid && is.finite(q)) p + q else NA_real_,
    diagnostic_p_plus_tie_rate_is_corrected_p = FALSE,
    returned_NA = returned_NA, returned_NaN = returned_NaN,
    returned_posInf = returned_posInf, returned_negInf = returned_negInf,
    returned_nonfinite = returned_nonfinite,
    predelete_NA_measured = NA_integer_, predelete_NaN_measured = NA_integer_,
    predelete_posInf_measured = NA_integer_, predelete_negInf_measured = NA_integer_,
    removed_NA_NaN_combined_source_inference = removed_combined_source_inference,
    predelete_status = predelete_status,
    exact_constant_boot_t = exact_constant,
    cp99_lower = cp[1], cp99_upper = cp[2], cp99_touches_alpha = cp_touches,
    cp99_status = cp_status, cp99_is_p_replacement = FALSE,
    mc_independent_sign_design_confirmed = mc_independent_signs_confirmed,
    diagnostic_hold = hold, hold_reasons = paste(reasons, collapse = ';'),
    local_status = status, production_eligibility_validated = FALSE,
    stringsAsFactors = FALSE
  )
}
