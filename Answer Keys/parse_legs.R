# Answer Keys/parse_legs.R
# Wagerzon prop-description parser and generic leg-based joint probability pricer.

#' Convert a Wagerzon fraction string to a numeric.
#'
#' Handles three formats Wagerzon uses interchangeably:
#'   "3"   → 3
#'   "2.5" → 2.5
#'   "2½"  → 2.5 (using the U+00BD unicode fraction)
#' Also handles ¼ (U+00BC → 0.25) and ¾ (U+00BE → 0.75) for future-proofing.
#'
#' Returns NA_real_ on unparseable input.
parse_fraction <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(NA_real_)
  frac_map <- c("¼" = 0.25, "½" = 0.5, "¾" = 0.75)
  # Try pure numeric first (handles "3", "2.5")
  n <- suppressWarnings(as.numeric(s))
  if (!is.na(n)) return(n)
  # Check for unicode fraction suffix
  last_char <- substr(s, nchar(s), nchar(s))
  if (last_char %in% names(frac_map)) {
    int_part <- suppressWarnings(as.numeric(substr(s, 1L, nchar(s) - 1L)))
    if (is.na(int_part)) int_part <- 0
    return(int_part + frac_map[[last_char]])
  }
  NA_real_
}
