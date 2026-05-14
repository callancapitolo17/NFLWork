# Answer Keys/MLB Dashboard/book_cell.R
# Render one grid cell for the bets-tab V8 card.
#
# Four visual states:
#   - empty:   book has no quote on this side -> dim em-dash
#   - exact:   book's line matches the bet's model line -> plain price
#   - alt:     book's closest line differs from model line -> amber background +
#              line tag above the price ("-1.5" for spreads, "O7.5" / "U5.5"
#              for totals)
#   - pick:    overrides exact/alt -> green border + green price (applied when
#              is_pick = TRUE; pick is exact by construction in expand_bets_to_book_prices)
#
# The book label (WZ, H88, etc.) is NOT rendered here — it sits in a column
# header above the grid. This keeps the cell compact.

#' Format a line value for the alt-line tag.
#' (Carried over from the old book_pill.R unchanged.)
.format_line_value <- function(x, signed = FALSE) {
  if (is.na(x)) return("")
  base <- if (x == round(x)) format(as.integer(x))
          else sub("\\.?0+$", "", sprintf("%.2f", x))
  if (signed && x > 0) paste0("+", base) else base
}

#' 2-way probit devig: given two American odds, return the no-vig
#' (fair) American odds for both sides.
#'
#' Mirrors `Tools.R::devig_american()` for the 2-way case using the
#' closed-form solution `c* = -(z1 + z2) / 2`. Inlined here to avoid
#' sourcing all of Tools.R into the dashboard. A parity test in
#' `Answer Keys/tests/test_devig_pair_matches_tools.R` guards against
#' drift across a table of representative inputs.
#'
#' @param odd1 American odds for side 1 (integer or numeric).
#' @param odd2 American odds for side 2.
#' @return list(fair1 = numeric, fair2 = numeric). Returns NA pair when
#'   either input is NA or zero.
.devig_american_pair <- function(odd1, odd2) {
  if (is.na(odd1) || is.na(odd2) || odd1 == 0 || odd2 == 0) {
    return(list(fair1 = NA_real_, fair2 = NA_real_))
  }
  # Implied probabilities from American odds
  p1 <- if (odd1 > 0) 100 / (odd1 + 100) else -odd1 / (-odd1 + 100)
  p2 <- if (odd2 > 0) 100 / (odd2 + 100) else -odd2 / (-odd2 + 100)
  # Probit z-shift: z' = z + c, with c chosen so p1' + p2' = 1
  eps <- 1e-9
  z1 <- qnorm(min(max(p1, eps), 1 - eps))
  z2 <- qnorm(min(max(p2, eps), 1 - eps))
  c_star <- -(z1 + z2) / 2
  q1 <- pnorm(z1 + c_star)
  q2 <- pnorm(z2 + c_star)
  # Convert devigged probabilities back to American odds
  to_amer <- function(p) {
    if (p >= 0.5) round(-100 * p / (1 - p))
    else          round( 100 * (1 - p) / p)
  }
  list(fair1 = to_amer(q1), fair2 = to_amer(q2))
}

#' Render one bets-tab grid cell.
#'
#' @param american_odds Integer odds (e.g., 125, -110). NA -> empty state.
#' @param opposite_american_odds Integer odds for the OTHER side at the same
#'   book (e.g., the Under price when this cell is the Over). Used to compute
#'   probit-devigged fair odds for the FAIR view of the toggle. NA -> no
#'   fair span emitted (cell shows raw only, behaves like legacy).
#' @param line_quoted Numeric line the book is showing on this side.
#' @param is_exact_line Boolean: TRUE when book's line matches the model line
#'   exactly; FALSE -> alt state.
#' @param is_pick TRUE if this is the picked book on the pick side; overrides
#'   exact/alt with the green pick state.
#' @param side_word "over" or "under" — only used for the O/U prefix on a
#'   mismatched totals line tag.
#' @param is_totals TRUE for totals markets (line tag gets O/U prefix);
#'   FALSE for spreads (signed line value, e.g. "-1.5").
#' @return HTML string for the cell. Contains both <span class="raw"> and
#'   (when devig is computable) <span class="fair">; CSS on the parent
#'   .price-grid container determines which is visible (toggle).
render_book_cell <- function(american_odds, line_quoted, is_exact_line,
                              is_pick = FALSE, side_word = "over",
                              is_totals = TRUE,
                              opposite_american_odds = NA_integer_) {
  # State 1: empty (no quote)
  if (is.na(american_odds)) {
    return('<div class="cell empty"><span class="raw">&mdash;</span><span class="fair">&mdash;</span></div>')
  }

  raw_str <- if (american_odds > 0) paste0("+", american_odds) else as.character(american_odds)

  # Compute devigged American for the FAIR span (if we have both sides).
  # Fallback to em-dash so every non-empty cell still emits <span class="fair">;
  # otherwise the Task 5 CSS toggle (.show-fair hides .raw) would render this
  # cell blank in FAIR view when the book quotes only one side.
  fair_html <- '<span class="fair">&mdash;</span>'  # fallback when no devig available
  if (!is.na(opposite_american_odds)) {
    pair <- .devig_american_pair(american_odds, opposite_american_odds)
    if (!is.na(pair$fair1)) {
      fair_str <- if (pair$fair1 > 0) paste0("+", as.integer(pair$fair1))
                  else as.character(as.integer(pair$fair1))
      fair_html <- sprintf('<span class="fair">%s</span>', fair_str)
    }
  }

  is_mismatched <- !isTRUE(is_exact_line)

  cell_class <- if (is_pick) "cell pick"
                else if (is_mismatched) "cell alt"
                else "cell exact"

  tag_html <- ""
  if (is_mismatched && !is_pick) {
    if (is_totals) {
      prefix <- if (side_word == "under") "U" else "O"
      tag_html <- sprintf('<span class="alt-line">%s%s</span>',
                          prefix, .format_line_value(line_quoted))
    } else {
      tag_html <- sprintf('<span class="alt-line">%s</span>',
                          .format_line_value(line_quoted, signed = TRUE))
    }
  }

  sprintf('<div class="%s">%s<span class="raw">%s</span>%s</div>',
          cell_class, tag_html, raw_str, fair_html)
}

# Backwards-compat shim: old code may still source book_pill.R via legacy paths.
# Provide render_book_pill as a thin wrapper so a stale source() doesn't crash
# during the transition. Remove after the V8 wiring lands and is stable.
render_book_pill <- function(book, american_odds, line_quoted, is_exact_line,
                              is_pick = FALSE, side = "over", is_totals = TRUE) {
  warning("render_book_pill is deprecated; use render_book_cell. The book ",
          "label is no longer part of the cell — it lives in the column header.",
          call. = FALSE)
  render_book_cell(american_odds, line_quoted, is_exact_line, is_pick,
                   side_word = side, is_totals = is_totals)
}
