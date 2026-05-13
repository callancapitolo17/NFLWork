# Answer Keys/MLB Dashboard/book_pill.R
# Render one book pill for the bets-tab odds-screen card.
#
# Three visual states (driven by precomputed is_exact_line flag):
#   - exact-line match    -> standard pill, just book + price
#   - mismatched          -> amber-tinted pill + amber line tag showing line_quoted
#   - no quote (NA odds)  -> muted dashed pill with em-dash
# Pick override (is_pick = TRUE) applies the green pick class on top of
# whichever state the row is in.
#
# Companion file: Answer Keys/books_strip.R (parlays-tab pill row)

#' Format a line value for the line tag (e.g., 5.0 -> "5", 5.5 -> "5.5").
#' When `signed = TRUE` a leading "+" is added for positive values, matching
#' the convention used by spread display ("+1.5" / "-1.5").
.format_line_value <- function(x, signed = FALSE) {
  if (is.na(x)) return("")
  base <- if (x == round(x)) format(as.integer(x))
          else sub("\\.?0+$", "", sprintf("%.2f", x))
  if (signed && x > 0) paste0("+", base) else base
}

#' Render one book pill.
#'
#' @param book Display label for the book (e.g., "WZ", "H88").
#' @param american_odds Integer odds (e.g., 125, -110). NA -> no-quote pill.
#' @param line_quoted The line this book is actually showing on this side.
#' @param is_exact_line Precomputed BOOLEAN: TRUE when book's line matches the
#'   model line exactly; FALSE triggers the amber mismatched state. Ignored
#'   when american_odds is NA.
#' @param is_pick TRUE if this is the pick book on the pick side; applies green.
#' @param side Either "over" (default for totals over / favorite spread) or
#'   "under" (totals under / dog spread). Drives O/U prefix on mismatched
#'   line tag (only when `is_totals = TRUE`).
#' @param is_totals TRUE when the bet is a totals market (O/U prefix added on
#'   mismatched line tag); FALSE for spread bets (signed line value shown
#'   instead, e.g. "-1.5"). Defaults to TRUE for backwards compatibility with
#'   the original totals-only call sites.
#' @return HTML string for the pill.
render_book_pill <- function(book, american_odds, line_quoted, is_exact_line,
                             is_pick = FALSE, side = "over", is_totals = TRUE) {
  # State 1: no quote
  if (is.na(american_odds)) {
    return(sprintf('<span class="pill muted"><span class="book">%s</span>&mdash;</span>',
                   htmltools::htmlEscape(book)))
  }

  price_str <- if (american_odds > 0) paste0("+", american_odds) else as.character(american_odds)

  is_mismatched <- !isTRUE(is_exact_line)

  classes <- "pill"
  if (is_pick) classes <- paste(classes, "pick")
  if (is_mismatched) classes <- paste(classes, "mismatched")

  tag_html <- ""
  if (is_mismatched) {
    if (is_totals) {
      prefix <- if (side == "under") "U" else "O"
      tag_html <- sprintf('<span class="line-tag">%s%s</span>',
                          prefix, .format_line_value(line_quoted))
    } else {
      # Spread: display the signed line value directly (e.g. "-1.5", "+0.5").
      tag_html <- sprintf('<span class="line-tag">%s</span>',
                          .format_line_value(line_quoted, signed = TRUE))
    }
  }

  sprintf('<span class="%s"><span class="book">%s</span>%s%s</span>',
          classes,
          htmltools::htmlEscape(book),
          tag_html,
          price_str)
}
