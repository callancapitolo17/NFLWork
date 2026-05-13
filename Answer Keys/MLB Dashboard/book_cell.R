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

#' Render one bets-tab grid cell.
#'
#' @param american_odds Integer odds (e.g., 125, -110). NA -> empty state.
#' @param line_quoted Numeric line the book is showing on this side.
#' @param is_exact_line Boolean: TRUE when book's line matches the model line
#'   exactly; FALSE -> alt state.
#' @param is_pick TRUE if this is the picked book on the pick side; overrides
#'   exact/alt with the green pick state.
#' @param side_word "over" or "under" — only used for the O/U prefix on a
#'   mismatched totals line tag.
#' @param is_totals TRUE for totals markets (line tag gets O/U prefix);
#'   FALSE for spreads (signed line value, e.g. "-1.5").
#' @return HTML string for the cell (a single <div class="cell ..."> ... </div>).
render_book_cell <- function(american_odds, line_quoted, is_exact_line,
                              is_pick = FALSE, side_word = "over",
                              is_totals = TRUE) {
  # State 1: empty (no quote)
  if (is.na(american_odds)) {
    return('<div class="cell empty"><span class="price">&mdash;</span></div>')
  }

  price_str <- if (american_odds > 0) paste0("+", american_odds) else as.character(american_odds)

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

  sprintf('<div class="%s">%s<span class="price">%s</span></div>',
          cell_class, tag_html, price_str)
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
