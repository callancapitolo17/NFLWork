# Answer Keys/books_strip.R
#
# Pure HTML renderer for the parlay-tab "books strip" — a single inline pill row
# showing the model + four book devigged fair probabilities + blended consensus.
# Used by `MLB Dashboard/mlb_dashboard.R` as a reactable cell renderer.
#
# All inputs are probabilities on [0, 1] or NA. Output is a character string.

render_books_strip <- function(model, dk, fd, px, nv, cons) {
  fmt_pill <- function(label, prob, css_class) {
    if (is.na(prob)) {
      sprintf('<span class="pill %s dim">%s &mdash;</span>', css_class, label)
    } else {
      sprintf('<span class="pill %s">%s %.1f</span>', css_class, label, prob * 100)
    }
  }

  paste0(
    '<span class="books-strip">',
    fmt_pill("M",    model, "model"),
    fmt_pill("DK",   dk,    "book"),
    fmt_pill("FD",   fd,    "book"),
    fmt_pill("PX",   px,    "book"),
    fmt_pill("NV",   nv,    "book"),
    fmt_pill("Cons", cons,  "cons"),
    '</span>'
  )
}
