# Answer Keys/books_strip.R
#
# Pure HTML renderer for the parlay-tab "books strip" — a single inline pill row
# showing the model + four book devigged fair probabilities + blended consensus.
# Used by `MLB Dashboard/mlb_dashboard.R` as a reactable cell renderer.
#
# All inputs are probabilities on [0, 1] or NA. Output is a character string.

render_books_strip <- function(model, dk, fd, px, nv, cons,
                               bmg = NA_real_, czr = NA_real_,
                               k_agree = NULL, n_agree = NULL) {
  fmt_pill <- function(label, prob, css_class) {
    if (is.na(prob)) {
      sprintf('<span class="pill %s dim">%s &mdash;</span>', css_class, label)
    } else {
      sprintf('<span class="pill %s">%s %.1f</span>', css_class, label, prob * 100)
    }
  }

  # Optional Agree k/n suffix pill — see compute_k_within() below.
  # Omit entirely when n_agree < 2 (no meaningful consensus possible) or
  # when k_agree is NA. The caller is responsible for computing k/n via
  # compute_k_within() over whichever voices are in scope.
  agree_pill <- if (!is.null(k_agree) && !is.null(n_agree) &&
                    !is.na(k_agree) && n_agree >= 2L) {
    sprintf('<span class="pill agree">Agree %d/%d</span>',
            as.integer(k_agree), as.integer(n_agree))
  } else {
    ""
  }

  paste0(
    '<span class="books-strip">',
    fmt_pill("M",    model, "model"),
    fmt_pill("DK",   dk,    "book"),
    fmt_pill("FD",   fd,    "book"),
    fmt_pill("PX",   px,    "book"),
    fmt_pill("NV",   nv,    "book"),
    fmt_pill("BMG",  bmg,   "book"),
    fmt_pill("CZR",  czr,   "book"),
    fmt_pill("Cons", cons,  "cons"),
    agree_pill,
    '</span>'
  )
}

# compute_k_within(): counts how many entries in `probs` fall within +/- `band`
# of the median of the non-NA entries. NAs are excluded from both the median
# and the count. Returns list(k, n) where n is the count of non-NA entries.
#
# When fewer than 2 non-NA entries exist, k is NA (no meaningful consensus
# check is possible with a single voice). The caller decides whether to
# render anything in that case.
#
# `probs` is on the [0, 1] scale (probabilities, not percentages). `band`
# is on the same scale, so 0.01 == 1 percentage point.
compute_k_within <- function(probs, band = 0.01) {
  probs <- probs[!is.na(probs)]
  n <- length(probs)
  if (n < 2) return(list(k = NA_integer_, n = n))
  med <- median(probs)
  list(k = sum(abs(probs - med) <= band), n = n)
}
