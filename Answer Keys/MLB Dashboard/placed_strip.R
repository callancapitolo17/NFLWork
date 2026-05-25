# Answer Keys/MLB Dashboard/placed_strip.R
#
# Pure helper that renders the "placed" hero strip for a single bet card
# after one or more Wagerzon placements. One <span class="placement-chip">
# per placement plus an optional dashed "+ another" button.
#
# Inputs:
#   chips           list of list(account=<chr>, risk=<num>, ticket=<chr>)
#   all_wz_accounts character vector of every WZ account label currently
#                   discovered by wagerzon_accounts.list_accounts()
#   bet_hash        character — the bet hash; emitted as data-bet-hash on
#                   the strip and on each chip so JS handlers can route
#                   per-bet (also doubles as the placeBet() routing key).
#   book            character — the bookmaker key of the pick. Retained for
#                   signature compatibility; NOT used to gate "+ another"
#                   (the pick book drifts as odds change). The "+ another"
#                   gate is account-based: it renders when any chip's account
#                   is in all_wz_accounts (i.e. the placements are on WZ).
#
# Returns a single HTML string (the entire <div class="hero-placed"> block).

.acct_short <- function(label) {
  # "Wagerzon" -> "WZ"; "WagerzonJ" -> "WZJ"; "WagerzonC" -> "WZC".
  # Falls back to label as-is for any non-Wagerzon book account.
  if (grepl("^Wagerzon", label)) {
    suffix <- sub("^Wagerzon", "", label)
    paste0("WZ", suffix)
  } else {
    label
  }
}

.chip_html <- function(chip, bet_hash) {
  sprintf(
    '<span class="placement-chip" data-account="%s" data-risk="%s" data-ticket="%s" data-bet-hash="%s"><span class="acct">%s</span>$%s <span class="ticket">#%s</span></span>',
    htmltools::htmlEscape(chip$account),
    as.character(round(as.numeric(chip$risk))),
    htmltools::htmlEscape(chip$ticket),
    htmltools::htmlEscape(bet_hash),
    htmltools::htmlEscape(.acct_short(chip$account)),
    as.character(round(as.numeric(chip$risk))),
    htmltools::htmlEscape(chip$ticket)
  )
}

render_placed_strip <- function(chips, all_wz_accounts, bet_hash, book) {
  stopifnot(is.list(chips), is.character(all_wz_accounts),
            is.character(bet_hash), is.character(book))

  chip_html <- paste(vapply(chips, .chip_html, character(1),
                            bet_hash = bet_hash),
                     collapse = "")

  # Gate "+ another" on whether any chip is on a WZ account, not on the
  # card's current pick book (book is retained for signature compatibility
  # but is no longer used for the WZ gate; pick-book drifts as odds change).
  placed_accounts <- vapply(chips, function(c) c$account, character(1))
  is_wz <- any(placed_accounts %in% all_wz_accounts)
  untouched <- setdiff(all_wz_accounts, placed_accounts)

  add_another_html <- ""
  if (is_wz) {
    disabled_attr <- if (length(untouched) == 0L) " disabled" else ""
    disabled_title <- if (length(untouched) == 0L) {
      ' title="All WZ accounts placed."'
    } else {
      ' title="Re-open the card to place on another WZ account"'
    }
    add_another_html <- sprintf(
      '<button type="button" class="add-another" data-bet-hash="%s" onclick="addAnother(this)"%s%s>+ another</button>',
      htmltools::htmlEscape(bet_hash), disabled_attr, disabled_title
    )
  }

  sprintf(
    '<div class="hero-placed" data-bet-hash="%s"><span class="placed-label">placed</span>%s%s</div>',
    htmltools::htmlEscape(bet_hash), chip_html, add_another_html
  )
}
