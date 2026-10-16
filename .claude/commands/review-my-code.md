# Pre-PR review

Run the `pre-pr-review` skill on $ARGUMENTS (default: the current branch against `main`). It runs the checks, tiers the change, reviews the diff against `REVIEW.md` with 4 agents plus a verifier on the other model, and reports Important / Nit / Pre-existing / Skipped, count first.
