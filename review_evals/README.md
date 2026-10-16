# review_evals

Real bugs that reached `main`, kept as a test set for the code reviewer. Each case is the diff that introduced the bug and the finding a reviewer should have made from it. The reviewer is scored by how many it catches; a change to `REVIEW.md`, the `pre-pr-review` skill, or the reviewer model is kept only if the score does not drop.

## Files

- `cases.json`: the cases. Fields per case: `id`, `title`, `area`, `bug_commit`, `diff_cmd` (reproduces the introducing diff), `fix_commit` (null while the bug is still on `main`), `expected_finding` (`file` + one-line `failure`), and `note` where the case needs a caveat.
- `score.py`: records one run. Appends to `results.duckdb::eval_runs` (gitignored) and prints every run's score.

## Run

Ask Claude to "score the reviewer" (the `pre-pr-review` skill, section "Score the reviewer"). It reviews each `diff_cmd` blind, grades the findings against `expected_finding`, then records:

```bash
python review_evals/score.py --model fable --caught A,B --missed C --note "what changed"
```

## Add a case

When a bug reaches `main` that the review missed: find the commit that introduced it (`git log -S'<buggy text>' -- <path>`), add a case with the exact `diff_cmd`, and add one line to `REVIEW.md` for the context the reviewer lacked. Product decisions that were reverted are not cases; only code defects a diff reviewer could catch.
