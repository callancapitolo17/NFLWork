---
name: pre-pr-review
description: >-
  The review a change gets before a PR opens, before review is re-requested,
  and before a feature branch merges. Runs the repo's checks and shows the
  output, tiers the change by what it touches, runs 4 reviewer agents plus a
  verifier on the other model against the repo's REVIEW.md, tags findings
  Important / Nit / Pre-existing, and caps rounds at 3. Load this when Callan
  says "review", "ready for a PR", "open the PR", "pre-merge review",
  "merge it", or "re-request review", and before any message that calls a PR
  ready to merge.
---

# Pre-PR review

This copy lives in the repo so cloud and project threads load it too; it is the same procedure as the user-level skill on Cal's Mac, plus the "Cloud and project threads" section at the end. This review is required by the repo's CLAUDE.md, so the subagents in steps 2 and 3 are required, not optional.

If the repo has a `REVIEW.md`, read it first: it says what failure looks like in each area, the diff base, which checks to run, which paths are High tier, and the reviewer focuses. Without one, use the repo's CLAUDE.md and the defaults below.

## 0. Run the checks and show the output

Run the checks the change touches (REVIEW.md's list, else the repo CLAUDE.md's lint, type check, and test commands) and paste the last lines. Never report a check as passing without its output. Name any touched file no check covers; silence is not a pass.

A failure called "also fails on main" needs proof: run the same check on the branch's merge-base commit and paste that output too.

## 1. Tier the change by what it touches, not who wrote it

- **High**: REVIEW.md's High tier list when it has one. Without one: code that moves money (orders, stakes, quotes, payments); patient or user data and its audit rows; auth, sessions, tokens; anything that changes what a patient or user hears or sees; CI workflows, deploy, schedulers.
- **Low**: docs or copy only, no code and no tests.
- **Medium**: everything else.

| Tier | Review |
|---|---|
| Low | Step 0 only. |
| Medium | Steps 2 to 4. |
| High | Steps 2 to 4, then `/security-review`, then Callan reads the diff and the step-0 output before the PR opens or the branch merges. |

## 2. Hunt

Launch 4 reviewer agents in one message, no more, on the other model: `model: "fable"` if you are Opus, `model: "opus"` if you are Fable. A fresh reviewer catches what the author misses. Each gets the diff (against REVIEW.md's diff base, else `git diff $(git merge-base origin/main HEAD)..HEAD`; with no commits on the branch yet, `git add -N .` first so new files show, then `git diff HEAD`), the repo's CLAUDE.md, REVIEW.md when there is one, and one focus. The focuses are REVIEW.md's when it lists them, else:

1. **Reuse**: for each new function, file, field, constant, or setting, find an existing one in this repo, the repo it talks to, or the library.
2. **Cut**: comments, docstrings, and code the change does not need.
3. **Gaps**: overlapping or interrupted calls, error paths, patient data reaching logs or tokens, missing audit rows.
4. **Tests and PR text**: tests that assert too little; PR description, ticket, or paired PR out of sync with the code.

Tell each agent: every finding cites `file:line` and the concrete failure (this input or state, this wrong result); at most 3 nits; report only, fix nothing.

## 3. Verify

One verifier agent, on the other model, gets every finding and the diff. For each it reads the cited code and returns CONFIRMED or REJECTED with a one-line reason, and runs `git blame` on the cited lines: older than the merge-base is **Pre-existing**.

## 4. Report, then fix

First line is the count: "2 to fix, 1 pre-existing, 3 skipped", or "No changes needed."

- **Important**: a real bug, or something a user, patient, or bet would see. Fix before the PR.
- **Nit**: fix if it is one line, otherwise skip.
- **Pre-existing**: not this branch's bug. List it and offer a separate issue; do not fix it here.
- **Skipped**: every rejected finding with its reason, so Callan can overrule.

Fix the class, not the line: search the diff for every other instance of the same pattern.

## 5. Rounds

- Round 1 reviews the whole branch. Later rounds review only the branch's own commits since the last round, `git log -p --no-merges --first-parent <last reviewed sha>..HEAD`, which leaves out commits merged in from main, and never re-raise a finding already reported or rejected.
- Stop after 3 rounds. If round 3 still has an Important finding, the approach is the problem: bring it to Callan before writing more fixes.
- Commits added after an approval (fixes, hand-resolved merges) get a diff-only round before re-requesting review or merging.

## 6. PR text

Fill in the `## Pre-merge review` section of `.github/pull_request_template.md` (tier, checks output, each reviewer's result, verifier, security review, findings and what was fixed), and add:

- **Risk**: the tier and the high-risk paths touched, so the reviewer reads those first.
- **Revert**: how to roll back, including any data, config, or scheduled job to undo.

A PR without that section filled in is not ready, and no message may call it ready.

## 7. Fix at the source

When a finding repeats from an earlier branch, or a bug reaches main that this review missed, name the context the reviewer lacked and draft one line for the repo's REVIEW.md (or CLAUDE.md when the repo has no REVIEW.md), citing the real incident. In a team repo it goes to the reviewer like any other change.

## Cloud and project threads

- Run the review before the first message that calls the PR ready, not after Cal asks for it (#249, #250 and #251 went to Cal with only `check.sh`, or a one-line self-check).
- A cloud clone has no local `main` ahead of `origin`: use `git diff $(git merge-base origin/main HEAD)..HEAD` as the diff base.
- `/security-review` needs `origin/HEAD`, which a fresh clone may lack: run `git remote set-head origin main` first (or point it at the merge-base) and run it before the PR is called ready, not after merge (#271 ran it late).
- High tier: before merging, or calling the PR ready to merge, get Cal's explicit word that he read the diff and the step-0 output. A bare "go" or "merge" is not that.
