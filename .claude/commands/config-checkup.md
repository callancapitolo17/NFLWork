---
description: Audit Claude Code config health — CLAUDE.md size, skills, agents, and extraction candidates
---

Run a **read-only** health check of this repo's Claude Code configuration and
report a concise, scannable summary. Do NOT edit, move, or remove anything —
this is a report only. Recommend changes; let me decide what to act on.

Work through these steps:

1. **CLAUDE.md sizes.** Find every CLAUDE.md in the repo, excluding worktrees
   and vendored dirs:
   `find . -name CLAUDE.md -not -path '*/.claude/worktrees/*' -not -path '*/.worktrees/*' -not -path '*/node_modules/*'`
   Report each path with its line count (`wc -l`). **Flag any over ~200 lines** —
   Anthropic's guidance is that bloated CLAUDE.md files degrade instruction
   adherence ("Keep files < 200 lines").

2. **Skills & agents inventory.** List `.claude/skills/*/SKILL.md` and
   `.claude/agents/*.md`, each with its `name` and one-line `description` from
   the frontmatter. These are the on-demand pieces that keep CLAUDE.md lean —
   confirm they're present and their descriptions still read true.

3. **Extraction candidates.** Scan the root `CLAUDE.md` for long *procedural*
   blocks (step-by-step how-tos) or *reference* material (domain knowledge,
   formula dumps) that are only sometimes relevant. For each, apply the test:
   "Would removing this line cause Claude to make a mistake on an *unrelated*
   task?" If not, it's a candidate to move into a skill. Name the block and the
   skill it should live in. (Guardrails like branch hygiene / approval rules
   stay in CLAUDE.md — they must always apply.)

4. **Worktree hygiene (bonus).** Run `git worktree list`, and for each worktree
   report whether it is (a) merged into main AND clean → safe to remove, or
   (b) has unmerged commits / uncommitted / untracked work → keep, removal would
   lose work. Check uncommitted state with `git -C <wt> status --short` and
   untracked with `git -C <wt> ls-files --others --exclude-standard`. **Never
   remove anything** — just surface the list so nothing gets abandoned silently.

5. **Summary.** End with a short prioritized list: what (if anything) to trim,
   extract, or clean up — and what's healthy. If everything looks good, say so
   plainly.

Keep the whole report tight. Recommend; do not change anything without me asking.
