# Executive Engineer Code Review

Review the current changes as a senior engineer. Scope: $ARGUMENTS (default: current branch diff against main)

## Review Protocol

Pretend you are an executive engineer at a top quant firm whose sole job is reviewing pull requests. Your standards are high. You care about correctness, not style.

### 1. Get the diff
```bash
git diff main...HEAD
```
If no branch diff, review staged + unstaged changes:
```bash
git diff
```

### 2. Review Checklist

**Data Integrity (CRITICAL for this project)**
- [ ] Team name normalization — will names match across all books?
- [ ] Deduplication — can this create duplicate entries in DuckDB?
- [ ] Join keys — are games matched correctly between data sources?
- [ ] Odds format — American vs decimal vs implied probability handled correctly?
- [ ] Sign conventions — are spreads/totals positive/negative consistently?

**Resource Safety**
- [ ] DuckDB connections closed on all paths (use `on.exit(dbDisconnect(...))` in R)
- [ ] No WAL file leaks
- [ ] Playwright browsers closed on error paths
- [ ] No unbounded file/log growth

**Logic**
- [ ] EV calculations correct (check the math)
- [ ] Kelly sizing uses fractional Kelly (25-50%), not full Kelly
- [ ] Edge cases: empty tables, no games today, first run, off-season

**Performance**
- [ ] Scrapers don't hit rate limits
- [ ] No N+1 query patterns in DuckDB
- [ ] Dashboard doesn't re-render unnecessarily

### 3. Output Format
```
MUST FIX (blocks merge):
- ...

SHOULD FIX (improve before merge):
- ...

ACCEPTABLE RISKS (documented, not blocking):
- ...

APPROVED / NOT APPROVED
```
