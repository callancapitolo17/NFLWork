# NFL Draft Portal — Improvements Wishlist (2027 Cycle)

**Purpose:** Living document of *additive* improvements to the NFL Draft
portal/dashboard identified during and after the 2026 draft. Every item
in this file is a layer on top of what the dashboard already does — not
a redesign or replacement of existing functionality. The current
dashboard works; these are the things that would take it to the next
level. Each entry captures *what* we want, *why* it's an edge or
workflow win, and *what we'd need to build it* — so we can pick this
back up cold next April without rediscovering the context.

## Design principles (from 2026 draft retrospective)

These two principles emerged from looking back at this year and should
shape how every item below is prioritized and built.

1. **Speed is the bottleneck, not edge.** The 2026 cycle confirmed there
   was plenty of soft, beatable line available across books — error
   prices, stale lines, mispriced over/unders. What stopped us from
   capturing more of it was *speed*: speed of getting information in
   (news, mocks, mispricings) and speed of getting bets out. Every item
   in this doc should be evaluated by how much it collapses the
   detect → decide → execute latency. If an improvement doesn't move the
   speed needle, it's lower priority than one that does.

2. **In 2026 the dashboard was a pre-draft tool only — making it live-
   capable is an open goal, not a closed door.** Lived experience this
   cycle: during the live draft the dashboard was *not used* because
   picks were happening too fast for the dashboard as it exists today
   to keep up. Where it delivered value was in the days leading up to
   the draft (research, finding edges), and where it could deliver
   more is after the fact (P&L, calibration, post-mortem).
   That said, the dashboard *can* be useful live — it just wasn't this
   year. Whether we get there is an active design question: do we
   invest in making the dashboard fast and ergonomic enough to be
   useful in a live picks-coming-off-the-board environment, or do we
   accept that live capture is machine-driven (automated execution
   item #2, loud alerts item #6, insider pipeline item #7) and let the
   dashboard be a research tool? Both are legitimate paths. Until we
   pick one, default the planning assumption to: "build automation for
   live, build the dashboard for everything else, and revisit live-
   dashboard utility as a separate decision."

## How to use this doc
- Add new items as numbered sections at the bottom. Don't edit closed items.
- Frame each item as an addition on top of the existing dashboard, not a
  rewrite of what's there.
- Each item should answer: Goal, Why it matters, Example/spec, Data required,
  Open questions.
- When an item ships, move it to the "Shipped" section at the bottom with a
  date and commit/PR reference.

---

## 1. Line Movement Tracking

**Status:** Open — captured 2026-04-23 (draft eve recap)

### Goal
Add a line-movement layer on top of the existing dashboard so it can show
*where the market is trending*, not just where it currently sits. The
current snapshot view stays as-is; this adds time-series tracking on top
so we can answer "which way is this moving and how fast" alongside
"what's the price right now."

### Why it matters
- **Edge:** The biggest +EV windows during draft week open while a line is
  mid-move and slow books haven't repriced. We can only see those windows
  if we're tracking direction and velocity, not just current price.
- **Read on the market:** Most of the actionable insight during draft week
  is in the deltas. Knowing a player's line tightened from "later" to "coin
  flip on pick 7" tells us something the static price doesn't.
- **Posting / sharing:** A clear movement view also gives us material to
  communicate the market read to others.

### Example of the kind of read this should enable
The text below is excerpted from a Rob Pizzola tweet during 2026 draft week
([source](https://x.com/robpizzola/status/2047428088649388366?s=20)),
written manually by reading multiple books across the day. The point isn't
to replicate this format exactly — it's to show the kind of trend-aware
narrative the dashboard should make easy to derive once we have proper
time-series of the lines:

> - Reese retook the favorite at pick 2 over Bailey. Full reversal from
>   earlier today when Bailey had pulled ahead. This one flipped hard.
> - Proctor reversed course. Was drifting later this morning. Now firming
>   earlier with real confidence. Market now has Proctor likely landing
>   before pick 17.
> - Tate has moved into the top 7 with conviction. His posted line tightened
>   from a slightly later number to essentially a coin flip on pick 7.
> - Cooper steamed earlier — line moved from low 20s into inside the top 22.
> - Peter Woods seemingly locked as DT1.
> - McCoy has barely moved. Still projects late R1 to early R2 corner.

In each of those bullets the underlying primitive is the same: *this line
was X at time T1 and is Y at time T2*. Build that primitive well and the
narrative reads write themselves.

### Foundation needed
Before any UI, get the data right:
- Confirm we are capturing dense enough snapshots of every market (across
  every book) to support useful lookback windows.
- Compare lines across books on a vig-free basis so movement comparisons
  are apples-to-apples.
- Stable player and market IDs across snapshots so a "delta" actually means
  the same line moved, not that a name got renamed.

### Open questions
- What lookback windows are most useful (last hour, last few hours, since
  open)? Probably worth supporting more than one.
- How do we want to surface this — chart per player, "movers" list, both?
- Player-pick vs position-rank vs round-bucket markets — do they all live
  in the same view or separate?
- Single-book moves vs cross-book consensus moves — same weight, or treat
  cross-book moves as the stronger signal?

### Before building
Pull a slice of this year's snapshot history and try to hand-write a
movement read off it. If the data isn't dense enough or doesn't join
cleanly across snapshots, the data layer is the first thing to fix —
not the UI.

---

## 2. Live Automated Betting on Kalshi During the Draft

**Status:** Open — captured 2026-04-23

### Goal
Add a live execution layer on top of the existing portal so we can fire
bets on Kalshi automatically as the draft unfolds, instead of clicking
through manually. Today the dashboard surfaces edges; a human still has to
act on them. During the live draft the windows are too short for that to
work reliably.

### Why it matters
- **Speed (lived experience from 2026):** During the live draft we
  literally could not click fast enough. The exact-pick markets moved
  too quickly to manually act on — by the time we identified an edge
  and clicked through, the price was gone. The over/under and bucket
  markets (top 10, first round, etc.) were the markets that *were*
  realistically attackable live, and even those were tight on time.
  This pattern reinforces item #3: automation should focus first on
  the bucket / OU markets where we already know the live edges exist.
- **Speed of execution was an explicit, user-stated bottleneck this
  cycle.** The 2026 retro identified that the soft edge was *there* —
  what limited us was not being fast enough to get the bets in. This
  item is the direct response. (See doc-level design principle #1.)
- **Consistency:** Pre-defining the rules ("if Player X is still on the
  board after pick N, hammer his under at Y price") removes in-the-moment
  hesitation and emotional decisions during the draft.
- **Coverage:** A human can only watch a handful of markets at once. An
  automated layer can watch every player simultaneously.

### Foundation needed
- Reliable Kalshi order-placement code (we likely have most of this from
  the MM bot work — verify how much is reusable for directional bets).
- A "playbook" definition format — what conditions trigger what bet at what
  price, with size and a max budget.
- Hard kill switch and per-bet/per-day limits before this ever touches
  real money.

### Open questions
- Which market types to automate first (over/unders are probably safer to
  automate than exact-pick markets given the strategy in item #3).
- How much of the existing Kalshi MM infrastructure (auth, order placement,
  position tracking) can be reused vs needs separate code.
- How do we handle Kalshi outages or rate limits mid-draft.
- Manual approval gate per bet, or fully autonomous within the playbook?

---

## 3. Strategy Framing: Where the Dashboard Should Bias Us

**Status:** Open — captured 2026-04-23

### Goal
This is less a dashboard feature and more an explicit strategy note that
should shape *how* the dashboard surfaces opportunities. The two distinct
betting modes we want the dashboard to support cleanly:

1. **Large swings on exact-position markets** — "who goes #6 overall",
   "who goes #7", "exact team-pick combos". Low hit rate, occasionally
   huge payouts. Lottery-ticket allocation.
2. **Consistent grind on bucket / over-under markets** — "over/unders on
   draft position", "first round yes/no", "top 5", "top 10". Higher hit
   rate, more capital deployed here, this is where the steady money comes
   from.

### Why it matters
Looking back at past results, the bread-and-butter winning has come from
the bucket / over-under markets, not the exact-position swings. The
dashboard today does not visibly distinguish those two modes, which makes
it easier to over-allocate to the long-shot exact picks because they look
exciting.

### Implication for dashboard
- Bucket and over-under markets should get first-class treatment in the
  UI, not be tucked under the exact-pick view.
- Sizing recommendations (or at least sizing *bands*) should reflect the
  mode of the market — e.g. exact-pick edges should default to small,
  bucket edges should default to standard sizing.

### Open questions
- Do we want the dashboard to actually enforce sizing differences, or
  just display the mode and let the user size?
- How do we cleanly tag every market as "exact" vs "bucket / OU"?
- Is there a third mode (props on team selection patterns,
  position-rank markets) that deserves its own bucket?

---

## 4. Mock Draft Consensus / Vegas Sheet as a Pricing Input

**Status:** Open — captured 2026-04-23

### Goal
Add an external-signal layer on top of the dashboard that pulls in
mock-draft consensus data (e.g. Grinding the Mocks–style aggregation) and
any shared "Vegas sheet" of insider mocks, and uses it as either a direct
pricing input or to drive a simulation we can compare against the market.

### Why it matters
- **Information edge:** Mock-draft consensus often leads the market — if
  we can quantify the consensus and compare it to the current market
  price, mismatches are a direct edge signal.
- **Simulation:** A draft-order simulator parameterized by mock consensus
  could let us derive our own implied probabilities for every market
  (every player's landing-pick distribution, every over/under, every
  team-pick combo) instead of only consuming book-derived implied probs.
- **Independent check:** Even when we don't act on it, having an
  internally-derived probability number to sanity-check against the
  market makes us less reliant on the consensus of books being right.

### Foundation needed
- Identify and pull from a reliable source of aggregated mock data
  (Grinding the Mocks is the obvious candidate; verify what's
  scrape-able and how often it updates).
- Decide what the "Vegas sheet" actually is in this context — is it a
  specific spreadsheet of insider mocks we already have access to, a
  paid product to subscribe to, or something we'd compile ourselves?
  (Clarify before next year.)
- A draft-order simulation framework (e.g. Monte Carlo over team
  preference distributions) if we want to convert mock data into market
  probabilities ourselves.

### Open questions
- Confirm what "Vegas sheet" refers to — likely a specific data source
  the user has in mind. Capture the source URL / contact when known.
- Is the goal to use this as a *pricing model* (we generate our own
  probs and bet anywhere we beat the book by X) or as a *signal*
  (flag markets where consensus and book disagree)?
- Refresh cadence — mock data evolves daily in the run-up to the draft;
  the pipeline needs to keep up.
- Weighting scheme — how do we trust newer mocks vs older ones, sharper
  analysts vs noisier ones?

---

## 5. Safer Trade Execution Interface for Kalshi

**Status:** Open — captured 2026-04-23

### Goal
Build a thin execution wrapper on top of Kalshi (either inside the
existing dashboard or as a separate trade panel) that makes misclicks
materially harder. Kalshi's native UI is too easy to act on the wrong
side of, and the cost of a single wrong click during live draft is
unacceptable.

### Why it matters
- **Lived experience from 2026:** A single buy-vs-sell misclick on
  Kalshi during the live draft cost thousands of dollars. The native
  interface doesn't visually differentiate or confirm enough for the
  speed at which we're acting. This was the single most expensive
  mistake of the night and it was a pure UX failure, not a bad read.
- **This risk does not go away with manual trading.** Even with the
  movement-tracking and pricing improvements above, every human-
  triggered order on Kalshi's UI is one slip away from another four-
  figure mistake. A safer interface caps the downside of being human.
- **Reinforces item #2.** The strongest case for automation isn't just
  speed — it's removing the human hand from the click entirely on the
  highest-value orders.

### Foundation needed
- Kalshi order-placement code (shared with item #2 — likely the same
  underlying client).
- A confirmation / size-limit / sanity-check layer:
  - Hard confirm gate above a configurable dollar threshold.
  - Visually unmistakable buy vs sell controls (color, position,
    label) — not the easy-to-confuse native layout.
  - Per-market max position to make over-sizing impossible.
  - Possibly a "review-then-fire" pattern where every order is staged
    in a list and committed in a separate explicit step.

### Open questions
- Build this as a panel inside our dashboard, or as a standalone Kalshi
  trade UI?
- What's the right dollar size for the hard confirm gate? (Per-trade
  vs per-market vs per-day cap.)
- Do we want this for *all* Kalshi trading or only during live draft?
- How much of this is solved by item #2 (just automate and remove the
  click) vs needing its own safer-UI track for trades we still want to
  fire by hand?

---

## 6. Loud Alerts for Drastically Mispriced / Error Lines

**Status:** Open — captured 2026-04-23

### Goal
Make obviously broken lines on any tracked book impossible to miss. Add
a high-priority alerting layer on top of the dashboard that surfaces
*catastrophically off-market* prices the moment they appear, with enough
context to act on them before the book corrects.

### Why it matters
- **Lived experience from 2026 (best win of the cycle):** The single
  biggest win of draft week was finding a line on Hoop88 that had been
  *incorrectly written* — a clear error, not a sharp disagreement. We
  caught it, hammered it, and it was a substantial profit. We caught
  this one. We almost certainly missed others on other books.
- **Highest-EV opportunities are book errors,** not consensus
  disagreements. When a book mistypes a line or fat-fingers a price,
  the EV per dollar can be enormous and the window before correction
  is short. We need to see these the second they appear.
- **Speed of info ingestion was a user-stated bottleneck.** Per
  design principle #1, getting information in faster is one of the
  two highest-leverage things we can build. Mispricing alerts are
  the cleanest application of that — direct from book → us → action,
  no intermediary.
- **Already partially in flight.** Existing plans
  (`2026-04-23-asymmetric-outlier-flags.md`,
  `2026-04-23-nfl-draft-staleness-filter.md`) address pieces of this.
  This item is the broader product framing those should ladder up to.

### Foundation needed
- Cross-book consensus per market (we have this) plus a robust
  "distance from consensus" metric that distinguishes:
  - **Error** — one book wildly out of line, others tightly clustered.
  - **Stale** — one book hasn't updated since news/movement on others.
  - **Sharp first-mover** — one book has moved and others will follow.
- A confidence layer so we're not crying wolf on every minor outlier.
  An alert is only useful if it's almost always real.
- An alert delivery mechanism that breaks through normal dashboard
  noise — push notification, dedicated panel, sound, etc.

### Open questions
- Threshold for "drastically mispriced" — absolute price gap, percent
  gap, vig-adjusted gap? Probably needs to be calibrated per market
  type (exact-pick markets are noisier than OUs).
- Should we differentiate alert *channels* by confidence (silent
  in-dashboard flag vs phone push vs everything-stops scream)?
- Are there book-specific signatures of error lines (round numbers,
  default prices, specific time-of-day patterns) that we can train
  on from this year's history?
- How do we keep a record of every detected mispricing, whether we
  acted or not, so we can backtest this detector next year?

---

## 7. Live Insider News Pipeline for Draft Night

**Status:** Open — captured 2026-04-23

### Goal
Build a real-time pipeline that ingests breaking insider posts (Schefter,
Rapoport, Pelissero, Garafolo, etc.) during the live draft and surfaces
"X going to Y at #N" calls the moment they're posted — fast enough to
act before the official pick is announced and before slow books reprice.
Distinct from item #4 (which is about *pre-draft* mock aggregation); this
item is specifically about the seconds-to-minutes window when an insider
breaks a pick.

### Why it matters
- **Speed of info ingestion was a user-stated bottleneck.** This is the
  other half of design principle #1 (the first half being execution
  speed, addressed by item #2). Insiders consistently break picks
  ahead of the official announcement and ahead of slow books
  repricing. That window is pure edge if we can ingest and act on it
  programmatically.
- **Live insider speed almost certainly demands automation.** Even if
  we eventually make the dashboard live-capable (see design principle
  #2 — open question), insider-post-to-pick windows are seconds long.
  This pipeline most likely only delivers full value if it feeds
  directly into automated alerts (item #6) or automated execution
  (item #2). A pure "look at this in the dashboard" UI is unlikely to
  solve the speed problem on this specific signal.
- **Source trust matters.** A Schefter pick is not the same signal as
  a random aggregator. Source-level confidence ranking has to be part
  of the product, not bolted on later.

### Foundation needed
- Reliable real-time ingest of a curated set of insider X/Twitter
  accounts (consider X API rate limits and pricing, scraping risk,
  Bluesky/Threads as backups).
- Source ranking / confidence — break out tier-1 (Schefter, Rapoport)
  from tier-2 from noise, with weights informed by historical
  accuracy if we can backtest it.
- Text parsing for the "Source: Player X to Team Y at pick N" pattern.
  Insider posts are templated enough that a simple parser plus an
  LLM fallback should cover most cases.
- A direct hand-off into either alerts (item #6) or execution
  (item #2). Without that, the pipeline is just a slower Twitter feed.

### Open questions
- X/Twitter API access — is the v2 API enough at scale, do we need
  paid tier, or do we scrape?
- Auto-act on insider posts vs human-in-the-loop confirmation? At
  Schefter speed, even a one-second human gate is probably too slow.
- How do we handle false alarms, bait posts, retweets that look like
  original calls, mistakes from tier-1 sources?
- Backtest harness using last cycle's insider post timeline so we can
  *quantify* the edge before going live with it.

---

## 8. Build for Portability Across Drafts (NBA, MLB, etc.)

**Status:** Open — captured 2026-04-23

### Goal
Design the next round of infrastructure (items #1–7) to be portable
across draft markets — NBA draft (June), MLB draft (July), NHL draft,
even derivative markets like awards or trade-deadline activity. The
same playbook (line-movement tracking, mispricing alerts, automated
execution, insider ingest) applies to anywhere there's a structured
event with discrete outcomes priced across multiple books.

### Why it matters
- **Same playbook, more shots on goal.** The NFL draft is one weekend
  a year. NBA draft, MLB draft, and NHL draft are all within a few
  months of it, and the same edge structure (book error, stale
  lines, mock consensus, insiders) exists in those markets too.
  Building NFL-only is leaving a multiple of the value on the table
  for the same engineering investment.
- **Most edge structure is shared.** The specifics differ (positions,
  teams, draft order) but movement tracking, mispricing alerts,
  automated execution, and insider ingestion are sport-agnostic. If
  we hard-code NFL terminology into every module, every other sport
  becomes a costly rewrite.
- **Near-term proof point.** NBA draft hits in late June / early
  July — a small, lower-stakes opportunity to test the post-NFL
  improvements before the next big NFL cycle.

### Foundation needed
- Sport-agnostic abstractions — player IDs, market types, team
  mappings, position taxonomies — instead of NFL-specific schema.
- A per-sport configuration layer where the playbook (which markets
  matter, which insiders to follow, which books cover) is data,
  not code.
- Inventory of what in current code is NFL-specific and what isn't,
  before adding any sports.

### Open questions
- Order of expansion: NBA draft (cheapest test) → MLB → NHL →
  awards futures? Or pick the sport with the softest existing
  market first?
- Do team-needs / position-rank markets generalize across sports,
  or do they need bespoke per-sport models?
- How much of item #4's "Vegas sheet" / mock consensus has equivalents
  in other sports, and where do we get those data sources?
- Master-schema-then-adapt vs rebuild-from-sport-neutral-base —
  which is the cheaper refactor given what we have today?

---

## Backlog (to be expanded)

Add additional improvement items below as numbered sections following the
same template (Goal / Why / Spec / Data / Open questions).

<!-- Template:
## N. <Short title>

**Status:** Open — captured YYYY-MM-DD

### Goal
### Why it matters
### Spec / example
### Data required
### Open questions
-->

---

## Shipped

<!-- When an item ships, move its section here with date + commit/PR ref:
- 2027-XX-XX — Item #N (commit abc1234 / PR #42)
-->
