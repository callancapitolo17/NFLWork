"""Record one run of the reviewer against review_evals/cases.json and print every run's score.

Inputs: --model (the reviewer model), --caught and --missed (comma-separated case ids
from cases.json), --note (what changed in the reviewer before this run).
Side effect: APPENDS one row per case to review_evals/results.duckdb::eval_runs
(run_time, reviewer_model, case_id, caught, note). Never updates or deletes.

    python review_evals/score.py --model fable --caught A,B,E --missed C --note "REVIEW.md: added money section"
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import duckdb

EVALS_DIR = Path(__file__).resolve().parent
CASES_PATH = EVALS_DIR / "cases.json"
RESULTS_DB_PATH = EVALS_DIR / "results.duckdb"


def parse_ids(text: str) -> list[str]:
    return [case_id.strip() for case_id in text.split(",") if case_id.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True)
    parser.add_argument("--caught", default="")
    parser.add_argument("--missed", default="")
    parser.add_argument("--note", default="")
    args = parser.parse_args()

    known_ids = {case["id"] for case in json.loads(CASES_PATH.read_text())["cases"]}
    caught = parse_ids(args.caught)
    missed = parse_ids(args.missed)
    unknown = (set(caught) | set(missed)) - known_ids
    if unknown:
        raise SystemExit(f"expected case ids from {sorted(known_ids)}, got unknown {sorted(unknown)}")
    unscored = known_ids - set(caught) - set(missed)
    if unscored:
        raise SystemExit(f"every case needs a verdict; missing {sorted(unscored)}")

    run_time = datetime.now(timezone.utc)
    rows = [(run_time, args.model, case_id, True, args.note) for case_id in caught]
    rows += [(run_time, args.model, case_id, False, args.note) for case_id in missed]

    with duckdb.connect(str(RESULTS_DB_PATH)) as con:
        con.execute(
            "CREATE TABLE IF NOT EXISTS eval_runs ("
            "run_time TIMESTAMPTZ, reviewer_model VARCHAR, case_id VARCHAR, caught BOOLEAN, note VARCHAR)"
        )
        con.executemany("INSERT INTO eval_runs (run_time, reviewer_model, case_id, caught, note) VALUES (?, ?, ?, ?, ?)", rows)
        history = con.execute(
            "SELECT run_time, reviewer_model, SUM(caught::INT) AS caught, COUNT(*) AS total, "
            "string_agg(CASE WHEN NOT caught THEN case_id END, ',' ORDER BY case_id) AS missed, any_value(note) AS note "
            "FROM eval_runs GROUP BY run_time, reviewer_model ORDER BY run_time"
        ).fetchall()

    for run_time, model, caught_count, total, missed_ids, note in history:
        print(f"{run_time:%Y-%m-%d %H:%M} {model:8} {caught_count}/{total} caught  missed: {missed_ids or '-'}  {note}")


if __name__ == "__main__":
    main()
