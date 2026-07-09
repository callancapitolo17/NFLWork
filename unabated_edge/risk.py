import datetime
from unabated_edge import config


def staleness_ok(generated_at, max_age_sec):
    now = datetime.datetime.now(datetime.timezone.utc)
    if generated_at.tzinfo is None:
        generated_at = generated_at.replace(tzinfo=datetime.timezone.utc)
    return (now - generated_at).total_seconds() <= max_age_sec


def tipoff_ok(commence_time, cancel_min, now=None):
    if commence_time is None:
        return True
    now = now or datetime.datetime.now(datetime.timezone.utc)
    return (commence_time - now).total_seconds() > cancel_min * 60


def kill_switch_ok():
    return not config.KILL_FILE.exists()
