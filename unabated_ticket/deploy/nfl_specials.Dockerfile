# The fecta pricer's runtime: Python 3.12 plus nfl_specials/requirements.txt.
# No code and no secrets are copied in: compose mounts the repo read-only at
# /app, so an image holds nothing private and a code update is a restart, not
# a rebuild. Wheels only (--only-binary): every requirement ships a manylinux
# aarch64 wheel, so the Ampere VM compiles nothing. Build context is the repo
# root; nfl_specials.Dockerfile.dockerignore keeps every file but
# nfl_specials/requirements.txt out of it.
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

COPY nfl_specials/requirements.txt /tmp/requirements.txt
RUN pip install --only-binary=:all: -r /tmp/requirements.txt && rm /tmp/requirements.txt

WORKDIR /app
CMD ["python", "-m", "nfl_specials.app"]
