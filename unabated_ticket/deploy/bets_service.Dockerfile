# The bets service's runtime: Python 3.12 plus bets_service/requirements.txt.
# No code and no secrets are copied in — compose mounts the repo read-only at
# /app and the credentials as files, so an image holds nothing private and a
# code update is a restart, not a rebuild. Wheels only (--only-binary): duckdb
# and cryptography ship manylinux aarch64 wheels, so the Ampere VM compiles
# nothing. Build context is unabated_ticket/; bets_service.Dockerfile.dockerignore
# keeps every file but requirements.txt out of it.
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

COPY bets_service/requirements.txt /tmp/requirements.txt
RUN pip install --only-binary=:all: -r /tmp/requirements.txt && rm /tmp/requirements.txt

WORKDIR /app
CMD ["python", "-m", "unabated_ticket.bets_service.service"]
