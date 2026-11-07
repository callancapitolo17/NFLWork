# The dry-run parlay quoter's runtime (#283): Python 3.12 plus
# kalshi_parlay_quoter/requirements.txt. Like the fecta pricer, no code and no
# secrets are copied in: compose mounts the repo read-only at /app. Wheels only.
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

COPY kalshi_parlay_quoter/requirements.txt /tmp/requirements.txt
RUN pip install --only-binary=:all: -r /tmp/requirements.txt && rm /tmp/requirements.txt

WORKDIR /app
CMD ["python", "-m", "kalshi_parlay_quoter.main"]
