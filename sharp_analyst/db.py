"""DuckDB storage layer for podcast transcript chunks and vector search."""

import hashlib
from datetime import datetime
from pathlib import Path

import duckdb

DB_PATH = Path(__file__).parent / "sharp_analyst.duckdb"

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS episodes (
    episode_id VARCHAR PRIMARY KEY,
    podcast_name VARCHAR,
    episode_title VARCHAR,
    episode_date DATE,
    vtt_path VARCHAR,
    chunk_count INTEGER,
    ingested_at TIMESTAMP
);

CREATE TABLE IF NOT EXISTS chunks (
    chunk_id VARCHAR PRIMARY KEY,
    episode_id VARCHAR,
    podcast_name VARCHAR,
    episode_title VARCHAR,
    episode_date DATE,
    speaker VARCHAR,
    start_sec FLOAT,
    end_sec FLOAT,
    text TEXT,
    embedding FLOAT[384],
    token_count INTEGER,
    ingested_at TIMESTAMP
);
"""


def get_connection(db_path=None, read_only=False):
    """Get a DuckDB connection."""
    path = str(db_path or DB_PATH)
    return duckdb.connect(path, read_only=read_only)


def init_database(db_path=None):
    """Create tables if they don't exist."""
    con = get_connection(db_path)
    try:
        con.execute(SCHEMA_SQL)
    finally:
        con.close()


def episode_exists(episode_id, db_path=None):
    """Check if an episode has already been ingested."""
    con = get_connection(db_path, read_only=True)
    try:
        result = con.execute(
            "SELECT 1 FROM episodes WHERE episode_id = ?", [episode_id]
        ).fetchone()
        return result is not None
    finally:
        con.close()


def make_episode_id(podcast_name, episode_title):
    """Generate a stable episode ID from podcast + title."""
    key = f"{podcast_name}::{episode_title}"
    return hashlib.sha256(key.encode()).hexdigest()[:16]


def make_chunk_id(podcast_name, episode_title, start_sec):
    """Generate a stable chunk ID."""
    key = f"{podcast_name}::{episode_title}::{start_sec:.1f}"
    return hashlib.sha256(key.encode()).hexdigest()[:16]


def insert_chunks(chunks, episode_meta, db_path=None):
    """Insert chunks and episode metadata into the database.

    Args:
        chunks: list of dicts with keys: chunk_id, speaker, start_sec, end_sec,
                text, embedding (list of floats), token_count
        episode_meta: dict with keys: episode_id, podcast_name, episode_title,
                      episode_date, vtt_path
    """
    con = get_connection(db_path)
    try:
        now = datetime.now()

        con.execute(
            """INSERT INTO episodes (episode_id, podcast_name, episode_title,
               episode_date, vtt_path, chunk_count, ingested_at)
               VALUES (?, ?, ?, ?, ?, ?, ?)""",
            [
                episode_meta["episode_id"],
                episode_meta["podcast_name"],
                episode_meta["episode_title"],
                episode_meta.get("episode_date"),
                episode_meta["vtt_path"],
                len(chunks),
                now,
            ],
        )

        for chunk in chunks:
            con.execute(
                """INSERT INTO chunks (chunk_id, episode_id, podcast_name,
                   episode_title, episode_date, speaker, start_sec, end_sec,
                   text, embedding, token_count, ingested_at)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                [
                    chunk["chunk_id"],
                    episode_meta["episode_id"],
                    episode_meta["podcast_name"],
                    episode_meta["episode_title"],
                    episode_meta.get("episode_date"),
                    chunk.get("speaker"),
                    chunk["start_sec"],
                    chunk["end_sec"],
                    chunk["text"],
                    chunk["embedding"],
                    chunk["token_count"],
                    now,
                ],
            )
    finally:
        con.close()


def search_chunks(query_embedding, top_k=8, db_path=None):
    """Search for the most similar chunks using cosine similarity.

    Args:
        query_embedding: list of 384 floats
        top_k: number of results to return
        db_path: optional database path

    Returns:
        list of dicts with chunk metadata and similarity score
    """
    con = get_connection(db_path, read_only=True)
    try:
        results = con.execute(
            """SELECT chunk_id, podcast_name, episode_title, episode_date,
                      speaker, start_sec, end_sec, text, token_count,
                      array_cosine_similarity(embedding, ?::FLOAT[384]) AS score
               FROM chunks
               ORDER BY score DESC
               LIMIT ?""",
            [query_embedding, top_k],
        ).fetchall()

        columns = [
            "chunk_id", "podcast_name", "episode_title", "episode_date",
            "speaker", "start_sec", "end_sec", "text", "token_count", "score",
        ]
        return [dict(zip(columns, row)) for row in results]
    finally:
        con.close()


def get_stats(db_path=None):
    """Get knowledge base statistics."""
    con = get_connection(db_path, read_only=True)
    try:
        chunk_count = con.execute("SELECT COUNT(*) FROM chunks").fetchone()[0]
        episode_count = con.execute("SELECT COUNT(*) FROM episodes").fetchone()[0]
        podcasts = con.execute(
            "SELECT DISTINCT podcast_name FROM episodes ORDER BY podcast_name"
        ).fetchall()
        podcast_names = [r[0] for r in podcasts]
        return {
            "chunks": chunk_count,
            "episodes": episode_count,
            "podcasts": podcast_names,
        }
    except duckdb.CatalogException:
        return {"chunks": 0, "episodes": 0, "podcasts": []}
    finally:
        con.close()
