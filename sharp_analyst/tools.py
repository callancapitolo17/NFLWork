"""MCP tool for searching podcast transcripts via the Agent SDK."""

from pathlib import Path

from claude_agent_sdk import tool, create_sdk_mcp_server

from . import db

DB_PATH = db.DB_PATH

# Lazy-loaded embedding model (loaded on first search, cached after)
_model = None


def _get_model():
    """Load sentence-transformers model (cached after first call)."""
    global _model
    if _model is None:
        from sentence_transformers import SentenceTransformer
        _model = SentenceTransformer("all-MiniLM-L6-v2")
    return _model


def _format_results(results):
    """Format search results for the agent."""
    if not results:
        return "No relevant transcript chunks found."

    parts = []
    for i, r in enumerate(results, 1):
        source = f"[{r['podcast_name']}, {r['episode_title']}"
        if r.get("episode_date"):
            source += f", {r['episode_date']}"
        source += "]"

        speaker_prefix = f"{r['speaker']}: " if r.get("speaker") else ""

        time_start = _format_time(r["start_sec"])
        time_end = _format_time(r["end_sec"])

        parts.append(
            f"--- Result {i} (score: {r['score']:.3f}) ---\n"
            f"Source: {source} ({time_start}-{time_end})\n"
            f"{speaker_prefix}{r['text']}"
        )

    return "\n\n".join(parts)


def _format_time(seconds):
    """Format seconds as MM:SS."""
    m, s = divmod(int(seconds), 60)
    return f"{m}:{s:02d}"


@tool(
    "search_podcasts",
    "Search the podcast transcript database for relevant clips about sports "
    "betting strategy, market analysis, sharp betting concepts, and specific "
    "sports/matchup analysis. Returns the most relevant transcript chunks with "
    "source attribution.",
    {"query": str},
)
async def search_podcasts(args):
    """Search podcast transcripts by semantic similarity."""
    query = args["query"]

    # Embed the query
    model = _get_model()
    query_embedding = model.encode(query).tolist()

    # Search DuckDB
    results = db.search_chunks(query_embedding, top_k=8, db_path=DB_PATH)

    # Format for the agent
    formatted = _format_results(results)

    return {"content": [{"type": "text", "text": formatted}]}


def create_server():
    """Create the MCP server with the search_podcasts tool."""
    return create_sdk_mcp_server("sharp-analyst", tools=[search_podcasts])
