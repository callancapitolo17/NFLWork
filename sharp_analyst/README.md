# Sharp Bettor Analyst

RAG-powered sharp bettor analyst that searches podcast transcripts to inform betting analysis. Runs through your Claude Code Max subscription — no additional API costs.

## Setup

```bash
cd sharp_analyst
pip install -r requirements.txt
```

First install downloads the `all-MiniLM-L6-v2` embedding model (~80MB). PyTorch is required by sentence-transformers (~2GB one-time install).

## Ingest Transcripts

Get VTT transcript files from Apple Podcasts, then:

```bash
# Ingest a directory of VTT files
python -m sharp_analyst.ingest --vtt-dir /path/to/transcripts/ --podcast-name "Circles Off"

# Ingest a single file
python -m sharp_analyst.ingest --vtt-file /path/to/episode.vtt --podcast-name "Bet The Process"
```

Already-ingested episodes are automatically skipped.

## Run the Analyst

```bash
python -m sharp_analyst.cli
```

Ask questions about betting strategy, market dynamics, modeling decisions, etc. The agent automatically searches podcast transcripts when relevant.

### Example Questions
- "How should I think about pricing alt lines in CBB?"
- "What's the best approach to same-game parlay correlation?"
- "Is my CLV tracking methodology sound? Here's how I'm doing it..."
- "How do sharps think about totals in conference tournament week?"

## How It Works

1. **Ingest**: VTT transcripts are parsed, chunked (~500 tokens each), and embedded locally using sentence-transformers
2. **Store**: Chunks + embeddings stored in DuckDB with cosine similarity search
3. **Query**: When you ask a question, the agent can search transcripts via a `search_podcasts` tool
4. **Analyze**: Claude synthesizes podcast knowledge with its own expertise to give informed analysis

## Data

- Embeddings: `all-MiniLM-L6-v2` (384 dimensions, runs locally)
- Storage: `sharp_analyst.duckdb` (auto-created, gitignored)
- Vector search: DuckDB `array_cosine_similarity()` brute-force scan
