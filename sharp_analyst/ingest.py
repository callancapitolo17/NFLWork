"""Ingest Apple Podcasts VTT transcripts into the vector database.

Usage:
    python -m sharp_analyst.ingest --vtt-dir /path/to/vtt/ --podcast-name "Podcast Name"
    python -m sharp_analyst.ingest --vtt-file /path/to/episode.vtt --podcast-name "Podcast Name"
"""

import argparse
import re
from pathlib import Path

from sentence_transformers import SentenceTransformer

from . import db


# Timestamp pattern: handles both HH:MM:SS.mmm and MM:SS.mmm
TIMESTAMP_RE = re.compile(
    r"(?:(\d+):)?(\d{2}):(\d{2})\.(\d{3})\s*-->\s*(?:(\d+):)?(\d{2}):(\d{2})\.(\d{3})"
)

# Speaker label pattern: "Speaker 1:", "John:", etc. at start of text
SPEAKER_RE = re.compile(r"^([A-Za-z][A-Za-z0-9 ]{0,30}):\s*")

# VTT formatting tags
TAG_RE = re.compile(r"<[^>]+>")

# Target chunk size in estimated tokens (~1.3 tokens per word)
TARGET_CHUNK_TOKENS = 500
OVERLAP_SECONDS = 30


def parse_timestamp(h, m, s, ms):
    """Convert timestamp components to seconds."""
    hours = int(h) if h else 0
    return hours * 3600 + int(m) * 60 + int(s) + int(ms) / 1000


def parse_vtt(vtt_path):
    """Parse a VTT file into a list of cues.

    Returns:
        list of dicts: [{start_sec, end_sec, text, speaker}, ...]
    """
    text = Path(vtt_path).read_text(encoding="utf-8", errors="replace")

    # Strip BOM and WEBVTT header
    text = text.lstrip("\ufeff")
    if text.startswith("WEBVTT"):
        text = text.split("\n", 1)[-1] if "\n" in text else ""

    cues = []
    lines = text.split("\n")
    i = 0

    while i < len(lines):
        line = lines[i].strip()

        # Look for timestamp line
        match = TIMESTAMP_RE.search(line)
        if match:
            start_sec = parse_timestamp(
                match.group(1), match.group(2), match.group(3), match.group(4)
            )
            end_sec = parse_timestamp(
                match.group(5), match.group(6), match.group(7), match.group(8)
            )

            # Collect text lines until blank line or next timestamp
            text_lines = []
            i += 1
            while i < len(lines):
                tl = lines[i].strip()
                if not tl or TIMESTAMP_RE.search(tl):
                    break
                # Clean formatting tags
                tl = TAG_RE.sub("", tl)
                if tl:
                    text_lines.append(tl)
                i += 1

            cue_text = " ".join(text_lines).strip()
            if not cue_text:
                continue

            # Extract speaker if present
            speaker = None
            speaker_match = SPEAKER_RE.match(cue_text)
            if speaker_match:
                speaker = speaker_match.group(1).strip()
                cue_text = cue_text[speaker_match.end():]

            cues.append({
                "start_sec": start_sec,
                "end_sec": end_sec,
                "text": cue_text,
                "speaker": speaker,
            })
        else:
            i += 1

    return cues


def estimate_tokens(text):
    """Rough token estimate: ~1.3 tokens per word."""
    return int(len(text.split()) * 1.3)


def chunk_cues(cues, podcast_name, episode_title):
    """Group cues into chunks of ~TARGET_CHUNK_TOKENS with overlap.

    Returns:
        list of dicts ready for embedding (minus the embedding field)
    """
    if not cues:
        return []

    chunks = []
    current_texts = []
    current_start = cues[0]["start_sec"]
    current_speaker = cues[0].get("speaker")
    current_tokens = 0

    for cue in cues:
        cue_tokens = estimate_tokens(cue["text"])

        # If adding this cue would exceed target, finalize current chunk
        if current_tokens + cue_tokens > TARGET_CHUNK_TOKENS and current_texts:
            chunk_text = " ".join(current_texts)
            chunk_end = cues[cues.index(cue) - 1]["end_sec"] if cue != cues[0] else cue["end_sec"]

            chunks.append({
                "chunk_id": db.make_chunk_id(podcast_name, episode_title, current_start),
                "start_sec": current_start,
                "end_sec": chunk_end,
                "text": chunk_text,
                "speaker": current_speaker,
                "token_count": estimate_tokens(chunk_text),
            })

            # Start new chunk with overlap: find cues within OVERLAP_SECONDS of the boundary
            overlap_start = chunk_end - OVERLAP_SECONDS
            overlap_texts = []
            overlap_tokens = 0
            new_start = cue["start_sec"]

            for prev_cue in cues:
                if prev_cue["start_sec"] >= overlap_start and prev_cue["start_sec"] < cue["start_sec"]:
                    overlap_texts.append(prev_cue["text"])
                    overlap_tokens += estimate_tokens(prev_cue["text"])
                    new_start = min(new_start, prev_cue["start_sec"])

            current_texts = overlap_texts
            current_start = new_start
            current_tokens = overlap_tokens
            current_speaker = cue.get("speaker") or current_speaker

        current_texts.append(cue["text"])
        current_tokens += cue_tokens
        if cue.get("speaker"):
            current_speaker = cue["speaker"]

    # Final chunk
    if current_texts:
        chunk_text = " ".join(current_texts)
        chunks.append({
            "chunk_id": db.make_chunk_id(podcast_name, episode_title, current_start),
            "start_sec": current_start,
            "end_sec": cues[-1]["end_sec"],
            "text": chunk_text,
            "speaker": current_speaker,
            "token_count": estimate_tokens(chunk_text),
        })

    return chunks


def embed_chunks(chunks, model):
    """Add embedding vectors to chunks using sentence-transformers.

    Args:
        chunks: list of chunk dicts (modified in place)
        model: SentenceTransformer model instance
    """
    if not chunks:
        return

    texts = [c["text"] for c in chunks]
    embeddings = model.encode(texts, show_progress_bar=False, batch_size=64)

    for chunk, emb in zip(chunks, embeddings):
        chunk["embedding"] = emb.tolist()


def ingest_vtt(vtt_path, podcast_name, model, episode_date=None, db_path=None):
    """Ingest a single VTT file.

    Returns:
        number of chunks created, or 0 if already ingested
    """
    vtt_path = Path(vtt_path)
    episode_title = vtt_path.stem  # filename without extension

    episode_id = db.make_episode_id(podcast_name, episode_title)
    if db.episode_exists(episode_id, db_path):
        return 0

    cues = parse_vtt(vtt_path)
    if not cues:
        print(f"  No cues found in {vtt_path.name}, skipping")
        return 0

    chunks = chunk_cues(cues, podcast_name, episode_title)
    if not chunks:
        return 0

    embed_chunks(chunks, model)

    episode_meta = {
        "episode_id": episode_id,
        "podcast_name": podcast_name,
        "episode_title": episode_title,
        "episode_date": episode_date,
        "vtt_path": str(vtt_path),
    }

    db.insert_chunks(chunks, episode_meta, db_path)
    return len(chunks)


def main():
    parser = argparse.ArgumentParser(description="Ingest podcast VTT transcripts")
    parser.add_argument("--vtt-dir", help="Directory containing VTT files")
    parser.add_argument("--vtt-file", help="Single VTT file to ingest")
    parser.add_argument("--podcast-name", required=True, help="Name of the podcast")
    parser.add_argument("--db", help="Database path (default: sharp_analyst.duckdb)")
    args = parser.parse_args()

    if not args.vtt_dir and not args.vtt_file:
        parser.error("Provide --vtt-dir or --vtt-file")

    db_path = args.db or db.DB_PATH
    db.init_database(db_path)

    print("Loading embedding model...")
    model = SentenceTransformer("all-MiniLM-L6-v2")

    # Collect VTT files
    if args.vtt_file:
        vtt_files = [Path(args.vtt_file)]
    else:
        vtt_dir = Path(args.vtt_dir)
        vtt_files = sorted(vtt_dir.glob("*.vtt"))

    if not vtt_files:
        print("No VTT files found")
        return

    print(f"Found {len(vtt_files)} VTT file(s)")

    total_chunks = 0
    ingested = 0

    for vtt_path in vtt_files:
        print(f"  Processing: {vtt_path.name}")
        n = ingest_vtt(vtt_path, args.podcast_name, model, db_path=db_path)
        if n > 0:
            print(f"    → {n} chunks")
            total_chunks += n
            ingested += 1
        else:
            print(f"    → already ingested, skipping")

    print(f"\nDone: {ingested} episodes, {total_chunks} chunks ingested")
    stats = db.get_stats(db_path)
    print(f"Knowledge base: {stats['episodes']} episodes, {stats['chunks']} chunks from {len(stats['podcasts'])} podcast(s)")


if __name__ == "__main__":
    main()
