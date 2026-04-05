"""Interactive CLI for the Sharp Bettor Analyst.

Usage:
    python -m sharp_analyst.cli
"""

import asyncio
import sys

from claude_agent_sdk import query, ClaudeAgentOptions, ResultMessage

from . import db
from .prompts import SHARP_BETTOR_PROMPT
from .tools import create_server


def print_stats():
    """Show knowledge base statistics."""
    stats = db.get_stats()
    if stats["chunks"] == 0:
        print("Knowledge base: empty")
        print("Run: python -m sharp_analyst.ingest --vtt-dir /path/to/vtt/ --podcast-name \"Name\"")
    else:
        podcasts = ", ".join(stats["podcasts"]) if stats["podcasts"] else "none"
        print(f"Knowledge base: {stats['episodes']} episodes, {stats['chunks']} chunks")
        print(f"Podcasts: {podcasts}")


async def chat(user_input, server, session_id=None):
    """Send a message and get the agent's response.

    Returns:
        (response_text, session_id) tuple
    """
    options = ClaudeAgentOptions(
        system_prompt=SHARP_BETTOR_PROMPT,
        mcp_servers={"sharp-analyst": server},
    )

    if session_id:
        options = ClaudeAgentOptions(resume=session_id)

    result_text = ""
    new_session_id = session_id

    async for message in query(prompt=user_input, options=options):
        if isinstance(message, ResultMessage):
            result_text = message.result
        elif hasattr(message, "type") and message.type == "system":
            if hasattr(message, "subtype") and message.subtype == "init":
                new_session_id = getattr(message, "data", {}).get("session_id")
                if not new_session_id:
                    new_session_id = getattr(message, "session_id", None)

    return result_text, new_session_id


async def repl():
    """Run the interactive REPL."""
    print("Sharp Bettor Analyst")
    print("=" * 40)
    print_stats()
    print("\nAsk me anything about betting strategy, market dynamics, or your models.")
    print("Type 'quit' to exit.\n")

    db.init_database()
    server = create_server()
    session_id = None

    while True:
        try:
            user_input = input("> ").strip()
        except (KeyboardInterrupt, EOFError):
            print("\nGoodbye.")
            break

        if not user_input:
            continue
        if user_input.lower() in ("quit", "exit", "q"):
            print("Goodbye.")
            break

        try:
            response, session_id = await chat(user_input, server, session_id)
            if response:
                print(f"\n{response}\n")
            else:
                print("\n(No response)\n")
        except Exception as e:
            print(f"\nError: {e}\n")
            session_id = None  # Reset session on error


def main():
    asyncio.run(repl())


if __name__ == "__main__":
    main()
