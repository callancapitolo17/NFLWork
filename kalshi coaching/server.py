#!/usr/bin/env python3
"""
Simple local server for the Kalshi Coaching Dashboard.
Serves the HTML and provides a /refresh endpoint to re-fetch data.

Run with: python server.py
Then open: http://localhost:5000
"""

import subprocess
import os
import mimetypes
from pathlib import Path
from flask import Flask, Response, jsonify

BASE_DIR = Path(__file__).parent.resolve()
app = Flask(__name__)


@app.route("/")
def index():
    """Serve the dashboard HTML."""
    html_path = BASE_DIR / "report.html"
    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()
    return Response(content, mimetype="text/html")


@app.route("/lib/<path:filename>")
def serve_lib(filename):
    """Serve static library files."""
    file_path = BASE_DIR / "lib" / filename
    if not file_path.exists():
        return "File not found", 404

    mimetype, _ = mimetypes.guess_type(str(file_path))
    with open(file_path, "rb") as f:
        content = f.read()
    return Response(content, mimetype=mimetype or "application/octet-stream")


@app.route("/refresh", methods=["POST"])
def refresh():
    """Re-run the scrapers and regenerate the report."""
    try:
        # Run Python scraper
        print("Running Python scraper...")
        result1 = subprocess.run(
            ["python3", str(BASE_DIR / "kalshi_coaching.py")],
            capture_output=True,
            text=True,
            cwd=BASE_DIR.parent
        )
        if result1.returncode != 0:
            return jsonify({
                "success": False,
                "error": f"Python scraper failed: {result1.stderr}"
            }), 500

        # Run R display script
        print("Running R display script...")
        result2 = subprocess.run(
            ["Rscript", str(BASE_DIR / "kalshi_coaching_display.R")],
            capture_output=True,
            text=True,
            cwd=BASE_DIR.parent
        )
        if result2.returncode != 0:
            return jsonify({
                "success": False,
                "error": f"R script failed: {result2.stderr}"
            }), 500

        print("Refresh complete!")
        return jsonify({"success": True, "message": "Data refreshed successfully"})

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


if __name__ == "__main__":
    print("\n" + "=" * 50)
    print("Kalshi Coaching Dashboard Server")
    print("=" * 50)
    print("\nOpen in browser: http://127.0.0.1:8080")
    print("Press Ctrl+C to stop\n")
    app.run(host="127.0.0.1", port=8080, debug=False)
