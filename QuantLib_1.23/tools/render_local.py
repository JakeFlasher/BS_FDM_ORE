#!/usr/bin/env python3
"""
Flatten a local directory into an XML document for LLM consumption,
and generate a summary HTML page showing file statistics and skip reasons.

No external dependencies required — uses only the Python standard library.

Usage examples:
    python renderlocal.py ./my-project
    python renderlocal.py ./my-project -o output.xml --html-out summary.html
    python renderlocal.py ./my-project --max-bytes 100000 --ignore node_modules __pycache__ .venv
    python renderlocal.py ./my-project --no-html          # XML only, skip HTML summary
"""

from __future__ import annotations

import argparse
import html
import os
import pathlib
import subprocess
import sys
import webbrowser
from dataclasses import dataclass
from typing import List


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MAX_DEFAULT_BYTES = 50 * 1024  # 50 KiB

BINARY_EXTENSIONS = {
    ".png", ".jpg", ".jpeg", ".gif", ".webp", ".bmp", ".svg", ".ico",
    ".pdf", ".zip", ".tar", ".gz", ".bz2", ".xz", ".7z", ".rar",
    ".mp3", ".mp4", ".mov", ".avi", ".mkv", ".wav", ".ogg", ".flac",
    ".ttf", ".otf", ".eot", ".woff", ".woff2",
    ".so", ".dll", ".dylib", ".class", ".jar", ".exe", ".bin",
}

DEFAULT_IGNORE = [".git", ".am"]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class RenderDecision:
    include: bool
    reason: str  # "ok" | "binary" | "too_large" | "ignored"


@dataclass
class FileInfo:
    path: pathlib.Path  # absolute path on disk
    rel: str            # path relative to root (forward-slash separated)
    size: int
    decision: RenderDecision


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def bytes_human(n: int) -> str:
    """Human-readable byte size with 1 decimal for KiB+."""
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    f = float(n)
    i = 0
    while f >= 1024.0 and i < len(units) - 1:
        f /= 1024.0
        i += 1
    if i == 0:
        return f"{int(f)} {units[i]}"
    return f"{f:.1f} {units[i]}"


def looks_binary(path: pathlib.Path) -> bool:
    """Heuristic: known extension, null bytes, or non-UTF-8."""
    if path.suffix.lower() in BINARY_EXTENSIONS:
        return True
    try:
        with path.open("rb") as f:
            chunk = f.read(8192)
        if b"\x00" in chunk:
            return True
        try:
            chunk.decode("utf-8")
        except UnicodeDecodeError:
            return True
        return False
    except Exception:
        return True


def read_text(path: pathlib.Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# File scanning
# ---------------------------------------------------------------------------

def decide_file(
    path: pathlib.Path,
    root: pathlib.Path,
    max_bytes: int,
    ignore_patterns: List[str],
) -> FileInfo:
    rel = str(path.relative_to(root)).replace(os.sep, "/")
    try:
        size = path.stat().st_size
    except FileNotFoundError:
        size = 0

    # Match ignore patterns against individual path components
    rel_wrapped = f"/{rel}/"
    for pat in ignore_patterns:
        if f"/{pat}/" in rel_wrapped or rel == pat:
            return FileInfo(path, rel, size, RenderDecision(False, "ignored"))

    if size > max_bytes:
        return FileInfo(path, rel, size, RenderDecision(False, "too_large"))
    if looks_binary(path):
        return FileInfo(path, rel, size, RenderDecision(False, "binary"))
    return FileInfo(path, rel, size, RenderDecision(True, "ok"))


def collect_files(
    root: pathlib.Path, max_bytes: int, ignore_patterns: List[str],
) -> List[FileInfo]:
    infos: List[FileInfo] = []
    for p in sorted(root.rglob("*")):
        if p.is_symlink():
            continue
        if p.is_file():
            infos.append(decide_file(p, root, max_bytes, ignore_patterns))
    return infos


# ---------------------------------------------------------------------------
# Directory tree
# ---------------------------------------------------------------------------

def _generate_tree_fallback(root: pathlib.Path, ignore_patterns: List[str]) -> str:
    lines: List[str] = []

    def should_skip(name: str) -> bool:
        return name in ignore_patterns

    def walk(dir_path: pathlib.Path, prefix: str = "") -> None:
        entries = [e for e in dir_path.iterdir() if not should_skip(e.name)]
        entries.sort(key=lambda e: (not e.is_dir(), e.name.lower()))
        for i, e in enumerate(entries):
            last = i == len(entries) - 1
            branch = "└── " if last else "├── "
            lines.append(prefix + branch + e.name)
            if e.is_dir():
                walk(e, prefix + ("    " if last else "│   "))

    lines.append(root.name)
    walk(root)
    return "\n".join(lines)


def generate_tree(root: pathlib.Path, ignore_patterns: List[str]) -> str:
    """Try the system `tree` command; fall back to a pure-Python version."""
    try:
        cmd = ["tree", "-a"]
        if ignore_patterns:
            cmd.extend(["-I", "|".join(ignore_patterns)])
        cmd.append(".")
        cp = subprocess.run(
            cmd, cwd=str(root), check=True, text=True, capture_output=True,
        )
        return cp.stdout
    except Exception:
        return _generate_tree_fallback(root, ignore_patterns)


# ---------------------------------------------------------------------------
# XML generation
# ---------------------------------------------------------------------------

def _cdata_wrap(text: str) -> str:
    """Wrap text in CDATA, safely splitting on any embedded ]]> sequences."""
    return "<![CDATA[" + text.replace("]]>", "]]]]><![CDATA[>") + "]]>"


def generate_xml(infos: List[FileInfo]) -> str:
    """Build a full XML document with one <document> per rendered file."""
    lines = ['<?xml version="1.0" encoding="UTF-8"?>', "<documents>"]
    rendered = [i for i in infos if i.decision.include]

    for index, info in enumerate(rendered, 1):
        lines.append(f'  <document index="{index}">')
        lines.append(f"    <source>{html.escape(info.rel)}</source>")
        try:
            text = read_text(info.path)
            lines.append(f"    <document_content>{_cdata_wrap(text)}</document_content>")
        except Exception as e:
            lines.append(
                f"    <document_content>[read error: {html.escape(str(e))}]</document_content>"
            )
        lines.append("  </document>")

    lines.append("</documents>")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Summary HTML generation
# ---------------------------------------------------------------------------

def build_summary_html(
    root: pathlib.Path,
    infos: List[FileInfo],
    xml_out_path: str,
    ignore_patterns: List[str],
    max_bytes: int,
) -> str:
    rendered     = [i for i in infos if i.decision.include]
    skip_binary  = [i for i in infos if i.decision.reason == "binary"]
    skip_large   = [i for i in infos if i.decision.reason == "too_large"]
    skip_ignored = [i for i in infos if i.decision.reason == "ignored"]
    total_skipped = len(skip_binary) + len(skip_large) + len(skip_ignored)
    rendered_size = sum(i.size for i in rendered)
    total_size    = sum(i.size for i in infos)

    tree_text = generate_tree(root, ignore_patterns)

    # ---- helper: build <tr> rows ----
    def _rendered_rows() -> str:
        return "\n".join(
            f'<tr><td><code>{html.escape(i.rel)}</code></td>'
            f'<td class="size">{bytes_human(i.size)}</td></tr>'
            for i in rendered
        )

    def _skip_rows(items: List[FileInfo], label: str) -> str:
        return "\n".join(
            f'<tr><td><code>{html.escape(i.rel)}</code></td>'
            f'<td class="size">{bytes_human(i.size)}</td>'
            f"<td>{label}</td></tr>"
            for i in items
        )

    all_skip_rows = (
        _skip_rows(skip_binary,  "Binary file")
        + _skip_rows(skip_large, f"Exceeds {bytes_human(max_bytes)}")
        + _skip_rows(skip_ignored, "Ignored pattern")
    )

    rendered_block = (
        "<p>No files were rendered.</p>" if not rendered else
        f'<details open><summary>{len(rendered)} files · {bytes_human(rendered_size)}</summary>'
        f'<table><thead><tr><th>File</th><th>Size</th></tr></thead>'
        f'<tbody>{_rendered_rows()}</tbody></table></details>'
    )

    skipped_block = (
        "<p>No files were skipped.</p>" if total_skipped == 0 else
        f'<details open><summary>{total_skipped} files skipped</summary>'
        f'<table><thead><tr><th>File</th><th>Size</th><th>Reason</th></tr></thead>'
        f'<tbody>{all_skip_rows}</tbody></table></details>'
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>Render Summary – {html.escape(root.name)}</title>
<style>
  *, *::before, *::after {{ box-sizing: border-box; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto,
                 Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji';
    margin: 0; padding: 0; line-height: 1.5; color: #24292f; background: #fff;
  }}
  .container {{ max-width: 960px; margin: 0 auto; padding: 1.5rem; }}
  h1 {{ font-size: 1.5rem; margin: 0 0 0.25rem; }}
  h2 {{
    font-size: 1.2rem; margin: 1.5rem 0 0.5rem;
    border-bottom: 1px solid #d0d7de; padding-bottom: 0.3rem;
  }}
  .subtitle {{ color: #656d76; font-size: 0.95rem; margin-bottom: 1rem; }}

  /* stat cards */
  .stats {{ display: flex; flex-wrap: wrap; gap: 0.75rem; margin: 1rem 0; }}
  .stat-card {{
    background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 6px;
    padding: 0.75rem 1rem; min-width: 140px; flex: 1;
  }}
  .stat-card .label {{
    font-size: 0.8rem; color: #656d76;
    text-transform: uppercase; letter-spacing: 0.05em;
  }}
  .stat-card .value {{ font-size: 1.4rem; font-weight: 600; }}
  .stat-card.good  .value {{ color: #1a7f37; }}
  .stat-card.bad   .value {{ color: #cf222e; }}

  pre {{
    background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 6px;
    padding: 1rem; overflow: auto; font-size: 0.85em; line-height: 1.4;
  }}
  code {{
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas,
                 'Liberation Mono', 'Courier New', monospace;
  }}

  table {{ width: 100%; border-collapse: collapse; font-size: 0.9rem; margin-bottom: 1rem; }}
  th, td {{ text-align: left; padding: 0.4rem 0.75rem; border-bottom: 1px solid #d0d7de; }}
  th {{ background: #f6f8fa; font-weight: 600; }}
  td.size {{ white-space: nowrap; color: #656d76; }}
  tr:hover {{ background: #f6f8fa; }}

  details {{ margin: 0.5rem 0; }}
  summary {{ cursor: pointer; font-weight: 600; padding: 0.3rem 0; }}
  summary:hover {{ color: #0969da; }}

  .note {{
    background: #ddf4ff; border: 1px solid #54aeff; border-radius: 6px;
    padding: 0.75rem 1rem; margin: 1rem 0; font-size: 0.9rem;
  }}
  .note code {{ background: #fff; padding: 0.1rem 0.3rem; border-radius: 3px; }}
</style>
</head>
<body>
<div class="container">

  <h1>📄 Render Summary</h1>
  <div class="subtitle">
    Directory: <strong>{html.escape(str(root.resolve()))}</strong>
  </div>

  <div class="note">
    XML output written to: <code>{html.escape(xml_out_path)}</code>
  </div>

  <div class="stats">
    <div class="stat-card">
      <div class="label">Total Files</div>
      <div class="value">{len(infos)}</div>
    </div>
    <div class="stat-card good">
      <div class="label">Rendered</div>
      <div class="value">{len(rendered)}</div>
    </div>
    <div class="stat-card bad">
      <div class="label">Skipped</div>
      <div class="value">{total_skipped}</div>
    </div>
    <div class="stat-card">
      <div class="label">Rendered Size</div>
      <div class="value">{bytes_human(rendered_size)}</div>
    </div>
    <div class="stat-card">
      <div class="label">Total Size</div>
      <div class="value">{bytes_human(total_size)}</div>
    </div>
  </div>

  <h2>📂 Directory Tree</h2>
  <pre>{html.escape(tree_text)}</pre>

  <h2>✅ Rendered Files ({len(rendered)})</h2>
  {rendered_block}

  <h2>⏭️ Skipped Files ({total_skipped})</h2>
  {skipped_block}

  <h2>⚙️ Settings</h2>
  <table>
    <tr><td><strong>Max file size</strong></td><td>{bytes_human(max_bytes)}</td></tr>
    <tr><td><strong>Ignore patterns</strong></td>
        <td><code>{html.escape(", ".join(ignore_patterns) or "(none)")}</code></td></tr>
  </table>

</div>
</body>
</html>"""


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Flatten a local directory into XML for LLM consumption, "
                    "with an optional HTML summary of included/skipped files.",
    )
    ap.add_argument(
        "directory",
        help="Path to the local directory to scan",
    )
    ap.add_argument(
        "-o", "--xml-out",
        help="Output XML file path (default: <dirname>.xml in current dir)",
    )
    ap.add_argument(
        "--html-out",
        help="Output summary HTML path (default: <dirname>_summary.html)",
    )
    ap.add_argument(
        "--no-html", action="store_true",
        help="Skip summary HTML generation",
    )
    ap.add_argument(
        "--max-bytes", type=int, default=MAX_DEFAULT_BYTES,
        help=f"Max file size in bytes to include (default: {MAX_DEFAULT_BYTES} = "
             f"{bytes_human(MAX_DEFAULT_BYTES)})",
    )
    ap.add_argument(
        "--ignore", nargs="*", default=[],
        help="Additional directory/file name patterns to ignore "
             "(e.g. node_modules __pycache__ .venv)",
    )
    ap.add_argument(
        "--no-open", action="store_true",
        help="Don't auto-open the summary HTML in a browser",
    )
    args = ap.parse_args()

    # --- Validate input directory ---
    root = pathlib.Path(args.directory).resolve()
    if not root.is_dir():
        print(f"❌ Error: '{args.directory}' is not a directory.", file=sys.stderr)
        return 1

    dir_name = root.name or "root"
    ignore_patterns = DEFAULT_IGNORE + (args.ignore or [])

    # --- Derive default output paths ---
    xml_out  = args.xml_out  or f"flattened/{dir_name}.xml"
    html_out = args.html_out or f"flattened/{dir_name}_summary.html"

    # --- Scan ---
    print(f"📂 Scanning {root} ...", file=sys.stderr)
    infos = collect_files(root, args.max_bytes, ignore_patterns)
    rendered_count = sum(1 for i in infos if i.decision.include)
    skipped_count  = len(infos) - rendered_count
    print(
        f"✓  Found {len(infos)} files "
        f"({rendered_count} renderable, {skipped_count} skipped)",
        file=sys.stderr,
    )

    # --- Generate XML ---
    print("📝 Generating XML ...", file=sys.stderr)
    xml_text = generate_xml(infos)
    xml_path = pathlib.Path(xml_out)
    xml_path.write_text(xml_text, encoding="utf-8")
    print(
        f"💾 Wrote XML: {xml_path.resolve()} ({bytes_human(xml_path.stat().st_size)})",
        file=sys.stderr,
    )

    # --- Generate summary HTML ---
    if not args.no_html:
        print("📊 Generating summary HTML ...", file=sys.stderr)
        html_text = build_summary_html(
            root, infos, str(xml_path.resolve()), ignore_patterns, args.max_bytes,
        )
        html_path = pathlib.Path(html_out)
        html_path.write_text(html_text, encoding="utf-8")
        print(
            f"💾 Wrote HTML: {html_path.resolve()} "
            f"({bytes_human(html_path.stat().st_size)})",
            file=sys.stderr,
        )
        if not args.no_open:
            print("🌐 Opening summary in browser ...", file=sys.stderr)
            webbrowser.open(f"file://{html_path.resolve()}")

    print("✅ Done.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
