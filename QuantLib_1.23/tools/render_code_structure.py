#!/usr/bin/env python3
"""
C++ Code Structure Discovery

Imports render_local.py and produces per-subdirectory XML files separated
by header (.hpp/.h) and source (.cpp/.cc), plus a top-level snapshot.

Output layout
─────────────
    <out>/
    ├── hpp/
    │   ├── cashflows.xml          ← all .hpp/.h under root/cashflows/  (recursive)
    │   ├── currencies.xml
    │   ├── experimental.xml
    │   └── …
    ├── cpp/
    │   ├── cashflows.xml          ← all .cpp/.cc under root/cashflows/ (recursive)
    │   ├── currencies.xml
    │   ├── experimental.xml
    │   └── …
    └── <root>_current.xml         ← top-level .hpp + .cpp only (non-recursive)

Prerequisites
─────────────
    Place render_local.py (the revised flattening script) next to this file.
    If your copy is named "render_local (1).py", rename it to "render_local.py".

Usage
─────
    python discover_structure.py ./ql
    python discover_structure.py ./ql -o ql_output
    python discover_structure.py ./ql --max-bytes 200000 --ignore build test
"""

from __future__ import annotations

import argparse
import fnmatch
import html as html_mod
import pathlib
import sys
import webbrowser
from typing import List, Set, Tuple

# ---------------------------------------------------------------------------
# Import the companion flattening module
# ---------------------------------------------------------------------------
try:
    from render_local import (                    # type: ignore[import-untyped]
        collect_files,
        generate_xml,
        FileInfo,
        bytes_human,
        DEFAULT_IGNORE,
        MAX_DEFAULT_BYTES,
    )
except ImportError:
    sys.exit(
        "❌  Cannot import render_local.\n"
        "    Place render_local.py next to this script (or on PYTHONPATH).\n"
        "    If the file is named 'render_local (1).py', rename it first."
    )


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

HEADER_EXTS: Set[str] = {".hpp", ".h", ".hxx", ".hh", ".h++"}
SOURCE_EXTS: Set[str] = {".cpp", ".c", ".cxx", ".cc", ".c++"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _filter_by_exts(infos: List[FileInfo], exts: Set[str]) -> List[FileInfo]:
    """Keep only FileInfo entries whose suffix is in *exts*."""
    return [i for i in infos if i.path.suffix.lower() in exts]


def _write_xml(infos: List[FileInfo], dest: pathlib.Path) -> int:
    """Write an XML document to *dest*; return number of rendered files."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(generate_xml(infos), encoding="utf-8")
    return sum(1 for i in infos if i.decision.include)


def _should_skip_dir(name: str, patterns: List[str]) -> bool:
    """True when a directory name matches any ignore pattern."""
    for pat in patterns:
        if name == pat:
            return True
        if ("*" in pat or "?" in pat) and fnmatch.fnmatch(name, pat):
            return True
    return False


# ---------------------------------------------------------------------------
# Summary HTML builder
# ---------------------------------------------------------------------------

def _build_summary_html(
    root_name: str,
    root_abs: str,
    out_abs: str,
    top_rendered: int,
    top_xml: str,
    rows: List[Tuple[str, int, str, int, str]],   # (subdir, n_hpp, hpp_sz, n_cpp, cpp_sz)
    total_hdr: int,
    total_src: int,
    total_hdr_xml: int,
    total_src_xml: int,
    ignore_patterns: List[str],
    max_bytes: int,
) -> str:
    """Return a self-contained HTML summary page."""

    def _e(s: str) -> str:
        return html_mod.escape(s)

    # ── table rows ──
    tbody = []
    for subdir, nh, hsz, ns, ssz in rows:
        hpp_cell = f"{nh} files ({hsz})" if nh else '<span class="muted">—</span>'
        cpp_cell = f"{ns} files ({ssz})" if ns else '<span class="muted">—</span>'
        tbody.append(
            f"<tr>"
            f'<td><code>{_e(subdir)}/</code></td>'
            f"<td>{hpp_cell}</td>"
            f"<td>{cpp_cell}</td>"
            f"<td>{nh + ns}</td>"
            f"</tr>"
        )
    tbody_html = "\n".join(tbody)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>Structure Discovery – {_e(root_name)}</title>
<style>
  *, *::before, *::after {{ box-sizing: border-box; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto,
                 Helvetica, Arial, sans-serif;
    margin: 0; padding: 0; line-height: 1.5; color: #24292f; background: #fff;
  }}
  .container {{ max-width: 1020px; margin: 0 auto; padding: 1.5rem; }}
  h1 {{ font-size: 1.5rem; margin: 0 0 0.25rem; }}
  h2 {{
    font-size: 1.15rem; margin: 1.5rem 0 0.5rem;
    border-bottom: 1px solid #d0d7de; padding-bottom: 0.25rem;
  }}
  .subtitle {{ color: #656d76; font-size: 0.92rem; margin-bottom: 1rem; }}
  .stats {{ display: flex; flex-wrap: wrap; gap: 0.75rem; margin: 1rem 0; }}
  .stat-card {{
    background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 6px;
    padding: 0.7rem 1rem; min-width: 130px; flex: 1;
  }}
  .stat-card .label {{
    font-size: 0.78rem; color: #656d76;
    text-transform: uppercase; letter-spacing: 0.04em;
  }}
  .stat-card .value {{ font-size: 1.35rem; font-weight: 600; }}
  .stat-card.hdr  .value {{ color: #0969da; }}
  .stat-card.src  .value {{ color: #1a7f37; }}
  .stat-card.top  .value {{ color: #8250df; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.9rem; margin-bottom: 1rem; }}
  th, td {{ text-align: left; padding: 0.4rem 0.75rem; border-bottom: 1px solid #d0d7de; }}
  th {{ background: #f6f8fa; font-weight: 600; }}
  td.r {{ text-align: right; }}
  tr:hover {{ background: #f6f8fa; }}
  code {{
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas,
                 'Liberation Mono', 'Courier New', monospace;
    font-size: 0.88em;
  }}
  .note {{
    background: #ddf4ff; border: 1px solid #54aeff; border-radius: 6px;
    padding: 0.75rem 1rem; margin: 1rem 0; font-size: 0.9rem;
  }}
  .note code {{ background: #fff; padding: 0.1rem 0.3rem; border-radius: 3px; }}
  .muted {{ color: #b0b8c0; }}
  .tree {{ margin: 1rem 0; }}
  .tree pre {{
    background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 6px;
    padding: 1rem; overflow: auto; font-size: 0.85em; line-height: 1.4;
  }}
</style>
</head>
<body>
<div class="container">
  <h1>🔍 C++ Structure Discovery</h1>
  <div class="subtitle">
    Source: <strong>{_e(root_abs)}</strong><br/>
    Output: <strong>{_e(out_abs)}</strong>
  </div>

  <div class="note">
    Header XMLs in <code>hpp/</code> &nbsp;·&nbsp;
    Source XMLs in <code>cpp/</code> &nbsp;·&nbsp;
    Top-level snapshot: <code>{_e(top_xml)}</code>
  </div>

  <div class="stats">
    <div class="stat-card hdr">
      <div class="label">Header XMLs</div>
      <div class="value">{total_hdr_xml}</div>
    </div>
    <div class="stat-card hdr">
      <div class="label">Header Files</div>
      <div class="value">{total_hdr}</div>
    </div>
    <div class="stat-card src">
      <div class="label">Source XMLs</div>
      <div class="value">{total_src_xml}</div>
    </div>
    <div class="stat-card src">
      <div class="label">Source Files</div>
      <div class="value">{total_src}</div>
    </div>
    <div class="stat-card top">
      <div class="label">Top-level</div>
      <div class="value">{top_rendered}</div>
    </div>
  </div>

  <h2>📂 Per-Subdirectory Breakdown</h2>
  <table>
    <thead>
      <tr>
        <th>Subdirectory</th>
        <th>Headers (.hpp/.h)</th>
        <th>Sources (.cpp/.cc)</th>
        <th>Total</th>
      </tr>
    </thead>
    <tbody>
      {tbody_html}
    </tbody>
  </table>

  <h2>⚙️ Settings</h2>
  <table>
    <tr><td><strong>Max file size</strong></td><td>{bytes_human(max_bytes)}</td></tr>
    <tr><td><strong>Ignore patterns</strong></td>
        <td><code>{_e(", ".join(ignore_patterns) or "(none)")}</code></td></tr>
    <tr><td><strong>Header extensions</strong></td>
        <td><code>{_e(", ".join(sorted(HEADER_EXTS)))}</code></td></tr>
    <tr><td><strong>Source extensions</strong></td>
        <td><code>{_e(", ".join(sorted(SOURCE_EXTS)))}</code></td></tr>
  </table>
</div>
</body>
</html>"""


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="C++ code structure discovery — produce separate hpp/cpp "
                    "XML files per subdirectory, plus a top-level snapshot.",
    )
    ap.add_argument(
        "directory",
        help="Root C++ source directory to scan (e.g. ./ql)",
    )
    ap.add_argument(
        "-o", "--output-dir", default=None,
        help="Output directory (default: <dirname>_structure/)",
    )
    ap.add_argument(
        "--max-bytes", type=int, default=MAX_DEFAULT_BYTES,
        help=f"Max file size in bytes "
             f"(default: {MAX_DEFAULT_BYTES} = {bytes_human(MAX_DEFAULT_BYTES)})",
    )
    ap.add_argument(
        "--ignore", nargs="*", default=[],
        help="Additional directory/file patterns to ignore (e.g. build test)",
    )
    ap.add_argument(
        "--no-html", action="store_true",
        help="Skip the summary HTML page",
    )
    ap.add_argument(
        "--no-open", action="store_true",
        help="Don't auto-open the summary HTML in a browser",
    )
    args = ap.parse_args()

    # ── validate ──────────────────────────────────────────────────────
    root = pathlib.Path(args.directory).resolve()
    if not root.is_dir():
        print(f"❌ '{args.directory}' is not a directory.", file=sys.stderr)
        return 1

    dir_name = root.name or "root"
    ignore = DEFAULT_IGNORE + (args.ignore or [])

    out = pathlib.Path(args.output_dir or f"{dir_name}_structure")
    hpp_dir = out / "hpp"
    cpp_dir = out / "cpp"

    # ── banner ────────────────────────────────────────────────────────
    print(f"📂 Root       : {root}", file=sys.stderr)
    print(f"📁 Output     : {out.resolve()}", file=sys.stderr)
    print(f"   ├── hpp/    (headers: {', '.join(sorted(HEADER_EXTS))})",
          file=sys.stderr)
    print(f"   ├── cpp/    (sources: {', '.join(sorted(SOURCE_EXTS))})",
          file=sys.stderr)
    print(f"   └── {dir_name}_current.xml  (top-level only)", file=sys.stderr)
    print(file=sys.stderr)

    total_hdr = 0
    total_src = 0
    total_hdr_xml = 0
    total_src_xml = 0
    html_rows: List[Tuple[str, int, str, int, str]] = []

    # ==================================================================
    # 1. Top-level (non-recursive) → <dir>_current.xml
    # ==================================================================
    top_infos = collect_files(root, args.max_bytes, ignore, no_recurse=True)
    top_relevant = _filter_by_exts(top_infos, HEADER_EXTS | SOURCE_EXTS)

    current_path = out / f"{dir_name}_current.xml"
    top_rendered = _write_xml(top_relevant, current_path)
    top_size = bytes_human(current_path.stat().st_size)

    print(
        f"  {'[top-level]':30s}  {top_rendered:4d} files  "
        f"→ {current_path.name} ({top_size})",
        file=sys.stderr,
    )

    # ==================================================================
    # 2. Per-subdirectory (recursive)
    # ==================================================================
    subdirs = sorted(
        d for d in root.iterdir()
        if d.is_dir() and not _should_skip_dir(d.name, ignore)
    )

    for sd in subdirs:
        sn = sd.name
        infos = collect_files(sd, args.max_bytes, ignore)

        hdr_infos = _filter_by_exts(infos, HEADER_EXTS)
        src_infos = _filter_by_exts(infos, SOURCE_EXTS)

        parts: List[str] = []
        nh, hsz, ns, ssz = 0, "", 0, ""

        if hdr_infos:
            nh = _write_xml(hdr_infos, hpp_dir / f"{sn}.xml")
            hsz = bytes_human((hpp_dir / f"{sn}.xml").stat().st_size)
            total_hdr += nh
            total_hdr_xml += 1
            parts.append(f"{nh} hpp ({hsz})")

        if src_infos:
            ns = _write_xml(src_infos, cpp_dir / f"{sn}.xml")
            ssz = bytes_human((cpp_dir / f"{sn}.xml").stat().st_size)
            total_src += ns
            total_src_xml += 1
            parts.append(f"{ns} cpp ({ssz})")

        label = ", ".join(parts) if parts else "(no hpp/cpp)"
        print(f"  {sn + '/':30s}  {label}", file=sys.stderr)
        html_rows.append((sn, nh, hsz, ns, ssz))

    # ==================================================================
    # 3. Summary
    # ==================================================================
    print(file=sys.stderr)
    print("═" * 56, file=sys.stderr)
    print(
        f"  Headers : {total_hdr:5d} files in {total_hdr_xml:3d} XMLs  (hpp/)",
        file=sys.stderr,
    )
    print(
        f"  Sources : {total_src:5d} files in {total_src_xml:3d} XMLs  (cpp/)",
        file=sys.stderr,
    )
    print(
        f"  Top-lvl : {top_rendered:5d} files"
        f"              ({current_path.name})",
        file=sys.stderr,
    )

    # ==================================================================
    # 4. Optional summary HTML
    # ==================================================================
    if not args.no_html:
        html_path = out / f"{dir_name}_summary.html"
        html_text = _build_summary_html(
            root_name=dir_name,
            root_abs=str(root),
            out_abs=str(out.resolve()),
            top_rendered=top_rendered,
            top_xml=current_path.name,
            rows=html_rows,
            total_hdr=total_hdr,
            total_src=total_src,
            total_hdr_xml=total_hdr_xml,
            total_src_xml=total_src_xml,
            ignore_patterns=ignore,
            max_bytes=args.max_bytes,
        )
        html_path.write_text(html_text, encoding="utf-8")
        print(
            f"  Summary : {html_path} ({bytes_human(html_path.stat().st_size)})",
            file=sys.stderr,
        )
        if not args.no_open:
            webbrowser.open(f"file://{html_path.resolve()}")

    print(f"✅ Done — {out.resolve()}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
