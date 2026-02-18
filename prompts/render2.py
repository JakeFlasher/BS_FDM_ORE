#!/usr/bin/env python3
"""
generate_prompts.py
===================
Assemble per-round LLM prompts from the multi-round QuantLib template.

For every round defined in the template the script:
  1. Inserts the shared preamble where [INSERT SHARED PREAMBLE …] appears.
  2. Locates actual QuantLib source files referenced by
     ``FILE: <path>`` / ``[INSERT FULL CONTENT …]`` and pastes their content.
  3. Replaces references to *previous-round* generated outputs with clearly
     marked attachment placeholders for the operator to fill in.

Usage
-----
    python generate_prompts.py  /path/to/QuantLib        # repo root (ql/ inside)
    python generate_prompts.py  /path/to/QuantLib/ql     # ql/ dir itself works too
    python generate_prompts.py  ../ql -t multi-round_prompt.md -o prompts -v
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

# ═══════════════════════════════════ CLI ═══════════════════════════════════


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate per-round prompts from the multi-round template.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "ql_dir",
        help="Path to the QuantLib source root (parent of ql/) "
        "or to the ql/ directory itself",
    )
    p.add_argument(
        "-t",
        "--template",
        default="multi-round_prompt.md",
        help="Path to the multi-round prompt template (default: %(default)s)",
    )
    p.add_argument(
        "-o",
        "--output-dir",
        default="generated_prompts",
        help="Output directory for the generated prompt files (default: %(default)s)",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    return p.parse_args()


# ═══════════════════════════════ File helpers ══════════════════════════════


def _read(path: str) -> str:
    """Read a text file, tolerating encoding quirks."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        return fh.read()


def _fix_known_typos(rel_path: str) -> str:
    """
    Fix specific typos found in the redo_the_paper_before_render_v2.md template.
    """
    # Fix directory typo in Round 4
    rel_path = rel_path.replace("finitediffere/", "finitedifferences/")
    
    # Fix extension typo in Round 2
    if rel_path.endswith("fdmbarrierprojectionconditipp"):
        rel_path = rel_path.replace(".conditipp", ".condition.cpp")
        
    # Fix extension typo in Round 1
    if rel_path.endswith("secondderivativeop.p"):
        rel_path = rel_path + "pp" # .p -> .cpp

    return rel_path


def _find_file(ql_dir: str, rel_path: str) -> str | None:
    """
    Resolve a template-style path such as ``ql/math/array.hpp``
    against the user-supplied *ql_dir*.  Tries three strategies so that
    pointing at either the repo root **or** the ``ql/`` directory works.
    """
    # Apply typo fixes before searching
    rel_path = _fix_known_typos(rel_path)

    base = str(Path(ql_dir).resolve())

    # Strategy 1: ql_dir is the repo root → <root>/ql/math/array.hpp
    c = os.path.join(base, rel_path)
    if os.path.isfile(c):
        return c

    # Strategy 2: ql_dir IS ql/ → strip leading ql/ from rel_path
    # Also handle cases where path starts with "FILE/" due to regex group capture issues
    clean_rel = rel_path
    if clean_rel.startswith(("ql/", "ql\\")):
        clean_rel = clean_rel[3:]
    
    c = os.path.join(base, clean_rel)
    if os.path.isfile(c):
        return c

    # Strategy 3: go one level up from ql_dir
    c = os.path.join(str(Path(base).parent), rel_path)
    if os.path.isfile(c):
        return c

    return None


def _find_by_name(ql_dir: str, filename: str) -> str | None:
    """Recursive search for *filename* under *ql_dir* and its parent."""
    for search_root in (str(Path(ql_dir).resolve()),
                        str(Path(ql_dir).resolve().parent)):
        for dirpath, _, filenames in os.walk(search_root):
            if filename in filenames:
                return os.path.join(dirpath, filename)
    return None


def _ql_relpath(abspath: str, ql_dir: str) -> str:
    """Best-effort ``ql/…`` style relative path for display."""
    a = str(Path(abspath).resolve())
    b = str(Path(ql_dir).resolve())
    try:
        rel = os.path.relpath(a, Path(b).parent)
    except ValueError:
        return os.path.basename(a)
    return rel.replace(os.sep, "/")


# ═══════════════════════════ Template parsing ═════════════════════════════

def _extract_preamble(doc: str) -> str:
    """Return the raw XML content of the SHARED PREAMBLE code block."""
    m = re.search(
        r"###\s+SHARED\s+PREAMBLE.*?```xml\s*\n(.*?)```", doc, re.DOTALL
    )
    if not m:
        sys.exit("ERROR: SHARED PREAMBLE section not found in template.")
    return m.group(1).rstrip("\n")

def _extract_rounds(doc: str) -> dict[int, dict]:
    """
    Return ``{round_num: {'title': …, 'content': …}}`` for every
    ``### ROUND N`` section found in the document.
    """
    rounds: dict[int, dict] = {}

    # The original template used --- horizontal rules between sections.
    # framework_0218.md uses ### headers instead.  Try --- first; if that
    # yields only a single blob, fall back to splitting on ### headers.
    sections = re.split(r"\n\s*-{3,}\s*\n", doc)
    if len(sections) <= 1:
        # Lookahead split keeps the ### header attached to its section
        sections = re.split(r"(?=^### )", doc, flags=re.MULTILINE)

    for section in sections:
        # Relaxed regex:
        # 1. Matches "ROUND <N>"
        # 2. Allows optional text like "(OPTIONAL)" before the dash
        # 3. Makes the dash and title optional (for "### ROUND 1")
        hdr = re.search(
            r"###\s+(?:Optional\s+)?ROUND\s+(\d+)(?:.*?)(?:[-—–]+\s*(.*))?",
            section,
            re.IGNORECASE
        )

        if not hdr:
            continue

        num = int(hdr.group(1))

        # Use the captured title if present, otherwise default to "Round N"
        title_text = hdr.group(2)
        title = title_text.strip() if title_text else f"Round {num}"

        xml_block = re.search(r"```xml\s*\n(.*?)```", section, re.DOTALL)
        if xml_block:
            rounds[num] = dict(title=title,
                               content=xml_block.group(1))
        else:
            print(f"  Warning: Round {num} — no ```xml block", file=sys.stderr)

    return rounds
# ═══════════════════ INSERT-directive classification ══════════════════════

# Matches [INSERT …] (may span multiple lines because of the [^\]]* class)
_RE_INSERT = re.compile(r"\[INSERT[^\]]*\]", re.DOTALL)

# Matches   FILE: <path>\n[INSERT …]
# MODIFIED: Allows "FILE/" or "FILE " or "FILE:" to handle template typos
_RE_FILE_INSERT = re.compile(
    r"(FILE[:\s/]+(\S+)[ \t]*\n)(\[INSERT[^\]]*\])", re.DOTALL
)

# Round-8 abbreviated source block: [... abbreviated — include … as source files ...]
_RE_ABBREVIATED = re.compile(
    r"\[\.\.\.\s*abbreviated\s*[-—–]+\s*include\s+(.*?)"
    r"\s+as\s+source\s+files\s*\.\.\.\s*\]",
    re.DOTALL | re.IGNORECASE,
)


def _is_source_insert(directive: str) -> bool:
    """True when the directive asks for actual source-file content."""
    return bool(re.search(r"INSERT\s+FULL\s+CONTENT", directive, re.IGNORECASE))


def _has_round_ref(directive: str) -> bool:
    """True when the directive references a previous round's output."""
    return bool(
        re.search(
            r"(?:THE\s+)?ROUND\s+\d+\s+OUTPUT|"
            r"from\s+Rounds?\s+\d|"
            r"all\s+\.hpp\s+files\s+from\s+Rounds?|"
            r"INSERT\s+(?:THE\s+)?ROUND\s+\d+",
            directive,
            re.IGNORECASE,
        )
    )


def _round_nums(text: str) -> list[int]:
    """Extract every round number mentioned in *text*."""
    nums: set[int] = set()
    for m in re.finditer(
        r"Rounds?\s+(\d+)(?:\s*[-/]\s*(\d+))?", text, re.IGNORECASE
    ):
        lo = int(m.group(1))
        hi = int(m.group(2)) if m.group(2) else lo
        nums.update(range(lo, hi + 1))
    return sorted(nums)


def _attachment_note(
    directive: str,
    ctx_round: int | None = None,
    file_hint: str | None = None,
) -> str:
    """
    Build a human-readable placeholder that tells the operator which
    previous-round output to paste here before sending the prompt.
    """
    inner = directive.strip("[] \n")
    nums = _round_nums(inner)

    # Try to pull a .hpp filename out of the directive text
    fm = re.search(r"INSERT\s+(?:THE\s+)?(\S+\.hpp)", inner)
    fname = fm.group(1) if fm else file_hint

    if not nums and ctx_round is not None:
        nums = [ctx_round]

    rs = ", ".join(str(n) for n in nums) if nums else "previous"

    if fname:
        return (
            f"[>>> ATTACH HERE: the {fname} file "
            f"generated in Round(s) {rs} <<<]"
        )
    return (
        f"[>>> ATTACH HERE: output file(s) "
        f"generated in Round(s) {rs} <<<]"
    )


# ══════════════════════ Per-round processing stages ═══════════════════════

# ── Stage 1: preamble insertion ──────────────────────────────────────────


def _stage_preamble(content: str, preamble: str) -> str:
    """Replace both flavours of preamble placeholder."""
    # Standard form used in Rounds 1-7
    content = content.replace("[INSERT SHARED PREAMBLE HERE]", preamble)

    # Round-8 form: modify the preamble to acknowledge the legacy framework
    r8 = re.compile(
        r"\[INSERT\s+SHARED\s+PREAMBLE,\s*replacing\s+the\s+line\s+about\s+"
        r"not\s+using\s+legacy\s+framework\s+with:\]",
        re.IGNORECASE,
    )
    if r8.search(content):
        modified_preamble = preamble.replace(
            "(not the deprecated legacy MixedScheme framework)",
            "(including the deprecated legacy MixedScheme/"
            "FiniteDifferenceModel framework — see note below)",
        )
        content = r8.sub(modified_preamble, content)
    return content


# ── Stage 2: FILE: <path> + [INSERT …] pairs ────────────────────────────


def _stage_file_inserts(content: str, ql_dir: str, verbose: bool) -> str:
    """Resolve every ``FILE: path\\n[INSERT …]`` block."""

    def _repl(m: re.Match) -> str:
        file_line: str = m.group(1)   # "FILE: <path>\n"
        file_path: str = m.group(2)   # the ql/… path
        directive: str = m.group(3)   # "[INSERT …]"
        
        # Clean up file path if it was captured via FILE/path
        if file_path.startswith("FILE/"):
            file_path = file_path.replace("FILE/", "")
            
        basename = os.path.basename(file_path)

        # ── Round-output reference? ─────────────────────────────────
        if _has_round_ref(directive):
            return file_line + _attachment_note(
                directive, file_hint=basename
            )

        # ── Real source file ────────────────────────────────────────
        found = _find_file(ql_dir, file_path)
        if found:
            if verbose:
                print(f"    ✓ {file_path}", file=sys.stderr)
            return file_line + _read(found)

        # ── Not found ───────────────────────────────────────────────
        optional = "if available" in directive.lower()
        tag = "INFO" if optional else "WARNING"
        extra = (
            " Reconstruction from header interface may be needed."
            if optional
            else " Ensure the QuantLib source directory path is correct."
        )
        print(f"    ✗ {file_path} — not found", file=sys.stderr)
        return file_line + f"[{tag}: source file not found at {file_path}.{extra}]"

    return _RE_FILE_INSERT.sub(_repl, content)


# ── Stage 3: Round-8 abbreviated source block ───────────────────────────


def _stage_abbreviated(content: str, ql_dir: str, verbose: bool) -> str:
    """Expand ``[… abbreviated — include … as source files …]``."""

    def _repl(m: re.Match) -> str:
        names = re.findall(r"(\w+\.hpp)", m.group(1))
        parts = ["<source_files>"]
        for nm in names:
            fp = _find_by_name(ql_dir, nm)
            if fp:
                rel = _ql_relpath(fp, ql_dir)
                parts.append(f"\nFILE: {rel}\n{_read(fp)}")
                if verbose:
                    print(f"    ✓ (search) {nm} → {rel}", file=sys.stderr)
            else:
                parts.append(f"\nFILE: {nm}\n[WARNING: {nm} not found]")
                print(f"    ✗ (search) {nm}", file=sys.stderr)
        parts.append("\n</source_files>")
        return "\n".join(parts)

    return _RE_ABBREVIATED.sub(_repl, content)


# ── Stage 4: standalone [INSERT …] (not preceded by FILE:) ──────────────


def _ctx_round(content: str, pos: int) -> int | None:
    """Look backwards from *pos* for a ``ROUND N OUTPUT`` context line."""
    window = content[max(0, pos - 400) : pos]
    hits = re.findall(r"ROUND\s+(\d+)\s+OUTPUT", window, re.IGNORECASE)
    return int(hits[-1]) if hits else None


def _resolve_test_infra(ql_dir: str, verbose: bool) -> str:
    """Try to locate QuantLib test-suite files and return their content."""
    candidates = [
        Path(ql_dir).resolve() / ".." / "test-suite",
        Path(ql_dir).resolve() / "test-suite",
        Path(ql_dir).resolve().parent / "test-suite",
        Path(ql_dir).resolve().parent.parent / "test-suite",
    ]
    td = next((str(d.resolve()) for d in candidates if d.is_dir()), None)
    if not td:
        return (
            "[Test-suite directory not found. Use Boost.Test "
            "BOOST_AUTO_TEST_SUITE pattern or check your QuantLib "
            "source tree layout.]"
        )

    parts: list[str] = []
    # Include key infrastructure headers
    for name in ("utilities.hpp", "quantlibtestsuite.hpp", "speedlevel.hpp"):
        fp = os.path.join(td, name)
        if os.path.isfile(fp):
            parts.append(f"FILE: test-suite/{name}\n{_read(fp)}")
            if verbose:
                print(f"    ✓ test-suite/{name}", file=sys.stderr)

    # Include one representative .cpp test for registration patterns
    try:
        cpp_files = sorted(f for f in os.listdir(td) if f.endswith(".cpp"))
    except OSError:
        cpp_files = []
    for name in cpp_files:
        if "option" in name.lower():
            fp = os.path.join(td, name)
            parts.append(
                f"FILE: test-suite/{name} (registration pattern reference)\n"
                + _read(fp)
            )
            if verbose:
                print(f"    ✓ test-suite/{name}", file=sys.stderr)
            break

    if parts:
        return "\n\n".join(parts)
    return (
        "[No usable test-suite files found. "
        "Use BOOST_AUTO_TEST_SUITE pattern.]"
    )


def _stage_standalone(content: str, ql_dir: str, verbose: bool) -> str:
    """Handle every remaining ``[INSERT …]`` directive."""
    pieces: list[str] = []
    prev_end = 0

    for m in _RE_INSERT.finditer(content):
        directive = m.group(0)
        inner = directive[1:-1]           # strip [ ]
        start, end = m.start(), m.end()
        pieces.append(content[prev_end:start])

        # (a) Preamble — should already be gone, but keep as safety net
        if "SHARED PREAMBLE" in inner:
            pieces.append(directive)

        # (b) <filename>.hpp HEADER  (standalone, e.g. in Round 5)
        elif re.search(r"\S+\.hpp\s+HEADER", inner, re.IGNORECASE):
            ctx = _ctx_round(content, start)
            pieces.append(_attachment_note(directive, ctx_round=ctx))

        # (c) <filename>.hpp from Round N  (standalone, e.g. in Round 6)
        elif re.search(r"from\s+Rounds?\s+\d", inner, re.IGNORECASE):
            pieces.append(_attachment_note(directive))

        # (d) all .hpp files from Rounds N-M  (Round 7)
        elif re.search(r"all\s+\.hpp\s+files", inner, re.IGNORECASE):
            pieces.append(_attachment_note(directive))

        # (e) test-suite / test patterns  (Round 7)
        elif "test-suite" in inner.lower() or "test pattern" in inner.lower():
            pieces.append(_resolve_test_infra(ql_dir, verbose))

        # (f) Generic ROUND N reference
        elif re.search(r"ROUND\s+\d+", inner, re.IGNORECASE):
            ctx = _ctx_round(content, start)
            pieces.append(_attachment_note(directive, ctx_round=ctx))
        
        # (g) Round 7 Golden Reference / Equation Chain inserts
        elif "Golden Reference" in inner or "Equation Chain" in inner:
            pieces.append(f"[>>> ATTACH HERE: {inner.strip()} <<<]")

        # (h) Unrecognised — leave as-is with a warning
        else:
            print(
                f"    ? unhandled INSERT: {directive[:72]}…", file=sys.stderr
            )
            pieces.append(directive)

        prev_end = end

    pieces.append(content[prev_end:])
    return "".join(pieces)


# ══════════════════════ Top-level round assembly ══════════════════════════


def process_round(
    num: int,
    content: str,
    preamble: str,
    ql_dir: str,
    verbose: bool,
) -> str:
    """Run all four processing stages on a single round's template text."""
    print(f"Round {num}:", file=sys.stderr)
    content = _stage_preamble(content, preamble)
    content = _stage_file_inserts(content, ql_dir, verbose)
    content = _stage_abbreviated(content, ql_dir, verbose)
    content = _stage_standalone(content, ql_dir, verbose)
    return content


# ══════════════════════════════ main ══════════════════════════════════════


def main() -> None:
    args = parse_args()

    if not os.path.isdir(args.ql_dir):
        sys.exit(f"ERROR: directory not found: {args.ql_dir}")
    if not os.path.isfile(args.template):
        sys.exit(f"ERROR: template file not found: {args.template}")

    # Quick sanity check — does this look like a QuantLib tree?
    if not (
        _find_file(args.ql_dir, "ql/qldefines.hpp")
        or _find_file(args.ql_dir, "ql/quantlib.hpp")
    ):
        print(
            "WARNING: ql/qldefines.hpp not found under the supplied path.\n"
            "         File resolution may fail.  Make sure you are pointing\n"
            "         at the QuantLib repo root or its ql/ subdirectory.\n",
            file=sys.stderr,
        )

    doc = _read(args.template)
    preamble = _extract_preamble(doc)
    rounds = _extract_rounds(doc)

    if not rounds:
        sys.exit("ERROR: no ROUND sections found in the template.")

    print(f"Parsed rounds: {sorted(rounds)}", file=sys.stderr)
    os.makedirs(args.output_dir, exist_ok=True)

    for num in sorted(rounds):
        info = rounds[num]
        body = process_round(
            num, info["content"], preamble, args.ql_dir, args.verbose
        )

        # Safety check: flag any [INSERT …] that slipped through
        leftover = _RE_INSERT.findall(body)
        if leftover:
            print(
                f"  ⚠  {len(leftover)} unresolved INSERT directive(s):",
                file=sys.stderr,
            )
            for item in leftover[:5]:
                print(f"      {item[:78]}…", file=sys.stderr)

        out_path = os.path.join(args.output_dir, f"round_{num}_prompt.txt")
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write(body)
        print(f"  → {out_path}", file=sys.stderr)

    print(
        f"\nDone — {len(rounds)} prompt(s) written to {args.output_dir}/",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
