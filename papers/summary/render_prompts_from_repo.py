#!/usr/bin/env python3
"""
generate_prompts.py
===================
Assemble per-round prompts for the QuantLib FDM improvement prompting
framework.

Reads the framework template document, resolves [INSERT ...] placeholders
against a local QuantLib v1.23 source tree, and writes one ready-to-send
prompt file per round into an output folder.

For each round the script:
  1. Concatenates the shared preamble with the round-specific template.
  2. Locates every referenced QuantLib source file on disk and inlines
     its full content.
  3. Replaces references to prior-round generated outputs with clear
     "please attach" notes.
  4. Saves the assembled prompt into <output-dir>/round_<N>_prompt.txt.

Usage:
    python generate_prompts.py <ql_source_dir>
    python generate_prompts.py ../ql/ -f prompting_framework.txt -o prompts
"""

import argparse
import os
import re
import sys


# ======================================================================
#  1.  FRAMEWORK DOCUMENT PARSING
# ======================================================================

def parse_framework(text):
    """
    Split the framework markdown into:
        preamble  – the shared system-context / math-context XML string
        rounds    – dict  {round_number: template_string}

    Each section is delimited by a ``### `` header; its payload is the
    content inside the *first* ```xml ... ``` fenced code block.
    """
    # Split at lines that begin a new ### section (lookahead keeps the
    # header attached to the block that follows).
    parts = re.split(r'\n(?=### )', text)

    preamble = None
    rounds = {}

    for part in parts:
        block = re.search(r'```xml\s*\n(.*?)```', part, re.DOTALL)
        if not block:
            continue

        header = part.lstrip('\n').split('\n', 1)[0]
        code = block.group(1)

        if 'SHARED PREAMBLE' in header:
            preamble = code
        else:
            m = re.search(r'ROUND\s+(\d+)', header)
            if m:
                rounds[int(m.group(1))] = code

    return preamble, rounds


# ======================================================================
#  2.  SOURCE-FILE RESOLUTION
# ======================================================================

def _try_read(path):
    """Return file contents as a string, or None on any failure."""
    try:
        if os.path.isfile(path):
            with open(path, encoding='utf-8', errors='replace') as fh:
                return fh.read()
    except OSError:
        pass
    return None


def resolve_file(ql_dir, rel_path):
    """
    Locate *rel_path* (e.g. ``ql/methods/.../foo.hpp``) on disk.

    The caller's *ql_dir* may point either to the QuantLib project root
    (which contains the ``ql/`` sub-tree) **or** directly to the ``ql/``
    directory itself.  We try both conventions plus an up-one-level probe.

    Returns the file content string or ``None``.
    """
    norm = os.path.normpath(ql_dir)
    candidates = [os.path.join(norm, rel_path)]

    # If rel_path starts with 'ql/', try stripping that prefix
    # (for when ql_dir already IS the ql/ directory).
    if rel_path.startswith('ql/') or rel_path.startswith('ql\\'):
        candidates.append(os.path.join(norm, rel_path[3:]))

    # Also probe one directory above ql_dir.
    parent = os.path.dirname(norm)
    if parent != norm:
        candidates.append(os.path.join(parent, rel_path))

    for c in candidates:
        content = _try_read(c)
        if content is not None:
            return content
    return None


def find_test_files(ql_dir):
    """
    Locate representative ``test-suite/`` files for Round 7.

    Returns a list of ``(relative_path, content)`` pairs, possibly empty.
    """
    norm = os.path.normpath(ql_dir)
    search_dirs = [
        os.path.join(norm, 'test-suite'),
        os.path.normpath(os.path.join(norm, '..', 'test-suite')),
        os.path.join(os.path.dirname(norm), 'test-suite'),
    ]
    tdir = None
    for d in search_dirs:
        if os.path.isdir(d):
            tdir = d
            break
    if tdir is None:
        return []

    found = []

    # Test-utility header
    for name in ('utilities.hpp', 'utilities.h'):
        c = _try_read(os.path.join(tdir, name))
        if c is not None:
            found.append(('test-suite/' + name, c))
            break

    # A representative test .cpp (prefer FDM-related)
    for name in ('fdmlinearop.cpp', 'europeanoption.cpp',
                 'americanoption.cpp'):
        c = _try_read(os.path.join(tdir, name))
        if c is not None:
            found.append(('test-suite/' + name, c))
            break

    return found


# ======================================================================
#  3.  [INSERT ...] PLACEHOLDER HELPERS
# ======================================================================

# We anchor to the start of a line (re.MULTILINE) so we never accidentally
# match the *example* mention of "[INSERT FULL CONTENT ...]" that appears
# inline inside the preamble's prose.
_INSERT_RE = re.compile(r'^\[INSERT[^\]]*\]', re.MULTILINE)


def _classify(tag):
    """
    Classify an ``[INSERT ...]`` tag.

    Returns one of ``'prior_round'``, ``'test_infra'``, or ``'source'``.
    """
    # ── prior-round output references ──
    if re.search(r'ROUND\s+\d', tag):
        return 'prior_round'
    if re.search(r'from\s+Round', tag, re.IGNORECASE):
        return 'prior_round'
    if re.search(r'Rounds\s+\d', tag):
        return 'prior_round'
    # [INSERT filename.hpp HEADER]  (but NOT "[INSERT FULL CONTENT …]")
    if (re.search(r'\w+\.hpp\s+HEADER\s*\]', tag)
            and 'FULL CONTENT' not in tag):
        return 'prior_round'

    # ── test infrastructure ──
    if 'test-suite' in tag.lower():
        return 'test_infra'

    # ── default: source file to be found on disk ──
    return 'source'


def _prior_round_note(tag):
    """
    Build a human-readable note telling the user to attach the relevant
    output from a previous round.
    """
    nums = re.findall(r'(?:Round|ROUND)\s+(\d+(?:[/\-]\d+)?)', tag)
    fm = re.search(r'(\w+\.(?:hpp|cpp))', tag)
    name = fm.group(1) if fm else None

    if 'all .hpp files' in tag.lower():
        return ('[Attach ALL .hpp output files generated from '
                'Rounds 1 through 6 here]')
    if name and nums:
        return ('[Attach the generated {} output from Round {} here]'
                .format(name, nums[0]))
    if nums:
        return ('[Attach the relevant output from Round {} here]'
                .format(nums[0]))
    if name:
        return ('[Attach the generated {} from its prior round here]'
                .format(name))
    return '[Attach the relevant prior-round output here]'


def _preceding_file_path(text, pos):
    """
    Return the ``FILE:`` path from the line immediately above *pos*,
    provided only whitespace separates the two.  Otherwise return ``None``.
    """
    before = text[:pos]
    hits = list(re.finditer(r'^FILE:\s*(.+)$', before, re.MULTILINE))
    if not hits:
        return None
    last = hits[-1]
    # Make sure nothing but whitespace sits between the FILE: line and
    # the INSERT tag.
    if before[last.end():].strip():
        return None
    return last.group(1).strip()


# ======================================================================
#  4.  TEMPLATE ASSEMBLY
# ======================================================================

def _handle_modified_preamble(template, preamble, match):
    """
    Round 8 special case.

    The INSERT reads:
        [INSERT SHARED PREAMBLE, replacing the line about
         not using legacy framework with:]
    and the next non-blank lines (up to the first blank line or XML tag)
    are the *replacement text* to splice into the preamble.

    Returns the updated template string.
    """
    start, end = match.start(), match.end()

    # ── gather replacement-text lines that follow the INSERT tag ──
    tail = template[end:]
    skip = 1 if tail.startswith('\n') else 0
    tail = tail[skip:]

    repl_lines = []
    eaten = 0
    for ln in tail.split('\n'):
        s = ln.strip()
        if s == '' or s.startswith('<'):
            break
        repl_lines.append(ln)
        eaten += len(ln) + 1           # +1 for the '\n' removed by split

    repl_text = '\n'.join(repl_lines)
    cut_end = end + skip + eaten        # first char we keep after the splice

    # ── patch the preamble ──
    mod = preamble.replace(
        '(not the deprecated legacy MixedScheme framework)', '')
    mod = re.sub(r'  +', ' ', mod)      # tidy up any resulting double-spaces

    marker = 'Fdm* framework.'
    idx = mod.find(marker)
    if idx >= 0:
        p = idx + len(marker)
        mod = mod[:p] + '\n\n' + repl_text + '\n' + mod[p:]
    else:
        # Fallback: simply prepend the replacement text.
        mod = repl_text + '\n\n' + mod

    return template[:start] + mod + template[cut_end:]


def assemble_prompt(template, preamble, ql_dir):
    """
    Build a complete, ready-to-send prompt from a round *template* by
    resolving every ``[INSERT ...]`` placeholder.

    Returns ``(prompt_string, [warnings])``.
    """
    warnings = []

    # ── Phase 1: preamble insertion ─────────────────────────────
    if '[INSERT SHARED PREAMBLE HERE]' in template:
        template = template.replace('[INSERT SHARED PREAMBLE HERE]',
                                    preamble)
    else:
        # Round 8 variant: preamble with inline modification
        mod_m = re.search(
            r'^\[INSERT SHARED PREAMBLE[^\]]*replacing[^\]]*\]',
            template, re.MULTILINE)
        if mod_m:
            template = _handle_modified_preamble(template, preamble, mod_m)

    # ── Phase 2: remaining [INSERT …] placeholders ──────────────
    # Process from *end* to *start* so that string-position arithmetic
    # stays valid after each replacement.
    for m in reversed(list(_INSERT_RE.finditer(template))):
        tag = m.group(0)
        start, end = m.start(), m.end()
        kind = _classify(tag)

        if kind == 'prior_round':
            repl = _prior_round_note(tag)

        elif kind == 'test_infra':
            tfiles = find_test_files(ql_dir)
            if tfiles:
                repl = '\n\n'.join(
                    'FILE: {}\n{}'.format(rp, ct) for rp, ct in tfiles)
            else:
                repl = ('[Test infrastructure files not found in source '
                        'tree. Please provide test-suite/utilities.hpp '
                        'and an example test .cpp file manually.]')
                warnings.append('test-suite/ directory not found')

        else:   # source file
            fpath = _preceding_file_path(template, start)
            if fpath:
                content = resolve_file(ql_dir, fpath)
                if content is not None:
                    repl = content
                else:
                    repl = '[SOURCE FILE NOT FOUND: {}]'.format(fpath)
                    warnings.append('File not found: {}'.format(fpath))
            else:
                # Cannot determine which file to inline — keep tag as-is.
                repl = tag
                warnings.append('No FILE: path for: {}'.format(
                    tag[:80] + ('...' if len(tag) > 80 else '')))

        template = template[:start] + repl + template[end:]

    return template, warnings


# ======================================================================
#  5.  ENTRY POINT
# ======================================================================

def main():
    ap = argparse.ArgumentParser(
        description=(
            'Assemble per-round prompts for the QuantLib FDM improvement '
            'prompting framework by inlining referenced source files.'),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples:\n'
            '  python generate_prompts.py ../ql/\n'
            '  python generate_prompts.py /src/QuantLib-1.23/ -o my_prompts\n'
            '  python generate_prompts.py ../ql/ -f prompting_framework.txt\n'))

    ap.add_argument(
        'ql_dir',
        help='Path to the QuantLib source tree (project root or ql/ subdir)')
    ap.add_argument(
        '-f', '--framework',
        default='prompting_framework.txt',
        help='Path to the framework document  [default: %(default)s]')
    ap.add_argument(
        '-o', '--output-dir',
        default='prompts',
        help='Output folder for generated prompt files  [default: %(default)s]')

    args = ap.parse_args()

    # ── validate inputs ──────────────────────────────────────────
    if not os.path.isfile(args.framework):
        sys.exit('Error: framework file not found: ' + args.framework)
    if not os.path.isdir(args.ql_dir):
        sys.exit('Error: directory not found: ' + args.ql_dir)

    # ── parse framework document ─────────────────────────────────
    with open(args.framework, encoding='utf-8') as fh:
        fw_text = fh.read()

    preamble, rounds = parse_framework(fw_text)
    if preamble is None:
        sys.exit('Error: shared preamble not found in the framework document.')
    if not rounds:
        sys.exit('Error: no round templates found in the framework document.')

    print('Parsed framework document')
    print('  Preamble : {:,} chars'.format(len(preamble)))
    print('  Rounds   : {}\n'.format(sorted(rounds.keys())))

    # ── generate one prompt file per round ───────────────────────
    os.makedirs(args.output_dir, exist_ok=True)
    all_warnings = {}

    for rn in sorted(rounds.keys()):
        prompt, ws = assemble_prompt(rounds[rn], preamble, args.ql_dir)

        out_path = os.path.join(
            args.output_dir, 'round_{}_prompt.txt'.format(rn))
        with open(out_path, 'w', encoding='utf-8') as fh:
            fh.write(prompt)

        size_kb = len(prompt.encode('utf-8')) / 1024
        flag = '  [!] {} warning(s)'.format(len(ws)) if ws else ''
        print('  Round {}: {}  ({:,.1f} KB){}'.format(
            rn, out_path, size_kb, flag))

        if ws:
            all_warnings[rn] = ws

    # ── summary ──────────────────────────────────────────────────
    print()
    if all_warnings:
        print('--- Warnings ---')
        for rn in sorted(all_warnings.keys()):
            for w in all_warnings[rn]:
                print('  Round {}: {}'.format(rn, w))
        print('\nSome placeholders could not be resolved (see above).')
        print('The output files contain [SOURCE FILE NOT FOUND: ...] '
              'markers at those locations.')
    else:
        print('All placeholders resolved successfully.')


if __name__ == '__main__':
    main()
