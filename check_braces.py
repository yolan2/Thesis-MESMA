#!/usr/bin/env python3
"""
check_braces.py

Simple utility to find unbalanced braces (and other delimiters) in source files (R, .R).
Ignores characters inside single-quoted and double-quoted strings and comments (# ...).

Usage:
  python check_braces.py file1.R file2.R
  python check_braces.py -r path/to/dir   # recurse and check .R files

Exit code 0 if all files balanced, 1 otherwise.
"""

import sys
import argparse
import os
from pathlib import Path

PAIRS = {'{': '}', '(': ')', '[': ']'}
OPEN = set(PAIRS.keys())
CLOSE = {v: k for k, v in PAIRS.items()}


def check_file(path):
    path = Path(path)
    errors = []

    try:
        text = path.read_text(encoding='utf-8')
    except Exception as e:
        errors.append((0, f'ERROR: Could not read file: {e}'))
        return errors

    stack = []  # list of tuples (char, line, col)

    in_single = False
    in_double = False
    in_backtick = False
    escaped = False

    lines = text.splitlines()
    for lineno, line in enumerate(lines, start=1):
        i = 0
        while i < len(line):
            ch = line[i]

            if escaped:
                escaped = False
                i += 1
                continue

            # handle escaping within strings
            if ch == "\\":
                # In R, backslash escapes are meaningful INSIDE strings only.
                # Don't treat a backslash outside a string as escaping the next char.
                if in_single or in_double:
                    escaped = True
                    i += 1
                    continue
                # otherwise treat as normal character
                i += 1
                continue

            # if currently in single-quoted string
            if in_single:
                if ch == "'":
                    in_single = False
                i += 1
                continue

            # if in double-quoted string
            if in_double:
                if ch == '"':
                    in_double = False
                i += 1
                continue

            # if in backtick (e.g., `[`), ignore content until next backtick
            if in_backtick:
                if ch == '`':
                    in_backtick = False
                i += 1
                continue

            # If comment starts, ignore rest of line
            if ch == '#':
                break

            # start of string
            if ch == "'":
                in_single = True
                i += 1
                continue
            if ch == '"':
                in_double = True
                i += 1
                continue
            # backtick starts: treat as quoted identifier `foo` (ignore content)
            if ch == '`':
                in_backtick = True
                i += 1
                continue

            # braces handling
            if ch in OPEN:
                stack.append((ch, lineno, i+1))
            elif ch in CLOSE:
                if stack and stack[-1][0] == CLOSE[ch]:
                    stack.pop()
                else:
                    # unmatched closing
                    # If debug flag is set, print current stack for context and exit quickly
                    if globals().get('_DEBUG_BRACES', False):
                        print('DEBUG: unmatched closing', ch, 'at', path, lineno, 'col', i+1)
                        print('DEBUG: current stack (top last):')
                        for s in stack[-10:]:
                            print('  ', s)
                        # show nearby lines
                        try:
                            txt = path.read_text(encoding='utf-8').splitlines()
                            start = max(0, lineno-6)
                            end = min(len(txt), lineno+3)
                            print('\nContext:')
                            for j in range(start, end):
                                marker = '->' if (j+1) == lineno else '  '
                                print(f"{marker} {j+1:5d}: {txt[j]}")
                        except Exception:
                            pass
                        sys.exit(1)
                    errors.append((lineno, f"Unmatched closing '{ch}' at col {i+1}"))
            i += 1

    # finished lines
    if in_single or in_double:
        errors.append((lineno, 'Unclosed string literal at EOF'))

    if stack:
        for ch, ln, col in stack:
            errors.append((ln, f"Unclosed opening '{ch}' at line {ln} col {col}"))

    return errors


def find_r_files(root):
    root = Path(root)
    for p in root.rglob('*.R'):
        yield p


def main():
    parser = argparse.ArgumentParser(description='Check unbalanced braces in R files')
    parser.add_argument('paths', nargs='*', help='Files or directories to check')
    parser.add_argument('-r', '--recursive', action='store_true', help='Recurse directories')
    args = parser.parse_args()

    targets = []
    if not args.paths:
        print('No files specified. Use -r to recurse a directory or pass file paths.')
        parser.print_help()
        sys.exit(2)

    for p in args.paths:
        ppath = Path(p)
        if ppath.is_file():
            targets.append(ppath)
        elif ppath.is_dir():
            if args.recursive:
                targets.extend(list(find_r_files(ppath)))
            else:
                # add .R files in that directory only
                for q in ppath.iterdir():
                    if q.is_file() and q.suffix.lower() == '.r':
                        targets.append(q)
        else:
            print(f'Warning: {p} not found')

    if not targets:
        print('No files to check')
        sys.exit(2)

    any_err = False
    for t in targets:
        print('Checking', t)
        errs = check_file(t)
        if not errs:
            print('  OK')
        else:
            any_err = True
            for ln, msg in errs:
                print(f'  {t}:{ln}: {msg}')
            # show context lines for first error
            first_ln = errs[0][0]
            try:
                lines = t.read_text(encoding='utf-8').splitlines()
                start = max(0, first_ln-3)
                end = min(len(lines), first_ln+2)
                print('  Context:')
                for i in range(start, end):
                    mark = '->' if (i+1) == first_ln else '  '
                    print(f'    {mark} {i+1:4d}: {lines[i]}')
            except Exception:
                pass
    sys.exit(1 if any_err else 0)

if __name__ == '__main__':
    main()
