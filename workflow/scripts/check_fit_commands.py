#!/usr/bin/env python3
"""Check that every path named in samples/*/fits/*/command.txt exists (no refinement is run).

Checks the script, --scan, --space and --instprm paths, and that the CIF of every phase key
(--main-key, --partners) is in the --space folder. Outputs under results/ are skipped.
Run from the repository root: python workflow/scripts/check_fit_commands.py
"""
import shlex
import sys
from pathlib import Path

PATH_FLAGS = {'--scan', '--space', '--instprm'}


def check(cmd_file):
    problems, n = [], 0
    for line in cmd_file.read_text().splitlines():
        if line.startswith('#') or 'python' not in line:
            continue
        tok = shlex.split(line)
        script = next(t for t in tok if t.endswith('.py'))
        paths = [script] + [tok[i + 1] for i, t in enumerate(tok[:-1]) if t in PATH_FLAGS]
        space = tok[tok.index('--space') + 1] if '--space' in tok else None
        keys = []
        for flag in ('--main-key', '--partners'):
            if flag in tok:
                i = tok.index(flag) + 1
                while i < len(tok) and not tok[i].startswith('--'):
                    keys.append(tok[i])
                    i += 1
        paths += [f'{space}/{k}.cif' for k in keys] if space else []
        for p in paths:
            n += 1
            if not Path(p).exists():
                problems.append(f'{cmd_file}: missing {p}')
    return n, problems


def main():
    files = sorted(Path('samples').glob('*/fits/*/command.txt'))
    total, problems = 0, []
    for f in files:
        n, p = check(f)
        total += n
        problems += p
    print('\n'.join(problems))
    print(f'{len(files)} command.txt files, {total} paths checked, {len(problems)} missing')
    sys.exit(bool(problems))


if __name__ == '__main__':
    main()
