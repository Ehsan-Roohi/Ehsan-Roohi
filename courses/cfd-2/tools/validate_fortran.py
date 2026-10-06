"""Compile every archived Fortran variant and run cases in isolated directories.

Original files are never edited. Compatibility edits are recorded separately.
A successful process exit is not a physical validation of a CFD solution.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
FLAGS = ['-std=legacy', '-fdec', '-ffixed-line-length-none',
         '-ffree-line-length-none', '-fallow-argument-mismatch', '-O0', '-g',
         '-fcheck=all', '-ffpe-trap=invalid,zero,overflow', '-fbacktrace',
         '-finit-real=snan', '-finit-integer=-999999']


def sanitize(text, compatibility=False):
    lines = text.splitlines()
    edits = []
    for i, line in enumerate(lines):
        if re.match(r'(?i)^\s*CALL\s+SYSTEM\s*\(\s*[\"\x27]del\s', line):
            lines[i] = '! Removed Windows cleanup from working copy.'
            edits.append({'line': i + 1, 'kind': 'remove_cleanup'})
        elif compatibility and re.match(r'(?i)^\s*USE\s+(MSIMSL|MSFLIB)\s*$', line):
            lines[i] = '! Removed unused legacy module import: ' + line.strip()
            edits.append({'line': i + 1, 'kind': 'remove_legacy_import'})
    if compatibility:
        implicit = next((i for i, line in enumerate(lines)
                         if re.match(r'(?i)^\s*implicit\b', line)), None)
        first_parameter = next((i for i, line in enumerate(lines)
                                if re.match(r'(?i)^\s*parameter\s*\(', line)), None)
        if implicit is not None and first_parameter is not None and implicit > first_parameter:
            statement = lines[implicit]
            lines[implicit] = '! IMPLICIT moved before PARAMETER for GNU Fortran.'
            lines.insert(first_parameter, statement)
            edits.append({'line': implicit + 1, 'kind': 'reorder_implicit',
                          'before_original_line': first_parameter + 1})
        # The internal VISCOUS routine assigns host I, an active caller DO index.
        # Give that routine its own integer I; preserve the caller's loop state.
        for i, line in enumerate(lines):
            if re.match(r'(?i)^\s*subroutine\s+VISCOUS\s*\(', line):
                lines.insert(i + 1, '      INTEGER :: I')
                edits.append({'kind': 'localize_viscous_loop_index',
                              'subroutine': 'VISCOUS', 'variable': 'I'})
                break
    return '\n'.join(lines) + '\n', edits


def run_process(command, directory, timeout, env, input_text=None):
    start = time.monotonic()
    try:
        p = subprocess.run(command, cwd=directory, env=env, input=input_text,
                           capture_output=True, text=True, errors='replace', timeout=timeout)
        return {'exit_code': p.returncode, 'timed_out': False,
                'seconds': round(time.monotonic() - start, 3),
                'stdout': p.stdout, 'stderr': p.stderr}
    except subprocess.TimeoutExpired as e:
        def decode(value):
            return value.decode(errors='replace') if isinstance(value, bytes) else value or ''
        return {'exit_code': None, 'timed_out': True,
                'seconds': round(time.monotonic() - start, 3),
                'stdout': decode(e.stdout), 'stderr': decode(e.stderr)}
    except OSError as e:
        return {'exit_code': None, 'timed_out': False, 'environment_error': str(e),
                'seconds': round(time.monotonic() - start, 3), 'stdout': '', 'stderr': str(e)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--compiler', required=True)
    p.add_argument('--work', default='work/validation-fortran')
    p.add_argument('--output', default='docs/validation/fortran-validation.json')
    p.add_argument('--timeout', type=float, default=30)
    p.add_argument('--filter', default='', help='Substring of archive-relative source path')
    a = p.parse_args()
    compiler = Path(a.compiler).resolve()
    work = (ROOT / a.work).resolve()
    if work.exists():
        p.error('Use a new work directory; existing outputs are not overwritten')
    work.mkdir(parents=True)
    env = dict(os.environ)
    env['PATH'] = str(compiler.parent) + os.pathsep + env.get('PATH', '')
    version = subprocess.check_output([str(compiler), '--version'], env=env, text=True).splitlines()[0]
    result = {'compiler': version, 'flags': FLAGS, 'scope': 'All selected archived sources',
              'run_timeout_seconds': a.timeout, 'sources': []}
    sources = sorted(f for f in (ROOT / 'projects').rglob('*')
                     if f.suffix.lower() in ('.for', '.f90') and a.filter in f.as_posix())
    for n, source in enumerate(sources, 1):
        rel = source.relative_to(ROOT).as_posix()
        directory = work / f'{n:02d}'
        directory.mkdir()
        raw = source.read_bytes()
        text = raw.decode('latin1')
        entry = {'source': rel, 'sha256': hashlib.sha256(raw).hexdigest(),
                 'format': 'free' if source.suffix.lower() == '.f90' else 'fixed'}
        suffix = '.f90' if entry['format'] == 'free' else '.for'
        working_source = directory / ('source' + suffix)
        exe = directory / ('run.exe' if os.name == 'nt' else 'run')
        for compatible in (False, True):
            prepared, edits = sanitize(text, compatibility=compatible)
            working_source.write_text(prepared, encoding='latin1')
            flags = FLAGS + ['-ffree-form' if entry['format'] == 'free' else '-ffixed-form']
            build = run_process([str(compiler), *flags, working_source.name, '-o', exe.name],
                                directory, 60, env)
            key = 'compatibility_build' if compatible else 'cleanup_only_build'
            entry[key] = {'edits': edits, **build}
            if build['exit_code'] == 0:
                entry['active_build'] = key
                break
        if not entry.get('active_build'):
            entry['status'] = 'compile_failed'
        else:
            required = sorted(set(re.findall(r"(?i)\bfile\s*=\s*['\"]([^'\"]+)['\"]", text)))
            # Only co-located archival input files are staged; never borrow an unrelated mesh.
            staged = []
            for f in source.parent.iterdir():
                if f.is_file() and f.suffix.lower() not in ('.for', '.f90', '.md'):
                    shutil.copyfile(f, directory / f.name)
                    staged.append(f.name)
            aliases = []
            for name in required:
                match = next((n for n in staged if n.casefold() == name.casefold()), None)
                if match and match != name and not (directory / name).exists():
                    shutil.copyfile(directory / match, directory / name)
                    aliases.append({'original': match, 'alias': name})
                    staged.append(name)
            entry['input_filename_aliases'] = aliases
            entry['staged_files'] = staged
            has_read = bool(re.search(r'(?i)\bread\s*\(', text))
            if '/variants/' in rel and has_read:
                entry['status'] = 'build_passed_inputs_not_co_located'
                entry['referenced_filenames'] = required
            elif source.name.lower() == 'contour.for':
                entry['status'] = 'build_passed_postprocessor_requires_flow_outputs'
            else:
                # Second0 asks for Roe flux; LAST GRID asks for SOR relaxation.
                stdin = '1.0\n' if source.name.upper() == 'LAST GRID.FOR' else '1\n'
                entry['stdin'] = stdin.strip()
                entry['run'] = run_process([str(exe)], directory, a.timeout, env, stdin)
                entry['status'] = ('execution_environment_error' if entry['run'].get('environment_error') else
                                   'timeout' if entry['run']['timed_out'] else
                                   'completed' if entry['run']['exit_code'] == 0 else 'runtime_failed')
                entry['outputs'] = [{'file': f.name, 'bytes': f.stat().st_size}
                                    for f in sorted(directory.iterdir()) if f.is_file()
                                    and f.name not in staged + [working_source.name, exe.name]]
        result['sources'].append(entry)
        out = (ROOT / a.output).resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf8')
        print(f'{n}/{len(sources)} {rel}: {entry["status"]}', flush=True)
    print(json.dumps({'statuses': {s: sum(x['status'] == s for x in result['sources'])
                                  for s in sorted({x['status'] for x in result['sources']})}}, indent=2))


if __name__ == '__main__':
    main()
