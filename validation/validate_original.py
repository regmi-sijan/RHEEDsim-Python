"""Execute original MATLAB source in Octave and compare all exported arrays.

Usage: python validation/validate_original.py --source ../RHEEDsim
Use --check-only to check the committed compressed reference fixture.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from core import Surface, simulate, local_LEEDgen

NAMES = ('GaN', 'coordinates1', 'coordinates2', 'sqrt3BL', 'sqrt3HD')
RTOL, ATOL = 1e-10, 1e-10


def expected_arrays(examples):
    values = {}
    for name in NAMES:
        surface = Surface.load(examples / f'{name}.txt')
        for angle in (0, 90):
            s, peaks, profile, streaks = simulate(surface, angle)
            prefix = f'{name}_angle{angle}'
            for suffix, array in [('reciprocal', surface.reciprocal), ('vector', s), ('peaks', peaks), ('profile', profile), ('streaks', streaks)]:
                values[f'{prefix}_{suffix}'] = array
        peaks, density = local_LEEDgen(3, np.linalg.norm(surface.reciprocal[0])/5, 32, surface.reciprocal, surface.coordinates)
        values[f'{name}_leed_peaks'] = peaks
        values[f'{name}_leed'] = density
    return values


def compare(reference, examples):
    current = expected_arrays(examples)
    if set(reference) != set(current):
        raise AssertionError('Reference fixture arrays do not match the required cases.')
    results = []
    for key, actual in current.items():
        original = reference[key]
        if original.shape != actual.shape:
            raise AssertionError(f'{key}: shape {actual.shape} != {original.shape}')
        np.testing.assert_allclose(actual, original, rtol=RTOL, atol=ATOL, err_msg=key)
        results.append({'array': key, 'shape': list(actual.shape), 'max_abs_error': float(np.max(abs(actual-original)))})
    return results


def quote(path):
    return "'" + str(Path(path).resolve()).replace("'", "''") + "'"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=HERE.parents[1] / 'RHEEDsim')
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    fixture = HERE / 'reference' / 'octave_original.npz'
    report_path = HERE / 'reference' / 'report.json'
    if args.check_only:
        with np.load(fixture, allow_pickle=False) as archive:
            results = compare(archive, HERE.parent/'examples')
        print(f'PASS: {len(results)} arrays match original-source reference outputs.')
        return
    octave = shutil.which('octave-cli') or shutil.which('octave')
    if not octave:
        parser.error('GNU Octave is required for generation. Install it or use --check-only with an existing fixture.')
    if not (args.source / 'local_IM.m').is_file():
        parser.error('Source directory does not contain the original MATLAB routines.')
    source_files = sorted(args.source.glob('*.m')) + [args.source/f'{name}.txt' for name in NAMES]
    for path in source_files:
        if not path.is_file():
            parser.error(f'Missing original input: {path}')
    with tempfile.TemporaryDirectory(prefix='rheedsim-original-') as tmp:
        expression = f"addpath({quote(HERE)}); export_original({quote(args.source)}, {quote(tmp)}, {quote(HERE/'octave_compat')});"
        subprocess.run([octave, '--quiet', '--no-gui', '--eval', expression], check=True)
        reference = {p.stem: np.loadtxt(p, delimiter=',') for p in Path(tmp).glob('*.txt')}
        results = compare(reference, HERE.parent/'examples')
        fixture.parent.mkdir(exist_ok=True)
        np.savez_compressed(fixture, **reference)
    report = {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'runtime': subprocess.check_output([octave, '--version'], text=True).splitlines()[0],
        'source_doi': '10.17632/dk7w4jwg7d.1',
        'python_version': sys.version,
        'numpy_version': np.__version__,
        'fixture_sha256': hashlib.sha256(fixture.read_bytes()).hexdigest(),
        'python_core_sha256': hashlib.sha256((HERE.parent/'core.py').read_bytes()).hexdigest(),
        'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
        'harness_sha256': {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'export_original.m', *sorted((HERE/'octave_compat').glob('*.m'))]},
        'rtol': RTOL, 'atol': ATOL, 'arrays': results,
        'notes': 'Original MATLAB numerical routines executed unmodified in GNU Octave. UI waitbar/close are replaced by no-ops; normpdf uses the standard analytic formula. This is not a MATLAB-runtime validation.'}
    report_path.write_text(json.dumps(report, indent=2)+'\n')
    print(f'PASS: {len(results)} arrays. Reference fixture: {fixture}')


if __name__ == '__main__':
    main()
