"""Regression against arrays exported from unmodified original MATLAB routines."""
import hashlib
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from validation.validate_original import compare


class OriginalReferenceTests(unittest.TestCase):
    def test_original_runtime_arrays(self):
        fixture = ROOT/'validation/reference/octave_original.npz'
        self.assertTrue(fixture.is_file(), 'Missing original-runtime fixture; regenerate with validation/validate_original.py')
        with np.load(fixture, allow_pickle=False) as arrays:
            results = compare(arrays, ROOT/'examples')
        self.assertEqual(len(results), 60)

    def test_example_provenance(self):
        report = json.loads((ROOT/'validation/reference/report.json').read_text())
        fixture = ROOT/'validation/reference/octave_original.npz'
        self.assertEqual(hashlib.sha256(fixture.read_bytes()).hexdigest(), report['fixture_sha256'])
        for path in (ROOT/'examples').glob('*.txt'):
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), report['source_sha256'][path.name])


if __name__ == '__main__':
    unittest.main()
