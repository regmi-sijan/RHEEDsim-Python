# Original-source validation

`validate_original.py` runs the original, unmodified `local_*.m` functions in GNU
Octave and compares their exported arrays to the Python implementation. It covers
all five supplied models, beam angles 0° and 90°, the complete 400×500 streak maps,
500-point line profiles, reciprocal vectors, peak lists, and 32×32 LEED maps with
radius 3. The smaller LEED grid keeps the original triple-loop code affordable.

## Reproduce the reference

Download and extract the authors' original archive:
https://data.mendeley.com/datasets/dk7w4jwg7d/1

Install GNU Octave, then from the Python repository run:

```sh
python3 validation/validate_original.py --source /path/to/original/RHEEDsim
```

The runner writes a compressed NumPy fixture and JSON report only after all arrays
pass. The report includes the Octave version, input/source SHA-256 hashes,
harness hashes, tolerances, dimensions and maximum absolute error for every array.
The original source is not bundled again or edited by the harness.

## Check without Octave

```sh
python3 validation/validate_original.py --check-only
python3 -m unittest discover -s tests -v
```

The saved arrays are reference data produced by the original code, not another
implementation of the formulas. Tolerances are `rtol=1e-10, atol=1e-10`; they allow
floating-point differences between runtimes and linear-algebra libraries.

## Scope and limitations

This validates execution of MATLAB **source in Octave**, not MathWorks MATLAB
itself. The GUI is not exercised by this numerical harness. Three narrowly scoped
compatibility functions make the original routines runnable headlessly:

- `waitbar`: returns a sentinel instead of drawing progress windows.
- `close`: accepts only that sentinel, replacing progress-window cleanup.
- `normpdf`: evaluates the standard analytic Gaussian PDF without a statistics
  package dependency.

No structure factors, grid coordinates, symmetry operations, peak selection,
Gaussian truncation, or LEED contributions are changed. The MATLAB export harness
can also be run from MATLAB; a separate MATLAB-specific report should be retained
when such a runtime becomes available. These tests establish numerical
compatibility for their stated cases, not experimental or physical accuracy.
