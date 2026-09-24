# RHEEDsim — Python port

This repository is a **Python translation of the original MATLAB RHEEDsim code**
by Kangkang Wang and Arthur R. Smith, CPC catalogue AEJC_v1_0 (2011).
It implements the same algorithms and functionality, including documented
historical numerical conventions, rather than introducing a new diffraction
model. The Python interface and implementation differ from MATLAB; see the
compatibility and verification notes below. This is an independent port, not an
official release or endorsement by the original authors.
Original archive: https://data.mendeley.com/datasets/dk7w4jwg7d/1
Paper: https://doi.org/10.1016/j.cpc.2011.04.023
Original software copyright © 2011 Kangkang Wang, all rights reserved; the archive
is distributed under the CPC licence. This port does not relicense that work.


## Please cite both original sources

If you use this Python port in research, teaching, presentations, or other work,
**please acknowledge the original authors and cite both the paper and the original
MATLAB code archive**. Citation of this GitHub repository alone does not replace
credit to the original work.

1. **Paper:** Kangkang Wang and Arthur R. Smith (2011). “Efficient kinematical
   simulation of reflection high-energy electron diffraction streak patterns for
   crystal surfaces.” *Computer Physics Communications*, **182**(10), 2208–2212.
   DOI: [10.1016/j.cpc.2011.04.023](https://doi.org/10.1016/j.cpc.2011.04.023).
2. **Original MATLAB code:** Kangkang Wang and Arthur R. Smith (2011).
   “Efficient kinematical simulation of reflection high-energy electron diffraction
   streak patterns for crystal surfaces.” *Mendeley Data*, Version 1;
   RHEEDsim, CPC catalogue **AEJC_v1_0**.
   DOI: [10.17632/dk7w4jwg7d.1](https://doi.org/10.17632/dk7w4jwg7d.1).

Suggested acknowledgement: “We used a Python translation of RHEEDsim, originally
developed in MATLAB by Wang and Smith (2011), based on the published paper and
original code archive cited above.” Also record the Python repository URL and
commit hash used for reproducibility. Machine-readable references are provided
in `CITATION.cff`, with BibTeX in `REFERENCES.bib`.

## Start

From this folder:

```sh
python3 -m pip install -r requirements.txt
python3 RHEEDsim.py
```

The desktop GUI requires Tk (tkinter), supplied by many Python distributions.
On macOS, run `sh run.command`. To enable double-click launching after download,
run `chmod +x run.command` first. The launcher
prefers `.venv/bin/python` if an environment has been created there. To create one:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python RHEEDsim.py
```

Use a Python installation that includes Tk for the desktop interface. Batch mode
does not require Tk. Virtual environments and generated outputs are not included
in this repository. `requirements-lock.txt` records the development environment;
`requirements.txt` lists the direct dependencies for installation.

No MATLAB, Octave or SciPy is needed at runtime. The original colormap has been
converted losslessly from the MAT file to `RHEEDcmap.npy`.

1. Click **Load model**, select `examples/sqrt3BL.txt` or `sqrt3HD.txt`.
2. Lattices are displayed automatically; set repeat counts and **Plot lattices**.
3. Set beam angle to 0 degrees; click **Line profile**, then **Streaks**.
4. Click **LEED** for the reciprocal-plane intensity map.
5. Use **Save … data** for ASCII arrays, **… new window** for separate plots.
   Matplotlib toolbars provide zoom, pan and image saving.

Angles are degrees counterclockwise from +x. Distances are Å and reciprocal
coordinates Å⁻¹. Broadening controls are Gaussian standard deviations (sigma),
not FWHM, despite an inaccurate comment in the original MATLAB helper.

## Batch mode

```sh
python3 RHEEDsim.py --input examples/sqrt3HD.txt --angle 0 --leed --output output/HD
python3 RHEEDsim.py --input examples/sqrt3BL.txt --angle 0 --output output/BL
```

Batch mode works without Tk or a display. Options include `--width`, `--radius`,
`--leed-width`, and `--pixels`; see `--help`. It saves peak tables, line profile,
400×500 streak matrix, optional LEED matrix and PNG plots. Existing output files
of the same name are overwritten, so use separate output folders for models.

## Use from Python

```python
from core import Surface, simulate, local_LEEDgen

surface = Surface.load('examples/sqrt3HD.txt')
s, peaks, profile, streaks = simulate(surface, angle=0)
leed_peaks, leed = local_LEEDgen(
    6, 0.2, 300, surface.reciprocal, surface.coordinates)
```

`profile` has columns [reciprocal distance, intensity]. `streaks` has the same
row orientation as MATLAB's exported GUI matrix (after `flipud`). Its plotting
y-axis is linspace(-10, 0, 400). `leed_peaks` has columns [h,k,Gx,Gy,intensity].
Lattice vectors are rows. Scattering coordinates are [form_factor,x,y,z]; input
files instead use 1-based species indices, which `Surface.load` replaces.

## Source mapping

| MATLAB | Python |
| --- | --- |
| RHEEDsim.m and RHEEDsim.fig | RHEEDsim.py: desktop controls, plots, export; core.Surface.load |
| local_recip.m | core.local_recip |
| local_distfinder.m | core.local_distfinder |
| local_IM.m | core.local_IM |
| local_gbroaden.m | core.local_gbroaden |
| local_IM2D.m | core.local_IM2D |
| local_LEEDgen.m | core.local_LEEDgen |
| lattice plotting callback | core.lattices |
| RHEEDcmap.mat | RHEEDcmap.npy |

Global MATLAB variables are replaced by explicit arguments. NumPy vectorization
speeds up intensities and broadening; LEED evaluates only each spot's footprint.

## Compatibility details

The goal is the original program's output, including these conventions:

- Reciprocal search uses the same traversal (first qualifying candidate, not a
  general shortest-vector solver), 100-index limit and absolute 0.001 tolerance.
- Orders extend to MATLAB round(8/spacing); positive half-ties round upward.
- A single atom gives unit line-peak intensities regardless of its form factor.
- Intensities are mirrored at every height, even for structures for which this
  would not be a general physical symmetry. No physics correction is introduced.
- Gaussian cutoff is four sigma. The original overlapping half-array mirror
  assignment is preserved, including its two central sample peculiarities.
  Only even broadening sample counts are supported: odd counts fail in MATLAB.
- Height samples are 0.025 through 10 (400 rows); the GUI flips rows and displays
  them against -10 through 0, preserving its historical coordinate mismatch.
- Exactly zero atomic heights use the original replicated-line shortcut.
- LEED includes the original extra unit-intensity origin entry, in addition to
  the computed origin peak. Its radial broadening is a 1D normal PDF of distance.
- LEED sums indices -50…50, samples coordinates `j*2*r/pixel-r`, but displays
  against linspace(-r,r,pixel), matching the original sampling/plot distinction.
- Provided examples are copied unchanged. The supplied HD file has z=0.2 Å;
  the paper discusses about 0.3 Å. No silent adjustment is made.

GUI layout and fonts are native Python/Matplotlib, not pixel-identical MATLAB
GUIDE controls. Shaded maps use Matplotlib Gouraud interpolation; rendered pixels
can differ from MATLAB while the intensity arrays match the translated algorithm.
Invalid inputs, failed reciprocal searches and cancelled dialogs are handled
explicitly. The original blank-pixel-field typo is corrected to default to 300.
Changing model/parameters recomputes results rather than reusing stale globals.
Line profile calculation also prepares the streak array for convenient use.

## Verification

```sh
python3 -m unittest discover -s tests -v
python3 tests/gui_smoke.py  # optional, requires a graphical display and Tk
```

Tests compare all five examples against separate scalar implementations of the
MATLAB equations, including sampled vertical rows; compare the full small LEED
map against a direct triple loop; and cover reciprocal duality, historical
mirroring, single-atom behavior, invalid models, and failed beam searches.
MATLAB/Octave is not installed here, so these are source-based numerical
compatibility checks, not a claimed execution comparison with MATLAB itself.
