# Sequential improvement plan

Each stage keeps the original calculation available and adds tests before the
next stage changes the scientific or experimental workflow.

1. **Original-source validation** — completed with GNU Octave 11.3.0: execute the authors' MATLAB
   routines in Octave, archive reference outputs and provenance, and compare
   Python arrays with those references. A MathWorks MATLAB run remains a
   separate validation target.
2. **Structure editor and 3D view** — editable species, factors and coordinates;
   model saving; inspect surface layers and height differences.
3. **Model comparisons and parameter sweeps** — common scales, side-by-side
   models, reproducible parameter grids and sensitivity plots.
4. **Experimental image comparison** — crop, background subtraction, reciprocal
   calibration, extracted profiles and aligned simulation overlays.
5. **Parameter fitting** — bounded fits, residuals, sensitivity and uncertainty;
   validate with synthetic recovery before experimental use.
6. **Optional physics extensions** — distinct from compatibility mode; validate
   independent diffraction orders, corrected sampling, curved Ewald geometry,
   scattering factors and disorder models individually. Multiple scattering is
   a separate research implementation, not part of a quick compatibility patch.
7. **Packaging and usability** — background jobs, cancellation, saved sessions,
   export provenance, automated checks and installation improvements.
