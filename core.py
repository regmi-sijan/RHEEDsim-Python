"""Numerical port of Wang and Smith's RHEEDsim (2011).

Functions correspond to the six original local_*.m routines. Arrays use rows
for lattice vectors and [form_factor, x, y, z] for scattering coordinates.
Historical numerical conventions are deliberately retained; see README.md.
"""
from dataclasses import dataclass
from pathlib import Path
import numpy as np


def local_recip(cell):
    cell = np.asarray(cell, dtype=float)
    if cell.shape not in ((2, 2), (3, 3)) or not np.isfinite(cell).all():
        raise ValueError('Cell must be a finite 2x2 or 3x3 matrix.')
    try:
        return 2 * np.pi * np.linalg.inv(cell).T
    except np.linalg.LinAlgError as exc:
        raise ValueError('Lattice vectors must be linearly independent.') from exc


@dataclass
class Surface:
    title: str
    cell: np.ndarray
    species: np.ndarray
    form_factors: np.ndarray
    coordinates: np.ndarray

    @property
    def reciprocal(self):
        return local_recip(self.cell)

    @classmethod
    def load(cls, path):
        lines = Path(path).read_text(encoding='utf-8-sig').splitlines()
        if len(lines) < 5:
            raise ValueError('Expected title, two cell vectors, form factors and atoms.')
        cell = np.array([list(map(float, s.split())) for s in lines[1:3]])
        local_recip(cell)
        if cell.shape != (2, 2):
            raise ValueError('Surface input requires two 2D cell vectors.')
        factors = np.array(list(map(float, lines[3].split())))
        atoms = np.array([list(map(float, s.split())) for s in lines[4:] if s.strip()])
        if atoms.ndim != 2 or atoms.shape[1] != 4 or not len(atoms):
            raise ValueError('Each atom row must contain: species x y z.')
        if not np.isfinite(atoms).all() or not np.isfinite(factors).all():
            raise ValueError('All input numbers must be finite.')
        types = atoms[:, 0].astype(int)
        if np.any(types != atoms[:, 0]) or np.any(types < 1) or np.any(types > len(factors)):
            raise ValueError('Species must be integer indices into the form-factor list (starting at 1).')
        atoms[:, 0] = factors[types - 1]
        return cls(lines[0], cell, types, factors, atoms)


def local_distfinder(x):
    x = np.asarray(x, dtype=float)
    if x.shape != (2, 3) or not np.isfinite(x).all():
        raise ValueError('Input must be finite 2x3: [b1x,b1y,angle]; [b2x,b2y,0].')
    beam = np.array([np.cos(np.deg2rad(x[0, 2])), np.sin(np.deg2rad(x[0, 2]))])
    for q in range(1, 101):
        for h in range(q + 1):
            for k in (q - h, h - q):
                g = np.array([h, k]) @ x[:, :2]
                if abs(beam @ g) <= 0.001:
                    return np.array([h, k, np.linalg.norm(g)])
    return np.zeros(3)


def _orders(s):
    if len(s) != 3 or not np.isfinite(s).all() or s[2] <= 0:
        raise ValueError('No matching reciprocal vector found within the 100-index search.')
    # MATLAB round for positive values (NumPy round uses ties-to-even).
    return int(np.floor(8 / s[2] + 0.5))


def local_IM(s, reciprocal, coordinates):
    n = _orders(s)
    r = np.asarray(coordinates)
    g = np.asarray(reciprocal).T @ s[:2]
    if len(r) == 1:
        positive = np.ones(n + 1)  # Original single-atom special case.
    else:
        positive = abs(r[:, 0] @ np.exp(1j * np.outer(r[:, 1:3] @ g, np.arange(n + 1)))) ** 2
    return np.column_stack((np.arange(-n, n + 1) * s[2], np.r_[positive[:0:-1], positive]))


def _kernel(positions, w, step=500):
    if not np.isfinite(w) or w <= 0:
        raise ValueError('Gaussian sigma must be positive and finite.')
    if not isinstance(step, (int, np.integer)) or step < 2 or step % 2:
        raise ValueError('Use an even sample count >= 2 (the original mirroring requires it).')
    x = np.linspace(positions[0], positions[-1], step)
    d = abs(x[:, None] - positions)
    kernel = np.exp(-0.5 * (d / w)**2) / (np.sqrt(2 * np.pi) * w)
    kernel[d > 4*w] = 0
    # Reproduce the overlapping assignment in local_gbroaden.m exactly.
    kernel[:step//2 + 1] = kernel[step//2 - 1:][::-1].copy()
    return x, kernel


def local_gbroaden(im, w, step=500):
    im = np.asarray(im)
    x, kernel = _kernel(im[:, 0], w, step)
    return np.column_stack((x, kernel @ im[:, 1]))


def local_IM2D(im, s, w, reciprocal, coordinates, progress=None):
    n = _orders(s)
    r = np.asarray(coordinates)
    g = np.asarray(reciprocal).T @ s[:2]
    phase_xy = np.outer(r[:, 1:3] @ g, np.arange(n + 1))
    intensities = np.empty((400, 2*n + 1))
    for row in range(400):
        phase = phase_xy + r[:, 3, None] * (row + 1) * 10 / 400
        positive = abs(r[:, 0] @ np.exp(1j * phase))**2
        intensities[row] = np.r_[positive[:0:-1], positive]
        if progress and (row + 1) % 10 == 0:
            progress((row + 1) / 400)
    _, kernel = _kernel(im[:, 0], w)
    return intensities @ kernel.T


def local_LEEDgen(r, wl, pixel, reciprocal, coordinates, progress=None):
    if not np.isfinite(r) or r <= 0 or not np.isfinite(wl) or wl <= 0:
        raise ValueError('LEED radius and Gaussian sigma must be positive and finite.')
    if not isinstance(pixel, (int, np.integer)) or pixel < 2:
        raise ValueError('Pixel count must be an integer >= 2.')
    atoms = np.asarray(coordinates)
    peaks = [[0., 0., 0., 0., 1.]]  # Original extra origin entry.
    for h in range(-50, 51):
        for k in range(-50, 51):
            g = np.array([h, k]) @ reciprocal
            if g @ g <= r*r:
                intensity = abs(atoms[:, 0] @ np.exp(1j * (atoms[:, 1:3] @ g)))**2
                peaks.append([h, k, *g, intensity])
    peaks = np.array(peaks)
    axis = np.arange(1, pixel + 1) * (2*r/pixel) - r
    density = np.zeros((pixel, pixel))
    for index, (_, _, gx, gy, intensity) in enumerate(peaks):
        # Restrict work to each Gaussian's footprint; same radial normal PDF.
        ix = np.flatnonzero(abs(axis - gx) <= 4*wl)
        iy = np.flatnonzero(abs(axis - gy) <= 4*wl)
        d = np.hypot(axis[iy, None] - gy, axis[None, ix] - gx)
        value = intensity * np.exp(-0.5*(d/wl)**2) / (np.sqrt(2*np.pi)*wl)
        value[d > 4*wl] = 0
        density[np.ix_(iy, ix)] += value
        if progress:
            progress((index + 1)/len(peaks))
    return peaks, density


def lattices(surface, nx=5, ny=5):
    if nx < 0 or ny < 0 or int(nx) != nx or int(ny) != ny:
        raise ValueError('Lattice repeat counts must be nonnegative integers.')
    indices = np.array([(h, k) for h in range(-nx, nx+1) for k in range(-ny, ny+1)])
    real = indices @ surface.cell
    reciprocal = indices @ surface.reciprocal
    basis = np.vstack([np.zeros((1, 2))] + [real + atom[1:3] for atom in surface.coordinates])
    return real, reciprocal, basis


def simulate(surface, angle=0., width=None):
    query = np.column_stack((surface.reciprocal, [angle, 0.]))
    s = local_distfinder(query)
    im = local_IM(s, surface.reciprocal, surface.coordinates)
    width = s[2]/5 if width is None else width
    fine = local_gbroaden(im, width)
    if np.all(surface.coordinates[:, 3] == 0):
        streaks = np.tile(fine[:, 1], (400, 1))
    else:
        streaks = local_IM2D(im, s, width, surface.reciprocal, surface.coordinates)
    return s, im, fine, streaks[::-1].copy()
