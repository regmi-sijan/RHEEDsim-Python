"""Independent scalar reference calculations following the MATLAB source."""
import sys
from pathlib import Path
import unittest
import tempfile
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core import Surface, local_recip, local_distfinder, local_IM, local_IM2D, local_gbroaden, local_LEEDgen, simulate

EXAMPLES = Path(__file__).resolve().parents[1] / 'examples'


def reference_broaden(im, sigma, count=500):
    result = np.column_stack((np.linspace(im[0, 0], im[-1, 0], count), np.zeros(count)))
    for i in range(count//2 - 1, count):
        for position, intensity in im:
            d = abs(position-result[i, 0])
            if d <= 4*sigma:
                result[i, 1] += intensity*np.exp(-d*d/(2*sigma*sigma))/(np.sqrt(2*np.pi)*sigma)
    result[:count//2+1, 1] = result[count//2-1:, 1][::-1].copy()
    return result


class CompatibilityTests(unittest.TestCase):
    def test_reciprocal_duality_and_negative_orientation(self):
        for cell in (np.array([[2., 1.], [0., -3.]]), np.array([[2., 1., 0.], [0., 3., 1.], [1., 0., 4.]])):
            np.testing.assert_allclose(cell @ local_recip(cell).T, 2*np.pi*np.eye(len(cell)), atol=1e-14)

    def test_all_examples_against_scalar_matlab_equations(self):
        for path in EXAMPLES.glob('*.txt'):
            surface = Surface.load(path)
            s, im, fine, streaks = simulate(surface)
            n = int(np.floor(8/s[2] + .5))
            g = surface.reciprocal.T @ s[:2]
            expected = []
            for j in range(n + 1):
                value = sum(f*np.exp(1j*(x*g[0]*j+y*g[1]*j)) for f, x, y, z in surface.coordinates)
                expected.append(1. if len(surface.coordinates) == 1 else abs(value)**2)
            np.testing.assert_allclose(im[:, 1], expected[:0:-1]+expected, atol=1e-12)
            np.testing.assert_allclose(fine, reference_broaden(im, s[2]/5), atol=1e-12)
            for row in (0, 199, 399):
                k = 400-row  # flipud performed by GUI
                positive = []
                for j in range(n+1):
                    value = sum(f*np.exp(1j*(x*g[0]*j+y*g[1]*j+z*k*10/400)) for f,x,y,z in surface.coordinates)
                    positive.append(abs(value)**2)
                row_im = np.column_stack((im[:, 0], positive[:0:-1]+positive))
                expected_row = fine[:, 1] if np.all(surface.coordinates[:, 3] == 0) else reference_broaden(row_im, s[2]/5)[:, 1]
                np.testing.assert_allclose(streaks[row], expected_row, atol=2e-12)
            self.assertEqual(streaks.shape, (400, 500))

    def test_broadening_preserves_overlapping_mirror_assignment(self):
        im = np.array([[-1., 2.], [0., 4.], [1., 9.]])
        np.testing.assert_allclose(local_gbroaden(im, .2, 12), reference_broaden(im, .2, 12), atol=1e-14)
        with self.assertRaises(ValueError):
            local_gbroaden(im, .2, 11)

    def test_leed_against_full_scalar_sum(self):
        cell = np.eye(2)*3
        reciprocal = local_recip(cell)
        atoms = np.array([[1., 0., 0., 0.], [.3, .7, .2, .1]])
        peaks, density = local_LEEDgen(3., .3, 16, reciprocal, atoms)
        expected_peaks = [[0., 0., 0., 0., 1.]]
        for h in range(-50, 51):
            for k in range(-50, 51):
                gx, gy = np.array([h, k]) @ reciprocal
                if gx*gx+gy*gy <= 9:
                    intensity = abs(sum(f*np.exp(1j*(x*gx+y*gy)) for f,x,y,z in atoms))**2
                    expected_peaks.append([h,k,gx,gy,intensity])
        np.testing.assert_allclose(peaks, expected_peaks, atol=1e-14)
        expected = np.zeros((16, 16))
        for l in range(1, 17):
            for j in range(1, 17):
                for h,k,gx,gy,intensity in expected_peaks:
                    d = np.sqrt((j*6/16-3-gx)**2+(l*6/16-3-gy)**2)
                    if d <= 1.2:
                        expected[l-1,j-1] += intensity*np.exp(-.5*(d/.3)**2)/(np.sqrt(2*np.pi)*.3)
        np.testing.assert_allclose(density, expected, atol=1e-14)
        self.assertEqual(np.sum(np.all(peaks[:, 2:4] == 0, axis=1)), 2)

    def test_search_failure_and_single_atom_special_case(self):
        reciprocal = np.eye(2)
        np.testing.assert_array_equal(local_distfinder([[1,0,0],[0,1,0]]), [0,1,1])
        s = local_distfinder([[1,0,13.123456],[0,1,0]])
        np.testing.assert_array_equal(s, [0,0,0])
        with self.assertRaises(ValueError):
            local_IM(s, reciprocal, [[1,0,0,0]])
        im = local_IM(np.array([0,1,1]), reciprocal, np.array([[7,0,0,0]]))
        np.testing.assert_array_equal(im[:, 1], 1)
        with self.assertRaises(ValueError):
            local_recip([[1,1],[2,2]])

    def test_invalid_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'bad.txt'
            path.write_text('Bad\n1 0\n0 1\n1\n2 0 0 0\n')
            with self.assertRaises(ValueError):
                Surface.load(path)


if __name__ == '__main__':
    unittest.main()
