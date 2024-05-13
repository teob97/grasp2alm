import unittest
import os
import numpy as np
import grasp2alm as g2a
import healpy as hp
from pathlib import Path
from scipy.special import ive
from math import factorial

class TestBeamCut(unittest.TestCase):
    def setUp(self):
        self.path = str(Path(__file__).parent / "beam_files" / "beam2alm.cut")
        self.lmax = 256
        pol = True
        nside = 512
        beam_fwhm = np.deg2rad(10)
        beam_sigma = beam_fwhm/(2*np.sqrt(2*np.log(2)))
        amplitude = 30.0

        vini: float = -180.0
        vinc: float = 0.1
        vnum = int(abs(vini)*2/vinc + 1)
        c = 0
        ncut = 20
        header_1: str = "Field data in cuts"
        header_2 = f"{vini} {vinc} {vnum} {c} 3 1 2"

        theta = np.linspace(vini, -vini, vnum)
        theta = np.deg2rad(theta)

        beam_co = self.gaussian_beam(amplitude, beam_sigma, theta)
        self.write2cut(header_1, header_2, vnum, ncut, beam_co)

        self.test_alm = g2a.grasp2alm(self.path, nside, interp_method='cubic', lmax=self.lmax, mmax=2, pol=pol)
        self.ideal_alm = hp.blm_gauss(beam_fwhm, lmax=self.lmax, pol=pol)

    def tearDown(self):
        if os.path.exists(self.path):
            os.remove(self.path)

    def write2cut(self, header_1, header_2, vnum, ncut, co):
        with open(self.path, 'w') as file:
            for n in range(ncut):
                file.write(header_1)
                file.write('\n')
                file.write(header_2)
                file.write('\n')
                for i in range(vnum):
                    co_i = np.emath.sqrt(co[i])
                    cx_i = 0.0
                    file.write(f"{np.real(co_i)} {np.imag(co_i)} {np.real(cx_i)} {np.imag(cx_i)}\n")

    def gaussian_beam(self, amplitude, sigma, theta):
        return amplitude * np.exp(- theta**2 / (2*sigma**2)) + np.nextafter(0,1)

    def test_cut_beam2alm_I(self):
        return np.allclose(self.test_alm[0], self.ideal_alm[0])

if __name__ == '__main__':
    unittest.main()