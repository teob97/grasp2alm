import unittest
import os
import numpy as np
import grasp2alm as g2a
from pathlib import Path
from scipy.special import ive

class TestBeamCut(unittest.TestCase):
    def setUp(self):
        self.path = str(Path(__file__).parent / "beam_files" / "beam2alm.cut")
        self.lmax = 512
        pol = True
        nside = 512
        beam_fwhm = np.deg2rad(5)
        beam_sigma = beam_fwhm/(2*np.sqrt(2*np.log(2)))
        amplitude = 1.0

        vini: float = -180.0
        vinc: float = 0.1
        vnum = int(abs(vini)*2/vinc + 1)
        c = 0
        ncut = 40
        header_1: str = "Field data in cuts"
        header_2 = f"{vini} {vinc} {vnum} {c} 3 1 2"

        theta = np.linspace(vini, -vini, vnum)
        theta = np.deg2rad(theta)

        beam_stokes = self.gaussian_beam(amplitude, beam_sigma, theta)
        beam_stokes = np.repeat([beam_stokes], 4, axis=0).T

        co, cx = self.stokes2cocross(beam_stokes)

        self.write2cut(header_1, header_2, vnum, ncut, co, cx)

        self.alm_to_test = g2a.grasp2alm(self.path, nside)
        self.ideal_window = [self.ideal_baeam_window(amplitude, beam_sigma, i) for i in np.array(range(self.lmax))]

    def write2cut(self, header_1, header_2, vnum, ncut, co, cx):
        with open(self.path, 'w') as file:
            for n in range(ncut):
                file.write(header_1)
                file.write('\n')
                file.write(header_2)
                file.write('\n')
                for i in range(vnum):
                    co_i = co[i]
                    cx_i = cx[i]
                    file.write(f"{np.real(co_i)} {np.imag(co_i)} {np.real(cx_i)} {np.imag(cx_i)}\n")

    def gaussian_beam(self, amplitude, sigma, theta):
        return amplitude * np.exp(- theta**2 / (2*sigma**2)) + np.nextafter(0,1)

    def stokes2cocross(self, beam):
        modco2 = 0.5 * (beam[:,0]+beam[:,1])
        modcx2 = 0.5 * (beam[:,0]-beam[:,1])
        co = np.emath.sqrt(modco2)
        cx = np.emath.sqrt(modcx2)
        return co, cx

    def ideal_baeam_window(self, amplitude, sigma, l):
        alpha = 1/sigma**2
        k = np.sqrt(np.pi/(2*alpha))
        return 4*np.pi*amplitude*k*ive(l+0.5, alpha)

    def test_cut_beam2alm_I(self):
        W_I = np.array([self.alm_to_test[0,i]*np.sqrt((4*np.pi)/(2*i+1)) for i in range(self.lmax)])
        assert np.allclose(W_I, self.ideal_window, atol=1e-4)

if __name__ == '__main__':
    unittest.main()