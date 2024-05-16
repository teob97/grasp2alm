import unittest
import os
from pathlib import Path
import numpy as np
import healpy as hp
import grasp2alm as g2a

class TestBeamCut(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.path: str = str(Path(__file__).parent / "beam_files" / "beam2alm.cut")
        cls.pol: bool = True
        cls.nside: int = 1024
        cls.lmax: int = 2*cls.nside
        beam_fwhm_deg: float = np.rad2deg(hp.nside2resol(cls.nside))*100
        cls.beam_fwhm: float = np.deg2rad(beam_fwhm_deg)
        beam_sigma: float = cls.beam_fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))
        amplitude: float = 1/(2*np.pi*beam_sigma*beam_sigma)

        vini: float = -beam_fwhm_deg*3
        vnum: int = 30001
        vinc: float = abs(vini)*2/vnum
        c: int = 0
        ncut: int = 40
        header_1: str = "Field data in cuts"
        header_2: str = f"{vini} {vinc} {vnum} {c} 3 1 2"

        theta = np.linspace(vini, -vini, vnum)
        theta = np.deg2rad(theta)

        beam_co = cls.gaussian_beam(amplitude, beam_sigma, theta)
        cls.write2cut(cls.path, header_1, header_2, vnum, ncut, beam_co)

    @classmethod
    def write2cut(cls, path:str, header_1:str, header_2:str, vnum:int, ncut:int, co):
        with open(path, 'w', encoding='utf-8') as file:
            for _ in range(ncut):
                file.write(header_1)
                file.write('\n')
                file.write(header_2)
                file.write('\n')
                for i in range(vnum):
                    co_i = np.emath.sqrt(co[i])
                    cx_i = 0.0
                    file.write(f"{np.real(co_i)} {np.imag(co_i)} {np.real(cx_i)} {np.imag(cx_i)}\n")

    @classmethod
    def gaussian_beam(cls, amplitude, sigma, theta):
        return amplitude * np.exp(- theta**2 / (2*sigma**2))

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.path):
            os.remove(cls.path)

    def setUp(self):
        self.test_alm = g2a.grasp2alm(
            self.path,
            self.nside,
            interp_method='linear',
            lmax=self.lmax,
            mmax=2,
            pol=self.pol
        )
        self.ideal_alm = self.ideal_alm_gauss(self.beam_fwhm, lmax=self.lmax, pol=self.pol)

    def ideal_alm_gauss(self, fwhm:float, lmax:int, pol:bool):
        mmax: int = 2
        ncomp: int = 3
        nval = hp.Alm.getsize(lmax, mmax)

        blm = np.zeros((ncomp, nval), dtype=np.complex128)
        sigmasq = fwhm * fwhm / (8 * np.log(2.0))

        for l in range(0, lmax + 1):
            blm[0, hp.Alm.getidx(lmax, l, 0)] = np.sqrt((2 * l + 1) / (4.0 * np.pi)) * np.exp(
                -0.5 * sigmasq * l * (l+1)
            )

        if pol:
            for l in range(2, lmax + 1):
                blm[1, hp.Alm.getidx(lmax, l, 2)] = np.sqrt(
                    (2 * l + 1) / (16 * np.pi)
                ) * np.exp(-0.5 * sigmasq * l * (l+1))
            blm[2] = 1j * blm[1]

        return blm

    def test_cut_beam2alm_I(self):
        index = hp.Alm.getidx(self.lmax, np.arange(self.lmax), 0)
        self.assertTrue(np.allclose(self.test_alm[0][index], self.ideal_alm[0][index], atol=1e-3))
    def test_cut_beam2alm_Q(self):
        index = hp.Alm.getidx(self.lmax, np.arange(self.lmax), 2)
        self.assertTrue(np.allclose(self.test_alm[1][index], self.ideal_alm[1][index], atol=1e-3))
    def test_cut_beam2alm_U(self):
        index = hp.Alm.getidx(self.lmax, np.arange(self.lmax), 2)
        self.assertTrue(np.allclose(self.test_alm[2][index], self.ideal_alm[2][index], atol=1e-3))

if __name__ == '__main__':
    unittest.main()
