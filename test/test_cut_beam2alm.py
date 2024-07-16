import unittest
import os
from pathlib import Path
import numpy as np
import healpy as hp
import grasp2alm as g2a

class TestBeamCut(unittest.TestCase):
    """
    Unit tests for the beam_cut module, specifically testing the functionality
    of converting a beam cut to alm coefficients using the grasp2alm function.

    Methods:
        setUp():
            Initializes the test environment. This includes setting up parameters
            for the beam and writing the beam cut file.

        tearDown():
            Cleans up the test environment by removing the generated beam cut file.

        test_cut_beam2alm():
            Tests the grasp2alm function by comparing the generated alm coefficients
            from the beam cut file to the ideal alm coefficients of a Gaussian beam.
            Asserts that the coefficients match within a specified tolerance.
    """

    def setUp(self):
        """
        Sets up the test environment by initializing parameters.
        This includes defining the path to the beam cut file, setting the polarization flag,
        nside parameter, lmax parameter, and calculating the beam FWHM and sigma.
        It also generates the beam profile and writes it to a file.
        """
        self.path: str = str(Path(__file__).parent / "beam_files" / "beam2alm.cut")

        fwhm_deg: float = 0.5
        beam_sigma: float = np.deg2rad(fwhm_deg)/(2.0*np.sqrt(2.0*np.log(2.0)))
        amplitude: float = 1/(2*np.pi*beam_sigma*beam_sigma)

        gauss = g2a.BeamGauss(amplitude, fwhm_deg)

        pol: bool = True
        nside: int = 1024
        self.lmax: int = 2*nside
        mmax: int = 2

        vini: float = -fwhm_deg*2
        vnum: int = 30001
        ncut: int = 40

        gauss.write2cut(self.path, vini, vnum, ncut)

        self.test_alm = g2a.grasp2alm(
            self.path,
            nside,
            interp_method='pchip',
            copol_axis="y",
            lmax=self.lmax,
            mmax=2,
            pol=pol
        )

        self.ideal_alm = gauss.get_alm(self.lmax, mmax, pol)

    def tearDown(self):
        """
        Cleans up the test environment by removing the generated beam cut file.
        """
        if os.path.exists(self.path):
            os.remove(self.path)

    def test_cut_beam2alm(self):
        """
        Tests the grasp2alm function by comparing the alm coefficients generated
        from the beam cut file to the ideal alm coefficients of a Gaussian beam.
        Asserts that the T, E, and B mode coefficients match within a specified tolerance.

        Asserts:
            True if the alm coefficients generated from the beam cut file
            are close to the ideal alm coefficients within a tolerance of 1e-3.
        """
        index_T = hp.Alm.getidx(self.lmax, np.arange(self.lmax), 0)
        self.assertTrue(np.allclose(
            self.test_alm[0][index_T],
            self.ideal_alm[0][index_T],
            atol=1e-3
            )
        )
        index_E = hp.Alm.getidx(self.lmax, np.arange(self.lmax), 2)
        self.assertTrue(np.allclose(
            self.test_alm[1][index_E],
            self.ideal_alm[1][index_E],
            atol=1e-3
            )
        )
        index_B = hp.Alm.getidx(self.lmax, np.arange(self.lmax), 2)
        self.assertTrue(np.allclose(
            self.test_alm[2][index_B],
            self.ideal_alm[2][index_B],
            atol=1e-3
            )
        )

if __name__ == '__main__':
    unittest.main()
