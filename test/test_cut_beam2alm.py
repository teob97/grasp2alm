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

        write2cut(path: str, header_1: str, header_2: str, vnum: int, ncut: int, co):
            Writes the beam data to a cut file in a specified format. Used to
            simulate the input file for the grasp2alm function.

        gaussian_beam(amplitude, sigma, theta):
            Generates a Gaussian beam profile based on the provided amplitude,
            sigma (beam width), and theta (angle) values.

        ideal_alm_gauss(fwhm: float, lmax: int, pol: bool):
            Generates the ideal alm coefficients for a Gaussian beam profile
            given the full width at half maximum (FWHM), maximum multipole moment (lmax),
            and a flag indicating whether polarization is included.

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
        self.pol: bool = True
        self.nside: int = 1024
        self.lmax: int = 2*self.nside
        beam_fwhm_deg: float = np.rad2deg(hp.nside2resol(self.nside))*100
        self.beam_fwhm: float = np.deg2rad(beam_fwhm_deg)
        beam_sigma: float = self.beam_fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))
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

        beam_co = self.gaussian_beam(amplitude, beam_sigma, theta)
        self.write2cut(self.path, header_1, header_2, vnum, ncut, beam_co)

        self.test_alm = g2a.grasp2alm(
            self.path,
            self.nside,
            interp_method='linear',
            lmax=self.lmax,
            mmax=2,
            pol=self.pol
        )
        self.ideal_alm = self.ideal_alm_gauss(self.beam_fwhm, lmax=self.lmax, pol=self.pol)

    def tearDown(self):
        """
        Cleans up the test environment by removing the generated beam cut file.
        """
        if os.path.exists(self.path):
            os.remove(self.path)

    def write2cut(self, path:str, header_1:str, header_2:str, vnum:int, ncut:int, co):
        """
        Writes beam data to a cut file.

        Args:
            path (str): Path to the output cut file.
            header_1 (str): First header line for the cut file.
            header_2 (str): Second header line for the cut file.
            vnum (int): Number of data points in the theta array.
            ncut (int): Number of phi cuts.
            co (array): Copolar component of the beam.

        Writes the formatted beam data to the specified path with the provided headers.
        Each cut contains vnum lines of data.
        """
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

    def gaussian_beam(self, amplitude:float, sigma:float, theta):
        """
        Generates a Gaussian beam profile.

        Args:
            amplitude (float): Amplitude of the Gaussian beam.
            sigma (float): Standard deviation of the Gaussian beam.
            theta (array): Array of theta values (angles) in radians.

        Returns:
            array: Array of complex values representing the Gaussian beam profile.
        """
        return amplitude * np.exp(- theta**2 / (2*sigma**2))

    def ideal_alm_gauss(self, fwhm:float, lmax:int, pol:bool):
        """
        Generates alm coefficients for an ideal Gaussian beam profile.

        Args:
            fwhm (float): Full width at half maximum of the Gaussian beam in radians.
            lmax (int): Maximum multipole moment.
            pol (bool): Flag indicating whether polarization is included.

        Returns:
            array: Array of complex alm coefficients.
        """
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
