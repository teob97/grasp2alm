import unittest
import os
from pathlib import Path
import numpy as np
import healpy as hp
import grasp2alm as g2a

class TestBeamGrid(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        
        #path grid file
        cls.path = str( Path(__file__).parent / "beam_files" / "beam2alm.grd" )

        #grid variables
        freq: float = 119.0
        frequnit: str = "GHz"

        ktype: int = 1

        nset: int = 1
        icomp: int = 3
        ncomp: int = 2
        igrid: int = 7

        ix: int = 0
        iy: int = 0

        xs: float = 0.0
        ys: float = 0.0
        xe: float = 360.0
        ye: float = 90.0

        nx: int = 1001
        ny: int = 1001
        klimit: int = 0

        amp: np.ndarray = None


        header: str = "VERSION: TICRA-EM-FIELD-V0.1\n" + \
                    "Field data in grid\n" + \
                    "SOURCE_FIELD_NAME: baffle_aperture.po\n" + \
                    "FREQUENCY_NAME: freq\n" + \
                    f"FREQUENCIES [{frequnit}]:\n" + \
                    f"{freq}\n" + \
                    "++++\n" + \
                    f"{ktype}\n" + \
                    f"{nset} {icomp} {ncomp} {igrid}\n" + \
                    f"{ix} {iy}\n" + \
                    f"{xs} {ys} {xe} {ye}\n" + \
                    f"{nx} {ny} {klimit}"
        
        #map and alm variables
        cls.pol: bool = True
        cls.nside: int = 1024
        cls.lmax: int = 2*cls.nside
        
        #Gauss beam
        fwhm_deg = np.rad2deg(hp.nside2resol(cls.nside))*100
        fwhm = np.deg2rad(fwhm_deg)
        sigma = fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))
        amplitude = 1/(2*np.pi*sigma**2)
        
        beam_co = cls.gaussian_beam(xs,ys,xe,ye,nx,ny,amplitude,sigma)
        cls.write2grid(cls.path,header,beam_co)
    
    @classmethod
    def gaussian_beam(cls,xs:float, ys:float, xe:float, ye:float, nx:int, ny:int, amplitude, sigma:float,):
        grid = np.deg2rad( np.meshgrid(np.linspace(xs,xe,nx,endpoint=False),np.linspace(ys,ye,ny,endpoint=False)) )
        gauss = amplitude * np.exp( -0.5*(grid[1]/sigma)**2 )
        return gauss
    
    @classmethod
    def write2grid(cls,path:str,header:str,beam):
        with open(path, 'w', encoding='utf-8') as file:
            file.write(header)
            file.write("\n")
            co = np.sqrt(beam.reshape(-1))
            cx = 0.
            for i in range(beam.size):
                file.write(f"{np.real(co[i])} {np.imag(co[i])} {np.real(cx)} {np.imag(cx)}\n")
    
    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.path):
            os.remove(cls.path)

if __name__ == '__main__':
    unittest.main()
