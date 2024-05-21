import unittest
import os
from pathlib import Path
import numpy as np
import healpy as hp
import grasp2alm as g2a

class TestBeamGrid(unittest.TestCase):

    def setUp(self):
        
        #path grid file
        self.path = str( Path(__file__).parent / "beam_files" / "beam2alm.grd" )

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

        nx: int = 61
        ny: int = 30000
        klimit: int = 0

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
        pol: bool = True
        nside: int = 1024
        lmax: int = 2*nside
        pol = True
        
        #Index of T, E, B
        index_0 = hp.Alm.getidx(lmax,np.arange(lmax),0) #(lmax,l,m)
        index_2 = hp.Alm.getidx(lmax,np.arange(lmax),2) #(lmax,l,m)
        self.indexes = { 'T':index_0, 'E':index_2, 'B':index_2 }
        
        #Gauss beam variables
        fwhm_deg = np.rad2deg(hp.nside2resol(nside))*100
        fwhm = np.deg2rad(fwhm_deg)
        sigma = fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))
        amplitude = 1/(2*np.pi*sigma**2)
        
        #Compute Gauss beam and write to grid file
        beam_co = self.gaussian_beam(xs,ys,xe,ye,nx,ny,amplitude,sigma)
        self.write2grid(self.path,header,beam_co)

        #Compute Alm
        self.test_alm = g2a.grasp2alm(
            self.path,
            nside,
            interp_method='linear',
            lmax=lmax,
            mmax=2,
            pol=pol
        )
        self.ideal_alm = self.ideal_alm_gauss(fwhm,lmax,pol)

    def tearDown(self):
        if os.path.exists(self.path):
            os.remove(self.path)

    def gaussian_beam(self,xs:float, ys:float, xe:float, ye:float, nx:int, ny:int, amplitude, sigma:float,):
        grid = np.deg2rad( np.meshgrid(np.linspace(xs,xe,nx,endpoint=False),np.linspace(ys,ye,ny,endpoint=False)) )
        gauss = amplitude * np.exp( -0.5*(grid[1]/sigma)**2 )
        return gauss
    

    def write2grid(self,path:str,header:str,beam):
        with open(path, 'w', encoding='utf-8') as file:
            file.write(header)
            file.write("\n")
            co = np.sqrt(beam.reshape(-1))
            cx = 0.
            for i in range(beam.size):
                file.write(f"{np.real(co[i])} {np.imag(co[i])} {np.real(cx)} {np.imag(cx)}\n")
    

    def ideal_alm_gauss(self,fwhm:float,lmax:int,pol:bool):

        mmax:int = 2
        ncomp:int = 3
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
                ) * np.exp(-0.5 * sigmasq * l * (l+1) )
            blm[2] = 1j * blm[1]

        return blm
    

    def test_grid_beam2alm(self):
        self.assertTrue( np.allclose(self.test_alm[0,self.indexes['T']],self.ideal_alm[0,self.indexes['T']],atol=1e-3)  ) 
        self.assertTrue( np.allclose(self.test_alm[1,self.indexes['E']],self.ideal_alm[1,self.indexes['E']],atol=1e-3)  ) 
        self.assertTrue( np.allclose(self.test_alm[2,self.indexes['B']],self.ideal_alm[2,self.indexes['B']],atol=1e-3)  ) 
    

if __name__ == '__main__':
    unittest.main()
