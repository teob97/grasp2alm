import unittest
import os
import numpy as np
import grasp2alm as g2a
import healpy as hp
from pathlib import Path
from scipy.special import ive
from math import factorial

class TestBeamGrid(unittest.TestCase):
    def setUp(self):
        self.path = str( Path(__file__).parent / "beam_files" / "beam2alm.grd" )

        #grid variables
        freq: float = 119.0
        frequnit: str = "GHz"

        ktype: int = 0
        
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
        
        nx: int = 0
        ny: int = 0
        klimit: int = 0

        amp: np.ndarray = None


        header: str = "VERSION: TICRA-EM-FIELD-V0.1\n" + \
                    "Field data in grid\n" + \
                    "SOURCE_FIELD_NAME: baffle_aperture.po\n" + \
                    "FREQUENCY_NAME: freq\n" + \
                    f"FREQUENCIES [{frequnit}]:\n" + \
                    f"{freq}\n" + \
                    "++++" + \
                    f"{ktype}" + \
                    f"{nset} {icomp} {ncomp} {igrid}" + \
                    f"{ix} {iy}" + \
                    f"{xs} {ys} {xe} {ye}" + \
                    f"{nx} {ny} {klimit}"
        
        
        



    
