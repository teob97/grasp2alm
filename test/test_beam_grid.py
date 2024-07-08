import unittest
import os
from pathlib import Path
import numpy as np
from grasp2alm import BeamGrid

class TestBeamGrid(unittest.TestCase):
    def setUp(self):
        self.path_to_test_grid = str(Path(__file__).parent / "beam_files" / "unit_test.grd")
    def tearDown(self):
        if os.path.exists(self.path_to_test_grid):
            os.remove(self.path_to_test_grid)
    def write_to_test_grid(self, txt:str):
        text_file = open(self.path_to_test_grid, "w", encoding='utf-8')
        text_file.write(txt)
        text_file.close()

    def test_input_extension_exception(self):
        with self.assertRaises(ValueError):
            BeamGrid("test.not_grid")

    def test_input_grid_format(self):
        txt_with_error = "Test header\n++++\n2"
        self.write_to_test_grid(txt_with_error)
        with self.assertRaises(ValueError):
            BeamGrid(self.path_to_test_grid)

    def test_input_beams_number(self):
        txt_with_error = \
            "Test header\n" + \
            "++++\n" + \
            "1\n" + \
            "2 3 2 7\n" + \
            "0 0\n\n" + \
            "0.0 0.0 360.0 90.0\n" + \
            "2 2 0\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1"
        self.write_to_test_grid(txt_with_error)
        with self.assertWarns(Warning):
            BeamGrid(self.path_to_test_grid)

    def test_input_beam_solid_angle(self):
        txt_with_error = \
            "Test header\n" + \
            "++++\n" + \
            "1\n" + \
            "1 3 2 7\n" + \
            "0 0\n" + \
            "0.0 0.0 340.0 80.0\n" + \
            "2 2 0\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1"        
        self.write_to_test_grid(txt_with_error)
        with self.assertWarns(Warning):
            BeamGrid(self.path_to_test_grid)

    def test_nan_exception(self):
        txt_with_error = \
            "Test header\n" + \
            "++++\n" + \
            "1\n" + \
            "1 3 2 7\n" + \
            "0 0\n" + \
            '0.0 0.0 360.0 90.0\n' + \
            "2 2 0\n" + \
            "1 1 1 1\n" + \
            "1 Nan 1 1\n" + \
            "1 1 1 1\n" + \
            "1 1 1 1"
        self.write_to_test_grid(txt_with_error)
        with self.assertRaises(ValueError):
            BeamGrid(self.path_to_test_grid)

    def test_grid_reading(self):
        txt = \
        "Test header\n" + \
        "VERSION: TICRA-EM-FIELD-V0.1\n" + \
        "FREQUENCY_NAME: freq\n" + \
        "FREQUENCIES [GHz]:\n" + \
        "119.0\n" + \
        "++++\n" + \
        "1\n" + \
        "1 3 2 7\n" + \
        "0 0\n" + \
        '0.0 0.0 360.0 90.0\n' + \
        "3 3 0\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1\n" + \
        "1 1 1 1"

        self.write_to_test_grid(txt)
        test_grid = BeamGrid(self.path_to_test_grid)

        assert test_grid.header == "Test header\n" + \
                                  "VERSION: TICRA-EM-FIELD-V0.1\n" + \
                                  "FREQUENCY_NAME: freq\n"
        assert test_grid.ktype == 1
        assert test_grid.nset == 1
        assert test_grid.klimit == 0
        assert test_grid.icomp == 3
        assert test_grid.ncomp == 2
        assert test_grid.igrid == 7
        assert test_grid.ix == 0
        assert test_grid.iy == 0
        assert test_grid.xs == 0.0
        assert test_grid.ys == 0.0
        assert test_grid.xe == 360.0
        assert test_grid.ye == 90.0
        assert test_grid.nx == 3
        assert test_grid.ny == 3
        assert test_grid.freq == 119.0
        assert test_grid.frequnit == "GHz"

        expected_amp = np.array([
            [[1.+1.j,1.+1.j,1.+1.j],
            [1.+1.j,1.+1.j,1.+1.j],
            [1.+1.j,1.+1.j,1.+1.j]],
  
            [[1.+1.j,1.+1.j,1.+1.j],
            [1.+1.j,1.+1.j,1.+1.j],
            [1.+1.j,1.+1.j,1.+1.j]],
        ])
        self.assertTrue(np.array_equal(test_grid.amp,expected_amp))

if __name__ == "__main__":
    unittest.main()
