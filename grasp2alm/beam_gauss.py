import numpy as np

class BeamGauss:
    """Class to sample a gaussian beam.

    Attributes:
        amplitude (float): Amplitude of the beam.
        fwhm_deg (float): Full width at half maximum (FWHM) of the beam.

    Methods:
        __init__(self, ): Initializes a BeamGauss object.
        write2cut(self, path:str, header_1:str, header_2:str, vnum:int, ncut:int, co): Writes the formatted beam data to the specified path with the provided headers. Each cut contains vnum lines of data.
        write2grid(self,path:str,header:str,beam): Writes the formatted beam data to the specified path with the provided headers. Each cut contains vnum lines of data.
    """

    def __init__(self, amplitude:float, fwhm_deg:float):
        """
        Initializes a BeamGauss object.

        Args:
            amplitude (float): Amplitude of the beam.
            fwhm (float): Full width at half maximum (FWHM) of the beam.

        Raises:
            ValueError: If amplitude or fwhm is not a number.
        """
        if not isinstance(amplitude, (int, float)):
            raise ValueError("Amplitude must be a number.")
        if not isinstance(fwhm_deg, (int, float)):
            raise ValueError("FWHM must be a number.")
        self.amplitude = amplitude
        self.fwhm_deg = fwhm_deg
    
    def gaussian_beam(self, theta:float):
        """
        Calculates the value of a Gaussian beam at a given angle.

        Args:
            theta (float): The angle at which to calculate the beam value.

        Returns:
            float: The value of the Gaussian beam at the given angle.
        """
        fwhm_rad = np.deg2rad(self.fwhm_deg)
        sigmasq = fwhm_rad * fwhm_rad / (8 * np.log(2.0))
        return self.amplitude * np.exp(-0.5 * theta**2 / sigmasq)

    def write2cut(self, path:str, vini:int ,vnum:int, ncut:int):
        """
        Writes the formatted beam data to the specified path with the provided headers. 
        Each cut contains vnum lines of data.

        Args:
            path (str): The path to write the beam data to.
            vnum (int): The number of lines in each cut.
            ncut (int): The number of cuts.
            co (array): The beam data.

        Raises:
            ValueError: If the path is not a string.
        """
        vinc: float = abs(vini)*2/vnum
        header_1: str = "Field data in cuts"

        theta = np.linspace(vini, -vini, vnum)
        theta = np.deg2rad(theta)
        phi = np.linspace(0, 180-180/ncut, ncut)
        beam = self.gaussian_beam(theta)

        if not isinstance(path, str):
            raise ValueError("Path must be a string.")
        with open(path, 'w', encoding='utf-8') as file:
            for j in range(ncut):
                file.write(header_1)
                file.write('\n')
                file.write(f"{vini} {vinc} {vnum} {phi[j]} 3 1 2")
                file.write('\n')
                for i in range(vnum):
                    co_i = np.emath.sqrt(beam[i])
                    cx_i = 0.0
                    file.write(f"{np.real(co_i)} {np.imag(co_i)} {np.real(cx_i)} {np.imag(cx_i)}\n")
    
    def write2thetaphigrid(self, path:str, xs:float, ys:float, xe:float, ye:float, nx:int, ny:int):
        """
        Writes the formatted beam data to the specified path with the provided headers. 
        Each cut contains vnum lines of data.

        Args:
            path (str): The path to write the beam data to.
            xs (float): The starting x-coordinate of the grid.
            ys (float): The starting y-coordinate of the grid.
            xe (float): The ending x-coordinate of the grid.
            ye (float): The ending y-coordinate of the grid.
            nx (int): The number of grid points in the x-direction.
            ny (int): The number of grid points in the y-direction.
        """
        header: str = "VERSION: TICRA-EM-FIELD-V0.1\n" + \
            "Field data in grid\n" + \
            "SOURCE_FIELD_NAME: baffle_aperture.po\n" + \
            "FREQUENCY_NAME: freq\n" + \
            "FREQUENCIES [GHz]:\n" + \
            "30\n" + \
            "++++\n" + \
            "1\n" + \
            "1 3 2 7\n" + \
            "0 0\n" + \
            f"{xs} {ys} {xe} {ye}\n" + \
            f"{nx} {ny} 0"
        phi = np.linspace(xs, xe, nx, endpoint=False)
        theta = np.linspace(ys, ye, ny, endpoint=False)
        grid = np.deg2rad(np.meshgrid(phi, theta))
        
        beam = self.gaussian_beam(grid[1])
        
        if not isinstance(path, str):
            raise ValueError("Path must be a string.")
        with open(path, 'w', encoding='utf-8') as file:
            file.write(header)
            file.write("\n")
            co = np.emath.sqrt(beam.flatten())
            for i in co:
                file.write(f"{np.real(i)} {np.imag(i)} 0.0 0.0\n")

        return grid