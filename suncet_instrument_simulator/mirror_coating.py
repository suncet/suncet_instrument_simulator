import os
import numpy as np
from pandas import read_fwf
import astropy.units as u

class MirrorCoating: 

    def __init__(self, config, target_wavelengths=np.array([171, 175,177, 180,185, 188, 194, 195, 202, 204, 211])*u.Angstrom):
        self.config = config
        self.wavelengths = target_wavelengths
        self.reflectivity = np.nan

        self.name = self.__get_coating_name()
        self.__interpolate()


    def __interpolate(self):
        filename = self.config.mirror_coating_reflectivity_filename
        source_data = self.load_mirror_data(filename)
        self.reflectivity = np.interp(self.wavelengths.value, source_data.wavelength.values, source_data.reflectivity.values)


    def __get_coating_name(self):
        return os.path.basename(self.config.mirror_coating_reflectivity_filename).split('_1')[0]    

    def load_mirror_data(self, filename):
        return read_fwf(os.getcwd() + '/' + filename, skiprows=18, header=0, names=['wavelength', 'reflectivity']) # wavelength should be in Angstroms, reflectivity as a fraction
    

if __name__ == "__main__":
    print("See test_mirror_coating.py for example of how to configure to run this code.")