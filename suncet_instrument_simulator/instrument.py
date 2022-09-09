import os
import numpy as np
from pandas import read_fwf
import astropy.units as u

class Hardware: 

    def __init__(self, config):
        self.config = config

        self.coating_name = self.__get_coating_name()
        self.wavelengths = self.__get_target_wavelengths()


    def __get_coating_name(self):
        return os.path.basename(self.config.mirror_coating_reflectivity_filename).split('_1')[0] 


    def __get_target_wavelengths(self): 
        if self.__radiance_maps_exist():
            headers = self.__load_radiance_map_headers()
            return self.__extract_wavelengths_from_headers()
        else: 
            return np.array([171, 175, 177, 180, 185, 188, 194, 195, 202, 204, 211])*u.Angstrom


    def __radiance_maps_exist(self):
        pass  # TODO implement radiance_maps_exist


    def __load_radiance_map_headers(self):
        pass  # TODO implement load_radiance_map_headers


    def __extract_wavelengths_from_headers(self):
        pass  # TODO implement extract_wavelengths_from_headers


    def interpolate_mirror_coating_reflectivity(self):
        filename = self.config.mirror_coating_reflectivity_filename
        source_data = self.__load_mirror_data(filename)
        self.reflectivity = np.interp(self.wavelengths.value, source_data.wavelength.values, source_data.reflectivity.values)
   

    def __load_mirror_data(self, filename):
        return read_fwf(os.getcwd() + '/' + filename, skiprows=18, header=0, names=['wavelength', 'reflectivity']) # wavelength should be in Angstroms, reflectivity as a fraction


    def extract_fov(self):
        pass # TODO implement extract FOV


class OnboardSoftware:
    def __init__(self, config):
        self.config = config


if __name__ == "__main__":
    print("See test_instrument.py for example of how to configure to run this code.")