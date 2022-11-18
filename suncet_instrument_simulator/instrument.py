import os
from glob import glob
import numpy as np
from pandas import read_fwf
import astropy.units as u
from astropy.io import fits
import sunpy.map
from suncet_instrument_simulator import config_parser # This is just for running as a script -- should delete when done testing
from scipy.io import readsav # TODO: remove this temporary hack once we have radiance map files in non-IDL-saveset format

class Hardware: 

    def __init__(self, config):
        self.config = config

        self.coating_name = self.__get_coating_name()
        self.wavelengths = self.__get_target_wavelengths()


    def __get_coating_name(self):
        return os.path.basename(self.config.mirror_coating_reflectivity_filename).split('_1')[0] 


    def __get_target_wavelengths(self):
        expected_filenames = os.getenv('suncet_data') + 'mhd/dimmest/rendered_euv_maps/euv_sim_300*.fits'
        radiance_maps_filenames = glob(expected_filenames)
        if len(radiance_maps_filenames) > 0: 
            return self.__extract_wavelengths_from_maps(radiance_maps_filenames)
        else:
            raise ValueError('Radiance maps not found. Make sure that your radiance maps are saved as {}'.format(expected_filenames))


    def __extract_wavelengths_from_maps(self, filenames):
        wavelengths = []
        wavelength_unit = []
        for filename in filenames:
            map = sunpy.map.Map(filename)
            wavelengths.append(map.wavelength.value)
            if len(wavelength_unit) > 0: 
                if map.wavelength.unit != wavelength_unit[0]:
                    raise ValueError('The wavelengths in the radiance maps do not match each other. This is dangerous so please fix it.')
            wavelength_unit.append(map.wavelength.unit)
        return wavelengths * wavelength_unit[0]


    def interpolate_mirror_coating_reflectivity(self):
        filename = self.config.mirror_coating_reflectivity_filename
        source_data = self.__load_mirror_data(filename)
        self.reflectivity = np.interp(self.wavelengths.value, source_data.wavelength.values, source_data.reflectivity.values)
   

    def __load_mirror_data(self, filename):
        return read_fwf(os.getenv('suncet_data') + '/mirror_reflectivity/' + filename, skiprows=18, header=0, names=['wavelength', 'reflectivity']) # wavelength should be in Angstroms, reflectivity as a fraction


    def extract_fov(self):
        
        pass # TODO: implement extract FOV


class OnboardSoftware:
    def __init__(self, config):
        self.config = config


if __name__ == "__main__":
    config_filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    config = config_parser.Config(config_filename)
    hardware = Hardware(config)
    print("See test_instrument.py for example of how to configure to run this code.")