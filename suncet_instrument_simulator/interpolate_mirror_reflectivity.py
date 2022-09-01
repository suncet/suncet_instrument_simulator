import os
import numpy as np
import pandas as pd
import astropy.units as u
from suncet_instrument_simulator import config_parser # For script testing only -- should delete this and run from the main wrapper script

class MirrorReflectivity: 

    def __init__(self, config, target_wavelengths=np.array([171, 175,177, 180,185, 188, 194, 195, 202, 204, 211])*u.Angstrom):
        self.df = pd.DataFrame()
        self.config = config
        self.target_wavelengths = target_wavelengths


    def interpolate(self):
        filename = self.config.mirror_coating_reflectivity_filename
        source_data = self.load_mirror_data(filename)
        return np.interp(self.target_wavelengths, source_data.wavelength.values, source_data.reflectivity.values)



    def load_mirror_data(self, filename):
        return pd.read_fwf(os.getcwd() + '/' + filename, skiprows=18, header=0, names=['wavelength', 'reflectivity']) # wavelength should be in Angstroms, reflectivity as a fraction
    

if __name__ == "__main__":
    reflectivity_filename = os.getenv('suncet_data') + 'mirror_reflectivity/B4C_Mo_Al_1-11000A.txt'
    config_filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    
    config = config_parser.Config(config_filename)
    mirror_reflectivity = MirrorReflectivity(config)
    mirror_reflectivity.interpolate()