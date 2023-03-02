import os
from glob import glob
import numpy as np
from pandas import read_fwf, read_csv
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
import sunpy.map
from suncet_instrument_simulator import config_parser # This is just for running as a script -- should delete when done testing
from scipy.io import readsav # TODO: remove this temporary hack once we have radiance map files in non-IDL-saveset format

class Hardware: 
    def __init__(self, config):
        self.config = config
        self.coating_name = self.__get_coating_name()


    def __get_coating_name(self):
        return os.path.basename(self.config.mirror_coating_reflectivity_filename).split('_1')[0] 


    def store_target_wavelengths(self, radiance_maps): # Will interpolate/calculate any subsequent wavelength-dependent quantities to this "target" wavelength array
        wavelengths = []
        wavelength_unit = []
        for map in radiance_maps:
            wavelengths.append(map.wavelength.value)
            if len(wavelength_unit) > 0: 
                if map.wavelength.unit != wavelength_unit[0]:
                    raise ValueError('The wavelengths in the radiance maps do not match each other. This is dangerous so please fix it.')
            wavelength_unit.append(map.wavelength.unit)
        self.wavelengths = wavelengths * wavelength_unit[0]


    def compute_effective_area(self):
        geometric_area = np.pi * ((self.config.entrance_aperture)/2)**2
        mirror_reflectivity = self.__interpolate_mirror_coating_reflectivity()
        filter_transmission = self.__interpolate_filter_transmission()
        
        self.effective_area = geometric_area * mirror_reflectivity * filter_transmission


    def __interpolate_mirror_coating_reflectivity(self):
        filename = self.config.mirror_coating_reflectivity_filename
        source_data = self.__load_mirror_data(filename)
        return np.interp(self.wavelengths.value, source_data.wavelength.values, source_data.reflectivity.values)
   

    def __load_mirror_data(self, filename):
        return read_fwf(os.getenv('suncet_data') + '/mirror_reflectivity/' + filename, skiprows=18, header=0, names=['wavelength', 'reflectivity']) # wavelength should be in Angstroms, reflectivity as a fraction


    def __interpolate_filter_transmission(self): 
        filter_entrance = self.__load_filter_data(self.config.filter_entrance_transmission_filename)
        filter_focal_plane = self.__load_filter_data(self.config.filter_focal_plane_transmission_filename)
        combined_transmission = filter_entrance
        combined_transmission['transmission'] *= filter_focal_plane['transmission']
        return np.interp(self.wavelengths.value, combined_transmission['wavelength [Å]'].values, combined_transmission['transmission'].values)


    def __load_filter_data(self, filename):
        df = read_csv(os.getenv('suncet_data') + '/filter_transmission/' + filename, skiprows=2, usecols=[i for i in range(2)], names=['wavelength [nm]', 'transmission'])
        df['wavelength [nm]'] *= 10
        df.rename(columns={'wavelength [nm]': 'wavelength [Å]'}, inplace=True)
        return df 


    def extract_fov(self, radiance_maps):
        fov_half_angles = self.config.fov / 2.0

        # TODO: Raise warning if instrument FOV > model FOV

        return sunpy.map.MapSequence([map.submap(top_right=SkyCoord(fov_half_angles[0], fov_half_angles[1], frame=map.coordinate_frame), 
                                                 bottom_left=SkyCoord(-fov_half_angles[0], -fov_half_angles[1], frame=map.coordinate_frame)) 
                                      for map in radiance_maps])
    

    def interpolate_spatial_resolution(self, radiance_maps):
        map_list = []
        for map in radiance_maps: 
            map_list.append(map.resample(self.config.image_dimensions, method='spline')) # TODO: Figure out why the CDELTs in map.fits_header['CDELT1'] don't update with the resample
        return sunpy.map.MapSequence(map_list)


    
    def __copmute_plate_scale(): 
        pass


class OnboardSoftware:
    def __init__(self, config):
        self.config = config


if __name__ == "__main__":
    config_filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    config = config_parser.Config(config_filename)
    hardware = Hardware(config)
    print("See test_instrument.py for example of how to configure to run this code.")