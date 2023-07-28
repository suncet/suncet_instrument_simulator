"""
This is the main wrapper for most/all(?) of the other instrument simulator related python files
"""
import os
from glob import glob
from datetime import datetime
import astropy.units as u
import sunpy.map
import pandas as pd
from suncet_instrument_simulator import config_parser, make_radiance_maps, instrument

class Simulator:
    def __init__(self, config_filename=os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'):
        self.config_filename = config_filename
        self.config = self.__read_config(config_filename)
        self.radiance_maps = 'not yet loaded'
        self.hardware = 'not yet loaded'
        self.onboard_software = 'not yet loaded'
        self.metadata = 'not yet loaded'

    def __read_config(self, config_filename):   
        return config_parser.Config(config_filename)


    def run(self):
        self.hardware = instrument.Hardware(self.config)
        self.onboard_software = instrument.OnboardSoftware(self.config)
        self.__sun_to_detector()
        self.__simulate_noise()
        self.__simulate_detector()
        self.__apply_camera_software()
        self.__calculate_snr()
        self.__complete_metadata()
        self.__output_files()
    

    def __sun_to_detector(self):
        if self.config.compute_new_radiance_maps:
            self.radiance_maps = make_radiance_maps.MakeRadianceMaps(self.config).run()
        else: 
            self.__load_radiance_maps()
        self.hardware.store_target_wavelengths(self.radiance_maps)

        self.hardware.compute_effective_area()
        self.radiance_maps = self.hardware.extract_fov(self.radiance_maps)
        self.radiance_maps = self.hardware.interpolate_spatial_resolution(self.radiance_maps)
        if self.config.apply_psf: 
            self.radiance_maps = self.hardware.apply_psf(self.radiance_maps)
        if self.config.apply_scattered_light_psf:
            self.radiance_maps = self.hardware.apply_scattered_light_psf(self.radiance_maps)
        self.radiance_maps = self.hardware.apply_effective_area(self.radiance_maps)
        self.radiance_maps = self.hardware.apply_exposure_times(self.radiance_maps)

    
    def __load_radiance_maps(self):
        filenames = glob(os.getenv('suncet_data') + 'mhd/dimmest/rendered_euv_maps/euv_sim_300*.fits')
        self.radiance_maps = sunpy.map.Map(filenames, sequence=True)


    def __simulate_noise(self):
        self.radiance_maps_pure = self.radiance_maps
        self.radiance_maps = self.hardware.apply_photon_shot_noise(self.radiance_maps)
        self.detector_images = self.hardware.convert_to_electrons(self.radiance_maps, apply_noise=True)
        self.detector_images_pure = self.hardware.convert_to_electrons(self.radiance_maps_pure, apply_noise=False)
        self.dark_frame = self.hardware.make_dark_frame()
        self.read_frame = self.hardware.make_read_frame()
        self.spike_frame = self.hardware.make_spike_frame()

        self.__combine_noise_sources()


    def __combine_noise_sources(self):
        self.noise_only = None
        pass # implement combine_noise_sources (self.noise_only, self.dark_frame, self.read_frame, and self.spike_frame)
    

    def __simulate_detector(self):
        self.detector_images = self.hardware.combine_signal_and_noise(self.detector_images, self.detector_images_pure, self.noise_only)
        self.detector_images = self.hardware.convert_to_dn(self.detector_images)
        self.detector_images = self.hardware.apply_screwy_pixels(self.detector_images, self.spike_frame)    


    def __apply_camera_software(self):
        if self.config.subtract_dark:
            self.onboard_processed_images = self.onboard_software.subtract_dark(self.detector_images, self.dark_frame)
        else: 
            self.onboard_processed_images = self.detector_images
        self.onboard_processed_images = self.onboard_software.apply_jitter(self.onboard_processed_images)
        self.onboard_processed_images = self.onboard_software.median_image_stack(self.onboard_processed_images)
        self.onboard_processed_images = self.onboard_software.bin_image(self.onboard_processed_images)


    def __calculate_snr(self):
        pass # TODO: implement calculate_snr

    
    def __complete_metadata(self):
        metadata_definition = self.__load_metadata_definition()
        fits_var_names = metadata_definition['FITS variable name'].tolist()
        fits_var_names = [name for name in fits_var_names if pd.notna(name)] # TODO Gets rid of the COMMENT lines but need to figure out a way to include them
        map = self.onboard_processed_images
        for fits_var_name in fits_var_names:
            if fits_var_name not in map.meta:
                typical_value = metadata_definition.loc[metadata_definition['FITS variable name'] == fits_var_name, 'typical value'].values[0]
                if not pd.isna(typical_value):
                    map.meta[fits_var_name] = typical_value
        
        # Populate metadata values
        map.meta['LEVEL'] = '0.5'
        map.meta['TIMESYS'] = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
        map.meta['IMAGEW'] = map.dimensions[0].value
        map.meta['IMAGEH'] = map.dimensions[1].value
        map.meta['NBIN'] = self.config.num_pixels_to_bin
        map.meta['DET_TEMP'] = self.config.detector_temperature.value
        
        # TODO: Figure out how to add comments after the value in fits output

        self.onboard_processed_images = map
        pass # TODO: implement complete_metadata
    

    def __load_metadata_definition(self):
        return pd.read_csv(os.getenv('suncet_data') + '/metadata/' + self.config.base_metadata_filename)


    def __output_files(self):
        self.__write_fits()
        self.__write_binary()
        self.__output_snr()
        pass # TODO: implement output_files
    

    def __write_fits(self):
        path = os.getenv('suncet_data') + '/synthetic/level0_raw/fits/'
        filename = os.path.splitext(os.path.basename(self.config_filename))[0] + '.fits' # TODO: will have to deal with unique filenames for different timestamps here
        map = self.__strip_units_for_fits_compatibility(self.onboard_processed_images)
        map.meta['FILENAME'] = filename
        
        map.save(path+filename, filetype='fits', overwrite=True)
        pass # TODO: implement write_fits()


    def __strip_units_for_fits_compatibility(self, map):
        meta = map.meta
        for key, value in meta.items():
            if isinstance(value, u.Quantity):
                value = value.value
            meta[key] = value
        return sunpy.map.Map(map.data, meta)


    def __write_binary(self):
        pass # TODO: implement write_binary() to mimic onboard (or downlinked?) storage


    def __output_snr(self):
        pass # TODO: implement output_snr()

if __name__ == "__main__":
    simulator = Simulator()
    simulator.run()

    # simulator.clean() maybe