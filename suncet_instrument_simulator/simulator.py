"""
This is the main wrapper for most/all(?) of the other instrument simulator related python files
"""
import os
from glob import glob
from datetime import datetime
import astropy.units as u
from astropy.io import fits
import sunpy.map
import pandas as pd
import ast
import numpy as np
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
        maps_by_index_and_wavelength = {}
        filenames = glob(os.getenv('suncet_data') + '/mhd/bright_fast/rendered_euv_maps/radiance_maps_04*.fits')
        for filename in filenames:
            index = os.path.basename(filename).split('_')[-1].replace('.fits', '')
            maps = sunpy.map.Map(filename)

            if index not in maps_by_index_and_wavelength:
                maps_by_index_and_wavelength[index] = {}

            for map in maps:
                wavelength = str(map.wavelength)
                maps_by_index_and_wavelength[index][wavelength] = map
            
        
        self.radiance_maps = maps_by_index_and_wavelength


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
        self.onboard_processed_images = self.onboard_software.filter_out_particle_hits(self.onboard_processed_images)
        self.onboard_processed_images = self.onboard_software.bin_image(self.onboard_processed_images)


    def __calculate_snr(self):
        pass # TODO: implement calculate_snr

    
    def __complete_metadata(self):
        metadata_definition = self.__load_metadata_definition()
        map = self.__strip_units_for_fits_compatibility(self.onboard_processed_images)

        tmp = fits.PrimaryHDU()
        header = tmp.header # A blank header with the required FITS keywords in the required order

        for i, row in metadata_definition.iterrows():
            if 'COMMENT' in str(row['Full variable name']): 
                header.set('COMMENT', value=row['Full variable name'].replace('COMMENT ', ''), after=len(header))  # line breaking comments for human readability
            elif row['FITS variable name'] not in header: 
                value = row['typical value']
                if pd.isna(value): 
                    value = 'N/A'
                else: 
                    try: 
                        value = ast.literal_eval(value)
                    except:
                        pass
                header.set(row['FITS variable name'], value=value, comment=row['Description'], after=len(header))

        # Populate metadata defined by the config file or resultant from the simulation
        header.set('NAXIS', value=3)
        header.set('NAXIS1', value=1) # TODO: this is the number of images in the stack through time.. but that should always be collapsed to 1 by this point in the simulator. Need to reconcile. 
        header.set('NAXIS2', value=map.dimensions[0].value)
        header.set('NAXIS3', value=map.dimensions[1].value)
        header.set('LEVEL', value='0.5')
        header.set('TIMESYS', value=datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ'))
        header.set('DATE-OBS', value=map.meta['DATE-OBS'])
        header.set('IMAGEW', value=map.dimensions[0].value)
        header.set('IMAGEH', value=map.dimensions[1].value)
        header.set('NBIN', value=(self.config.num_pixels_to_bin[0] * self.config.num_pixels_to_bin[1]))
        header.set('NBIN1', value=self.config.num_pixels_to_bin[0])
        header.set('NBIN2', value=self.config.num_pixels_to_bin[1])
        header.set('DET_TEMP', value=self.config.detector_temperature.value)
        header.set('EXPTIME', map.meta['EXPTIME'])

        hdu = fits.PrimaryHDU(map.data, header) # FIXME: the line above setting nan to None or '' doesn't work either. Need to figure out what FITS allows for empty. 
        hdul = fits.HDUList(hdu)

        self.fits = hdul
        pass
    

    def __load_metadata_definition(self):
        return pd.read_csv(os.getenv('suncet_data') + '/metadata/' + self.config.base_metadata_filename)


    def __output_files(self):
        self.__write_fits()
        self.__write_binary()
        self.__output_snr()
        pass # TODO: implement output_files
    

    def __write_fits(self):
        path = os.getenv('suncet_data') + '/synthetic/level0_raw/fits/'
        fits = self.fits
        filename = os.path.splitext(os.path.basename(self.config_filename))[0] + '_OBS_' + fits[0].header['DATE-OBS'] + '.fits'
        fits[0].header.set('FILENAME', value=filename)
        
        fits.writeto(path+filename, overwrite=True)
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


# Convenience function primarily for debugging -- provide the data for the image, e.g., map.data or fits.data and this will plot with good scalings
def __display_data(data): 
    import matplotlib.pyplot as plt
    plt.imshow(np.log10(np.clip(data, 0.1, None)))


if __name__ == "__main__":
    simulator = Simulator()
    simulator.run()

    # simulator.clean() maybe