"""
This is the main wrapper for most/all(?) of the other instrument simulator related python files
"""
import os
from glob import glob
import sunpy.map
from suncet_instrument_simulator import config_parser, make_radiance_maps, instrument

class Simulator:
    def __init__(self, config_filename=os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'):
        self.config = self.__read_config(config_filename)
        self.radiance_maps = 'not yet loaded'
        self.hardware = 'not yet loaded'
        self.onboard_software = 'not yet loaded'

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
        self.__output_files()
    

    def __sun_to_detector(self):
        if self.config.compute_new_radiance_maps:
            self.radiance_maps = make_radiance_maps.MakeRadianceMaps(os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini').run()
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
        pass # implement combine_noise_sources (self.noise_only, self.dark_frame, self.read_frame, and self.spike_frame)
    

    def __simulate_detector(self):
        self.detector_images = self.hardware.combine_signal_and_noise(self.detector_images, self.pure_signal, self.noise_only)
        self.detector_images = self.hardware.convert_to_dn(self.detector_images)
        self.detector_images = self.hardware.apply_screwy_pixels(self.detector_images, self.spike_frame)    


    def __apply_camera_software(self):
        if self.config.subtract_dark:
            self.onboard_processed_images = self.onboard_software.subtract_dark(self.detector_images, self.dark_frame)
        else: 
            self.onboard_processed_images = self.detector_images
        self.split_images = self.onboard_software.separate_images(self.onboard_processed_images)
        self.split_images = self.onboard_software.apply_jitter(self.split_images)
        self.onboard_processed_images = self.onboard_software.median_image_stack(self.split_images)
        self.onboard_processed_images - self.onboard_software.bin_image(self.onboard_processed_images)


    def __calculate_snr(self):
        pass # implement calculate_snr


    def __output_files(self):
        pass # implement output_files


if __name__ == "__main__":
    simulator = Simulator()
    simulator.run()

    # simulator.clean() maybe