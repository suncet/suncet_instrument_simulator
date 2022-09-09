"""
This is the main wrapper for most/all(?) of the other instrument simulator related python files
"""
import os
from suncet_instrument_simulator import config_parser, make_radiance_maps, instrument

class Simulator:
    def __init__(self, config_filename=os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'):
        self.config = self.__load_config(config_filename)
        self.radiance_maps = 'not yet loaded'
        self.hardware = 'not yet loaded'
        self.onboard_software = 'not yet loaded'

    def __load_config(self, config_filename):   
        self.config = config_parser.Config(config_filename)


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
            make_radiance_maps()
        self.__load_radiance_maps()

        self.hardware.interpolate_mirror_coating_reflectivity()
        self.hardware.extract_fov()

    
    def __load_radiance_maps(self):
        pass  # TODO implement load_radiance_maps


    def __simulate_noise(self):
        pass # TODO implement simulate_noise
    

    def __simulate_detector(self):
        pass # TODO implement simulate_detector


    def __apply_camera_software(self):
        pass # TODO implement apply_camera_software


    def __calculate_snr(self):
        pass # implement calculate_snr


    def __output_files(self):
        pass # implement output_files


if __name__ == "__main__":
    simulator = Simulator()
    simulator.run()