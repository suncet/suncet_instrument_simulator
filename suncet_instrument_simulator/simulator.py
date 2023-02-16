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

        self.hardware.compute_effective_area()
        self.radiance_maps = self.hardware.extract_fov(self.radiance_maps)
        self.radiance_maps = self.hardware.interpolate_spatial_resolution(self.radiance_maps)
        pass

    
    def __load_radiance_maps(self):
        filenames = glob(os.getenv('suncet_data') + 'mhd/dimmest/rendered_euv_maps/euv_sim_300*.fits')
        self.radiance_maps = sunpy.map.Map(filenames, sequence=True)


    def __simulate_noise(self):
        pass # TODO: implement simulate_noise
    

    def __simulate_detector(self):
        pass # TODO: implement simulate_detector


    def __apply_camera_software(self):
        pass # TODO: implement apply_camera_software


    def __calculate_snr(self):
        pass # implement calculate_snr


    def __output_files(self):
        pass # implement output_files


if __name__ == "__main__":
    simulator = Simulator()
    simulator.run()

    # simulator.clean() maybe