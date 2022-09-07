"""
This is the main wrapper for most/all(?) of the other instrument simulator related python files
"""
import os
from suncet_instrument_simulator import config_parser, make_radiance_maps, instrument as instr

class Simulator:
    def __init__(self, config_filename=os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'):
        self.config = self.__load_config(config_filename)
        self.radiance_maps = 'not yet loaded'

    def __load_config(self, config_filename):   
        self.config = config_parser.Config(config_filename)


    def run(self): 
        if self.config.compute_new_radiance_maps:
            make_radiance_maps()
        self.__load_radiance_maps()

        instrument = instr.Instrument(self.config)
        instrument.interpolate_mirror_coating_reflectivity()
        instrument.extract_fov()

    
    def __load_radiance_maps(self):
        pass  # TODO implement load_radiance_maps


if __name__ == "__main__":
    simulator = Simulator()
    simulator.run()