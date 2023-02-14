import os
from suncet_instrument_simulator import config_parser

class GoSunCET: 

    def __init__(self, config): 
        self.config = config


if __name__ == "__main__":
    config_filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    config = config_parser.Config(config_filename)
    hardware = GoSunCET(config)
    print("SunCET is cool.")
    