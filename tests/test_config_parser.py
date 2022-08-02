import os
from suncet_instrument_simulator import config_parser

def test_config_parser():
    filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    config = config_parser.Config(filename)

    num_parameters = 0
    with open(filename) as file:
        for line in file:
            if line.strip() and line[0] != '[' and line[0] != '#': 
                num_parameters += 1

    assert len(vars(config)) != 0
    assert len(vars(config)) == num_parameters
