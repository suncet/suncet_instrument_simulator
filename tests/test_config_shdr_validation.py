import configparser
import os
import tempfile

import pytest

from suncet_instrument_simulator import config_parser


def _write_config(filter_hits, num_short, num_long):
    config = configparser.ConfigParser()
    config.read(
        os.getcwd() + '/suncet_instrument_simulator/config_files/config_mars_telecom_network.ini'
    )
    config['behavior']['filter_out_particle_hits'] = str(filter_hits)
    config['shdr']['num_short_exposures_to_stack'] = str(num_short)
    config['shdr']['num_long_exposures_to_stack'] = str(num_long)

    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.ini', delete=False)
    config.write(temp_file)
    temp_file.close()
    return temp_file.name


def test_filter_true_with_single_member_stack_raises():
    filename = _write_config(True, 1, 4)
    try:
        with pytest.raises(ValueError, match='filter_out_particle_hits requires at least two integrations'):
            config_parser.Config(filename)
    finally:
        os.remove(filename)
