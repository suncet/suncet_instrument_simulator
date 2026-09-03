import os
import pytest

from suncet_instrument_simulator import config_parser

DERIVED_CONFIG_ATTRIBUTES = {
    'short_integration_window',
    'long_integration_window',
    'observation_window',
    'output_stride_model_steps',
}


CONFIG_FILENAMES = (
    'config_default.ini',
    'config_default_no_particle_filter.ini',
    'config_mars_telecom_network.ini',
    'config_three_viewpoint_45deg.ini',
    'config_three_viewpoint_halo.ini',
    'config_three_viewpoint_limb.ini',
)


@pytest.mark.parametrize('config_filename', CONFIG_FILENAMES)
def test_config_parser(config_filename):
    filename = os.path.join(
        os.getcwd(), 'suncet_instrument_simulator', 'config_files', config_filename)
    config = config_parser.Config(filename)

    num_parameters = 0
    with open(filename) as file:
        for line in file:
            if line.strip() and line[0] != '[' and line[0] != '#': 
                num_parameters += 1

    assert len(vars(config)) != 0
    assert len(vars(config)) == num_parameters + len(DERIVED_CONFIG_ATTRIBUTES)
    assert config.output_stride_model_steps >= 1
    assert hasattr(config, 'em_map_directory_name')
    assert hasattr(config, 'euv_radiance_map_directory_name')


@pytest.mark.parametrize('viewpoint', ('45deg', 'halo', 'limb'))
def test_three_viewpoint_config_paths(viewpoint):
    filename = os.path.join(
        os.getcwd(),
        'suncet_instrument_simulator',
        'config_files',
        'config_three_viewpoint_{}.ini'.format(viewpoint),
    )
    config = config_parser.Config(filename)
    expected_model_directory = '/mhd/three_viewpoint_2000_kmps/{}/'.format(viewpoint)

    assert config.model_directory_name == expected_model_directory
    assert config.model_data_folder == expected_model_directory
    assert config.em_map_directory_name == '/em_maps/'
    assert config.euv_radiance_map_directory_name == '/rendered_euv_maps/'
