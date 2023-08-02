import os
import sys
import astropy.units as u
from suncet_instrument_simulator import config_parser, instrument, make_radiance_maps
import setup_minimum_required_folders_files

root_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
print(root_directory)
print(os.path.exists(root_directory + '/setup_minimum_required_folders_files.py'))
sys.path.insert(0, root_directory)

def test_instrument():
    if os.getenv('suncet_data') == None:
        os.environ['suncet_data'] = './'
        setup_minimum_required_folders_files.run()
    
    hardware = setup_instrument_hardware()
    run_mirror_coating_tests(hardware)
  

def setup_instrument_hardware():
    config_filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    config = config_parser.Config(config_filename)
    hardware = instrument.Hardware(config)
    radiance_maps = make_radiance_maps.MakeRadianceMaps(config).run()
    hardware.store_target_wavelengths(radiance_maps)
    hardware.compute_effective_area()
    return hardware


def run_mirror_coating_tests(hardware):
    assert hardware.coating_name == 'B4C_Mo_Al'
    assert hardware.wavelengths.unit == u.Angstrom
    assert hardware.effective_area.unit == u.cm**2


if __name__ == "__main__":
    test_instrument()
