import os
import sys
from glob import glob
import astropy.units as u
import sunpy.map
from suncet_instrument_simulator import config_parser, instrument, make_radiance_maps

root_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
print(root_directory)
print(os.path.exists(root_directory + '/setup_minimum_required_folders_files.py'))
sys.path.insert(0, root_directory)
import setup_minimum_required_folders_files

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
    make_radiance_maps.MakeRadianceMaps(config).run()
    filenames = glob(os.getenv('suncet_data') + '/mhd/bright_fast/rendered_euv_maps/radiance_maps_200.fits')
        
    maps_by_index_and_wavelength = {}
    for filename in filenames:
        index = os.path.basename(filename).split('_')[-1].replace('.fits', '')
        maps = sunpy.map.Map(filename)

        if index not in maps_by_index_and_wavelength:
            maps_by_index_and_wavelength[index] = {}

        for map in maps:
            wavelength = str(map.wavelength)
            maps_by_index_and_wavelength[index][wavelength] = map
    
    radiance_maps = maps_by_index_and_wavelength

    hardware.store_target_wavelengths(radiance_maps)
    hardware.compute_effective_area()
    return hardware


def run_mirror_coating_tests(hardware):
    assert hardware.wavelengths.unit == u.Angstrom
    assert hardware.effective_area.unit == u.cm**2


if __name__ == "__main__":
    test_instrument()
