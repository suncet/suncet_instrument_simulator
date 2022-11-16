import os
from pathlib import Path
import urllib.request
import ssl
import numpy as np
import astropy.units as u
from urllib.parse import urlparse
from suncet_instrument_simulator import config_parser, instrument

tmp_file_urls = ["https://www.dropbox.com/s/bctrdr7de28m99o/B4C_Mo_Al_1-11000A.txt?dl=1", 
                 "https://www.dropbox.com/s/f51fep2nu1vr7ai/euv_sim_300_171A.fits?dl=1", 
                 "https://www.dropbox.com/s/5tkfvphczs0a929/euv_sim_300_193A.fits?dl=1"] # dl=1 is important

def test_instrument():
    setup_environment()
    hardware = setup_instrument_hardware()
    run_mirror_coating_tests(hardware)
    delete_tmp_file()


def setup_environment(): 
    download_test_data()
    
    if os.getenv('suncet_data') == None: 
        os.environ['suncet_data'] = './'
        path = Path('./mhd/dimmest/rendered_euv_maps')
        path.mkdir(parents=True, exist_ok=True)
        for url in tmp_file_urls: 
            filename = url.rsplit('/', 1)[-1].rsplit('?', 1)[0]
            Path('./' + filename).rename(path / filename)


def setup_instrument_hardware():
    config_filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    config = config_parser.Config(config_filename)
    hardware = instrument.Hardware(config)
    hardware.interpolate_mirror_coating_reflectivity()
    return hardware


def download_test_data():
    ssl._create_default_https_context = ssl._create_unverified_context

    for url in tmp_file_urls:
        thisurl = urllib.request.urlopen(url)  
        data = thisurl.read()
        thisurl.close()
        filename = get_filename_from_url(url)
        with open(filename, "wb") as f:
            f.write(data)


def get_filename_from_url(url):
    parsed_url = urlparse(url)
    return os.path.basename(parsed_url.path)


def run_mirror_coating_tests(hardware):
    assert hardware.coating_name == 'B4C_Mo_Al'
    assert hardware.wavelengths.unit == u.Angstrom
    assert isinstance(hardware.reflectivity, np.ndarray)


def delete_tmp_file():
    for url in tmp_file_urls:
        filename = get_filename_from_url(url)
        os.remove(filename)


if __name__ == "__main__":
    test_instrument()