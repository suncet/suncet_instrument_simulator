import os
import urllib.request
import ssl
import numpy as np
import astropy.units as u
from suncet_instrument_simulator import config_parser, instrument

tmp_filename = 'B4C_Mo_Al_1-11000A.txt'


def test_instrument():
    hardware = setup_instrument_hardware()
    run_mirror_coating_tests(hardware)
    delete_tmp_file()


def setup_instrument_hardware():
    download_mirror_coating_data_file()
    config_filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    config = config_parser.Config(config_filename)
    hardware = instrument.Hardware(config)
    hardware.interpolate_mirror_coating_reflectivity()
    return hardware


def download_mirror_coating_data_file():
    ssl._create_default_https_context = ssl._create_unverified_context
    url = urllib.request.urlopen("https://www.dropbox.com/s/bctrdr7de28m99o/B4C_Mo_Al_1-11000A.txt?dl=1")  # dl=1 is important
    data = url.read()
    url.close()
    with open(tmp_filename, "wb") as f:
        f.write(data)


def run_mirror_coating_tests(hardware):
    assert hardware.coating_name == 'B4C_Mo_Al'
    assert hardware.wavelengths.unit == u.Angstrom
    assert isinstance(hardware.reflectivity, np.ndarray)


def delete_tmp_file():
    os.remove(tmp_filename)


if __name__ == "__main__":
    test_instrument()