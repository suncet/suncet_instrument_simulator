import os
import urllib.request
import ssl
import numpy as np
import astropy.units as u
from suncet_instrument_simulator import config_parser, mirror_coating as mc

tmp_filename = 'B4C_Mo_Al_1-11000A.txt'


def test_mirror_coating():
    mirror_coating = setup_mirror_coating()

    assert mirror_coating.name == 'B4C_Mo_Al'
    assert mirror_coating.wavelengths.unit == u.Angstrom
    assert isinstance(mirror_coating.reflectivity, np.ndarray)

    delete_tmp_file()


def setup_mirror_coating():
    download_mirror_coating_data_file()
    config_filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    config = config_parser.Config(config_filename)
    return mc.MirrorCoating(config)


def download_mirror_coating_data_file():
    ssl._create_default_https_context = ssl._create_unverified_context
    url = urllib.request.urlopen("https://www.dropbox.com/s/bctrdr7de28m99o/B4C_Mo_Al_1-11000A.txt?dl=1")  # dl=1 is important
    data = url.read()
    url.close()
    with open(tmp_filename, "wb") as f:
        f.write(data)


def delete_tmp_file():
    os.remove(tmp_filename)


if __name__ == "__main__":
    test_mirror_coating()