from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest

from suncet_instrument_simulator.make_radiance_maps import MakeRadianceMaps


def _make_header(model_directory_name, frame_index=0):
    config = SimpleNamespace(
        model_directory_name=model_directory_name,
        model_timestep=10 * u.second,
    )
    radiance_maps = MakeRadianceMaps(config)
    radiance_maps.em_maps_filename = 'em_map_{:03d}.sav'.format(frame_index)
    radiance_maps.raw_radiance = np.empty((1, 1024, 1024))
    return radiance_maps._MakeRadianceMaps__make_header_template()


@pytest.mark.parametrize('model_name', ('bright_fast', 'dimmest', 'bright_slow'))
def test_legacy_models_calculate_plate_scale_from_fov(model_name):
    header = _make_header('/mhd/{}/'.format(model_name), frame_index=2)
    expected_plate_scale = 2 * 5.6 * 960 / 1024

    assert header['DATE-OBS'] == '2023-02-14T17:00:20.000'
    assert header['CDELT1'] == pytest.approx(expected_plate_scale)
    assert header['CDELT2'] == pytest.approx(expected_plate_scale)
    assert header['CDELT1'] == pytest.approx(10.5)


def test_three_viewpoint_header_uses_model_fov_and_frame_time():
    header = _make_header('/mhd/three_viewpoint_2000_kmps/', frame_index=2)
    expected_plate_scale = 2 * 6.4 * 960 / 1024

    assert header['DATE-OBS'] == '2012-03-08T20:00:20.000'
    assert header['CDELT1'] == pytest.approx(expected_plate_scale)
    assert header['CDELT2'] == pytest.approx(expected_plate_scale)
    assert header['CDELT1'] == pytest.approx(12.0)
