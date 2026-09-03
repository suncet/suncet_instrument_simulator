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

    assert header['DATE-OBS'] == '2011-02-15T17:00:20.000'
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


def test_three_viewpoint_uses_standard_model_subdirectories(tmp_path, monkeypatch):
    model_directory = tmp_path / 'mhd' / 'three_viewpoint_2000_kmps' / 'halo'
    em_map_directory = model_directory / 'em_maps'
    em_map_directory.mkdir(parents=True)
    em_map = em_map_directory / 'em_map_002.sav'
    em_map.touch()
    monkeypatch.setenv('suncet_data', str(tmp_path))

    config = SimpleNamespace(
        model_directory_name='/mhd/three_viewpoint_2000_kmps/halo/',
        em_map_directory_name='/em_maps/',
        euv_radiance_map_directory_name='/rendered_euv_maps/',
        timesteps_to_process=[2, 2, 1],
    )
    radiance_maps = MakeRadianceMaps(config)
    files = radiance_maps._MakeRadianceMaps__parse_filenames()
    radiance_maps.em_maps_filename = str(em_map)

    assert list(files) == [str(em_map)]
    assert radiance_maps._MakeRadianceMaps__make_outgoing_filename() == str(
        model_directory / 'rendered_euv_maps' / 'radiance_maps_002.fits'
    )
