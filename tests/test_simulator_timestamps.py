from types import SimpleNamespace

import astropy.units as u
from astropy.io import fits
from astropy.time import Time
import numpy as np
import pandas as pd
import pytest
import sunpy.map

from suncet_instrument_simulator import instrument, stack_schedule
from suncet_instrument_simulator.simulator import Simulator, _set_observation_times_from_map


def test_observation_timestamps_use_radiance_map_time_and_duration():
    header = fits.Header()
    header['TIMESYS'] = 'UTC'
    header['DATE-BEG'] = '2023-01-14T17:00:00.000'
    header['DATE-OBS'] = '2023-01-14T17:00:00.000'
    header['DATE-END'] = '2023-01-14T17:00:59.650'
    header['DATE'] = '2024-01-14T17:00:00.000'
    map_header = fits.Header()
    map_header['DATE-OBS'] = '2012-03-08T20:08:00.000'

    _set_observation_times_from_map(header, map_header=map_header, duration=15)

    assert header['DATE-BEG'] == '2012-03-08T20:08:00.000'
    assert header['DATE-OBS'] == '2012-03-08T20:08:00.000'
    assert header['DATE-END'] == '2012-03-08T20:08:15.000'
    assert header['DATE'] == '2024-01-14T17:00:00.000'


def test_missing_radiance_map_observation_timestamp_is_rejected():
    header = fits.Header()
    map_header = fits.Header()

    with pytest.raises(ValueError, match='DATE-OBS'):
        _set_observation_times_from_map(header, map_header=map_header, duration=15)


def test_observation_timestamps_convert_source_time_scale_to_utc():
    header = fits.Header({'TIMESYS': 'UTC'})
    map_header = fits.Header({
        'DATE-OBS': '2012-03-08T20:00:34.000', 'TIMESYS': 'TAI',
    })

    _set_observation_times_from_map(header, map_header, duration=15 * u.s)

    assert header['TIMESYS'] == 'UTC'
    assert header['DATE-OBS'] == '2012-03-08T20:00:00.000'
    assert header['DATE-END'] == '2012-03-08T20:00:15.000'


def test_missing_source_timesys_defaults_to_utc_not_metadata_template():
    header = fits.Header({'TIMESYS': 'TAI'})
    map_header = fits.Header({'DATE-OBS': '2012-03-08T20:00:00.000'})

    _set_observation_times_from_map(header, map_header, duration=15)

    assert header['TIMESYS'] == 'UTC'
    assert header['DATE-OBS'] == '2012-03-08T20:00:00.000'


@pytest.mark.parametrize('start_seconds,stack_count,expected_start,expected_end', [
    (0., 4, '2012-03-08T20:00:00.000', '2012-03-08T20:01:00.000'),
    (15., 1, '2012-03-08T20:00:15.000', '2012-03-08T20:00:30.000'),
    (480., 4, '2012-03-08T20:08:00.000', '2012-03-08T20:09:00.000'),
])
def test_stack_to_written_fits_keeps_observation_time(
        monkeypatch, tmp_path, start_seconds, stack_count, expected_start, expected_end):
    simulator = Simulator.__new__(Simulator)
    simulator.config = SimpleNamespace(
        exposure_time_short=0.0035 * u.s, exposure_time_long=15 * u.s,
        model_timestep=10 * u.s, num_short_exposures_to_stack=stack_count,
        num_long_exposures_to_stack=stack_count, num_pixels_to_bin=[1, 1],
        detector_temperature=-10 * u.deg_C, observation_window=15 * stack_count * u.s,
        inner_fov_circle_radius=1 * u.dimensionless_unscaled,
    )
    simulator.stack_schedule = stack_schedule.build_stack_schedule_at_time(
        start_seconds, simulator.config)
    simulator.radiance_by_model_index = {}
    for index in simulator.stack_schedule.unique_model_indices:
        metadata = {
            'CTYPE1': 'HPLN-TAN', 'CTYPE2': 'HPLT-TAN',
            'CUNIT1': 'arcsec', 'CUNIT2': 'arcsec',
            'CDELT1': 12., 'CDELT2': 12., 'EXPTIME': 1.,
            'DATE-OBS': (Time('2012-03-08T20:00:00') + index * 10 * u.s).isot,
        }
        simulator.radiance_by_model_index[index] = {
            '170 Angstrom': sunpy.map.Map(np.full((2, 2), index + 1.), metadata),
        }

    # Exercise real stack creation, exposure application, stack collapse,
    # composite metadata, and FITS output without expensive optics or noise.
    simulator.hardware = instrument.Hardware.__new__(instrument.Hardware)
    simulator.hardware.config = simulator.config
    monkeypatch.setattr(simulator.hardware, 'store_target_wavelengths', lambda maps: None)
    monkeypatch.setattr(simulator.hardware, 'compute_effective_area', lambda: None)
    monkeypatch.setattr(simulator, '_Simulator__process_radiance_through_optics', lambda maps: maps)
    simulator._Simulator__sun_to_detector()
    detector_images = {
        exposure: {index: wavelengths['170 Angstrom'] for index, wavelengths in members.items()}
        for exposure, members in simulator.radiance_maps.items()
    }
    software = instrument.OnboardSoftware(simulator.config)
    monkeypatch.setattr(
        software, '_OnboardSoftware__get_solar_disk_center_and_radius_in_pixels',
        lambda image: (np.array([0.5, 0.5]) * u.pix, 1 * u.pix))
    if stack_count > 1:
        detector_images = software.filter_out_particle_hits(detector_images)
    collapsed = software.collapse_exposure_stacks(detector_images)
    simulator.onboard_processed_images = software.create_composite(collapsed)
    metadata_definition = pd.DataFrame([
        {'Field Name': keyword, 'FITS variable name': keyword, 'typical value': value}
        for keyword, value in {
            'TIMESYS': 'UTC', 'DATE-OBS': '2023-01-14T17:00:00.000',
            'DATE-BEG': '2023-01-14T17:00:00.000', 'DATE-END': '2023-01-14T17:00:59.650',
        }.items()
    ])
    monkeypatch.setattr(simulator, '_Simulator__load_metadata_definition', lambda: metadata_definition)
    simulator._Simulator__complete_metadata()

    output_dir = tmp_path / 'synthetic' / 'level0' / 'fits'
    output_dir.mkdir(parents=True)
    monkeypatch.setenv('suncet_data', str(tmp_path))
    simulator.config_filename = 'config_three_viewpoint_limb.ini'
    simulator.current_output_index = 8
    simulator._Simulator__write_fits()

    filename = 'config_three_viewpoint_limb_OBS_{}_008.fits'.format(expected_start)
    with fits.open(output_dir / filename, checksum=True) as hdul:
        header = hdul[0].header
        assert header['TIMESYS'] == 'UTC'
        assert header['DATE-OBS'] == expected_start
        assert header['DATE-BEG'] == expected_start
        assert header['DATE-END'] == expected_end
        assert header['FILENAME'] == filename
        assert hdul[0].verify_checksum() == 1
