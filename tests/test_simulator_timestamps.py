from astropy.io import fits

from suncet_instrument_simulator.simulator import _set_observation_times


def test_observation_timestamps_use_instrument_start_and_duration():
    header = fits.Header()
    header['TIMESYS'] = 'UTC'
    header['DATE-BEG'] = '2023-01-14T17:00:00.000'
    header['DATE-OBS'] = '2023-01-14T17:00:00.000'
    header['DATE-END'] = '2023-01-14T17:00:59.650'
    header['DATE'] = '2024-01-14T17:00:00.000'

    _set_observation_times(header, start_seconds=3600, duration=15)

    assert header['DATE-BEG'] == '2023-01-14T18:00:00.000'
    assert header['DATE-OBS'] == '2023-01-14T18:00:00.000'
    assert header['DATE-END'] == '2023-01-14T18:00:15.000'
    assert header['DATE'] == '2024-01-14T17:00:00.000'


def test_missing_optional_observation_timestamp_is_ignored():
    header = fits.Header()
    header['TIMESYS'] = 'UTC'
    header['DATE-OBS'] = '2023-01-14T17:00:00.000'

    _set_observation_times(header, start_seconds=15, duration=15)

    assert header['DATE-OBS'] == '2023-01-14T17:00:15.000'
    assert 'DATE-BEG' not in header
    assert 'DATE-END' not in header
