from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time

from repair_simulator_output_timestamps import apply_repair, plan_repairs, validate_file


CONFIG = (Path(__file__).resolve().parents[1] / 'suncet_instrument_simulator'
          / 'config_files' / 'config_three_viewpoint_45deg.ini')


def make_files(tmp_path):
    radiance = (tmp_path / 'mhd' / 'three_viewpoint_2000_kmps' / '45deg'
                / 'rendered_euv_maps' / 'radiance_maps_048.fits')
    radiance.parent.mkdir(parents=True)
    fits.PrimaryHDU(header=fits.Header({'DATE-OBS': '2012-03-08T20:08:00.000'})).writeto(radiance)
    output = (tmp_path / 'synthetic' / 'level0' / 'fits'
              / 'config_three_viewpoint_45deg_OBS_2023-01-14T17:08:00.000_008.fits')
    output.parent.mkdir(parents=True)
    fits.PrimaryHDU(np.arange(12, dtype=np.uint16).reshape(3, 4), header=fits.Header({
        'DATE-OBS': '2023-01-14T17:08:00.000', 'FILENAME': output.name,
    })).writeto(output, checksum=True)
    return output


def test_repair_preserves_pixels_and_checksums_and_is_idempotent(tmp_path):
    source = make_files(tmp_path)
    original_bytes = source.read_bytes()
    plan, = plan_repairs(CONFIG, tmp_path)
    assert plan.action == 'repair'
    assert plan.fields['DATE-OBS'] == '2012-03-08T20:08:00.000'
    assert plan.fields['DATE-END'] == '2012-03-08T20:09:00.000'
    assert source.read_bytes() == original_bytes  # Audit never changes data.
    backup = tmp_path / 'backups'
    apply_repair(plan, backup)
    assert not source.exists()
    assert (backup / source.name).read_bytes() == original_bytes
    assert validate_file(plan.destination, plan.fields)[0]
    with fits.open(plan.destination) as hdul:
        np.testing.assert_array_equal(hdul[0].data, np.arange(12).reshape(3, 4))
        assert hdul[0].verify_checksum() == 1
        assert hdul[0].verify_datasum() == 1
    assert [item.action for item in plan_repairs(CONFIG, tmp_path)] == ['correct']


def test_existing_correct_output_is_preserved_and_old_duplicate_archived(tmp_path):
    source = make_files(tmp_path)
    plan, = plan_repairs(CONFIG, tmp_path)
    fits.PrimaryHDU(np.ones((3, 4), dtype=np.uint16), fits.Header(plan.fields)).writeto(
        plan.destination, checksum=True)
    correct_bytes = plan.destination.read_bytes()
    original_bytes = source.read_bytes()
    duplicate = next(item for item in plan_repairs(CONFIG, tmp_path) if item.source == source)
    assert duplicate.action == 'archive_duplicate'
    backup = tmp_path / 'backups'
    apply_repair(duplicate, backup)
    assert not source.exists()
    assert plan.destination.read_bytes() == correct_bytes
    assert (backup / source.name).read_bytes() == original_bytes


def test_radiance_hdus_use_their_own_time_scales(tmp_path):
    make_files(tmp_path)
    radiance = (tmp_path / 'mhd' / 'three_viewpoint_2000_kmps' / '45deg'
                / 'rendered_euv_maps' / 'radiance_maps_048.fits')
    instant = Time('2012-03-08T20:08:00.000', scale='utc')
    fits.HDUList([
        fits.PrimaryHDU(header=fits.Header({'DATE-OBS': instant.isot, 'TIMESYS': 'UTC'})),
        fits.ImageHDU(header=fits.Header({'DATE-OBS': instant.tai.isot, 'TIMESYS': 'TAI'})),
    ]).writeto(radiance, overwrite=True)
    plan, = plan_repairs(CONFIG, tmp_path)
    assert plan.fields['DATE-OBS'] == instant.isot
    assert plan.fields['TIMESYS'] == 'UTC'
