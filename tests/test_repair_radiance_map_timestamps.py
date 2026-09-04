from astropy.io import fits
import numpy as np
import pytest

from repair_radiance_map_timestamps import (
    BACKUP_SUFFIX, expected_date_obs, find_radiance_maps, main, repair_file,
    restore_backup,
)


def test_repair_file_dry_run_and_write(tmp_path):
    output_directory = (
        tmp_path
        / 'mhd'
        / 'three_viewpoint_2000_kmps'
        / '45deg'
        / 'rendered_euv_maps'
    )
    output_directory.mkdir(parents=True)
    path = output_directory / 'radiance_maps_002.fits'
    old_date = '2023-02-14T17:00:00.000'
    hdul = fits.HDUList([
        fits.PrimaryHDU(data=np.arange(35, dtype=np.int16).reshape(5, 7),
                        header=fits.Header({'DATE-OBS': old_date})),
        fits.ImageHDU(data=np.arange(12, dtype=np.float64).reshape(3, 4),
                      header=fits.Header({'DATE-OBS': old_date})),
    ])
    hdul.writeto(path, checksum=True)
    original_bytes = path.read_bytes()
    backup_path = path.with_name(path.name + BACKUP_SUFFIX)
    with fits.open(path) as original:
        locations = [hdu.fileinfo() for hdu in original]

    dry_run = repair_file(path, cadence_seconds=10, write=False)
    assert dry_run.needs_update
    assert dry_run.expected_date_obs == '2012-03-08T20:00:20.000'
    assert path.read_bytes() == original_bytes
    assert not backup_path.exists()
    with fits.open(path) as unchanged:
        assert {hdu.header['DATE-OBS'] for hdu in unchanged} == {old_date}

    repair_file(path, cadence_seconds=10, write=True)
    assert backup_path.exists()
    repaired_bytes = path.read_bytes()
    assert len(repaired_bytes) == len(original_bytes)
    for info in locations:
        start, size = info['datLoc'], info['datSpan']
        assert repaired_bytes[start:start + size] == original_bytes[start:start + size]
    with fits.open(path) as repaired:
        assert {hdu.header['DATE-OBS'] for hdu in repaired} == {
            '2012-03-08T20:00:20.000'
        }
        assert [hdu.fileinfo() for hdu in repaired] == [
            dict(info, file=repaired._file) for info in locations
        ]
        assert all(hdu.verify_checksum() == 1 for hdu in repaired)
        assert all(hdu.verify_datasum() == 1 for hdu in repaired)

    backup_bytes = backup_path.read_bytes()
    assert not repair_file(path, write=True).needs_update
    assert path.read_bytes() == repaired_bytes
    assert backup_path.read_bytes() == backup_bytes

    assert restore_backup(backup_path) == path
    assert path.read_bytes() == original_bytes
    repair_file(path, write=True)
    assert path.read_bytes() == repaired_bytes
    assert backup_path.read_bytes() == backup_bytes


def test_existing_datasum_avoids_reading_image_data(tmp_path, monkeypatch):
    directory = tmp_path / 'bright_fast'
    directory.mkdir()
    path = directory / 'radiance_maps_001.fits'
    fits.PrimaryHDU(np.arange(15, dtype=np.int16),
                    fits.Header({'DATE-OBS': '2023-01-14T17:00:00.000'})).writeto(
                        path, checksum=True)

    def unexpected_data_read(*args, **kwargs):
        pytest.fail('Repair tried to checksum image data despite an existing DATASUM')

    monkeypatch.setattr(fits.PrimaryHDU, '_calculate_datasum', unexpected_data_read)
    repair_file(path, write=True)
    monkeypatch.undo()
    with fits.open(path) as repaired:
        assert repaired[0].header['DATE-OBS'] == '2011-02-15T17:00:10.000'
        assert repaired[0].verify_checksum() == 1


def test_checksum_without_datasum_is_updated(tmp_path):
    directory = tmp_path / 'dimmest'
    directory.mkdir()
    path = directory / 'radiance_maps_005.fits'
    hdu = fits.PrimaryHDU(np.arange(15, dtype=np.float32),
                          fits.Header({'DATE-OBS': '2023-01-14T17:00:00.000'}))
    hdu.add_checksum(override_datasum=True)
    hdu.writeto(path)
    repair_file(path, write=True)
    with fits.open(path) as repaired:
        assert repaired[0].verify_checksum() == 1
        assert repaired[0].verify_datasum() == 1


@pytest.mark.parametrize('cadence', [0, -1, float('nan'), float('inf')])
def test_invalid_cadence_is_rejected(tmp_path, cadence):
    with pytest.raises(ValueError, match='Cadence must be finite'):
        expected_date_obs(tmp_path / 'bright_fast' / 'radiance_maps_001.fits', cadence)
    assert main([str(tmp_path), '--cadence-seconds', str(cadence)]) == 1


def test_nonstandard_names_are_ignored_when_scanning(tmp_path):
    standard = tmp_path / 'radiance_maps_001.fits'
    standard.touch()
    duplicate = tmp_path / 'radiance_maps_121(4aZG).fits'
    duplicate.touch()
    prefixed = tmp_path / 'copy_radiance_maps_001.fits'
    prefixed.touch()
    with pytest.warns(UserWarning, match='Ignoring nonstandard'):
        assert find_radiance_maps([tmp_path, standard]) == [standard]
    for invalid in (duplicate, prefixed):
        with pytest.raises(ValueError):
            expected_date_obs(invalid, 10)
        with pytest.raises(ValueError):
            find_radiance_maps([invalid])


def test_restore_refuses_unrelated_header_edits(tmp_path):
    directory = tmp_path / 'bright_slow'
    directory.mkdir()
    path = directory / 'radiance_maps_000.fits'
    fits.PrimaryHDU(header=fits.Header({'DATE-OBS': '2023-01-14T17:00:00.000'})).writeto(path)
    repair_file(path, write=True)
    with fits.open(path, mode='update') as hdul:
        hdul[0].header['OBSERVER'] = 'subsequent edit'
    changed_bytes = path.read_bytes()
    with pytest.raises(ValueError, match='Header differs'):
        restore_backup(path.with_name(path.name + BACKUP_SUFFIX))
    assert path.read_bytes() == changed_bytes
