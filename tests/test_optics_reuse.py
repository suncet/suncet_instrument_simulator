"""Compare cached scenes with separately processed integration stacks."""
from copy import deepcopy
from pathlib import Path

import astropy.units as u
from astropy.io import fits
from astropy.time import Time
import numpy as np
import pytest
import sunpy.map

from suncet_instrument_simulator import config_parser, instrument, stack_schedule
from suncet_instrument_simulator.simulator import Simulator

CONFIG = Path(__file__).resolve().parents[1] / 'suncet_instrument_simulator/config_files/config_default.ini'


def make_simulator(start, counts=(9, 4), available=None):
    sim = Simulator.__new__(Simulator)
    sim.config = config_parser.Config(CONFIG)
    sim.config.image_dimensions = [24, 16] * u.pix
    sim.config.fov = [144, 96] * u.arcsec
    sim.config.plate_scale = 6 * u.arcsec / u.pix
    sim.config.num_short_exposures_to_stack, sim.config.num_long_exposures_to_stack = counts
    sim.stack_schedule = stack_schedule.build_stack_schedule_at_time(start, sim.config, available)
    indices = sim.stack_schedule.unique_model_indices
    sim.radiance_by_model_index = {}
    yy, xx = np.indices((12, 16))
    for index in indices:
        wavelengths = {}
        for wavelength in range(170, 174):
            metadata = {
                'CTYPE1': 'HPLN-TAN', 'CTYPE2': 'HPLT-TAN',
                'CUNIT1': 'arcsec', 'CUNIT2': 'arcsec',
                'CRPIX1': 8.5, 'CRPIX2': 6.5, 'CRVAL1': 0., 'CRVAL2': 0.,
                'CDELT1': 12., 'CDELT2': 12., 'EXPTIME': 1.,
                'DATE-OBS': (Time('2012-03-08T20:00:00') + index * 10 * u.s).isot,
                'TIMESYS': 'UTC', 'WAVELNTH': wavelength, 'WAVEUNIT': 'Angstrom',
                'BUNIT': 'ph / (s sr Angstrom)', 'DSUN_OBS': 1.496e11,
                'HGLN_OBS': 0., 'HGLT_OBS': 0., 'RSUN_REF': 6.96e8,
            }
            data = (1 + index / 100 + wavelength / 170) * (2 + np.cos(xx) + yy / 3)
            wavelengths[f'{wavelength}.0 Angstrom'] = sunpy.map.Map(data, metadata)
        sim.radiance_by_model_index[index] = wavelengths
    sim.hardware = instrument.Hardware.__new__(instrument.Hardware)
    sim.hardware.config = sim.config
    kernel = np.arange(1, 26, dtype=float).reshape(5, 5)
    kernel /= kernel.sum()
    sim.hardware.mesh_diffraction_psf = fits.HDUList([fits.ImageHDU(kernel) for _ in range(4)])
    sim.hardware.mirror_scatter_psf = np.ones((3, 3)) / 100

    def effective_area():
        sim.hardware.effective_area = {
            str(w): 0.1 * u.cm**2 for w in sim.hardware.wavelengths
        }
    sim.hardware.compute_effective_area = effective_area
    return sim


def independently_process(sim):
    result = {}
    for exposure, members, exposure_time in (
        ('short exposure', sim.stack_schedule.short_members, sim.config.exposure_time_short),
        ('long exposure', sim.stack_schedule.long_members, sim.config.exposure_time_long),
    ):
        scenes = stack_schedule.build_radiance_by_stack_member(
            members, sim.radiance_by_model_index,
            start_seconds=sim.stack_schedule.start_seconds, exposure_time=exposure_time,
            model_timestep=sim.config.model_timestep)
        sim.hardware.store_target_wavelengths(scenes)
        sim.hardware.compute_effective_area()
        scenes = sim._Simulator__process_radiance_through_optics(scenes)
        result[exposure] = sim.hardware.apply_exposure_times_for_stack(scenes, exposure_time)
    return result


@pytest.mark.parametrize('start,counts,available,expected_passes', [
    (0., (9, 4), None, 5),
    (15., (9, 4), None, 5),
    (15., (1, 1), None, 2),
    (3600., (9, 4), [360], 3),
])
def test_reused_optics_matches_independent_members(start, counts, available, expected_passes):
    reference = make_simulator(start, counts, available)
    expected = independently_process(reference)
    sim = make_simulator(start, counts, available)
    before = {
        (i, w): (m.data.copy(), deepcopy(m.meta))
        for i, waves in sim.radiance_by_model_index.items() for w, m in waves.items()
    }
    calls = []
    process = sim._Simulator__process_radiance_through_optics

    def counted(members):
        calls.append(len(members))
        return process(members)
    sim._Simulator__process_radiance_through_optics = counted
    sim._Simulator__sun_to_detector()
    assert calls == [1] * expected_passes
    for exposure, members in expected.items():
        for index, waves in members.items():
            for wavelength, expected_map in waves.items():
                actual = sim.radiance_maps[exposure][index][wavelength]
                np.testing.assert_array_equal(actual.data, expected_map.data)
                assert actual.meta == expected_map.meta
    for (index, wavelength), (data, metadata) in before.items():
        source = sim.radiance_by_model_index[index][wavelength]
        np.testing.assert_array_equal(source.data, data)
        assert source.meta == metadata
    # Metadata mutated by later exposure, electron, and jitter stages must not alias.
    first = sim.radiance_maps['short exposure'][0]['170.0 Angstrom']
    last = sim.radiance_maps['long exposure'][counts[1] - 1]['170.0 Angstrom']
    first.meta['crpix1'] = -123
    assert last.meta['crpix1'] != -123
    if counts[0] > 1:
        second = sim.radiance_maps['short exposure'][1]['170.0 Angstrom']
        assert np.shares_memory(first.data, second.data)
        assert not first.data.flags.writeable
        with pytest.raises(ValueError, match='read-only'):
            first.data[0, 0] = -1


def test_member_key_preserves_order_and_exact_weights():
    contribution = stack_schedule.MapContribution
    assert stack_schedule.radiance_member_key([contribution(1, .5), contribution(2, .5)]) != (
        stack_schedule.radiance_member_key([contribution(2, .5), contribution(1, .5)]))
    assert stack_schedule.radiance_member_key([contribution(1, .5), contribution(2, .5)]) != (
        stack_schedule.radiance_member_key([contribution(1, .5 + 1e-15), contribution(2, .5)]))


@pytest.mark.parametrize('retain', [False, True])
def test_pure_reference_is_optional_without_changing_noisy_path(retain):
    sim = make_simulator(0., (1, 1))
    sim.retain_pure_reference = retain
    sim._Simulator__sun_to_detector()
    assert (sim.radiance_maps_pure is not None) == retain
    assert sim.detector_images_pure is None
    expected_radiance = sim.radiance_maps
    calls = []
    sim.hardware.apply_photon_shot_noise = lambda images: images

    def convert(images, apply_noise):
        calls.append(apply_noise)
        return images
    sim.hardware.convert_to_electrons = convert
    for name in ['make_dark_frame', 'make_read_frame', 'make_hot_pixel_mask', 'make_dead_pixel_mask']:
        setattr(sim.hardware, name, lambda: None)
    sim.hardware.make_spike_masks = lambda images: None
    sim._Simulator__simulate_noise()
    assert calls == ([True, False] if retain else [True])
    assert sim.detector_images is expected_radiance
    assert (sim.detector_images_pure is not None) == retain


@pytest.mark.parametrize('dtype', [np.float64, np.uint16, np.uint32])
def test_in_memory_fits_header_preserves_saved_values(tmp_path, dtype):
    metadata = {'CTYPE1': 'HPLN-TAN', 'CTYPE2': 'HPLT-TAN',
                'CUNIT1': 'arcsec', 'CUNIT2': 'arcsec', 'CDELT1': 6., 'CDELT2': 6.,
                'DATE-OBS': '2012-03-08T20:00:00.035', 'EXPTIME': .035}
    image = sunpy.map.Map(np.arange(12).reshape(3, 4).astype(dtype), metadata)
    filename = tmp_path / 'legacy.fits'
    image.save(filename)
    with fits.open(filename) as hdul:
        expected = dict(hdul[0].header)
    actual = Simulator.__new__(Simulator)._Simulator__convert_sunpy_meta_to_fits_header(image)
    assert dict(actual) == expected
