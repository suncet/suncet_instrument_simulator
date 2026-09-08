"""Numerical and allocation regressions for the bounded instrument kernels."""
from types import SimpleNamespace

import astropy.units as u
from astropy.convolution import convolve_fft
import numpy as np
import pytest
from scipy.integrate import simpson
import sunpy.map

from suncet_instrument_simulator import instrument


def _map(data, wavelength=170., unit=u.ph / (u.Angstrom * u.pix**2)):
    return sunpy.map.Map(np.asanyarray(data), {
        'CTYPE1': 'HPLN-TAN', 'CTYPE2': 'HPLT-TAN',
        'CUNIT1': 'arcsec', 'CUNIT2': 'arcsec',
        'CDELT1': 4.8, 'CDELT2': 4.8,
        'DATE-OBS': '2012-03-08T20:00:00.035', 'TIMESYS': 'UTC',
        'EXPTIME': 0.035, 'WAVELNTH': wavelength, 'WAVEUNIT': 'Angstrom',
        'BUNIT': unit.to_string(),
    })


@pytest.mark.parametrize('scene', ['constant', 'impulse', 'structured'])
def test_finite_convolution_matches_original_and_owns_compact_output(scene, monkeypatch):
    data = np.ones((9, 12))
    if scene == 'impulse':
        data[:] = 0
        data[1, 3] = 1e8
    elif scene == 'structured':
        data = np.geomspace(1e-5, 1e8, data.size).reshape(data.shape)
    y, x = np.mgrid[-5:6, -5:6]
    kernel = np.exp(-(x*x + y*y) / 12.)
    kernel *= 0.03 / kernel.sum()
    expected = convolve_fft(data, kernel, boundary='wrap', normalize_kernel=False)
    calls = []

    def record(*args, **kwargs):
        calls.append(kwargs)
        return convolve_fft(*args, **kwargs)

    monkeypatch.setattr(instrument, 'convolve_fft', record)
    actual = instrument._convolve_psf(data, kernel)

    np.testing.assert_allclose(actual, expected, rtol=3e-14, atol=1e-8)
    assert calls[0]['nan_treatment'] == 'fill'
    assert calls[0]['fftn'].keywords['workers'] == 2
    assert calls[0]['ifftn'].keywords['workers'] == 2
    assert actual.shape == data.shape
    assert actual.dtype == np.float64
    assert actual.flags.owndata and actual.flags.c_contiguous
    assert actual.base is None
    assert actual.nbytes == data.size * 8


@pytest.mark.parametrize('special', ['nan', 'inf', 'mask', 'kernel_nan'])
def test_nonfinite_convolution_preserves_original_fallback(special, monkeypatch):
    data = np.arange(30., dtype=float).reshape(5, 6)
    kernel = np.ones((3, 3)) / 9.
    if special == 'mask':
        data = np.ma.array(data, mask=data == 7)
    elif special == 'kernel_nan':
        kernel[0, 0] = np.nan
    else:
        data[1, 1] = np.nan if special == 'nan' else np.inf
    expected = convolve_fft(data, kernel, boundary='wrap', normalize_kernel=False)
    calls = []

    def record(*args, **kwargs):
        calls.append(kwargs)
        return convolve_fft(*args, **kwargs)

    monkeypatch.setattr(instrument, 'convolve_fft', record)
    actual = instrument._convolve_psf(data, kernel)

    np.testing.assert_array_equal(actual, expected)
    assert calls == [{'boundary': 'wrap', 'normalize_kernel': False}]
    assert actual.flags.owndata


@pytest.mark.parametrize('kernel', [np.zeros((3, 3)), np.eye(3) * 1e-10])
def test_nonnormalizable_kernel_keeps_original_error(kernel):
    with pytest.raises(ValueError, match='Cannot interpolate NaNs'):
        instrument._convolve_psf(np.ones((5, 6)), kernel)


def test_hardware_preserves_psf_crops_scatter_fraction_and_metadata():
    hardware = instrument.Hardware.__new__(instrument.Hardware)
    kernel = np.arange(1., 145.).reshape(12, 12)
    kernel /= kernel.sum()
    hardware.mesh_diffraction_psf = [SimpleNamespace(data=kernel)]
    hardware.mirror_scatter_psf = kernel[:-1, :-1] * 0.01
    source = _map(np.arange(108., dtype=float).reshape(9, 12))

    diffracted = hardware.apply_diffraction_psf({0: {'170': source}})[0]['170']
    expected_diffraction = convolve_fft(
        source.data, kernel[:-1, :-1], boundary='wrap', normalize_kernel=False)
    np.testing.assert_allclose(diffracted.data, expected_diffraction, rtol=3e-14)
    assert diffracted.data.flags.owndata
    assert diffracted.meta == source.meta

    scattered = hardware.apply_mirror_scattered_light_psf({0: {'170': source}})[0]['170']
    expected_scatter = source.data * (1 - hardware.mirror_scatter_psf.sum()) + convolve_fft(
        source.data, hardware.mirror_scatter_psf, boundary='wrap', normalize_kernel=False)
    np.testing.assert_allclose(scattered.data, expected_scatter, rtol=3e-14)
    assert scattered.meta == source.meta


def _hardware(wavelengths, full_well=1e30):
    hardware = instrument.Hardware.__new__(instrument.Hardware)
    hardware.wavelengths = np.asarray(wavelengths) * u.Angstrom
    hardware.config = SimpleNamespace(
        fano_factor=0.119, pixel_full_well=full_well * u.electron / u.pix**2)
    return hardware


def _bands(hardware, mixed_units=False, shape=(5, 6)):
    result = {}
    for i, wavelength in enumerate(hardware.wavelengths.value):
        unit = u.ph / (u.Angstrom * u.pix**2)
        if mixed_units and i % 2:
            unit = u.ph / (u.nm * u.pix**2)
        data = (np.arange(np.prod(shape)).reshape(shape) + 0.25) * (i + 1)
        result[str(wavelength)] = _map(data, wavelength, unit)
    return result


def _legacy_electron_map(hardware, bands):
    quantum_yield = hardware._Hardware__compute_quantum_yields()
    first = next(iter(bands.values()))
    values = [m.data * m.unit * quantum_yield[i] * (1 * u.ct / u.electron)
              for i, m in enumerate(bands.values())]
    integrated = simpson(np.stack(values, axis=-1), x=hardware.wavelengths, axis=-1)
    quantity = (np.zeros_like(first.data) * first.unit * quantum_yield.unit
                * (1 * u.ct / u.electron) * hardware.wavelengths.unit)
    quantity += integrated * quantity.unit
    metadata = first.meta.copy()
    metadata['bunit'] = quantity.unit.to_string()
    return sunpy.map.Map(quantity, metadata)


@pytest.mark.parametrize('count', [1, 2, 5, 40])
@pytest.mark.parametrize('mixed_units', [False, True])
def test_streaming_electrons_match_simpson_policy_and_units(count, mixed_units):
    wavelengths = 170 + np.cumsum(np.linspace(0.4, 1.7, count))
    hardware = _hardware(wavelengths)
    bands = _bands(hardware, mixed_units=mixed_units)
    expected = _legacy_electron_map(hardware, bands)
    inputs = {'short exposure': {0: bands}, 'long exposure': {}}

    actual = hardware.convert_to_electrons(inputs, apply_noise=False)['short exposure'][0]

    np.testing.assert_allclose(actual.data, expected.data, rtol=3e-14, atol=1e-10)
    assert actual.unit == expected.unit
    assert actual.meta['date-obs'] == '2012-03-08T20:00:00.035'
    assert actual.meta['timesys'] == 'UTC'
    assert actual.meta['exptime'] == 0.035
    assert actual.meta['wavelnth'] == pytest.approx(hardware.wavelengths.value.mean())


def test_streaming_preserves_fano_draw_order_and_full_well():
    hardware = _hardware([170., 171., 173.], full_well=800)
    inputs = {exposure: {i: _bands(hardware) for i in range(2)}
              for exposure in ('short exposure', 'long exposure')}
    np.random.seed(7812)
    expected = []
    for members in inputs.values():
        for bands in members.values():
            image = _legacy_electron_map(hardware, bands)
            image = hardware._Hardware__apply_fano_noise(image)
            expected.append(hardware._Hardware__clip_at_full_well(image).data)
    expected_rng = np.random.get_state()

    np.random.seed(7812)
    actual = hardware.convert_to_electrons(inputs, apply_noise=True)
    actual_rng = np.random.get_state()

    for image, reference in zip(
            [m for members in actual.values() for m in members.values()], expected):
        np.testing.assert_array_equal(image.data, reference)
        assert image.data.max() == 800
    assert actual_rng[0] == expected_rng[0]
    np.testing.assert_array_equal(actual_rng[1], expected_rng[1])
    assert actual_rng[2:] == expected_rng[2:]


def test_streaming_rejects_missing_bands_and_broadcastable_shapes():
    hardware = _hardware([170., 171.])
    with pytest.raises(ValueError, match='count'):
        hardware.convert_to_electrons(
            {'short exposure': {0: {'170': _map(np.ones((5, 6)))}}, 'long exposure': {}},
            apply_noise=False)
    bands = _bands(hardware)
    bands['171.0'] = _map(np.ones((1, 6)), wavelength=171.)
    with pytest.raises(ValueError, match='matching shapes'):
        hardware.convert_to_electrons(
            {'short exposure': {0: bands}, 'long exposure': {}}, apply_noise=False)


class _MapStub:
    def __init__(self, data, meta=None):
        self.data = np.asanyarray(data)
        self.meta = {} if meta is None else meta


@pytest.mark.parametrize('mixed_dtype', [False, True])
def test_running_particle_filter_matches_uint32_stack_and_preserves_inputs(monkeypatch, mixed_dtype):
    arrays = [
        np.array([[4294967295, 100], [65536, 0]], dtype=np.uint32),
        np.array([[4294967294, 60000], [65537, 1]], dtype=np.uint32),
        np.array([[60000, 60000], [65538, 2]], dtype=np.uint32),
    ]
    if mixed_dtype:
        arrays[1] = arrays[1].astype(np.float64) - 0.5
    original = [array.copy() for array in arrays]
    stacked = np.stack(arrays, axis=-1).astype(np.uint32)
    expected = np.sum(stacked, axis=-1, dtype=np.uint32) - np.max(stacked, axis=-1)
    inputs = {exposure: {i: _MapStub(array, {'member': i})
                         for i, array in reversed(list(enumerate(arrays)))}
              for exposure in ('short exposure', 'long exposure')}
    monkeypatch.setattr(instrument.sunpy.map, 'Map', _MapStub)

    actual = instrument.OnboardSoftware(SimpleNamespace()).filter_out_particle_hits(inputs)

    for image in actual.values():
        np.testing.assert_array_equal(image.data, expected)
        assert image.data.dtype == np.uint32
        assert image.meta == {'member': 0}
    for array, reference in zip(arrays, original):
        np.testing.assert_array_equal(array, reference)
