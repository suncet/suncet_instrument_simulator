from types import SimpleNamespace
import warnings

import numpy as np
import pytest
import xarray as xr

from suncet_instrument_simulator import make_radiance_maps


def _legacy_radiance(em_map, emissivity, native_binsize, binsize):
    """Independent reference retaining the original per-pixel rounding stages."""
    native = np.empty((emissivity.shape[1], *em_map.shape[1:]), dtype=np.float32)
    for x in range(em_map.shape[1]):
        for y in range(em_map.shape[2]):
            native[:, x, y] = np.matmul(em_map[:, x, y], emissivity)
    ratio = round(binsize / native_binsize, 1)
    result = np.empty((int(native.shape[0] / ratio), *em_map.shape[1:]), dtype=np.float32)
    for n in range(result.shape[0]):
        result[n] = np.sum(
            native[int(n * ratio):int(n * ratio + (ratio - 1))] * native_binsize,
            axis=0,
        )
    return result


@pytest.mark.parametrize('dtype', (np.float32, np.float64))
@pytest.mark.parametrize('emissivity_dtype', (np.float32, np.float64))
@pytest.mark.parametrize('binsize', (0.1, 0.35, 1.0, 5.0))
@pytest.mark.parametrize('noncontiguous', (False, True))
def test_batched_radiance_retains_legacy_rounding_stages(monkeypatch, dtype, emissivity_dtype, binsize, noncontiguous):
    rng = np.random.default_rng(417)
    em_map = (10.0 ** rng.uniform(20, 29, size=(12, 9, 14))).astype(dtype)
    if noncontiguous:
        em_map = em_map[:, :, ::2]
    emissivity = (10.0 ** rng.uniform(-32, -18, size=(12, 40))).astype(emissivity_dtype)
    radiance = make_radiance_maps.MakeRadianceMaps(SimpleNamespace())
    radiance.em_map = em_map
    radiance.logt_axis = np.arange(12)
    radiance.native_wave_axis = np.arange(40)
    radiance.native_binsize = np.float64(0.1)
    radiance.binsize = binsize
    radiance.emiss = xr.Dataset(
        {'total': (('logte', 'wave'), emissivity)},
        coords={'logte': radiance.logt_axis, 'wave': radiance.native_wave_axis},
    )
    # Force several chunk boundaries and a partial final chunk.
    monkeypatch.setattr(make_radiance_maps, '_RADIANCE_PIXEL_CHUNK_SIZE', 11)

    expected = _legacy_radiance(em_map, emissivity, radiance.native_binsize, binsize)
    with np.errstate(all='raise'):
        actual = radiance._MakeRadianceMaps__compute_radiance()

    assert actual.dtype == np.float32
    assert np.all(np.isfinite(actual))
    if dtype == emissivity_dtype == np.float32 and actual.size:
        # Single-precision BLAS can round a batched reduction differently.
        # Bound the effect for physical (nonnegative) emission/emissivity data.
        np.testing.assert_array_max_ulp(actual, expected, maxulp=2)
    else:
        np.testing.assert_array_equal(actual, expected)


def test_batched_radiance_keeps_native_rounding_before_binning(monkeypatch):
    rng = np.random.default_rng(293)
    em_map = rng.uniform(0, 1, size=(12, 2, 3)).astype(np.float32)
    emissivity = rng.uniform(0, 1, size=(12, 20))
    radiance = make_radiance_maps.MakeRadianceMaps(SimpleNamespace())
    radiance.em_map = em_map
    radiance.logt_axis = np.arange(12)
    radiance.native_wave_axis = np.arange(20)
    radiance.native_binsize = np.float64(0.1)
    radiance.binsize = 1.0
    radiance.emiss = xr.Dataset(
        {'total': (('logte', 'wave'), emissivity)},
        coords={'logte': radiance.logt_axis, 'wave': radiance.native_wave_axis},
    )
    monkeypatch.setattr(make_radiance_maps, '_RADIANCE_PIXEL_CHUNK_SIZE', 4)

    expected = _legacy_radiance(em_map, emissivity, radiance.native_binsize, radiance.binsize)
    pre_binned = np.stack([emissivity[:, n * 10:n * 10 + 9].sum(axis=1) * 0.1 for n in range(2)])
    moved_rounding = (pre_binned @ em_map.reshape(12, -1)).reshape(expected.shape).astype(np.float32)
    assert np.any(expected != moved_rounding)
    np.testing.assert_array_equal(radiance._MakeRadianceMaps__compute_radiance(), expected)


def test_single_precision_emissivity_has_bounded_last_bit_roundoff():
    rng = np.random.default_rng(417)
    em_map = rng.uniform(0, 1, (12, 9, 14)).astype(np.float32)
    emissivity = rng.uniform(0, 1, (12, 600)).astype(np.float32)
    radiance = make_radiance_maps.MakeRadianceMaps(SimpleNamespace())
    radiance.em_map = em_map
    radiance.logt_axis = np.arange(12)
    radiance.native_wave_axis = np.arange(600)
    radiance.native_binsize = np.float64(0.1)
    radiance.binsize = 1.0
    radiance.emiss = xr.Dataset(
        {'total': (('logte', 'wave'), emissivity)},
        coords={'logte': radiance.logt_axis, 'wave': radiance.native_wave_axis},
    )
    # The legacy matmul may report BLAS status flags even for these inputs in
    # [0, 1]. Capture only that reference calculation and verify its results;
    # the optimized path must run with floating-point errors raised.
    with warnings.catch_warnings(record=True) as legacy_warnings:
        expected = _legacy_radiance(em_map, emissivity, radiance.native_binsize, radiance.binsize)
    assert np.all(np.isfinite(expected))
    expected_messages = {
        'divide by zero encountered in matmul',
        'overflow encountered in matmul',
        'invalid value encountered in matmul',
    }
    assert all(w.category is RuntimeWarning and str(w.message) in expected_messages for w in legacy_warnings)
    with np.errstate(all='raise'):
        actual = radiance._MakeRadianceMaps__compute_radiance()

    np.testing.assert_array_max_ulp(actual, expected, maxulp=2)


@pytest.mark.parametrize('emissivity_value', (1e30, 1e300))
def test_batched_radiance_still_reports_actual_overflow(emissivity_value):
    radiance = make_radiance_maps.MakeRadianceMaps(SimpleNamespace())
    radiance.em_map = np.full((12, 1, 1), 1e20)
    radiance.logt_axis = np.arange(12)
    radiance.native_wave_axis = np.arange(20)
    radiance.native_binsize = np.float64(0.1)
    radiance.binsize = 1.0
    radiance.emiss = xr.Dataset(
        {'total': (('logte', 'wave'), np.full((12, 20), emissivity_value))},
        coords={'logte': radiance.logt_axis, 'wave': radiance.native_wave_axis},
    )

    # Exercise both float32 cast overflow and overflow of the float64 product.
    with np.errstate(over='raise', invalid='ignore', divide='ignore'):
        with pytest.raises(FloatingPointError, match='overflow'):
            radiance._MakeRadianceMaps__compute_radiance()
