import numpy as np
import pytest

from suncet_instrument_simulator.image_statistics import local_rms


def _original_rms(data, window):
    height, width = data.shape
    before = int(np.floor(window / 2))
    after = int(np.ceil(window / 2))
    result = np.zeros(data.shape)
    for x in range(width):
        for y in range(height):
            values = data[max(0, y-before):min(height, y+after+1),
                          max(0, x-before):min(width, x+after+1)]
            result[y, x] = np.sqrt(np.mean([value ** 2 for value in values]))
    return result


@pytest.mark.parametrize('window', [0, 1, 2, 3, 4, 9])
@pytest.mark.parametrize('dtype', [np.uint16, np.uint32, np.int32, np.float16, np.float32, np.float64])
def test_local_rms_preserves_window_edges_dtype_and_overflow(window, dtype):
    rng = np.random.default_rng(162)
    if np.issubdtype(dtype, np.integer):
        data = rng.integers(0, 65000, size=(7, 11)).astype(dtype)
    else:
        data = rng.normal(size=(7, 11)).astype(dtype)
    with np.errstate(invalid='ignore'):
        expected = _original_rms(data, window)
        actual = local_rms(data, window)
    if np.issubdtype(dtype, np.floating):
        tolerance = 4 * np.finfo(dtype).eps
        np.testing.assert_allclose(actual, expected, rtol=tolerance, atol=0., equal_nan=True)
        if dtype == np.float32:
            np.testing.assert_array_max_ulp(actual.astype(dtype), expected.astype(dtype), maxulp=4)
    else:
        np.testing.assert_array_equal(actual, expected)
    assert actual.dtype == np.float64


@pytest.mark.parametrize('dtype', [np.float32, np.float64])
def test_nonfinite_values_only_affect_overlapping_windows(dtype):
    data = np.ones((10, 12), dtype=dtype)
    data[1, 2] = np.nan
    data[7, 10] = np.inf
    expected = _original_rms(data, 3)
    np.testing.assert_array_equal(local_rms(data, 3), expected)


def test_zero_windows_and_small_signal_survive_large_dynamic_range():
    data = np.ones((12, 15))
    data[0, 0] = 1e100
    data[8:, 8:] = 0.
    expected = _original_rms(data, 3)
    np.testing.assert_allclose(local_rms(data, 3), expected, rtol=1e-15, atol=0.)


def test_float32_reduction_overflow_is_preserved():
    data = np.full((5, 7), np.sqrt(np.finfo(np.float32).max / 2), dtype=np.float32)
    with np.errstate(over='ignore'):
        np.testing.assert_array_equal(local_rms(data, 3), _original_rms(data, 3))


def test_single_pixel_and_noncontiguous_input():
    np.testing.assert_array_equal(local_rms(np.array([[3.]]), 9), [[3.]])
    data = np.arange(100, dtype=np.uint32).reshape(10, 10)[::2, ::2]
    np.testing.assert_array_equal(local_rms(data, 3), _original_rms(data, 3))
