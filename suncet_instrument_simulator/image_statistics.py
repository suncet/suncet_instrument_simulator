"""Local image statistics used by the simulator's optional SNR diagnostic."""

import numpy as np
from scipy.ndimage import correlate


def _local_rms_reference(data, window_min, window_max):
    """Retain NumPy's original reduction behavior for exceptional dtypes/ranges."""
    height, width = data.shape
    result = np.zeros((height, width))
    for x in range(width):
        for y in range(height):
            window = data[
                max(0, y - window_min):min(height, y + window_max + 1),
                max(0, x - window_min):min(width, x + window_max + 1),
            ]
            result[y, x] = np.sqrt(np.mean([row ** 2 for row in window]))
    return result


def local_rms(data, window):
    """Calculate the simulator's existing local RMS without per-pixel Python.

    ``window`` retains the historical inclusive bounds: 3 selects offsets -1
    through +2 on each axis, and therefore a 4-by-4 neighborhood. Edge windows
    include only available pixels. Squaring retains the input dtype (including
    integer overflow), and the mean/square root retain NumPy's dtype policy
    before storage in the original float64 result buffer. Floating reduction
    order can differ by rounding error.
    """
    data = np.asanyarray(data)
    if data.ndim != 2:
        raise ValueError('local_rms requires a two-dimensional image')
    window_min = int(np.floor(window / 2))
    window_max = int(np.ceil(window / 2))
    if window_min < 0:
        raise ValueError('window must be nonnegative')
    if data.size == 0:
        return np.zeros(data.shape)
    if data.dtype.kind not in 'buif':
        return _local_rms_reference(data, window_min, window_max)

    size = window_min + window_max + 1
    squared = data ** 2
    mean_dtype = np.mean(np.zeros(1, dtype=squared.dtype)).dtype
    # NumPy accumulates float32 means in float32, so even finite squared values
    # can overflow during summation. Keep that behavior instead of silently
    # changing it through the float64 accumulator used by ndimage.
    if squared.dtype.kind == 'f' and squared.dtype.itemsize <= 4:
        accumulation_dtype = np.result_type(squared.dtype, np.float32)
        if np.any(squared > np.finfo(accumulation_dtype).max / (size * size)):
            return _local_rms_reference(data, window_min, window_max)

    # A direct box correlation avoids rolling-sum cancellation at large dynamic
    # ranges and confines NaN/Inf effects to the same windows as the old loop.
    sums = correlate(
        squared.astype(np.float64, copy=False), np.ones((size, size)),
        mode='constant', cval=0., origin=window_min - size // 2,
    )
    height, width = data.shape
    y = np.arange(height)
    x = np.arange(width)
    rows = np.minimum(height, y + window_max + 1) - np.maximum(0, y - window_min)
    cols = np.minimum(width, x + window_max + 1) - np.maximum(0, x - window_min)
    means = (sums / (rows[:, None] * cols[None, :])).astype(mean_dtype, copy=False)
    return np.sqrt(means).astype(np.float64, copy=False)
