from types import SimpleNamespace

import astropy.units as u
import numpy as np

from suncet_instrument_simulator import instrument


class _MapStub:
    def __init__(self, data, meta=None):
        self.data = np.asanyarray(data)
        self.meta = {} if meta is None else meta

    def __mul__(self, factor):
        value = factor.value if hasattr(factor, 'value') else factor
        return _MapStub(self.data * value, self.meta)


def test_unsigned_16_bit_conversion_wraps_instead_of_clipping():
    values = np.array([-1, 0, 65535, 65536, 65537, 131071, 131072])

    actual = instrument._wrap_to_unsigned(values, 16)

    np.testing.assert_array_equal(actual, [0, 0, 65535, 0, 1, 65535, 0])
    assert actual.dtype == np.uint16


def test_detector_dn_conversion_wraps_at_16_bits(monkeypatch):
    monkeypatch.setattr(instrument.sunpy.map, 'Map', _MapStub)
    hardware = instrument.Hardware.__new__(instrument.Hardware)
    hardware.config = SimpleNamespace(
        detector_gain=1 * u.dN / u.electron,
        readout_bits=16 * u.bit,
    )
    values = np.array([[-1, 65535, 65536, 65537]], dtype=float)
    detector_images = {
        'short exposure': {0: _MapStub(values)},
        'long exposure': {0: _MapStub(values)},
    }

    converted = hardware.convert_to_dn(detector_images)

    np.testing.assert_array_equal(
        converted['short exposure'][0].data,
        [[0, 65535, 0, 1]],
    )


def test_particle_filter_uses_32_bit_buffer_without_16_bit_saturation(monkeypatch):
    monkeypatch.setattr(instrument.sunpy.map, 'Map', _MapStub)
    software = instrument.OnboardSoftware(SimpleNamespace())
    images = {
        exposure: {
            0: _MapStub([[60000]]),
            1: _MapStub([[60000]]),
            2: _MapStub([[60000]]),
        }
        for exposure in ('short exposure', 'long exposure')
    }

    filtered = software.filter_out_particle_hits(images)

    assert filtered['short exposure'].data[0, 0] == 120000
    assert filtered['short exposure'].data.dtype == np.uint32


def test_post_filter_right_shift_wraps_when_written_to_16_bits(monkeypatch):
    monkeypatch.setattr(instrument.sunpy.map, 'Map', _MapStub)
    software = instrument.OnboardSoftware(
        SimpleNamespace(
            num_shift_bits_32_to_16=2,
            readout_bits=16 * u.bit,
        )
    )
    image = _MapStub([[262140, 262144, 262148]])

    shifted = software.bit_shift_data(image)

    np.testing.assert_array_equal(shifted.data, [[65535, 0, 1]])
    assert shifted.data.dtype == np.uint16
