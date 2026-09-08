import astropy.units as u
import numpy as np
import pytest
from scipy.signal import fftconvolve

from suncet_instrument_simulator.Diffraction import suncet_psf


def _legacy_component(meshinfo, angles, orders, size, focal_plane=False):
    """Reference calculation over every pixel, including subnormal tails."""
    result = np.zeros((size, size), dtype=float)
    x = np.outer(np.ones(size), np.arange(size) + 0.5)
    y = np.outer(np.arange(size) + 0.5, np.ones(size))
    width = meshinfo['width'].value
    spacing = meshinfo['spacing_fp'] if focal_plane else meshinfo['spacing_e']
    mesh_ratio = (meshinfo['mesh_pitch'] / meshinfo['mesh_width']).decompose().value
    spacing_x = spacing * np.cos(angles)
    spacing_y = spacing * np.sin(angles)
    for order in orders:
        if order == 0:
            continue
        intensity = np.sinc(order / mesh_ratio) ** 2
        for dx, dy in zip(spacing_x.value, spacing_y.value):
            xc = x - (0.5 * size + dx * order + 0.5)
            yc = y - (0.5 * size + dy * order + 0.5)
            result += np.exp(-width * xc * xc - width * yc * yc) * intensity
    core = np.exp(-width * (x - 0.5 * size - 0.5) ** 2
                  - width * (y - 0.5 * size - 0.5) ** 2)
    area = meshinfo['Area']
    return (1 - area) * result / result.sum() + area * core / core.sum()


@pytest.mark.parametrize('wavelength,lpi', ((170, 20), (185, 65), (210, 20)))
@pytest.mark.parametrize('angles', ((0, 90), (17, 103)))
@pytest.mark.parametrize('focal_plane', (False, True))
def test_local_support_preserves_all_representable_gaussian_values(wavelength, lpi, angles, focal_plane):
    angles = np.array(angles) * u.deg
    meshinfo = suncet_psf.filter_mesh_parameters(wavelength, lpi=lpi)
    orders = np.array([-1500, -70, -11, -1, 0, 1, 7, 39, 1500])
    expected = _legacy_component(meshinfo, angles, orders, 96, focal_plane)
    actual = suncet_psf._psf(meshinfo, angles, orders, [96, 96], focal_plane, use_gpu=False)

    np.testing.assert_array_equal(actual, expected)


def test_local_support_retains_subnormal_tails():
    meshinfo = suncet_psf.filter_mesh_parameters(185)
    meshinfo['width'] = 1.0 * u.pixel
    angles = np.array([0, 90]) * u.deg
    orders = np.array([-1, 1])
    expected = _legacy_component(meshinfo, angles, orders, 96)
    actual = suncet_psf._psf(meshinfo, angles, orders, [96, 96], use_gpu=False)

    assert np.any((expected > 0) & (expected < np.finfo(np.float64).tiny))
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize('mesh_pitch', (0.0, np.nan))
def test_nonfinite_intensity_retains_dense_nan_propagation(mesh_pitch):
    meshinfo = suncet_psf.filter_mesh_parameters(185)
    meshinfo['mesh_pitch'] = mesh_pitch * u.um
    angles = np.array([0, 90]) * u.deg
    orders = np.array([-4000, 4000])
    with np.errstate(all='ignore'):
        expected = _legacy_component(meshinfo, angles, orders, 32)
        actual = suncet_psf._psf(meshinfo, angles, orders, [32, 32], use_gpu=False)

    assert np.all(np.isnan(expected))
    np.testing.assert_array_equal(actual, expected)


def test_composite_psf_energy_wings_and_convolved_image_are_unchanged():
    size = 96
    orders = np.arange(-20, 21)
    entrance_angles = np.array([0, 90]) * u.deg
    focal_angles = np.array([13, 103]) * u.deg
    meshinfo = suncet_psf.filter_mesh_parameters(185, angle_arm=entrance_angles,
                                               angles_focal_plane=focal_angles)
    entrance = _legacy_component(meshinfo, entrance_angles, orders, size)
    focal = _legacy_component(meshinfo, focal_angles, orders, size, focal_plane=True)
    expected = abs(np.fft.fft2(np.fft.fft2(focal) * np.fft.fft2(entrance)))
    expected = np.roll(np.roll(expected, size // 2, axis=1), size // 2, axis=0) / size**2
    actual = suncet_psf.psf(185, diffraction_orders=orders, angle_arm=entrance_angles,
                           angles_focal_plane=focal_angles, output_size=[size, size], use_gpu=False)
    np.testing.assert_array_equal(actual, expected)

    y, x = np.indices((size, size))
    radius = np.hypot(x - size // 2, y - size // 2)
    for r in (1, 3, 10, 30):
        np.testing.assert_array_equal(actual[radius <= r].sum(), expected[radius <= r].sum())
    np.testing.assert_array_equal(actual[radius > 30].sum(), expected[radius > 30].sum())
    scene = np.exp(-radius / 10)
    scene[size // 3, size // 2] += 10
    np.testing.assert_array_equal(fftconvolve(scene, actual, mode='same'),
                                  fftconvolve(scene, expected, mode='same'))
