import importlib

from astropy.io import fits
import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
import pytest

from suncet_instrument_simulator import make_movie


def _legacy_png(data, filename, *, vmin, vmax, cmap):
    """The previous savefig path is the pixel-level compatibility reference."""
    height, width = data.shape
    fig = plt.figure(frameon=False)
    fig.set_size_inches(width / fig.dpi, height / fig.dpi)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(data, vmin=vmin, vmax=vmax, cmap=cmap, aspect='auto')
    fig.savefig(filename, dpi=fig.dpi, transparent=True)
    plt.close(fig)
    return imageio.imread(filename)


@pytest.mark.parametrize('scale', [None, 'log', 'sqrt', '1/4', '1/3', '1/8'])
def test_canvas_matches_legacy_png_pixels_and_optional_png(tmp_path, scale):
    data = np.tile(np.array([[0., .1, 1., 100.], [np.nan, np.inf, 4095., 16.]]), (8, 9))
    transforms = {
        None: lambda x: x,
        'log': lambda x: np.log10(np.clip(x, .1, None)),
        'sqrt': np.sqrt,
        '1/4': lambda x: x ** (1/4),
        '1/3': lambda x: x ** (1/3),
        '1/8': lambda x: np.clip(x ** (1/8), .5, None),
    }
    expected = _legacy_png(
        transforms[scale](data), tmp_path / 'reference.png',
        vmin=.08, vmax=21., cmap='inferno',
    )
    actual = make_movie.plot_scaled_image(data, scale=scale)
    saved = make_movie.plot_scaled_image(data, tmp_path / 'optional.png', scale=scale)

    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(saved, expected)
    np.testing.assert_array_equal(imageio.imread(tmp_path / 'optional.png'), expected)
    assert actual.shape == (*data.shape, 4)
    assert actual.dtype == np.uint8


def test_difference_canvas_matches_legacy_unsigned_arithmetic(tmp_path):
    data = np.tile(np.array([[0, 1, 65535], [100, 30000, 60000]], dtype=np.uint16), (8, 9))
    prior = np.full_like(data, 100)
    expected = _legacy_png(
        data - prior, tmp_path / 'reference.png',
        vmin=-10000, vmax=10000, cmap='gray',
    )

    np.testing.assert_array_equal(make_movie.plot_difference_image(data, prior), expected)


@pytest.mark.parametrize('difference', [False, True])
@pytest.mark.parametrize('save_png', [False, True])
def test_movie_streams_frames_without_png_reads_and_reads_each_fits_once(
        monkeypatch, tmp_path, difference, save_png):
    files = []
    data_frames = []
    for index in range(3):
        data = np.arange(16 * 20, dtype=float).reshape(16, 20) * (index + 1)
        data[2, 3] = -1
        data_frames.append(data)
        filename = tmp_path / f'config_default_OBS_2023-01-14T17:00:00.000_{index:03}.fits'
        fits.PrimaryHDU(data).writeto(filename)
        files.append(filename)

    corrected = [make_movie.replace_negative_values(data.copy()) for data in data_frames]
    expected = []
    for index, data in enumerate(corrected):
        # Preserve the historical ~do_difference behavior: radial filtering is
        # also applied to the current difference frame, but not its predecessor.
        radial = make_movie.apply_radial_filter(data, 300)
        if difference:
            if index == 0:
                continue
            display = radial - corrected[index - 1]
            limits = {'vmin': -10000, 'vmax': 10000, 'cmap': 'gray'}
        else:
            display = radial ** (1/4)
            limits = {'vmin': .08, 'vmax': 21., 'cmap': 'inferno'}
        expected.append(_legacy_png(display, tmp_path / f'reference_{index}.png', **limits))

    frames = []
    reads = []
    original_open = fits.open

    def tracked_open(filename, *args, **kwargs):
        reads.append(filename)
        return original_open(filename, *args, **kwargs)

    class Writer:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def append_data(self, frame):
            frames.append(frame.copy())

    def get_writer(filename, fps):
        assert filename == tmp_path / 'movie.mp4'
        assert fps == 20
        return Writer()

    def reject_png_read(*args, **kwargs):
        pytest.fail('Movie generation should not reread intermediate PNG files')

    monkeypatch.setattr(fits, 'open', tracked_open)
    monkeypatch.setattr(make_movie.imageio, 'get_writer', get_writer)
    monkeypatch.setattr(make_movie.imageio, 'imread', reject_png_read)
    png_directory = tmp_path / 'png'
    png_directory.mkdir()

    make_movie.make_movie(
        files, tmp_path / 'movie.mp4', do_difference=difference,
        png_directory=png_directory if save_png else None,
    )

    assert reads == files
    assert len(frames) == len(expected)
    for frame, reference in zip(frames, expected):
        np.testing.assert_array_equal(frame, reference)
    assert len(list(png_directory.glob('*.png'))) == (len(expected) if save_png else 0)


def test_import_does_not_require_data_root_or_start_a_movie(monkeypatch):
    monkeypatch.delenv('suncet_data', raising=False)

    def reject_writer(*args, **kwargs):
        pytest.fail('Importing movie helpers must not run movie generation')

    monkeypatch.setattr(make_movie.imageio, 'get_writer', reject_writer)
    importlib.reload(make_movie)


def test_asinh_matches_tracker_stretch_and_keeps_fixed_limits(tmp_path):
    samples = [np.arange(1000, dtype=float), np.arange(1000, dtype=float) * 2]
    low, high, width = make_movie._asinh_limits(samples)
    np.testing.assert_allclose([low, high], [14.985, 1494.0045])
    assert width == pytest.approx((high - low) * .03)
    limits = (low, high, width)
    for data in [np.ones((16, 20)) * 30, np.ones((16, 20)) * 300]:
        expected = _legacy_png(
            np.arcsinh((data - low) / width), tmp_path / 'asinh.png',
            vmin=0, vmax=np.arcsinh((high - low) / width), cmap='inferno',
        )
        np.testing.assert_array_equal(
            make_movie.plot_scaled_image(data, scale='asinh', asinh_limits=limits),
            expected,
        )
