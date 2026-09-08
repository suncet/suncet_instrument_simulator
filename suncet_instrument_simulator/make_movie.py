import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy as np
from astropy.io import fits
import os
from glob import glob
import imageio.v2 as imageio
import datetime

def apply_radial_filter(data, sigma):
    # Determine the center of the image
    xc, yc = data.shape[1] / 2, data.shape[0] / 2

    # Create a radial distance array
    Y, X = np.ogrid[:data.shape[0], :data.shape[1]]
    r = np.sqrt((X - xc)**2 + (Y - yc)**2)

    # Define a radial filter function, e.g., a Gaussian
    # Adjust sigma to control the spread of the Gaussian
    radial_filter = np.exp(-(r**2 / (2. * sigma**2)))

    # Apply the filter
    filtered_data = data * radial_filter

    return filtered_data


def replace_negative_values(data): 
    neg_indices = np.where(data < 0)

    for i, j in zip(*neg_indices):
        # Get neighboring indices
        neighbors = data[max(0, i-1):i+2, max(0, j-1):j+2]

        # Calculate mean of neighbors, excluding the negative value itself
        mean_val = np.mean(neighbors[neighbors >= 0])

        # Replace the negative value with the mean
        data[i, j] = mean_val

    return data


def _render_frame(data, *, vmin, vmax, cmap, output_filename=None):
    """Render the same borderless RGBA pixels used by the PNG movie workflow."""
    height, width = data.shape[:2]
    fig = plt.figure(frameon=False)
    fig.set_size_inches(width / fig.dpi, height / fig.dpi)
    canvas = FigureCanvas(fig)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    # Match savefig(transparent=True), including masked/nonfinite image pixels.
    ax.patch.set_facecolor('none')
    ax.patch.set_edgecolor('none')
    fig.add_axes(ax)
    try:
        ax.imshow(data, vmin=vmin, vmax=vmax, cmap=cmap, aspect='auto')
        canvas.draw()
        frame = np.asarray(canvas.buffer_rgba()).copy()
        if output_filename is not None:
            imageio.imwrite(output_filename, frame)
        return frame
    finally:
        plt.close(fig)


def plot_difference_image(data, data_prior, output_filename=None):
    return _render_frame(
        data - data_prior, vmin=-10000, vmax=10000, cmap='gray',
        output_filename=output_filename,
    )


def _asinh_limits(images):
    """Match CME tracker movie.py: median sampled percentiles, 3% softening."""
    lows, highs = [], []
    for data in images:
        finite = np.asarray(data, dtype=np.float64)
        finite = finite[np.isfinite(finite)]
        if finite.size:
            low, high = np.percentile(finite, [1.0, 99.7])
            lows.append(float(low))
            highs.append(float(high))
    if not lows:
        raise ValueError('Movie input contains no finite image pixels.')
    low, high = float(np.median(lows)), float(np.median(highs))
    if not high > low:
        high = low + max(abs(low), 1.0) * np.finfo(np.float64).eps
    width = max((high - low) * .03, np.finfo(np.float64).eps)
    return low, high, width


def plot_scaled_image(data, output_filename=None, scale=None, asinh_limits=None):
    if scale == 'asinh':
        low, high, width = asinh_limits or _asinh_limits([data])
        return _render_frame(
            np.arcsinh((np.asarray(data, dtype=np.float64) - low) / width),
            vmin=0, vmax=float(np.arcsinh((high - low) / width)),
            cmap='inferno', output_filename=output_filename,
        )
    scale_funcs = {
        'log': lambda x: np.log10(np.clip(x, a_min=0.1, a_max=None)),
        'sqrt': np.sqrt,
        '1/4': lambda x: x**(1/4),
        '1/3': lambda x: x**(1/3),
        '1/8': lambda x: np.clip(x**(1/8), a_min=0.5, a_max=None),
    }
    scale_func = scale_funcs.get(scale, lambda x: x)  # Default to no scaling if not found

    return _render_frame(
        scale_func(data), vmin=0.08, vmax=21.0, cmap='inferno',
        output_filename=output_filename,
    )


# Configure script here
filenames = 'config_default_OBS_*.fits'
do_difference = False
scale = '1/4'  # Also supports 'asinh', matching the CME tracker display stretch.
# Enable only when individual PNGs are also wanted; movies stream from memory.
save_png_frames = False


def make_movie(fits_files, movie_filename, *, do_difference=False, png_directory=None, scale='1/4'):
    """Stream FITS images into a movie, optionally retaining individual PNGs."""
    fits_files = list(fits_files)
    asinh_limits = None
    if scale == 'asinh' and not do_difference:
        indices = np.unique(np.linspace(0, len(fits_files) - 1, min(12, len(fits_files))).round().astype(int))

        def sample_images():
            for index in indices:
                with fits.open(fits_files[index]) as hdul:
                    yield apply_radial_filter(replace_negative_values(hdul[0].data), 300)

        asinh_limits = _asinh_limits(sample_images())
    data_prior = None
    with imageio.get_writer(movie_filename, fps=20) as writer:
        for file in fits_files:
            with fits.open(file) as hdul:
                unfiltered_data = replace_negative_values(hdul[0].data)
                data = unfiltered_data
                # Preserve the historical Boolean/difference behavior. Correcting
                # this expression changes image values and is a separate change.
                if ~do_difference:
                    data = apply_radial_filter(data, 300)

            output_filename = None
            if png_directory is not None:
                # Keep the existing PNG naming convention.
                basename = str(file).split('.')[1]
                suffix = '_difference.png' if do_difference else '.png'
                output_filename = os.path.join(png_directory, basename + suffix)

            if do_difference:
                if data_prior is not None:
                    frame = plot_difference_image(data, data_prior, output_filename)
                    writer.append_data(frame)
                # The previous frame was historically reread without its radial
                # filter. Retain that same array instead of reading it a second time.
                data_prior = unfiltered_data
            else:
                frame = plot_scaled_image(data, output_filename, scale=scale, asinh_limits=asinh_limits)
                writer.append_data(frame)


def main():
    data_root = os.getenv('suncet_data')
    path = data_root + '/synthetic/level0/fits/'
    output_directory = data_root + '/synthetic/images and movies/'
    movie_filename = output_directory + 'synthetic_suncet_movie_' + datetime.datetime.now().strftime('%Y-%m-%d')
    movie_filename += '_difference.mp4' if do_difference else '.mp4'
    make_movie(
        sorted(glob(path + filenames)), movie_filename,
        do_difference=do_difference, scale=scale,
        png_directory=output_directory if save_png_frames else None,
    )


if __name__ == '__main__':
    main()
