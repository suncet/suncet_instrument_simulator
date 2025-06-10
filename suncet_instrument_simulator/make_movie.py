import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy as np
from astropy.io import fits
import os
from glob import glob
import imageio
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


def plot_difference_image(data,data_prior, output_filename):
    height, width = data.shape[:2]
    
    # Matplotlib for some reason makes saving an image to disk with the exact right dimensions and no white space non-trivial, hence all this elaborate setup
    fig = plt.figure(frameon=False)
    fig.set_size_inches(width / fig.dpi, height / fig.dpi)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    diff = data - data_prior
    ax.imshow(diff, vmin=-10000, vmax=10000, cmap='gray', aspect='auto')
    plt.savefig(output_filename, dpi=fig.dpi, transparent=True)
    plt.close()


def plot_scaled_image(data, output_filename, scale=None):
    height, width = data.shape[:2]
    
    # Matplotlib for some reason makes saving an image to disk with the exact right dimensions and no white space non-trivial, hence all this elaborate setup
    fig = plt.figure(frameon=False)
    fig.set_size_inches(width / fig.dpi, height / fig.dpi)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    
    scale_funcs = {
        'log': lambda x: np.log10(np.clip(x, a_min=0.1, a_max=None)),
        'sqrt': np.sqrt,
        '1/4': lambda x: x**(1/4),
        '1/3': lambda x: x**(1/3),
        '1/8': lambda x: np.clip(x**(1/8), a_min=0.5, a_max=None),
    }
    scale_func = scale_funcs.get(scale, lambda x: x)  # Default to no scaling if not found

    ax.imshow(scale_func(data), vmin=0.08, vmax=21.0, cmap='inferno', aspect='auto')
    plt.savefig(output_filename, dpi=fig.dpi, transparent=True)
    plt.close()


# Configure script here
path = os.getenv('suncet_data') + '/synthetic/level0/fits/'
filenames = 'config_default_OBS_2023-02-14T17:00:00.000_*.fits'
do_difference = False


fits_files = sorted(glob(path + filenames))
image_files = []

for i, file in enumerate(fits_files):
    with fits.open(file) as hdul:
        data = hdul[0].data
        data = replace_negative_values(data)

        if ~do_difference:
            data = apply_radial_filter(data, 300)
    
    if do_difference: 
        if i > 0: 
            with fits.open(fits_files[i-1]) as hdul: 
                data_prior = hdul[0].data
                data_prior = replace_negative_values(data_prior)
        
    output_filename = os.getenv('suncet_data')+ '/synthetic/images and movies/' + f"{file.split('.')[1]}"
    if do_difference: 
        if i > 0: 
            output_filename += f"_difference.png"
            plot_difference_image(data, data_prior, output_filename)
            image_files.append(output_filename)
    else: 
        output_filename += f".png"
        plot_scaled_image(data, output_filename, scale='1/4')
        image_files.append(output_filename)
    
    

movie_filename = os.getenv('suncet_data')+ '/synthetic/images and movies/synthetic_suncet_movie_' + datetime.datetime.now().strftime('%Y-%m-%d')
if do_difference: 
    movie_filename += f"_difference.mp4"
else: 
    movie_filename += f".mp4"

with imageio.get_writer(movie_filename, fps=20) as writer:
    for filename in image_files:
        image = imageio.imread(filename)
        writer.append_data(image)
