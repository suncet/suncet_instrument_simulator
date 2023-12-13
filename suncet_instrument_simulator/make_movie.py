import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import os
from glob import glob
import imageio

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


def plot_scaled_image(data, title, output_filename, scale=None):
    plt.figure(figsize=(24, 24))
    scale_funcs = {
        'log': lambda x: np.log10(np.clip(x, 0.1, None)),
        'sqrt': np.sqrt,
        '1/4': lambda x: x**(1/4),
        '1/8': lambda x: np.clip(x**(1/8), 0.5, None),
    }
    scale_func = scale_funcs.get(scale, lambda x: x)  # Default to no scaling if not found
    
    plt.imshow(scale_func(data), cmap='inferno')
    plt.axis('off')
    plt.savefig(output_filename, bbox_inches='tight', pad_inches=0)
    plt.close()

path = os.getenv('suncet_data') + '/synthetic/level0_raw/fits/'
filenames = 'config_default_OBS_2023-02-14T17:00:00.000_*.fits'
fits_files = sorted(glob(path + filenames))
image_files = []

for file in fits_files:
    with fits.open(file) as hdul:
        data = hdul[0].data
        filtered_data = apply_radial_filter(data, 300)

    output_filename = os.getenv('suncet_data')+ '/synthetic/images and movies/' + f"{file.split('.')[1]}.png"
    plot_scaled_image(filtered_data, '', output_filename, scale='1/4')
    image_files.append(output_filename)

with imageio.get_writer(os.getenv('suncet_data')+ '/synthetic/images and movies/movie.mp4', fps=10) as writer:
    for filename in image_files:
        image = imageio.imread(filename)
        writer.append_data(image)
