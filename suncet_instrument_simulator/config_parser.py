import os
import configparser
import json
import astropy.units as u


class Config:

    def __init__(self, filename):
        config = configparser.ConfigParser()
        print('Reading configuration file: {}'.format(filename))
        if os.path.isfile(filename): 
            config.read(filename)
        else: 
            raise Exception('Config file does not exist or path is incorrect.')

        # behavior
        self.compute_new_radiance_maps = config['behavior'].getboolean('compute_new_radiance_maps')
        self.apply_psf = config['behavior'].getboolean('apply_psf')
        self.apply_scattered_light_psf = config['behavior'].getboolean('apply_scattered_light_psf')
        self.subtract_dark = config['behavior'].getboolean('subtract_dark_onboard')
        self.timesteps_to_process = json.loads(config.get('behavior', 'timesteps_to_process'))
        self.model_directory_name = config['behavior']['model_directory_name']
        self.em_directory_name = config['behavior']['em_directory_name']
        self.map_directory_name = config['behavior']['map_directory_name']

        # limits
        self.wavelength_limits = json.loads(config.get('limits', 'wavelength_limits')) * u.Angstrom
        self.wavelength_bin = config['limits'].getfloat('wavelength_bin') * u.Angstrom

        # telescope
        self.entrance_aperture = config['telescope'].getfloat('entrance_aperture') * u.cm
        self.primary_mirror_truncation = config['telescope'].getfloat('primary_mirror_truncation') * u.cm**2
        self.secondary_mirror_obscuration_fraction = config['telescope'].getfloat('secondary_mirror_obscuration') 
        self.fov = json.loads(config.get('telescope', 'fov')) * u.deg
        self.mirror_coating_reflectivity_filename = config['telescope']['mirror_coating_reflectivity_filename']
        self.filter_entrance_transmission_filename = config['telescope']['filter_entrance_transmission_filename']
        self.filter_focal_plane_transmission_filename = config['telescope']['filter_focal_plane_transmission_filename']

        # detector
        self.detector_physical_area = config['detector'].getfloat('detector_physical_area') * u.cm**2
        self.image_dimensions = json.loads(config.get('detector', 'image_dimensions')) * u.pix
        self.plate_scale = config['detector'].getfloat('plate_scale') * u.arcsec / u.pix
        self.quantum_efficiency_filename = config['detector']['quantum_efficiency_filename']
        self.read_noise = config['detector'].getfloat('read_noise') * u.electron
        self.pixel_full_well = config['detector'].getfloat('pixel_full_well') * u.electron
        self.num_pixels_to_bin = json.loads(config.get('detector', 'num_pixels_to_bin'))
        self.readout_bits = config['detector'].getint('readout_bits') * u.bit
        self.detector_temperature = config['detector'].getfloat('detector_temperature') * u.deg_C
        self.detector_gain = config['detector'].getfloat('detector_gain') * u.dn / u.electron
        self.fraction_dead_pixels = config['detector'].getfloat('fraction_dead_pixels')
        self.fraction_hot_pixels = config['detector'].getfloat('fraction_hot_pixels')
        self.spike_rate = config['detector'].getfloat('spike_rate') / u.second / u.cm**2
        self.fano_factor = config['detector'].getfloat('fano_factor') * u.dimensionless_unscaled

        # shdr
        self.exposure_time_short = config['shdr'].getfloat('exposure_time_short') * u.second
        self.exposure_time_long = config['shdr'].getfloat('exposure_time_long') * u.second
        self.num_short_exposures_to_stack = config['shdr'].getint('num_short_exposures_to_stack')
        self.num_long_exposures_to_stack = config['shdr'].getint('num_long_exposures_to_stack')
        self.num_inner_rows_for_short_exposure = config['shdr'].getint('num_inner_rows_for_short_exposure')
        self.inner_fov_circle_radius = config['shdr'].getfloat('inner_fov_circle_radius') * u.solRad
        self.jitter = config['shdr'].getfloat('jitter') * u.arcsec / u.second

        # structure
        self.base_metadata_filename = config['structure']['base_metadata_filename']
        self.model_data_folder = config['structure']['model_data_folder']



if __name__ == "__main__":
    filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    config = Config(filename)
    print(config.apply_psf)  # Just an example to see something return when running as a script