from suncet_instrument_simulator import config_parser

import os
import xarray as xr
import numpy as np
import scipy.io
from glob import glob
import sunpy.map
from astropy.io import fits
import re
import warnings
import astropy.units as u

class MakeRadianceMaps:
    def __init__(self, config = None):
        if config is None:
            config_filename = os.getcwd() + '/config_files/config_default.ini'
            self.config = self.__read_config(config_filename)
        else:
            self.config = config


    def run(self):
        self.emiss = self.__get_emissivity()
        self.files = self.__parse_filenames()

        if len(self.files) == 0:
            warnings.warn('No files to process in specified range. Stopping.')
            return

        for self.em_maps_filename in self.files:
            self.em_map = self.__read_em_map()
            self.logt_axis, self.native_wave_axis, self.wave_axis = self.__parameter_setup()
            self.raw_radiance = self.__compute_radiance()
            self.map_seq = self.__construct_maps()
            self.__save_map_sequence()

    def __read_config(self, config_filename):
        return config_parser.Config(config_filename)


    def __get_emissivity(self):
        expected_filename = os.getenv('suncet_data') + '/ancillary/emissivity/aia_V9_fullemiss.nc'
        emissivity_filename = glob(expected_filename)
        if len(emissivity_filename) > 0:
            return xr.load_dataset(emissivity_filename[0])
        else:
            raise ValueError('Emissivity file not found. Make sure that your emissivity file is saved as {}'.format(expected_filename))

    def __parse_filenames(self):
        expected_em_filenames = os.getenv('suncet_data') + self.config.model_directory_name + '/' + self.config.em_directory_name + 'em_map_*.sav'
        em_maps_filename = glob(expected_em_filenames)
        file_number_strings = [re.match('.*(\d\d\d).*', filename).group(1) for filename in em_maps_filename]
        file_numbers = np.array(list(map(int, file_number_strings)))
        n_files_to_process = np.floor(
            (self.config.timesteps_to_process[1] + 1 - self.config.timesteps_to_process[0]) / self.config.timesteps_to_process[2])

        files_to_process = []
        for n, number in enumerate(range(self.config.timesteps_to_process[0], self.config.timesteps_to_process[1] + 1, self.config.timesteps_to_process[2])):
            result = np.where(file_numbers == number)[0]
            if result.size != 0:
                files_to_process.append(result[0])
        if len(files_to_process) != 0:
            filenames = np.take(em_maps_filename, files_to_process)
        else:
            filenames = []
        return filenames

    def __read_em_map(self):
        em_data = scipy.io.readsav(self.em_maps_filename)
        return em_data.em_maps_plus * 10**26. # EM Maps have units of 10**26 cm^-5

    def __parameter_setup(self):
        wavelength_units = self.config.wavelength_limits.unit
        limits = self.config.wavelength_limits.value
        self.binsize = self.config.wavelength_bin.value
        self.native_binsize = np.round(np.median(self.emiss.wave.values[1:] - self.emiss.wave.values[0:-1]), 1)

        # This stuff is fixed for all of Meng's simulations, but would be better to make configurable
        lgTmin = 5.45  # minimum for lgT axis for inversion
        dlgT = 0.1  # width of lgT bin
        nlgT = 12  # number of lgT bins

        # LogT Axis, Native Wavelength Axis for Emissivity, Desired Wavelength Axis for Radiances
        return np.arange(lgTmin, lgTmin + ((nlgT - 1) * dlgT), dlgT), \
               np.arange(limits[0], limits[1], self.native_binsize) * wavelength_units, np.arange(limits[0], limits[1], self.binsize) * wavelength_units


    def __compute_radiance(self):
        # setup necessary variables
        em_map_dim = np.shape(self.em_map)
        native_wave_dims = np.shape(self.native_wave_axis)[0]
        suncet_emissivity = self.emiss.total.sel(logte=self.logt_axis, wave=self.native_wave_axis, method='nearest')

        # generate the array to receive full radiance cube
        emiss_wave_array_full = np.empty([native_wave_dims, em_map_dim[1], em_map_dim[2]], dtype=np.float32)

        # compute radiance from emissivity
        for x in range(em_map_dim[1]):
            for y in range(em_map_dim[2]):
                emiss_wave_array_full[:, x, y] = np.matmul(self.em_map[:, x, y], suncet_emissivity.values)

        # how much are we binning in spectral space?
        bin_ratio = round(self.binsize/self.native_binsize, 1)

        # rebin the data with proper accounting (i.e. accounting for the 1/Å in native binsize)
        emiss_wave_array = np.empty([int(native_wave_dims/bin_ratio), em_map_dim[1], em_map_dim[2]], dtype=np.float32)
        for n in range(int(native_wave_dims/bin_ratio)):
            emiss_wave_array[n, :, :] = np.sum(emiss_wave_array_full[int(n * bin_ratio):int(n * bin_ratio + (bin_ratio - 1)), :, :]
                                               * self.native_binsize, axis=0)
        return emiss_wave_array


    def __make_header_template(self):
        # TODO: figure out what format SunPy maps want the solar radius keyword
        header = {}
        header['DATE-OBS'] = '2023-02-14T17:00:00.000'
        header['CTYPE1'] = 'HPLN-TAN'
        header['CTYPE2'] = 'HPLT-TAN'
        header['CUNIT1'] = 'arcsec  '
        header['CUNIT2'] = 'arcsec  '
        header['CRVAL1'] = 0.0
        header['CRVAL2'] = 0.0
        header['LONPOLE'] = 180.0
        header['CRPIX1'] = 512.5
        header['CRPIX2'] = 512.5
        header['CDELT1'] = 10.5
        header['CDELT2'] = 10.5
        header['CROTA2'] = 0.0
        header['PC1_1'] = 1.0
        header['PC1_2'] = 0.0
        header['PC2_1'] = 0.0
        header['PC2_2'] = 1.0
        header['WAVELNTH'] = 170.0 # NOTE: FITS doesn't have a keyword to store the units for wavelength, it's usually listed in the header comment instead
        header['WAVEUNIT'] = 'Angstrom'
        header['BUNIT'] = 'photon/(cm2 s sr Angstrom)'
        header['WCSNAME'] = 'Helioprojective-cartesian'
        header['HGLT_OBS'] = 0.0
        header['HGLN_OBS'] = 0.0
        header['DSUN_OBS'] = 149597870691
        header['RSUN'] = 970.
        header['TELESCOP'] = 'SunCET'
        header['INSTRUME'] = 'SunCET'
        header['EXPTIME'] = 1.0
        return header
    

    def __construct_maps(self):
        header_template = self.__make_header_template()
        map_dims = np.shape(self.raw_radiance)

        # for each wavelength bin, make a header, and turn the data/header into a map
        # compile maps into a Map Sequence
        for n in range(map_dims[0]):
            header = header_template
            header['WAVELNTH'] = self.wave_axis[n].value # TODO: Make this more elegant
            map_out = sunpy.map.Map(self.raw_radiance[n, :, :], header)
            if n == 0:
                map_seq = sunpy.map.MapSequence(map_out)
            else:
                map_seq.maps.append(map_out)
        return map_seq

    def __solve_solar_params(self, header = None):
        date = header['DATE-OBS']
        solar_distance = int(sunpy.coordinates.sun.earth_distance(time = date).to(u.m).value)
        solar_radius = sunpy.coordinates.sun.angular_radius(t=date)
        return (solar_distance, solar_radius)

    def __save_map_sequence(self):
        map_file_out = self.__make_outgoing_filename()
        for n, map in enumerate(self.map_seq):
            solar_params = self.__solve_solar_params(header=map.fits_header)
            if n == 0:
                hdu = fits.PrimaryHDU(map.data, map.fits_header)
                hdu.header['DSUN_OBS'] = solar_params[0]
                hdu.header['RSUN'] = solar_params[1].value
                hdul = fits.HDUList(hdu)
            else:
                hdu = fits.ImageHDU(map.data, map.fits_header)
                hdu.header['DSUN_OBS'] = solar_params[0]
                hdu.header['RSUN'] = solar_params[1].value
                hdul.append(hdu)
        print('Saving Map: ' + map_file_out)
        try:
            hdul.writeto(map_file_out)
        except OSError:
            error_msg = "Couldn't write {}. Might be because one with the same name already existed.".format(map_file_out)
            warnings.warn(error_msg)


    def __make_outgoing_filename(self):
        file_number_string = re.match('.*(\d\d\d).*', self.em_maps_filename).group(1)
        print(self.em_maps_filename)
        map_file_out = os.getenv('suncet_data') + self.config.model_directory_name + '/' + self.config.map_directory_name + '/radiance_maps_' + file_number_string + '.fits'
        return map_file_out

if __name__ == "__main__":
    radiance_map = MakeRadianceMaps()
    radiance_map.run()

