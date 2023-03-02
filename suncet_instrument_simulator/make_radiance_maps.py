from suncet_instrument_simulator import config_parser

import os
import xarray as xr
import numpy as np
import scipy.io
from glob import glob
import sunpy.map
from astropy.io import fits

class MakeRadianceMaps:
    def __init__(self, config_filename=os.getcwd() + '/config_files/config_default.ini'):
        self.config = self.__read_config(config_filename)


    def run(self, save=False):
        self.emiss = self.__get_emissivity()
        self.em_map = self.__read_em_map()
        self.logt_axis, self.native_wave_axis, self.wave_axis = self.__parameter_setup()
        self.raw_radiance = self.__compute_radiance()
        self.map_seq = self.__construct_maps()
        if save:
            self.__save_map_sequence()
        return self.map_seq


    def __read_config(self, config_filename):
        return config_parser.Config(config_filename)


    def __get_emissivity(self):
        expected_filename = os.getenv('suncet_data') + '/ancillary/emissivity/aia_V9_fullemiss.nc'
        emissivity_filename = glob(expected_filename)
        if len(emissivity_filename) > 0:
            return xr.load_dataset(emissivity_filename[0])
        else:
            raise ValueError('Emissivity file not found. Make sure that your emissivity file is saved as {}'.format(expected_filename))


    def __read_em_map(self):
        # TODO: Make this work on more than one file
        expected_em_filenames = os.getenv('suncet_data') + '/mhd/bright_fast/em_maps/em_map_240.sav'
        em_maps_filenames = glob(expected_em_filenames)
        if len(em_maps_filenames) > 0:
            em_data = scipy.io.readsav(em_maps_filenames[0])
            return em_data.em_maps_plus * 10**26. # EM Maps have units of 10**26 cm^-5
        else:
            raise ValueError('EM maps not found. Make sure that your EM maps are saved as {}'.format(expected_em_filenames))


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
        header['BUNIT'] = 'photon/(cm2 s sr Angstrom)'
        header['WCSNAME'] = 'Helioprojective-cartesian'
        header['HGLT_OBS'] = 0.0
        header['HGLN_OBS'] = 0.0
        header['DSUN_OBS'] = 149597870691
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
            header['WAVELNTH'] = self.wave_axis[n].value
            map_out = sunpy.map.Map(self.raw_radiance[n, :, :], header)
            if n == 0:
                map_seq = sunpy.map.MapSequence(map_out)
            else:
                map_seq.maps.append(map_out)
        return map_seq


    def __save_map_sequence(self):
        for n, map in enumerate(self.map_seq):
            if n == 0:
                hdu = fits.PrimaryHDU(map.data, map.fits_header)
                hdul = fits.HDUList(hdu)
            else:
                hdu = fits.ImageHDU(map.data, map.fits_header)
                hdul.append(hdu)
        # TODO: Make this work on more than one file
        hdul.writeto('SunCET_MapSeq.fits')


if __name__ == "__main__":
    radiance_map = MakeRadianceMaps()
    radiance_map.run(save = True)

