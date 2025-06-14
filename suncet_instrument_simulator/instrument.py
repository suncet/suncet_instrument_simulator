import os
import copy
import warnings
from glob import glob
import numpy as np
from scipy.integrate import simpson
from scipy.ndimage import shift, gaussian_filter
from pandas import read_fwf, read_csv
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy import constants as const
import astropy.io.fits as fits
from astropy.convolution import convolve_fft
import sunpy.map
from pillow_jpls import Image
from sunpy.coordinates import frames
from suncet_instrument_simulator import config_parser  # This is just for running as a script -- should delete when done testing


class Hardware:
    def __init__(self, config):
        self.config = config
        self.coating_name_mirror1 = self.__get_coating_name(1)
        self.coating_name_mirror2 = self.__get_coating_name(2)
        if self.config.apply_mirror_scattered_light_psf:
            self.mirror_scatter_psf = self.__load_mirror_scatter_psf()
        if self.config.apply_mesh_diffraction:
            self.mesh_diffraction_psf = self.__load_mesh_diffraction_psf()

    def __get_coating_name(self, mirror_number):
        if mirror_number == 1:
            return os.path.basename(self.config.mirror1_coating_reflectivity_filename).split('_1')[0]
        elif mirror_number == 2:
            return os.path.basename(self.config.mirror2_coating_reflectivity_filename).split('_1')[0]
        else:
            raise ValueError('Need to provide a mirror number to get the corresponding coating name.')

    def __load_mirror_scatter_psf(self):
        filename = self.config.mirror_scatter_filename
        full_path = os.getenv('suncet_data') + '/mirror_scatter/' + filename
        psf = fits.open(full_path)[0].data
        psf_crop = psf[0:-1, 0:-1]  # because convolve_fft requires us to have an odd number of rows
        return psf_crop

    def __load_mesh_diffraction_psf(self):
        filename = self.config.filter_diffraction_filename
        full_path = os.getenv('suncet_data') + '/filter_mesh_diffraction/' + filename
        diffraction_psf = fits.open(full_path)
        return diffraction_psf

    def store_target_wavelengths(self, radiance_maps):  # Will interpolate/calculate any subsequent wavelength-dependent quantities to this "target" wavelength array
        if self.__is_nested_time_and_wavelength(radiance_maps):
            radiance_maps = next(iter(radiance_maps.values()))
        values = [u.Quantity(key).value for key in radiance_maps.keys()]
        unit = u.Quantity(list(radiance_maps.keys())[0]).unit
        self.wavelengths = u.Quantity(values, unit=unit)

    def __is_nested_time_and_wavelength(self, d):
        if not isinstance(d, dict):
            return False
        return any(isinstance(value, dict) for value in d.values())

    def compute_effective_area(self):
        geometric_area = np.pi * ((self.config.entrance_aperture) / 2) ** 2
        geometric_area *= (1 - self.config.secondary_mirror_obscuration_fraction)
        mirror_reflectivity = self.__interpolate_mirror_coating_reflectivity()  # This already multiplies the reflectivities together for each mirror
        filter_transmission = self.__interpolate_filter_transmission()  # This already multiplies the transmission together for each filter
        quantum_efficiency = self.__interpolate_quantum_efficiency()

        wavelength_strings = [f"{wavelength.value} {wavelength.unit}" for wavelength in self.wavelengths]
        effective_area = geometric_area * mirror_reflectivity * filter_transmission * quantum_efficiency
        self.effective_area = {wavelength_str: a_eff for wavelength_str, a_eff in zip(wavelength_strings, effective_area)}

    def __interpolate_mirror_coating_reflectivity(self):
        filename_mirror1 = self.config.mirror1_coating_reflectivity_filename
        filename_mirror2 = self.config.mirror2_coating_reflectivity_filename

        source_data1 = self.__load_mirror_data(filename_mirror1)
        source_data2 = self.__load_mirror_data(filename_mirror2)

        interpolated_data1 = np.interp(self.wavelengths.value, source_data1.wavelength.values, source_data1.reflectivity.values)
        interpolated_data2 = np.interp(self.wavelengths.value, source_data2.wavelength.values, source_data2.reflectivity.values)

        return interpolated_data1 * interpolated_data2
   

    # This is the function where we encapsulate all the messiness of filenames and what they mean
    def __load_mirror_data(self, filename):
        if 'measurement' in filename:
            data = read_csv(os.getenv('suncet_data') + '/mirror_reflectivity/' + filename, header=0)
        elif 'final' in filename:
            data = read_csv(os.getenv('suncet_data') + '/mirror_reflectivity/2024-03-21 rigaku measurements final/' + filename)
        else:
            data = read_fwf(os.getenv('suncet_data') + '/mirror_reflectivity/' + filename, skiprows=18, header=0, names=['wavelength', 'reflectivity'])

        if 'nm' in data.columns[0]:
            data[data.columns[0]] *= 10  # Converts nm to Angstrom
        elif '[A]' or 'Angstrom' in data.columns[0]:
            warnings.warn('Found indication that coating wavelength data is in angstroms so making the ASSUMPTION that this is actually true.')
        else:
            raise ValueError('Do not know what units the coating wavelength is in.')
        data.columns = ['wavelength', 'reflectivity']

        return data  # wavelength should be in Angstroms, reflectivity as a fraction

    def __interpolate_filter_transmission(self):
        filter_entrance = self.__load_filter_data(self.config.filter_entrance_transmission_filename)
        filter_focal_plane = self.__load_filter_data(self.config.filter_focal_plane_transmission_filename)
        combined_transmission = filter_entrance
        combined_transmission['transmission'] *= filter_focal_plane['transmission']
        return np.interp(self.wavelengths.value, combined_transmission['wavelength [Å]'].values, combined_transmission['transmission'].values)

    def __load_filter_data(self, filename):
        df = read_csv(os.getenv('suncet_data') + '/filter_transmission/' + filename, skiprows=2, usecols=[i for i in range(2)], names=['wavelength [nm]', 'transmission'])
        df['wavelength [nm]'] *= 10
        df.rename(columns={'wavelength [nm]': 'wavelength [Å]'}, inplace=True)
        return df

    def __interpolate_quantum_efficiency(self):
        quantum_efficiency = self.__load_quantum_efficiency_curve()
        return np.interp(self.wavelengths.value, quantum_efficiency['wavelength [Å]'].values, quantum_efficiency['qe'].values)

    def __load_quantum_efficiency_curve(self):
        df = read_csv(os.getenv('suncet_data') + '/quantum_efficiency/' + self.config.quantum_efficiency_filename, skiprows=1)
        df.columns = ['wavelength [Å]', 'qe']
        return df

    def extract_fov(self, radiance_maps):
        fov_half_angles = self.config.fov / 2.0

        # TODO: Raise warning if instrument FOV > model FOV

        radiance_maps_new_fov = {}
        for timestep, radiance_maps_by_wavelength in radiance_maps.items():
            radiance_maps_new_fov[timestep] = {}
            for wavelength, map in radiance_maps_by_wavelength.items():
                radiance_maps_new_fov[timestep][wavelength] = map.submap(
                    top_right=SkyCoord(fov_half_angles[0], fov_half_angles[1], frame=map.coordinate_frame),
                    bottom_left=SkyCoord(-fov_half_angles[0], -fov_half_angles[1], frame=map.coordinate_frame)
                )
        radiance_maps.update(radiance_maps_new_fov)
        return radiance_maps
    

    def interpolate_spatial_resolution(self, radiance_maps):
        radiance_maps_new_resolution = {}
        for timestep, radiance_maps_by_wavelength in radiance_maps.items():
            radiance_maps_new_resolution[timestep] = {}
            for wavelength, map in radiance_maps_by_wavelength.items():
                resampled_map = map.resample(self.config.image_dimensions, method='linear') * ((self.config.plate_scale / map.scale[0]).value) ** 2
                resampled_map.data[resampled_map.data < 0] = 0  # Clip any negative numbers that result from the interpolation above -- bottom out at 0
                resampled_map.meta['cdelt1'] = self.config.plate_scale.value  # Note 1: this is only needed because sunpy (v4.0.1) resample updates dimensions but NOT plate scale
                resampled_map.meta['cdelt2'] = self.config.plate_scale.value  # Note 2: there is also risk here because a) the user must be responsible in the config file to ensure the image_dimensions and plate_scale are compatible, and b) the units they input for plate_scale must be the same as those already in the map
                radiance_maps_new_resolution[timestep][wavelength] = resampled_map
        return radiance_maps_new_resolution
    
    
    def convert_steradians_to_pixels(self, radiance_maps):
        plate_scale_rad = self.config.plate_scale.to(u.rad / u.pix)
        pixel_area_steradians = (plate_scale_rad ** 2).to(u.steradian / (u.pix ** 2))
        return self.__apply_function_to_leaves(radiance_maps, lambda x: x * pixel_area_steradians)

    def apply_diffraction_psf(self, radiance_maps):
        return self.__apply_function_to_leaves(radiance_maps, self.__convolve_diffraction_psf)

    def __convolve_diffraction_psf(self, map):
        # todo: fix the bad hack to find the correct diffraction pattern
        target_psf = self.mesh_diffraction_psf[int(map.meta['wavelnth']) - 170].data
        target_psf_crop = target_psf[0:-1, 0:-1]
        convolved_data = convolve_fft(map.data, target_psf_crop, boundary='wrap', normalize_kernel=False)
        return sunpy.map.Map(convolved_data, map.meta)

    def apply_mirror_scattered_light_psf(self, radiance_maps):
        return self.__apply_function_to_leaves(radiance_maps, self.__convolve_mirror_scatter)

    def __convolve_mirror_scatter(self, map):
        background = convolve_fft(map.data, self.mirror_scatter_psf, boundary='wrap', normalize_kernel=False)
        eta = np.sum(self.mirror_scatter_psf)
        convolved_data = map.data * (1 - eta) + background
        return sunpy.map.Map(convolved_data, map.meta)

    def apply_effective_area(self, radiance_maps):
        '''
        The effective area of SunCET is the product of the detector quantum 
        efficiency and the collecting area.  This module multiplies input 
        Sunpy radiance map sequence by the instrument effective area

        Parameters
        ----------
        radiance_maps : [sunpy.map.MapSequence]
            A sunpy radiance map sequence
        effective_area : [self.effective_area]
            Instrument exposure array

        Returns
        ---
        radiance_maps_long : ['sunpy.map.MapSequence']
            A sunpy radiance map sequence normalised to long exposure times

        TODO
        ----
        TODO:apply size test for maps and self.effective_area

        '''
        radiance_maps_scaled_by_effective_area = {}
        for timestep, radiance_maps_by_wavelength in radiance_maps.items():
            radiance_maps_scaled_by_effective_area[timestep] = {}
            for wavelength, map in radiance_maps_by_wavelength.items():
                radiance_maps_scaled_by_effective_area[timestep][wavelength] = map * self.effective_area[wavelength]
        return radiance_maps_scaled_by_effective_area

    def apply_exposure_times(self, radiance_maps):
        '''
        An input Sunpy radiance map sequence is normalized by exposure time 
        from a config file.  Two radiance map sequences are output, the first
        with long exposure times, the second with short exposure times.  The 
        output map sequences have the same dimensions as the input sequence.

        Parameters
        ----------
        radiance_maps : [sunpy.map.MapSequence]
            A sunpy radiance map sequence
        exposure_time : [self.config.exposure_time_long]
            A config exposure time
        exposure_time : [self.config.exposure_time_short]
            A config exposure time
 
        
        Returns
        ---
        radiance_maps_long : ['sunpy.map.MapSequence']
            A sunpy radiance map sequence normalised to long exposure times
        radiance_maps_short : ['sunpy.map.MapSequence']
            A sunpy radiance map sequence normalised to short exposure times

        TODO
        ---
        TODO: Test to see if the first return is long, and the secong short
        '''
        radiance_maps_short = {}
        radiance_maps_long = {}
        for timestep, radiance_maps_by_wavelength in radiance_maps.items():
            radiance_maps_short[timestep] = {}
            radiance_maps_long[timestep] = {}
            for wavelength, map in radiance_maps_by_wavelength.items():
                map_short = self.__set_exposure_time(map, self.config.exposure_time_short)
                radiance_maps_short[timestep][wavelength] = map_short * self.config.exposure_time_short

                map_long = self.__set_exposure_time(map, self.config.exposure_time_long)
                radiance_maps_long[timestep][wavelength] = map_long * self.config.exposure_time_long

        return {"short exposure": radiance_maps_short, "long exposure": radiance_maps_long}
    
    
    def __set_exposure_time(self, sunpy_map, new_exposure_time):
        if "EXPTIME" in sunpy_map.meta:
            sunpy_map.meta["EXPTIME"] = new_exposure_time
        elif "EXPOSURE" in sunpy_map.meta:
            sunpy_map.meta["EXPOSURE"] = new_exposure_time
        else:
            raise KeyError("The header does not contain an exposure time keyword (EXPTIME or EXPOSURE).")
        return sunpy_map

    def apply_photon_shot_noise(self, radiance_maps):
        return self.__apply_function_to_leaves(radiance_maps, self.__apply_poisson)
    
    
    def __apply_function_to_leaves(self, orig_dict, func, **kwargs):
        new_dict = {}
        for key, value in orig_dict.items():
            if isinstance(value, dict):
                new_dict[key] = self.__apply_function_to_leaves(value, func, **kwargs)
            else:
                new_dict[key] = func(value, **kwargs)
        return new_dict

    def __apply_poisson(self, map):
        noisy_photons = np.random.poisson(lam=map.data)
        return sunpy.map.Map(noisy_photons, map.meta)

    def convert_to_electrons(self, radiance_maps, apply_noise=True):
        quantum_yield = self.__compute_quantum_yields()
        quantum_yield_units_hacked = 1 * u.count/u.electron

        for exposure_type in ['short exposure', 'long exposure']:
            for timestep in radiance_maps[exposure_type]:
                first_map = self.__get_any_map(radiance_maps[exposure_type][timestep])
                summed_electrons = []

                for i, map in enumerate(radiance_maps[exposure_type][timestep].values()):
                    summed_electrons.append(map.data * map.unit * quantum_yield[i] * quantum_yield_units_hacked) # [ct / (Angstrom pix2)] but TODO: Should really be in electrons, not counts. Update once this issue has been resolved https://github.com/sunpy/sunpy/issues/6823

                # Perform numerical integration across wavelengths; Assuming self.wavelengths is in the correct order and units
                stacked_data = np.stack(summed_electrons, axis=-1)
                integrated_data = simpson(stacked_data, x=self.wavelengths, axis=-1)  # Warning: there isn't presently (2023-11-29) an integration method that propagates astropy units so we have to do it carefully ourselves

                # Store the integrated data in the same format as the original maps
                detector_image = np.zeros_like(first_map.data) * first_map.unit * quantum_yield.unit * quantum_yield_units_hacked * self.wavelengths.unit # This is where we're manually accounting for the unit change due to integration over wavelength, by multiplying by the wavelength unit
                detector_image += integrated_data * detector_image.unit

                updated_meta = self.__update_map_meta(first_map.meta, detector_image)
                detector_image = sunpy.map.Map(detector_image, updated_meta)

                if apply_noise:
                    detector_image = self.__apply_fano_noise(detector_image)

                detector_image = self.__clip_at_full_well(detector_image)

                radiance_maps[exposure_type][timestep] = detector_image

        return radiance_maps
    
    
    def __get_any_map(self, nested_dict):
        for key, value in nested_dict.items():
            if isinstance(value, dict):
                return self.__get_any_map(value)
            else:
                return value

    def __update_map_meta(self, meta, detector_image):
        meta['wavelnth'] = self.wavelengths.mean().value
        meta['waveunit'] = self.wavelengths.mean().unit.to_string()
        meta['bunit'] = detector_image.unit.to_string()
        return meta

    def __compute_quantum_yields(self):
        photoelectron_in_silicon = 3.63 * u.eV * u.ph / u.electron  # the magic number 3.63 here the typical energy [eV] to release a photoelectron in silicon
        return (const.h * const.c / self.wavelengths.to(u.m)).to(u.eV) / photoelectron_in_silicon

    def __apply_fano_noise(self, detector_image):
        sigma = np.sqrt(detector_image.data * self.config.fano_factor)
        gaussian_noise = np.random.normal(0, sigma, size=detector_image.data.shape)  # TODO: Fano isn't _really_ Gaussian in shape. It's like Poisson but narrowed by the Fano factor... but that's not quite right either because it can be negative (i.e., slightly fewer than the typical number of photoelectrons generated). The distribution is tiny though, so the precise distribution doesn't matter much. Gaussian gets the job done.
        noisy_data = np.clip(detector_image.data + gaussian_noise, a_min=0, a_max=None)  # negative electrons is unphysical
        noisy_data = np.rint(noisy_data).astype(int)  # can only have integer numbers of electrons

        return sunpy.map.Map(noisy_data, detector_image.meta)


    def __clip_at_full_well(self, detector_image):
        # TODO: uncomment this block once we can actually store things in electrons instead of counts (see GitHub issue cited above)
        #if detector_image.unit != self.config.pixel_full_well.unit: 
        #    raise u.UnitsError('The units for the detector image does not match the units of the pixel full well.') 
        if detector_image.unit == (u.ct/u.pix**2): 
            print('Just a reminder that we have to use counts (u.ct) until sunpy allows us to store with electrons (u.electron)')
            mask = detector_image.data > self.config.pixel_full_well.value
            detector_image.data[mask] = self.config.pixel_full_well
        
        return detector_image
    
    def __clip_at_bit_depth(self, detector_image):
        # Clip values above maximum bit depth
        mask_high = detector_image.data > 2**self.config.readout_bits.value - 1
        detector_image.data[mask_high] = 2**self.config.readout_bits.value - 1
        
        # Clip negative values to 0
        mask_low = detector_image.data < 0
        detector_image.data[mask_low] = 0
        return detector_image
    

    def make_dark_frame(self):
        if self.config.detector_temperature.unit != u.deg_C:
            raise u.UnitsError("The detector temperature must be in units of Celsius. Either convert in the config or edit the math in this function accordingly.")

        dark_current_mean = 20 * 2 ** ((self.config.detector_temperature.value - 20) / 5.5)
        dark_current_std = 12 * 2 ** ((self.config.detector_temperature.value - 20) / 5.5)
        dim_x, dim_y = int(self.config.image_dimensions[0].value), int(self.config.image_dimensions[1].value)

        dark_frame = np.random.normal(dark_current_mean, dark_current_std, (dim_y, dim_x))
        dark_frame = np.clip(dark_frame, 0, None)  # Can't have negative values
        dark_frame *= u.electron / u.second / u.pix ** 2

        dark_frame_short = dark_frame * self.config.exposure_time_short
        dark_frame_long = dark_frame * self.config.exposure_time_long

        self.dark_frame = {"short exposure": dark_frame_short, "long exposure": dark_frame_long}

    def make_read_frame(self):
        if self.config.read_noise.unit != (u.electron / u.pixel ** 2):
            raise u.UnitsError("The read noise must be in units of electrons/pixel^2.")

        dim_x, dim_y = int(self.config.image_dimensions[0].value), int(self.config.image_dimensions[1].value)

        read_frame = np.random.normal(0, self.config.read_noise.value, (dim_y, dim_x))  # Note: read noise is allowed to be negative (https://www.photometrics.com/learn/advanced-imaging/pattern-noise-dsnu-and-prnu#:~:text=Noise%20in%20scientific%20cameras%20is,from%20a%20software%20point%20of)
        read_frame *= (u.electron / u.pix ** 2)
        self.read_frame = read_frame

    def make_spike_masks(self, detector_images):
        spike_masks = {}

        for exposure_type, images in detector_images.items():
            # Determine the number of spikes based on exposure time
            if exposure_type == 'short exposure':
                number_spikes = int(self.config.detector_physical_area * self.config.spike_rate * self.config.exposure_time_short)
            elif exposure_type == 'long exposure':
                number_spikes = int(self.config.detector_physical_area * self.config.spike_rate * self.config.exposure_time_long)

            spike_masks[exposure_type] = {}

            # For each timestep nested in this exposure type, generate a unique spike mask
            for timestep in images:
                spike_masks[exposure_type][timestep] = self.__populate_mask_with_artifacts(number_spikes)

        # Store the generated spike masks
        self.spike_masks = spike_masks
    
    
    def __populate_mask_with_artifacts(self, number_artifacts):
        dim_x, dim_y = int(self.config.image_dimensions[0].value), int(self.config.image_dimensions[1].value)
        mask = np.zeros((dim_y, dim_x), dtype=bool)  # this is the way numpy and matplotlib interpret horizontal/vertical. sunpy understands it.

        x_coords = np.random.randint(0, self.config.image_dimensions[0].value, number_artifacts)
        y_coords = np.random.randint(0, self.config.image_dimensions[1].value, number_artifacts)
        mask[y_coords, x_coords] = True
        return mask

    def make_hot_pixel_mask(self):
        number_hot_pixels = int(self.config.fraction_hot_pixels * self.config.image_dimensions[0].value * self.config.image_dimensions[1].value)
        self.hot_pixel_mask = self.__populate_mask_with_artifacts(number_hot_pixels)

    def make_dead_pixel_mask(self):
        number_dead_pixels = int(self.config.fraction_dead_pixels * self.config.image_dimensions[0].value * self.config.image_dimensions[1].value)
        self.dead_pixel_mask = self.__populate_mask_with_artifacts(number_dead_pixels)

    def combine_signal_and_noise(self, detector_images):
        detector_images = self.__add_dark_frame_per_exposure_time(detector_images)
        detector_images = self.__add_read_frame(detector_images)
        detector_images = self.__add_spikes(detector_images)
        detector_images = self.__add_hot_pixels(detector_images)
        detector_images = self.__add_dead_pixels(detector_images)
        detector_images = self.__apply_function_to_leaves(detector_images, self.__clip_at_full_well)

        return detector_images

    def __add_dark_frame_per_exposure_time(self, detector_images):
        short = self.__apply_function_to_leaves(detector_images['short exposure'], self.__add_dark_frame, exposure='short exposure')
        long = self.__apply_function_to_leaves(detector_images['long exposure'], self.__add_dark_frame, exposure='long exposure')
        return {"short exposure": short, "long exposure": long}

    def __add_dark_frame(self, map, exposure='short exposure'):
        # TODO: uncomment this block once we can actually store things in electrons instead of counts (see GitHub issue cited above)
        # if map.unit != self.dark_frame.unit:
        #    raise u.UnitsError('The units for the detector image does not match the units of dark frame.')
        data = map.data + self.dark_frame[exposure].value
        return sunpy.map.Map(data, map.meta)

    def __add_read_frame(self, detector_images):
        read_frame_units_hacked = self.read_frame * u.count / u.electron  # TODO: Won't need this once we can store things in electrons
        return self.__apply_function_to_leaves(detector_images, lambda x: x + read_frame_units_hacked)

    def __add_spikes(self, detector_images):
        for exposure_type in detector_images:
            for timestep in detector_images[exposure_type]:
                map = detector_images[exposure_type][timestep]
                mask = self.spike_masks[exposure_type][timestep]

                image_data = map.data
                image_data[mask] = self.config.pixel_full_well
                detector_images[exposure_type][timestep] = sunpy.map.Map(image_data, map.meta)

        return detector_images

    def __add_hot_pixels(self, detector_images):
        hot_pixel_values = self.__generate_hot_pixel_values()
        return self.__apply_function_to_leaves(detector_images, self.__add_hot_pixels_to_single_map, hot_pixel_values=hot_pixel_values)

    def __generate_hot_pixel_values(self):
        dim_x, dim_y = int(self.config.image_dimensions[0].value), int(self.config.image_dimensions[1].value)

        # The magic numbers here are based on Dan's experience with the intensity range and distribution of hot pixels
        # The resultant distribution should range from 0 to about 158 or so
        dist = np.random.normal(loc=15, scale=5, size=(dim_y, dim_x))
        hot_pixel_values = (np.clip(dist, 25, None) - 25) * 10

        return hot_pixel_values * (u.count / u.pix ** 2)  # TODO: change units to u.electron once sunpy supports (see GitHub issue above)

    def __add_hot_pixels_to_single_map(self, map, hot_pixel_values=None):
        data = map.data.copy()
        data[self.hot_pixel_mask] += hot_pixel_values[self.hot_pixel_mask].value  # FIXME: this underrepresents the hot pixels because it requires that the randomly generated mask is true at pixel x, y and also that the separately randomly generated hot pixel value at x, y is ≠0.
        return sunpy.map.Map(data, map.meta)

    def __add_dead_pixels(self, detector_images):
        def set_dead_pixels_in_map(map):
            data = map.data.copy()
            data[self.dead_pixel_mask] = 0
            return sunpy.map.Map(data, map.meta)

        return self.__apply_function_to_leaves(detector_images, set_dead_pixels_in_map)

    def convert_to_dn(self, detector_images):
        gain_factor = self.config.detector_gain * u.electron / u.count  # TODO: fix these units once we can actually store things in electrons instead of counts (see GitHub issue cited above)
        detector_images = self.__apply_function_to_leaves(detector_images, lambda x: x * gain_factor)
        detector_images = self.__apply_function_to_leaves(detector_images, self.__clip_at_bit_depth)
        detector_images = self.__apply_function_to_leaves(detector_images, lambda x: sunpy.map.Map(x.data.astype(np.uint16), x.meta))
        return detector_images
    
        # TODO: Clip to DN max (Alan set the gain so that pixel full well [33k] electrons is 90% of the ADC dynamic range); therefore, if the image has already been clipped to full well, this clipping function shouldn't ever do anything
        # return sunpy.map.MapSequence([map *  self.config.detector_gain # TODO: This is all that'll be needed once TODO above is addressed
        #                              for map in detector_images])


class OnboardSoftware:
    def __init__(self, config):
        self.config = config

    def subtract_dark(self, detector_images):
        '''
        A dark frame map, drawn from a random distribution on the basis of a 
        config item for the average dark current in the detector as a function
        of temperature, is subtracted from each detector image. A map of
        detector images is returned
        

        Parameters
        ----------
        detector_images : [sunpy.map.MapSequence]
            A sunpy detector images map sequence

        Returns
        ---
        detector_images : [sunpy.map.MapSequence]
            A dark current subtracted sunpy detector images map sequence 

        TODO
        ---
        TODO: Test to see if detector images are same dimension as the dark frame
        '''
        return detector_images
        # return sunpy.map.MapSequence([map - self.dark_frame / map.meta['EXPOSURE'] for map in detector_images]) # TODO: Won't work until we have something sensible returned for make_dark_frame()

    def apply_jitter(self, onboard_processed_images):  # Note really onboard software, but this is the place in the logic it belongs
        jitter_arcsec_per_timestep = self.config.jitter * self.config.model_timestep  # Units work easily here. This is going to be used for a single, coherent, instantaneous movement affecting the position of every feature in the image but not their sharpness.

        for exposure_type, image_dict in onboard_processed_images.items():
            if exposure_type == 'short exposure':
                exposure_time = self.config.exposure_time_short
            else:
                exposure_time = self.config.exposure_time_long

            blur_amount = self.config.jitter * np.sqrt(exposure_time)  # The units here don't work but the sqrt is how we account for this being a continuous accumulation of the random walk. Features maintain their position, but lose their definition.

            # Apply jitter (shift and blur) to each image
            for index, image in image_dict.items():
                # Convert jitter from arcsec to pixels
                scale = image.scale[0]  # Assuming uniform scale in x and y
                jitter_pixels = jitter_arcsec_per_timestep / scale
                blur_pixels = blur_amount / scale

                # Generate random jitter offsets for shift
                angle = np.random.uniform(0, 2*np.pi)
                jitter_magnitude = np.random.normal(0, jitter_pixels.value)
                dx = jitter_magnitude * np.cos(angle)
                dy = jitter_magnitude * np.sin(angle)

                # Shift the image data to simulate the instantaneous start/stop boundary between timesteps (cadence steps)
                shifted_data = shift(image.data, shift=[dy, dx], mode='nearest')

                # Apply Gaussian blur to simulate motion during exposure
                blurred_data = gaussian_filter(shifted_data, sigma=blur_pixels.value)

                # Make sure that the shifting and blurring don't blow up in edge cases by clipping to the input -- verified that the mean is nearly identical
                blurred_data_clipped = np.clip(blurred_data, a_min=np.min(image.data), a_max=np.max(image.data))

                # Create a new SunPy map with the jittered and blurred data
                onboard_processed_images[exposure_type][index] = sunpy.map.Map(blurred_data_clipped, image.meta)

        return onboard_processed_images

    
    def filter_out_particle_hits(self, onboard_processed_images):        
        for exposure_type in ['short exposure', 'long exposure']:
            map_indices = list(onboard_processed_images[exposure_type].keys())

            if len(map_indices) != 3:
                raise ValueError("Need exactly 3 images to apply particle filtering per Alan's CSIE FPGA logic.")

            sorted_indices = sorted(map_indices)[:3]
            maps = [onboard_processed_images[exposure_type][index] for index in sorted_indices]
            stack = np.stack([map_obj.data for map_obj in maps], axis=-1)

            stack = stack.astype(np.int32) # Note: onboard firmware actually uses a 24-bit buffer for the following sum but numpy only has options of 16 or 32. This is fine in most cases for simulation because even 24 bit is really, really big. 

            # Remove max value along the stack axis and sum the remaining values
            max_removed = np.sort(stack, axis=-1)[:, :, :-1]  # Sort and remove last (max) value
            summed_data = np.sum(max_removed, axis=-1)  # Sum the remaining values

            # Shift the data to the right by the number of bits specified in the config to get big numbers to still fit in 16 bits
            shift_factor = 2 ** self.config.num_shift_bits_24_to_16
            shifted_data = (summed_data / shift_factor).astype(np.uint16)

            new_map = sunpy.map.Map(shifted_data, maps[0].meta)
            onboard_processed_images[exposure_type] = new_map

        return onboard_processed_images
    
    
    def create_composite(self, onboard_processed_images):
        map1 = onboard_processed_images['short exposure']
        map2 = onboard_processed_images['long exposure']
        solar_disk_center, solar_disk_radius = self.__get_solar_disk_center_and_radius_in_pixels(map1)

        # Create a boolean mask where True values represent the solar ~disk area
        y_grid, x_grid = np.mgrid[:map1.data.shape[0], :map1.data.shape[1]]
        disk_mask = np.sqrt(((x_grid - solar_disk_center[0].value)**2 + (y_grid - solar_disk_center[1].value)**2)) <= (self.config.inner_fov_circle_radius.value * solar_disk_radius.value)

        # Combine the two maps using the boolean mask
        combined_data = np.where(disk_mask, map1.data, map2.data)
        return sunpy.map.Map(combined_data, map1.meta)


    def __get_solar_disk_center_and_radius_in_pixels(self, map): 
        solar_disk_center = np.array([
            map.world_to_pixel(map.center.transform_to(frames.Helioprojective(observer=map.observer_coordinate))).x.value,
            map.world_to_pixel(map.center.transform_to(frames.Helioprojective(observer=map.observer_coordinate))).y.value
        ]) * u.pix
        solar_disk_radius = (map.rsun_obs / map.scale[0]).decompose()

        return solar_disk_center, solar_disk_radius
    
    
    def create_image_histogram(self, onboard_processed_images):
        data = onboard_processed_images.data
        histogram = np.histogram(data, bins=self.config.num_bins, range=(0, np.max(data)))
        return histogram
    

    def bin_image(self, onboard_processed_images, xbin=None, ybin=None):
        '''
        An input Sunpy radiance map sequence is binned by the factors 
        xbin and ybin.

        Parameters
        ----------
        radiance_maps : [sunpy.map.MapSequence]
            A sunpy radiance map sequence
        xbin : [int]
            A factor by which the x dimension of the image is binned
        ybin : [int]
            A factor by which the x dimension of the image is binned
        
        Returns
        ---
        detector_images : [sunpy.map.MapSequence]
            A binned sunpy detector images map sequence 

        TODO
        ---
        check if input is a sunpy map
        '''
        # check if bin dimensions added
        if xbin == None:
            xbin = self.config.num_pixels_to_bin[0]
        if ybin == None:
            ybin = self.config.num_pixels_to_bin[1]

        # See input array shape
        array_shape = np.shape(onboard_processed_images.data)

        # Check if a number is a factor of input array dimensions.
        if array_shape[1] % xbin != 0:
            raise ValueError('xbin value should be an integer value and a factor of the array x dimension')
        elif array_shape[0] % ybin != 0:
            raise ValueError('ybin value should be an integer value and a factor of the array y dimension')
        else:    
            new_x_dim = int(array_shape[1]/xbin)
            new_y_dim = int(array_shape[0]/ybin)
            
            new_dimensions = [new_x_dim, new_y_dim] * u.pixel
            
            onboard_processed_images = onboard_processed_images.resample(new_dimensions, method='nearest')
            onboard_processed_images *= (xbin * ybin) # Conserve energy. DN come from electrons. resample doesn't account for the fact that the total number of electrons (DN) recorded across the detector is the same regardless of how you bin them in software
            
            # Flight firmware uses 24 bits when doing summing to avoid rollover that would otherwise be possible at 16 bits
            max_24bit = 2**24 - 1
            data_24_bit = np.clip(onboard_processed_images.data, 0, max_24bit).astype(np.uint32)
            
            new_map = sunpy.map.Map(data_24_bit, onboard_processed_images.meta)
            new_map.meta['cdelt1'] = self.config.plate_scale.value * xbin # Note 1: this is only needed because sunpy (v4.0.1) resample updates dimensions but NOT plate scale
            new_map.meta['cdelt2'] = self.config.plate_scale.value * ybin # Note 2: there is also risk here because a) the user must be responsible in the config file to ensure the image_dimensions and plate_scale are compatible, and b) the units they input for plate_scale must be the same as those already in the map

        return new_map
    

    def compress_image(self, onboard_processed_images): 
        normalized_data = (onboard_processed_images.data - np.min(onboard_processed_images.data)) / (np.max(onboard_processed_images.data) - np.min(onboard_processed_images.data))
        max_value = 2**self.config.readout_bits.value - 1
        
        if self.config.readout_bits <= (8 * u.bit):
            dtype = np.uint8
            mode = 'L'
        elif self.config.readout_bits <= (16 * u.bit): # if readout_bits is, e.g., 12 or 14, we've scaling accordingly and just won't use the excess bits in these standard numpy containers (8, 16, 32)
            dtype = np.uint16
            mode = 'I;16'
        elif self.config.readout_bits <= (32 * u.bit):
            dtype = np.uint32
            mode = 'I'
        else:
            raise ValueError("Unsupported bit depth")
        
        scaled_data = (normalized_data * max_value).astype(dtype)
        image = Image.fromarray(scaled_data, mode=mode)
        return sunpy.map.Map(np.array(image), onboard_processed_images.meta)
    

    def __display_data(data): 
        import matplotlib.pyplot as plt
        plt.imshow(np.log10(np.clip(data, 0.1, None)))


if __name__ == "__main__":
    config_filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    config = config_parser.Config(config_filename)
    hardware = Hardware(config)
    print("See test_instrument.py for example of how to configure to run this code.")