"""
This is the main wrapper for most/all(?) of the other instrument simulator related python files
"""
import os
from glob import glob
import astropy.units as u
from astropy.io import fits
from astropy.time import Time, TimeDelta
import sunpy.map
import pandas as pd
import ast
import numpy as np
from suncet_instrument_simulator import config_parser, make_radiance_maps, instrument, stack_schedule
from suncet_instrument_simulator.image_statistics import local_rms


def _set_observation_times_from_map(header, map_header, duration):
    """Set output timestamps from the already-timed input radiance map."""
    duration_seconds = (
        float(duration.to_value(u.s)) if hasattr(duration, 'to_value') else float(duration)
    )
    duration_offset = TimeDelta(duration_seconds, format='sec')
    timestamp = map_header.get('DATE-OBS')
    if timestamp in (None, '', 'N/A'):
        raise ValueError('Processed radiance map does not contain a valid DATE-OBS.')

    time_scale = str(map_header.get('TIMESYS', 'UTC')).strip().lower()
    observation_start = Time(timestamp, format='isot', scale=time_scale, precision=3).utc
    header['TIMESYS'] = 'UTC'
    header['DATE-BEG'] = observation_start.isot
    header['DATE-OBS'] = observation_start.isot
    header['DATE-END'] = (observation_start + duration_offset).isot

    return header


class Simulator:
    def __init__(self, config_filename=os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'):
        self.config_filename = config_filename
        self.config = self.__read_config(config_filename)
        self.radiance_maps = 'not yet loaded'
        self.hardware = 'not yet loaded'
        self.onboard_software = 'not yet loaded'
        self.metadata = 'not yet loaded'


    def __read_config(self, config_filename):   
        return config_parser.Config(config_filename)


    def run(self, observation_indices=None, *, retain_pure_reference=False):
        """Simulate observations, optionally retaining the noise-free diagnostic."""
        self.retain_pure_reference = retain_pure_reference
        self.hardware = instrument.Hardware(self.config)
        self.onboard_software = instrument.OnboardSoftware(self.config)

        observations = self.__get_observations()
        if observation_indices is not None:
            requested = {int(index) for index in observation_indices}
            available = {observation.output_index for observation in observations}
            unknown = sorted(requested - available)
            if unknown:
                raise ValueError('Unknown observation indices requested: {}'.format(unknown))
            observations = [
                observation for observation in observations
                if observation.output_index in requested
            ]

        for observation in observations:
            self.current_observation = observation
            self.current_output_index = observation.output_index
            self.current_timestep = str(self.current_output_index).zfill(3)
            self.__sun_emission()
            if self.radiance_maps_found:
                self.__sun_to_detector()
                self.__simulate_noise()
                self.__simulate_detector()
                self.__apply_camera_software()
                #self.__calculate_snr()
                #self.__plot_snr() # FIXME: Remove this. It's just for debugging or reference. 
                self.__complete_metadata()
                self.__output_files()
    

    def __get_observations(self):
        return stack_schedule.build_observation_sequence(self.config)


    def __sun_emission(self): 
        if self.config.compute_new_radiance_maps:
            self.radiance_maps = make_radiance_maps.MakeRadianceMaps(self.config).run()
            self.config.compute_new_radiance_maps = False # Only need to do it the one time, not every loop
        self.__load_radiance_maps()


    def __load_radiance_maps(self):
        self.radiance_maps_found = True
        start_seconds = self.current_observation.start_seconds
        self.stack_schedule = stack_schedule.build_stack_schedule_at_time(
            start_seconds, self.config)
        required_indices = self.stack_schedule.unique_model_indices
        filenames = self.__get_radiance_map_filenames(required_indices)
        available_indices = sorted(
            int(os.path.basename(filename).split('_')[-1].replace('.fits', ''))
            for filename in filenames
        )
        available_prefix = required_indices[:len(available_indices)]
        if not available_indices or available_indices != available_prefix:
            print(
                'Need {} radiance maps for output {} at {:.6g} seconds but only found {}.'.format(
                    len(required_indices), self.current_output_index, start_seconds, len(filenames)
                )
            )
            self.radiance_maps_found = False
            return
        if len(available_indices) < len(required_indices):
            print(
                'Missing trailing radiance maps {} for output {} at {:.6g} seconds; '
                'padding with model index {}.'.format(
                    required_indices[len(available_indices):],
                    self.current_output_index,
                    start_seconds,
                    available_indices[-1],
                )
            )

        self.radiance_by_model_index = {}
        for filename in filenames:
            index = int(os.path.basename(filename).split('_')[-1].replace('.fits', ''))
            maps = sunpy.map.Map(filename)
            maps_by_wavelength = {}
            for map in maps:
                wavelength = str(map.wavelength)
                if 'RSUN_REF' not in map.meta:
                    map.meta['RSUN_REF'] = 6.96e8  # meters but can't use astropy units here because its a FITS header standard
                maps_by_wavelength[wavelength] = map
            self.radiance_by_model_index[index] = maps_by_wavelength

        available_indices = sorted(self.radiance_by_model_index.keys())
        self.stack_schedule = stack_schedule.build_stack_schedule_at_time(
            start_seconds, self.config, available_indices=available_indices
        )


    def __get_radiance_map_filenames(self, model_indices):
        base_directory = os.path.join(
            os.getenv('suncet_data'),
            self.config.model_data_folder.strip('/'),
            self.config.euv_radiance_map_directory_name.strip('/'),
        )
        filenames = []
        for model_index in model_indices:
            file_pattern = 'radiance_maps_' + str(model_index).zfill(3) + '.fits'
            filenames.extend(glob(os.path.join(base_directory, file_pattern)))
        return sorted(filenames)


    def __process_radiance_through_optics(self, radiance_maps_by_member):
        radiance_maps = self.hardware.extract_fov(radiance_maps_by_member)
        radiance_maps = self.hardware.interpolate_spatial_resolution(radiance_maps)
        radiance_maps = self.hardware.convert_steradians_to_pixels(radiance_maps)
        if self.config.apply_mesh_diffraction:
            radiance_maps = self.hardware.apply_diffraction_psf(radiance_maps)
        if self.config.apply_mirror_scattered_light_psf:
            radiance_maps = self.hardware.apply_mirror_scattered_light_psf(radiance_maps)
        radiance_maps = self.hardware.apply_effective_area(radiance_maps)
        return radiance_maps


    def __sun_to_detector(self):
        # Cache only deterministic optical data within this observation. Noise
        # still runs for every integration, in the original order. Processing
        # one distinct member at a time also bounds intermediate FFT storage.
        optical_cache = {}
        exposed_cache = {}
        self.radiance_maps = {}
        for exposure, members, exposure_time in (
            ('short exposure', self.stack_schedule.short_members, self.config.exposure_time_short),
            ('long exposure', self.stack_schedule.long_members, self.config.exposure_time_long),
        ):
            radiance_by_member = {}
            for member_index, contributions in enumerate(members):
                key = stack_schedule.radiance_member_key(contributions)
                start_seconds = (self.stack_schedule.start_seconds
                                 + member_index * exposure_time.to_value(u.s))
                if key not in optical_cache:
                    member = stack_schedule.combine_radiance_for_member(
                        contributions, self.radiance_by_model_index,
                        start_seconds=start_seconds, model_timestep=self.config.model_timestep)
                    if not optical_cache:
                        self.hardware.store_target_wavelengths(member)
                        self.hardware.compute_effective_area()
                    optical_cache[key] = self.__process_radiance_through_optics({0: member})[0]

                exposure_key = (key, exposure_time.to_value(u.s))
                if exposure_key not in exposed_cache:
                    independent_metadata = {
                        wavelength: sunpy.map.Map(image.data, image.meta.copy())
                        for wavelength, image in optical_cache[key].items()
                    }
                    exposed = self.hardware.apply_exposure_times_for_stack(
                        {0: independent_metadata}, exposure_time)[0]
                    # Expected photon counts are deterministic too. Sharing
                    # their immutable arrays avoids retaining nine copies of
                    # the short scene while drawing independent Poisson noise.
                    for image in exposed.values():
                        image.data.setflags(write=False)
                    exposed_cache[exposure_key] = exposed

                # Preserve the cached crop/resample WCS, but give each exposure
                # its own timestamp and metadata: later stages mutate metadata.
                reference_index = key[0][0]
                copied_member = {}
                for wavelength, cached_map in exposed_cache[exposure_key].items():
                    metadata = cached_map.meta.copy()
                    reference = self.radiance_by_model_index[reference_index][wavelength]
                    metadata.update(stack_schedule.integration_time_metadata(
                        reference, reference_index, start_seconds, self.config.model_timestep))
                    copied_member[wavelength] = sunpy.map.Map(cached_map.data, metadata)
                radiance_by_member[member_index] = copied_member
            self.radiance_maps[exposure] = radiance_by_member

        self.radiance_maps_pure = None
        self.detector_images_pure = None
        if getattr(self, 'retain_pure_reference', False):
            self.radiance_maps_pure = {
                exposure: dict(members) for exposure, members in self.radiance_maps.items()
            }


    def __simulate_noise(self):
        self.radiance_maps = self.hardware.apply_photon_shot_noise(self.radiance_maps)
        self.detector_images = self.hardware.convert_to_electrons(self.radiance_maps, apply_noise=True)
        if self.radiance_maps_pure is not None:
            self.detector_images_pure = self.hardware.convert_to_electrons(self.radiance_maps_pure, apply_noise=False)
        self.hardware.make_dark_frame()
        self.hardware.make_read_frame()
        self.hardware.make_spike_masks(self.detector_images)
        self.hardware.make_hot_pixel_mask()
        self.hardware.make_dead_pixel_mask()
    

    def __simulate_detector(self):
        self.detector_images = self.hardware.combine_signal_and_noise(self.detector_images)
        self.detector_images = self.hardware.convert_to_dn(self.detector_images)


    def __apply_camera_software(self):
        if self.config.subtract_dark:
            self.onboard_processed_images = self.onboard_software.subtract_dark(self.detector_images)
        else: 
            self.onboard_processed_images = self.detector_images
        self.onboard_processed_images = self.onboard_software.apply_jitter(self.onboard_processed_images)
        if self.config.filter_out_particle_hits:
            self.onboard_processed_images = self.onboard_software.filter_out_particle_hits(self.onboard_processed_images)
        # Always collapse stacks so create_composite receives single maps even when
        # particle-hit filtering is disabled.
        self.onboard_processed_images = self.onboard_software.collapse_exposure_stacks(
            self.onboard_processed_images,
            collapse_method='first'
        )
        self.onboard_processed_images = self.onboard_software.create_composite(self.onboard_processed_images)
        self.image_histogram = self.onboard_software.create_image_histogram(self.onboard_processed_images)
        self.onboard_processed_images = self.onboard_software.bin_image(self.onboard_processed_images)
        self.onboard_processed_images = self.onboard_software.bit_shift_data(self.onboard_processed_images)
        if self.config.compress_image:
            self.onboard_processed_images = self.onboard_software.compress_image(self.onboard_processed_images)


    def __calculate_snr(self):
        if self.detector_images_pure is None:
            raise ValueError('SNR requires run(retain_pure_reference=True).')
        
        # generate no-noise image with compatible parameters to compare to simulated image
        composite_images_pure = self.hardware.convert_to_dn(self.detector_images_pure)
        if self.config.filter_out_particle_hits:
            composite_images_pure = self.onboard_software.filter_out_particle_hits(composite_images_pure)
        else:
            composite_images_pure = self.onboard_software.collapse_exposure_stacks(
                composite_images_pure, collapse_method='first'
            )
        composite_image_pure = self.onboard_software.create_composite(composite_images_pure)
        composite_image_pure_binned = self.onboard_software.bin_image(composite_image_pure)

        # generate pure noise image
        noise_image = self.onboard_processed_images.data - composite_image_pure_binned.data

        local_std = local_rms(noise_image, self.config.SNR_window.value)

        # deal with 0 noise pixels that would blow up the SNR
        zero_noise_indices = np.where(local_std == 0)
        if zero_noise_indices[0].size > 0:  # Only modify if there are zero noise indices
            local_std[zero_noise_indices] = 1
            data = composite_image_pure_binned.data
            data[zero_noise_indices] = float('inf')
        else:
            data = composite_image_pure_binned.data

        self.snr_image = data/local_std


    def __plot_snr(self):
        import matplotlib.pyplot as plt
        from scipy.ndimage import uniform_filter

        snr_map = self.snr_image

        # Deal with infinities
        neutral_value = np.nanmedian(snr_map[np.isfinite(snr_map)])  # For example, the median of finite values
        snr_map_no_inf = np.where(np.isinf(snr_map), neutral_value, snr_map)


        # Smooth the SNR map
        window_size = 20 # Note: IDL does "snr_smooth = smooth(rebin_pure_image/local_rms, 20, /edge_truncate)" and uses that for plotting the contours and pulling the +3.5 Rs SNR
        extended_array = np.pad(snr_map_no_inf, pad_width=window_size//2, mode='edge')
        smoothed_extended = uniform_filter(extended_array, size=window_size, mode='constant')
        start = window_size // 2
        end_offset = window_size - start
        smoothed_snr_map = smoothed_extended[start:-end_offset, start:-end_offset]

        # Create a figure and a set of subplots
        fig, ax = plt.subplots(figsize=(10, 10))

        # Display the image
        im = ax.imshow(self.onboard_processed_images.data, cmap='gray', origin='lower')

        # Overlaying the SNR contours
        contour_levels = [10, 40]
        ax.contour(smoothed_snr_map, levels=[10], colors='red', linewidths=2)
        ax.contour(smoothed_snr_map, levels=[40], colors='dodgerblue', linewidths=2)

        # Add the horizontal line for the trace
        vertical_center = round(smoothed_snr_map.shape[0] // 2.5)
        ax.axhline(y=vertical_center, color='limegreen', linestyle='--', linewidth=2)

        # Set the title and labels
        ax.set_title('Signal + Noise Image with SNR Contours')
        ax.set_xlabel('X Pixel')
        ax.set_ylabel('Y Pixel')

        # Adjusting the ticks for solar radii
        height, width = np.shape(snr_map)
        center_x, center_y = width // 2, height // 2
        scale_factor = 100  # 100 pixels per solar radius

        x_ticks = np.arange(0, width, scale_factor)
        x_labels = (x_ticks - center_x) / scale_factor
        y_ticks = np.arange(0, height, scale_factor)
        y_labels = (y_ticks - center_y) / scale_factor

        ax.set_xticks(x_ticks)
        ax.set_xticklabels([f"{x:.1f}" for x in x_labels])
        ax.set_yticks(y_ticks)
        ax.set_yticklabels([f"{y:.1f}" for y in y_labels])

        ax.set_xlabel("Solar Radii")
        ax.set_ylabel("Solar Radii")

        # Plot the horizontal trace
        horizontal_trace = smoothed_snr_map[vertical_center, :]

        plt.figure(figsize=(10, 6))
        plt.plot(horizontal_trace, color='limegreen')

        # Add vertical line at 3.5 Rs
        pixel_position = center_x + int(3.5 * scale_factor)
        plt.axvline(x=pixel_position, color='black', linestyle='--', linewidth=1)

        # Adjust the x-axis to represent solar radii
        num_pixels = len(horizontal_trace)
        x_ticks = np.arange(0, num_pixels, scale_factor)
        x_labels = (x_ticks - num_pixels // 2) / scale_factor
        plt.xticks(x_ticks, [f"{x:.1f}" for x in x_labels])

        # Set labels and title
        plt.title('SNR Horizontal Trace at Vertical Center, MHD model {}, frame {}'.format(self.config.model_directory_name, self.current_timestep))
        plt.xlabel('Solar Radii')
        plt.ylabel('SNR')
        plt.ylim(0, np.max(contour_levels))
        plt.grid(True)
        plt.show()

        # Print SNR value at +3.5 Rs
        snr_at_3_5 = smoothed_snr_map[vertical_center, pixel_position]
        print('SNR at +3.5 Rs = {:.0f}'.format(snr_at_3_5))

        pass


    def __metadata_field_name(self, row):
        if pd.isna(row.get('Field Name')):
            return ''
        return str(row['Field Name'])

    def __metadata_fits_keyword(self, row):
        if pd.isna(row.get('FITS variable name')):
            return None
        keyword = str(row['FITS variable name']).strip()
        if not keyword:
            return None
        return keyword

    def __metadata_description(self, row):
        if pd.isna(row.get('Description')):
            return ''
        return str(row['Description'])

    def __complete_metadata(self):
        metadata_definition = self.__load_metadata_definition()
        map = self.__strip_units_for_fits_compatibility(self.onboard_processed_images)

        header = fits.Header()
        old_header = self.__convert_sunpy_meta_to_fits_header(map)

        for _, row in metadata_definition.iterrows():
            field_name = self.__metadata_field_name(row)
            if 'COMMENT' in field_name:
                comment_text = field_name.replace('COMMENT ', '').strip()
                if comment_text:
                    header.set('COMMENT', value=comment_text, after=len(header))
                continue

            fits_keyword = self.__metadata_fits_keyword(row)
            if fits_keyword is None:
                continue

            if fits_keyword in header:
                continue

            value = row.get('typical value')
            if pd.isna(value):
                value = 'N/A'
            else:
                try:
                    value = ast.literal_eval(str(value))
                except (ValueError, SyntaxError):
                    value = str(value)

            try:
                header.set(
                    fits_keyword,
                    value=value,
                    comment=self.__metadata_description(row),
                    after=len(header),
                )
            except (ValueError, KeyError) as error:
                print('Skipping metadata row with FITS keyword {!r}: {}'.format(fits_keyword, error))



        # Populate metadata defined by the config file or resultant from the simulation
        header.set('LEVEL', value='0.5')
        header.set('TIMESYS', value='UTC')
        header.set('NBIN', value=(self.config.num_pixels_to_bin[0] * self.config.num_pixels_to_bin[1]))
        header.set('NBIN1', value=self.config.num_pixels_to_bin[0])
        header.set('NBIN2', value=self.config.num_pixels_to_bin[1])
        header.set('DET_TEMP', value=self.config.detector_temperature.value)
        _set_observation_times_from_map(
            header,
            map_header=old_header,
            duration=self.config.observation_window,
        )

        # Most values still come from the metadata definition above. Observation
        # timestamps are deliberately sourced from the simulated radiance maps.
        # header.set('EXPTIME', map.meta['EXPTIME'])
        # header.set('RSUN_REF', map.rsun_meters.value)
        # header.set('IMAGEW', value=map.dimensions[0].value)
        # header.set('IMAGEH', value=map.dimensions[1].value)

        # previously commented out
        # header.set('HISTOGRAM_SHORT', value=self.image_histogram[0]) # FIXME: arrays can't be stored in the header. For flight they'd be in the metadata and on the ground in the hdf5 file
        # header.set('HISTOGRAM_LONG', value=self.image_histogram[1])

        hdu = fits.PrimaryHDU(map.data, header)
        hdul = fits.HDUList(hdu)

        self.fits = hdul
    

    def __load_metadata_definition(self):
        return pd.read_csv(os.getenv('suncet_data') + '/metadata/' + self.config.base_metadata_filename)
    
    
    def __strip_units_for_fits_compatibility(self, map):
        meta = map.meta
        for key, value in meta.items():
            if isinstance(value, u.Quantity):
                value = value.value
            meta[key] = value
        return sunpy.map.Map(map.data, meta)
    

    def __convert_sunpy_meta_to_fits_header(self, map):
        # PrimaryHDU supplies structural and unsigned-data scaling keywords,
        # just as Map.save does, without writing the image merely for its header.
        return fits.PrimaryHDU(map.data, header=map.fits_header).header.copy()


    def __output_files(self):
        self.__write_fits()
        self.__write_binary()
        self.__output_snr()
        pass # TODO: implement output_files
    

    def __write_fits(self):
        path = os.getenv('suncet_data') + '/synthetic/level0/fits/'
        output_index = str(self.current_output_index).zfill(3)
        filename = os.path.splitext(os.path.basename(self.config_filename))[0] + '_OBS_' + self.fits[0].header['DATE-OBS'] + '_' + output_index + '.fits'
        self.fits[0].header.set('FILENAME', value=filename)
        
        self.fits.writeto(path+filename, overwrite=True, checksum=True)
        print('Wrote file: {}'.format(path+filename))


    def __write_binary(self):
        pass # TODO: implement write_binary() to mimic onboard (or downlinked?) storage


    def __output_snr(self):

        pass # TODO: implement output_snr()


# Convenience function primarily for debugging -- provide the data for the image, e.g., map.data or fits.data and this will plot with good scalings
def __display_data(data): 
    import matplotlib.pyplot as plt
    plt.imshow(np.log10(np.clip(data, 0.1, None)))


if __name__ == "__main__":
    simulator = Simulator()
    simulator.run()
