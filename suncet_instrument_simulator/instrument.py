import os
import warnings
from glob import glob
import numpy as np
import numpy.ma as ma
from pandas import read_fwf, read_csv
import astropy.units as u
from astropy.units import UnitsError
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import constants as const
import sunpy.map
from sunpy.coordinates import frames
from suncet_instrument_simulator import config_parser # This is just for running as a script -- should delete when done testing
from scipy.io import readsav # TODO: remove this temporary hack once we have radiance map files in non-IDL-saveset format

class Hardware: 
    def __init__(self, config):
        self.config = config
        self.coating_name = self.__get_coating_name()


    def __get_coating_name(self):
        return os.path.basename(self.config.mirror_coating_reflectivity_filename).split('_1')[0] 


    def store_target_wavelengths(self, radiance_maps): # Will interpolate/calculate any subsequent wavelength-dependent quantities to this "target" wavelength array
        wavelengths = []
        wavelength_unit = []
        for map in radiance_maps:
            wavelengths.append(map.wavelength.value)
            if len(wavelength_unit) > 0: 
                if map.wavelength.unit != wavelength_unit[0]:
                    raise ValueError('The wavelengths in the radiance maps do not match each other. This is dangerous so please fix it.')
            wavelength_unit.append(map.wavelength.unit)
        self.wavelengths = wavelengths * wavelength_unit[0]


    def compute_effective_area(self):
        geometric_area = np.pi * ((self.config.entrance_aperture)/2)**2
        mirror_reflectivity = self.__interpolate_mirror_coating_reflectivity()
        filter_transmission = self.__interpolate_filter_transmission()
        quantum_efficiency = self.__interpolate_quantum_efficiency()
        
        self.effective_area = geometric_area * mirror_reflectivity * filter_transmission * quantum_efficiency


    def __interpolate_mirror_coating_reflectivity(self):
        filename = self.config.mirror_coating_reflectivity_filename
        source_data = self.__load_mirror_data(filename)
        return np.interp(self.wavelengths.value, source_data.wavelength.values, source_data.reflectivity.values)
   

    def __load_mirror_data(self, filename):
        return read_fwf(os.getenv('suncet_data') + '/mirror_reflectivity/' + filename, skiprows=18, header=0, names=['wavelength', 'reflectivity']) # wavelength should be in Angstroms, reflectivity as a fraction


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
        df = read_csv(os.getenv('suncet_data') + 'quantum_efficiency/' + self.config.quantum_efficiency_filename, skiprows=1)
        df.columns = ['wavelength [Å]', 'qe']
        return df


    def extract_fov(self, radiance_maps):
        fov_half_angles = self.config.fov / 2.0

        # TODO: Raise warning if instrument FOV > model FOV

        return sunpy.map.MapSequence([map.submap(top_right=SkyCoord(fov_half_angles[0], fov_half_angles[1], frame=map.coordinate_frame), 
                                                 bottom_left=SkyCoord(-fov_half_angles[0], -fov_half_angles[1], frame=map.coordinate_frame)) 
                                      for map in radiance_maps])
    

    def interpolate_spatial_resolution(self, radiance_maps):
        map_list = []
        for map in radiance_maps:
            resampled_map = map.resample(self.config.image_dimensions, method='linear')  * ((self.config.plate_scale / map.scale[0]).value)**2
            resampled_map.data[resampled_map.data < 0] = 0 # Clip any negative numbers that result from the interpolation above -- bottom out at 0
            resampled_map.meta['cdelt1'] = self.config.plate_scale.value # Note 1: this is only needed because sunpy (v4.0.1) resample updates dimensions but NOT plate scale
            resampled_map.meta['cdelt2'] = self.config.plate_scale.value # Note 2: there is also risk here because a) the user must be responsible in the config file to ensure the image_dimensions and plate_scale are compatible, and b) the units they input for plate_scale must be the same as those already in the map
            map_list.append(resampled_map)
        return sunpy.map.MapSequence(map_list)


    def __compute_plate_scale(): 
        pass # TODO: is this needed? Even if not, would it be useful to include some computed quantities as instance variables just for reference?


    def apply_psf(self, radiance_maps):
        return radiance_maps
        pass # TODO: implement apply_psf


    def apply_scattered_light_psf(self, radiance_maps):
        return radiance_maps
        pass # TODO: implement apply_scattered_light_psf

    
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
        maps = []
        for i, map in enumerate(radiance_maps):
            maps.append(map * self.effective_area[i])
        return sunpy.map.MapSequence(maps)
        

    
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
        long_exposure_map_list = []
        short_exposure_map_list = []

        for map in radiance_maps: 
            map = self.__set_exposure_time(map, self.config.exposure_time_long)
            long_exposure_map_list.append(map * self.config.exposure_time_long)

            map = self.__set_exposure_time(map, self.config.exposure_time_short)
            short_exposure_map_list.append(map * self.config.exposure_time_short)
        
        return sunpy.map.MapSequence(long_exposure_map_list, short_exposure_map_list)
    

    def __set_exposure_time(self, sunpy_map, new_exposure_time):
        if "EXPTIME" in sunpy_map.meta:
            sunpy_map.meta["EXPTIME"] = new_exposure_time
        elif "EXPOSURE" in sunpy_map.meta:
            sunpy_map.meta["EXPOSURE"] = new_exposure_time
        else:
            raise KeyError("The header does not contain an exposure time keyword (EXPTIME or EXPOSURE).")
        return sunpy_map


    def apply_photon_shot_noise(self, radiance_maps):
        map_list = []
        for map in radiance_maps:
            noisy_photons = np.random.poisson(lam=map.data)
            map_list.append(sunpy.map.Map(noisy_photons, map.meta))
        
        return sunpy.map.MapSequence(map_list)


    def convert_to_electrons(self, radiance_maps, apply_noise=True):
        radiance_maps_long, radiance_maps_short = self.__filter_mapsequence_by_exposure(radiance_maps)

        detector_image_long = self.__convert_to_electrons(radiance_maps_long, apply_noise=apply_noise)
        detector_image_short = self.__convert_to_electrons(radiance_maps_short, apply_noise=apply_noise)
        return sunpy.map.MapSequence(detector_image_long, detector_image_short)


    def __filter_mapsequence_by_exposure(self, mapsequence):
        short_exposures = [m for m in mapsequence if m.exposure_time.value == self.config.exposure_time_short.value]
        long_exposures = [m for m in mapsequence if m.exposure_time.value == self.config.exposure_time_long.value]
        return sunpy.map.MapSequence(long_exposures), sunpy.map.MapSequence(short_exposures)
    

    def __convert_to_electrons(self, radiance_maps, apply_noise=True):
        quantum_yield = self.__compute_quantum_yields()
        detector_image = radiance_maps[0].data * 0 # empty array
        for i, map in enumerate(radiance_maps): 
            #detector_images.append(map * quantum_yield[i]) # TODO: Update once this issue has been resolved https://github.com/sunpy/sunpy/issues/6823
            quantum_yield_units_hacked = quantum_yield[i] * u.count/u.electron
            detector_image = detector_image + (map.data * map.unit) * quantum_yield_units_hacked # TODO: Should really be in electrons, not counts
        detector_image = self.__update_map_meta(detector_image, radiance_maps[0].meta, quantum_yield_units_hacked) # TODO: Update to "quantum_yield" once above issue resolved

        if apply_noise:
            detector_image = self.__apply_fano_noise(detector_image)

        detector_image = self.__clip_at_full_well(detector_image)
        return detector_image
    

    def __update_map_meta(self, detector_image, meta, quantum_yield):
        meta['wavelnth'] = self.wavelengths.mean().value
        meta['waveunit'] = self.wavelengths.mean().unit.to_string()
        
        map = sunpy.map.Map(detector_image, meta)

        map = map * (1 * quantum_yield.unit)
        return map


    def __compute_quantum_yields(self):
        photoelectron_in_silicon = 3.63 * u.eV * u.ph / u.electron # the magic number 3.63 here the typical energy [eV] to release a photoelectron in silicon
        return (const.h * const.c / self.wavelengths.to(u.m)).to(u.eV) / photoelectron_in_silicon 


    def __apply_fano_noise(self, detector_image):
        fano_factor = 0.1190 # Based on material, silicon in this case, and fairly insensitive to wavelength. This number comes from Rodrigues+ 2021 (https://doi.org/10.1016/j.nima.2021.165511)
        header = detector_image.meta
        detector_image = np.random.poisson(lam=detector_image.data) * np.sqrt(fano_factor)
        return sunpy.map.Map(detector_image, header)


    def __clip_at_full_well(self, detector_image):
        if detector_image.unit != self.config.pixel_full_well.unit: 
            # raise UnitsError('The units for the detector image does not match the units of the pixel full well.') # TODO: uncomment this once earlier functions are written to ensure the units actually are the same
            return detector_image
        mask = detector_image.data > self.config.pixel_full_well.value
        detector_image.data[mask] = self.config.pixel_full_well
        
        return detector_image
    
    def make_dark_frame(self):
        pass # TODO: implement make_dark_frame


    def make_read_frame(self):
        pass # TODO: implement make_read_frame

    
    def make_spike_frame(self):
        pass # TODO: implement make_spike_frame
    

    def combine_signal_and_noise(self, detector_images, pure_signal, noise_only):
        return detector_images
        pass # TODO: implement combine_signal_and_noise

    
    def convert_to_dn(self, detector_images):
        return sunpy.map.MapSequence([map * (self.config.detector_gain * u.dN/u.dn * u.electron/u.count) # TODO: Remove all these unit gymnastics once sunpy issue has been resolved https://github.com/sunpy/sunpy/issues/6823
                                      for map in detector_images])
    
        # TODO: Clip to DN max (Alan set the gain so that pixel full well [33k] electrons is 90% of the ADC dynamic range); therefore, if the image has already been clipped to full well, this clipping function shouldn't ever do anything
        #return sunpy.map.MapSequence([map *  self.config.detector_gain # TODO: This is all that'll be needed once TODO above is addressed
        #                              for map in detector_images])

    
    def apply_screwy_pixels(self, detector_images, spike_frame):
        return detector_images
        pass # TODO: implement apply_screwy_pixels (spikes, dead pixels, hot pixels; ensure that none of these are above the full well)


class OnboardSoftware:
    def __init__(self, config):
        self.config = config
    

    def subtract_dark(self, detector_images, dark_frame):
        '''
        A dark frame map, drawn from a random distribution on the basis of a 
        config item for the average dark current in the detector as a function
        of temperature, is subtracted from each detector image. A map of
        detector images is returned
        

        Parameters
        ----------
        detector_images : [sunpy.map.MapSequence]
            A sunpy detector images map sequence
        dark_frame : [self.dark_frame]
            A config dark_frame map

        Returns
        ---
        detector_images : [sunpy.map.MapSequence]
            A dark current subtracted sunpy detector images map sequence 

        TODO
        ---
        TODO: Test to see if detector images are same dimension as the dark frame
        '''
        return detector_images
        #return sunpy.map.MapSequence([map - dark_frame / map.meta['EXPOSURE'] for map in detector_images]) # TODO: Won't work until we have something sensible returned for make_dark_frame()


    def apply_jitter(self, onboard_processed_images): # Note really onboard software, but this is the place in the logic it belongs
        return onboard_processed_images
        pass # TODO: implement apply_jitter (will be different for the short exposure central disk and long exposure off-disk)

    
    def median_image_stack(self, onboard_processed_images):
        # TODO: handle the median stacking
        
        map1 = onboard_processed_images[0]
        map2 = onboard_processed_images[1]
        composite_map = self.__composite_maps(map1, map2)
        
        return composite_map # TODO: figure out if this needs to be an array through time or not
        pass # TODO: impelement median_image_stack (returns a composite images merging the on- and off-disk, and spanning time up to the duration corresponding to how many images to stack (e.g., four 15-second exposures stacked for median will result in a 1-minute composite))
    
    
    def __composite_maps(self, map1, map2):
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
        #if isinstance(onboard_processed_images, 'sunpy.map.mapbase.GenericMap'):
        #    print("nope")
        #if type(onboard_processed_images) != 'sunpy.map.mapbase.GenericMap':
        #    print(type(onboard_processed_images))
        #    warnings.warn('bin_image module expected sunpy.map')

        # check if bin dimensions added
        if xbin == None:
            xbin = self.config.num_pixels_to_bin[0]
        if ybin == None:
            ybin = self.config.num_pixels_to_bin[1]

        # See input array shape
        array_shape=np.shape(onboard_processed_images.data)

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
            onboard_processed_images.meta['cdelt1'] = self.config.plate_scale.value * xbin # Note 1: this is only needed because sunpy (v4.0.1) resample updates dimensions but NOT plate scale
            onboard_processed_images.meta['cdelt2'] = self.config.plate_scale.value * ybin # Note 2: there is also risk here because a) the user must be responsible in the config file to ensure the image_dimensions and plate_scale are compatible, and b) the units they input for plate_scale must be the same as those already in the map

        return onboard_processed_images




if __name__ == "__main__":
    config_filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    config = config_parser.Config(config_filename)
    hardware = Hardware(config)
    print("See test_instrument.py for example of how to configure to run this code.")