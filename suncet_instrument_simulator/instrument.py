import os
import warnings
from glob import glob
import numpy as np
from pandas import read_fwf, read_csv
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import constants as const
import sunpy.map
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
            map_list.append(map.resample(self.config.image_dimensions, method='nearest')) # TODO: Figure out why the CDELTs in map.fits_header['CDELT1'] don't update with the resample
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
        return sunpy.map.MapSequence([map*self.effective_area for map in radiance_maps])
        

    
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
            map.fits_header['EXPTIME']=self.config.exposure_time_long
            long_exposure_map_list.append(map*self.config.exposure_time_long)

            map.fits_header['EXPTIME']=self.config.exposure_time_short
            short_exposure_map_list.append(map*self.config.exposure_time_short)
        
        return sunpy.map.MapSequence(long_exposure_map_list), sunpy.map.MapSequence(short_exposure_map_list)
    

    def apply_photon_shot_noise(self, radiance_maps):
        map_list = []
        for map in radiance_maps:
            noisy_photons = np.random.poisson(lam=map.data)
            map_list.append(sunpy.map.Map(noisy_photons, map.meta))
        
        return sunpy.map.MapSequence(map_list)


    def convert_to_electrons(self, radiance_maps, apply_noise=True): # TODO: What to do with photon_shot_noise?
        quantum_yield = self.__compute_quantum_yields()
        detector_image = radiance_maps[0].data * 0 # empty array
        for i, map in enumerate(radiance_maps): 
            #detector_images.append(map * quantum_yield[i]) # TODO: Update once this issue has been resolved https://github.com/sunpy/sunpy/issues/6823
            quantum_yield_units_hacked = quantum_yield[i] * u.count/u.electron
            detector_image = detector_image + (map.data * map.unit) * quantum_yield_units_hacked # TODO: Should really be in electrons, not counts
        detector_image = sunpy.map.Map(detector_image, radiance_maps[0].meta) # TODO: Update any of the meta? Wavelength?

        if apply_noise:
            detector_image = self.__apply_fano_noise(detector_image)

        detector_image = self.__clip_at_full_well(detector_image)
        return detector_image
    
    
    def __compute_quantum_yields(self):
        photoelectron_in_silicon = 3.63 * u.eV * u.ph / u.electron # the magic number 3.63 here the typical energy [eV] to release a photoelectron in silicon
        return (const.h * const.c / self.wavelengths.to(u.m)).to(u.eV) / photoelectron_in_silicon 


    def __apply_fano_noise(self, detector_image):
        fano_factor = 0.1190 # Based on material, silicon in this case, and fairly insensitive to wavelength. This number comes from Rodrigues+ 2021 (https://doi.org/10.1016/j.nima.2021.165511)
        header = detector_image.meta
        detector_image = np.random.poisson(lam=detector_image.data) * np.sqrt(fano_factor)
        return sunpy.map.Map(detector_image, header)


    def __clip_at_full_well(self, detector_image):
        mask = detector_image.data > self.config.pixel_full_well.value
        detector_image.data[mask] = self.config.pixel_full_well
        
        if detector_image.unit != self.config.pixel_full_well.unit: 
            unit_mismatch = True
        if unit_mismatch: 
            warnings.warn('The units for the detector image does not match the units of the pixel full well.')
        
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
        return sunpy.map.MapSequence([map *  (self.config.detector_gain * u.dN/u.dn * u.electron/u.count) # TODO: Remove all these unit gymnastics once sunpy issue has been resolved https://github.com/sunpy/sunpy/issues/6823
                                      for map in detector_images])
    
        # TODO: Clip to DN max (Alan set the gain so that pixel full well [33k] electrons is 90% of the ADC dynamic range); therefore, if the image has already been clipped to full well, this clipping function shouldn't ever do anything
        #return sunpy.map.MapSequence([map *  self.config.detector_gain # TODO: This is all that'll be needed once TODO above is addressed
        #                              for map in detector_images])

    
    def apply_screwy_pixels(self, detector_images, spike_frame):
        pass # TODO: implement apply_screwy_pixels (spikes, dead pixels, hot pixels; ensure that none of these are above the full well)


class OnboardSoftware:
    def __init__(self, config):
        self.config = config
    

    def subtract_dark(self, detector_images, dark_frame):
        pass # TOOD: implement subtract_dark

    
    def separate_images(self, onboard_processed_images):
        pass # TODO: implement separate_images


    def apply_jitter(self, split_images): # Note really onboard software, but this is the place in the logic it belongs
        pass # TODO: implement apply_jitter (will be different for the short exposure central disk and long exposure off-disk)

    
    def median_image_stack(self, split_images):
        pass # TODO: impelement median_image_stack (returns a composite images merging the on- and off-disk, and spanning time up to the duration corresponding to how many images to stack (e.g., four 15-second exposures stacked for median will result in a 1-minute composite))
    
    
    def bin_image(self, onboard_processed_images):
        pass # TODO: implement bin_image

if __name__ == "__main__":
    config_filename = os.getcwd() + '/suncet_instrument_simulator/config_files/config_default.ini'
    config = config_parser.Config(config_filename)
    hardware = Hardware(config)
    print("See test_instrument.py for example of how to configure to run this code.")