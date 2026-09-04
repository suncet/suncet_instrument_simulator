# SunCET Instrument Simulator

## Simulator timestamps and existing files

Radiance-map `DATE-OBS` is the model epoch plus the filename's frame index times
`model_timestep`. The three-viewpoint epoch is `2012-03-08T20:00:00.000` UTC;
the bright-fast, bright-slow, and dimmest epoch is `2011-02-15T17:00:00.000` UTC.
Instrument integrations inherit the source map's time, with an offset only when
their start differs from that model frame. Final `DATE-BEG` and `DATE-OBS` come
from the first integration; `DATE-END` adds the full observation window. The
metadata CSV's typical timestamps are not used.

These one-off utilities default to inspection only. Run them in the simulator's
Python environment, with `suncet_data` set to your data directory:

```bash
python repair_radiance_map_timestamps.py "$suncet_data/mhd"
python repair_simulator_output_timestamps.py suncet_instrument_simulator/config_files/config_three_viewpoint_limb.ini
```

Append `--write` to apply repairs after inspection. The radiance utility assumes
a 10-second cadence unless `--cadence-seconds` is supplied, edits every wavelength
HDU, and saves compressed original headers beside each repaired file. To undo
a repair, pass its `.timestamp-headers.json.gz` backup with `--restore`.
Nonstandard filenames are reported and excluded from the recursive scan.

The output utility uses each supplied config's schedule and the corresponding
cached radiance timestamps. It repairs time cards and filenames, verifies image
data/checksums, and backs up originals under
`$suncet_data/synthetic/level0/timestamp_repair_backups/`. If a correct output
already occupies the destination, the older duplicate is archived there.
Use the configuration settings from the original run when repairing outputs.

## Make realistic Level 0.5 products

`make_level0_5.py` wraps one or more synthetic FITS images in the same FITS-header
and companion JSON metadata layout currently written by
`suncet_processing_pipeline.make_level0_5_slowly_built_up`:

```bash
python make_level0_5.py synthetic/level0/fits/*.fits \
  --image-id-start 10000 \
  --suffix=-synthetic
```

Output defaults to `$suncet_data/synthetic/level0_5`. Use `--output-dir` only when
an explicit alternate destination is needed.

Image dimensions and `NBIN1`/`NBIN2` are used to populate the CSIE metadata;
`EXPTIME` is also used when present. For telemetry values that cannot be produced
by the simulator, the bundled CSIE v1.1.2 template contains representative flight
values. To use housekeeping values from a particular real exposure, pass its
pipeline-produced JSON file:

```bash
python make_level0_5.py synthetic.fits \
  --metadata-template image_4112_meta.json --image-id-start 10000
```

`--row-bin`, `--col-bin`, `--exposure-ms`, and `--encoding` can override inferred
values. Inputs must contain a 2-D numeric image whose values fit in unsigned 16 bits.

![Python 3.11](https://img.shields.io/badge/python-3.11-blue)
[![Tests](https://github.com/suncet/suncet_instrument_simulator/actions/workflows/tests.yml/badge.svg)](https://github.com/suncet/suncet_instrument_simulator/actions/workflows/tests.yml)

[![Example of synthetic SunCET observations generated with this code](https://img.youtube.com/vi/GFjgdv2Y7gY/0.jpg)](https://www.youtube.com/watch?v=GFjgdv2Y7gY)


## Overview
This  code ingests the output of solar simulations, "follows the photons" through the telescope to the detector accounting for physical effects, then follows the electron-hole pairs generated in the detector, and finally applies flight-like onboard processing to make the data look like how it is stored onboard. We published more detail about this process in [Mason et al. 2022](http://dx.doi.org/10.3847/1538-4357/ac33a1). 
Another code repository includes our data processing pipeline, which takes the data from this stage up through all of our data product levels. 

## Inputs
* `em_map_*.sav`: These emission measure (EM) maps were provided by Dr. Meng Jin. These data are the result of magnetohydrodynamic (MHD) simulations of coronal mass ejections (CMEs) on the sun, following the procedure described in [Jin et al. 2017](https://doi.org/10.3847/1538-4357/834/2/173), which uses the Alfvén Wave Solar Model ([AWSoM; van der Holst et al. 2014](https://doi.org/10.1088/0004-637X/782/2/81)) within the Space Weather Modeling Framework ([Tóth et al. 2012](https://doi.org/10.1016/j.jcp.2011.02.006)). The files are IDL savesets. Each file corresponds to a step in time. Within each file is a set of 1024x1024 maps at 12 different temperatures, from logT = 5.4-5.5, 5.5-5.6, ... , 6.5-6.6. The emission measure is in units of cm<sup>-5</sup>.
* `config_*.ini`: This file contains the configuration for the synthetic instrument. The base file `config.ini` is an example, with each `_*.ini` file being a particular configuration summarized in the filename. The file configures things like aperture size, image dimensions, short/long exposure times, and the number of images to stack for taking a median. 

## Outputs
Note: the below is very notional as of 2022-06-13. 
* Rendered extreme ultraviolet (EUV) maps: This is an intermediate output that only needs to be generated once. These files are the result of converting the `em_map_*.sav` files into EUV light with units of erg/cm<sup>2</sup>/s/pixel. 
* `signal.?`: This file contains the clean signal passed through the full instrument, for all times passed in. 
* `noise.?`: Same as above, except for the pure noise. 
* `synthetic_data.?`: Same as above, except this is the realsitic combination of the signal and the noise. 
* `snr.?`: Same as above, except this contains the signal to noise ratio (SNR) data for generating contour maps. 
* `synthetic_data_*.png`: For each step in time, an image is rendered and saved to disk for quick reference. The image scaling, filtering, and color table are defined by `config.ini`. 
* `snr_*.png`: For each step in time, a set of SNR contours (at levels defined by `config.ini`) are output for quick reference. 
* `synthetic_data.mov`: The individual images rendered into a movie. 
