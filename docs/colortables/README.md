# SunCET color table study

Open index.html locally. There are 40 distinct tables, each rendered at 1000 x 750 pixels with two shared stretches. The gallery starts with log10.

Input: `config_default_OBS_2023-01-14T17:00:00.000_300_level1_v2.0.0_provisional.fits`

The input is the pipeline's provisional exposure-normalized Level 1 frame (DN/s). Both stretches use 43.2777 DN/s to the image maximum, with no radial filter. Effective exposures are recorded in manifest.json; this is not full Level 1 calibration. Images retain the existing orientation and borderless renderer. Log10 and asinh (200 DN/s softening) match local experiment options 1 and 7. The directory key 'current' is retained for compatibility and now means log10.

Standard mission tables come from SunPy. Shared mission tables are combined into one option: AIA171/SUVI171/EUI174, AIA193/SUVI195, AIA335/SUVI284, AIA304/SUVI304/EUI304, AIA131/SUVI131, and AIA94/SUVI94. EUI-inspired amber is custom. Instrument names describe color tables, not a change in the simulated passband.

Source: https://docs.sunpy.org/en/stable/_modules/sunpy/visualization/colormaps/cm.html

SunCET brand anchors use the PDF's printed hex values (some printed RGB numbers differ). Poster anchors are exact frequently occurring RGB pixels from SunCET Poster v3.png. Custom tables are interpolated in sRGB and are aesthetic candidates, not asserted to be perceptually uniform.

The manifest records color anchors and rendering settings; luts/ contains reusable 256-row RGB CSVs. Load a CSV with matplotlib.colors.ListedColormap(np.loadtxt(path, delimiter=',', skiprows=1)/255). The CSVs reproduce the exported 8-bit color mapping. No production defaults were changed.
