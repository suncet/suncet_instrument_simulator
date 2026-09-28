# SunCET color table study

Open index.html locally. There are 40 distinct tables, each rendered at 1000 x 750 pixels with two shared stretches. The gallery starts with the current PNG stretch.

Input: `config_default_OBS_2023-01-14T17:00:00.000_300.fits`

The current reference is pixel-identical to make_movie.plot_scaled_image(scale='1/4'). All images use the existing negative-value handling, sigma=300 radial filter, orientation, and borderless renderer. The optional asinh images use the current renderer's single-frame 1st/99.7th-percentile limits and 3% softening, fixed across all tables. These are not movie-wide sampled limits. The dark ring comes from the input exposure composite.

Standard mission tables come from SunPy. Shared mission tables are combined into one option: AIA171/SUVI171/EUI174, AIA193/SUVI195, AIA335/SUVI284, AIA304/SUVI304/EUI304, AIA131/SUVI131, and AIA94/SUVI94. EUI-inspired amber is custom. Instrument names describe color tables, not a change in the simulated passband.

Source: https://docs.sunpy.org/en/stable/_modules/sunpy/visualization/colormaps/cm.html

SunCET brand anchors use the PDF's printed hex values (some printed RGB numbers differ). Poster anchors are exact frequently occurring RGB pixels from SunCET Poster v3.png. Custom tables are interpolated in sRGB and are aesthetic candidates, not asserted to be perceptually uniform.

The manifest records color anchors and rendering settings; luts/ contains reusable 256-row RGB CSVs. Load a CSV with matplotlib.colors.ListedColormap(np.loadtxt(path, delimiter=',', skiprows=1)/255). The CSVs reproduce the exported 8-bit color mapping. No production defaults were changed.
