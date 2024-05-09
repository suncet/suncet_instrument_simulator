"""
Calculate the point spread function (PSF) for the SunCET telescopes.
Kluged together from AIApy PSF generation function
"""
import numpy as np

import astropy.units as u

try:
    import cupy
    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False

__all__ = ["psf", "filter_mesh_parameters", "_psf"]

def filter_mesh_parameters(lambda_value, angle_arm=None, angles_focal_plane=None, lpi=None):
    """
    Geometric parameters for meshes in SunCET filters used to calculate the point
    spread function.

    Returns
    -------
    meshinfo : `dict`
        Dictionary with filter mesh information for each channel. Each channel
        entry then contains another dictionary with the following keys
        describing filter mesh properties of that channel
        (see Table 2 of [1]_):

        * `angle_arm`: Angles of the four entrance filter arms
        * `error_angle_arm`: Error in angle of the four entrance filter arms
        * `spacing_e`: Distance between diffraction spikes from entrance filter
        * `spacing_fp`: Distance between diffraction spikes from focal plane filter
        * `mesh_pitch`: Pitch of the mesh
        * `mesh_width`: Width of the mesh
        * `width`: Width applied to the Gaussian such that *after*
          convolution we have the proper width
          (:math:`4/3` at :math:`1/e` of max)
        * `Area`: Fractional filter transmitting area outside of grid

    References
    ----------
    .. [1] `Grigis, P., Su, Y., Weber M., et al., 2012,
            AIA PSF Characterization and Deconvolution
            <https://sohoftp.nascom.nasa.gov/solarsoft/sdo/aia/idl/psf/DOC/psfreport.pdf>`_

    See Also
    --------
    psf : Calculate the composite point spread function
    """
    #todo: Get key parameters from config instead of hard coding
    if lpi is None:
         lines_per_inch = 20.0 # units not supported
    else:
        lines_per_inch = lpi
    line_pitch_inches = 1.0/lines_per_inch # still no units
    inch_per_micron = 3.93701E-5
    line_pitch_micron = line_pitch_inches/inch_per_micron * u.um

    line_width = 65.0 * u.um

    lambda_ang = lambda_value * u.angstrom
    lambda_micron = lambda_ang.to(u.micron)

    A = (1 - line_width / line_pitch_micron) ** 2.

    theta_rad = 1 * lambda_micron / line_pitch_micron * u.radian
    theta_arcsec = theta_rad.to(u.arcsec)

    px_scale = 4.8 * u.arcsec/u.pixel
    focal_len = 300 * u.mm
    focal_plane_filter_dist = 3.5 * u.mm

    psf_80pct_arcsec = 25. * u.arcsec # PSF 80% encircled energy width, using typical value from Alan's analysis
    psf_80pct_px = psf_80pct_arcsec / px_scale  # get 80% encircled in pixels
    psf_sigma = psf_80pct_px  * 0.6178878 # Compute sigma from 80% encircled value

    if angle_arm is None:
        angle_arm = [0, 90] * u.deg
    if angles_focal_plane is None:
        angles_focal_plane = [90, 0] * u.deg

    return {
            "wavelength": lambda_ang,
            "angle_arm": angle_arm,
            "angles_focal_plane": angles_focal_plane,
            "error_angle_arm": [0.02, 0.02] * u.deg,
            "spacing_e": theta_arcsec/px_scale,
            "mesh_pitch": line_pitch_micron,
            "mesh_width": line_width,
            "spacing_fp": theta_arcsec/px_scale * focal_plane_filter_dist/focal_len,
            "width": psf_sigma,
            "CDELT": [px_scale * u.pixel, px_scale * u.pixel],
            "Area": A
    }


def psf(channel: u.angstrom, use_preflightcore=False, diffraction_orders=None, angle_arm=None, angles_focal_plane=None,
        lpi=None, output_size=None, use_gpu=True):
    r"""
    Calculate the composite PSF for a given channel, including diffraction and
    core effects.

    .. note:: This function has been adapted from
              `aia_calc_psf.pro <https://sohoftp.nascom.nasa.gov/solarsoft/sdo/aia/idl/psf/PRO/aia_calc_psf.pro>`_.

    .. note:: If the `~cupy` package is installed
              and your machine has an NVIDIA GPU, the PSF calculation will
              automatically be accelerated with CUDA. This can lead to
              several orders of magnitude in performance increase compared to
              pure `numpy` on a CPU.

    The point spread function (PSF) can be modeled as a 2D Gaussian function
    of the radial distance :math:`r` from the center,

    .. math::

        I(r, \theta) = I_0 \exp\left(\frac{-r^2}{2\sigma^2}\right)

    where,

    - :math:`I_0` : the intensity of a diffraction spike
    - :math:`r` : the radial distance from the center
    - :math:`\theta = m\lambda/d`
    - :math:`m` : diffraction order
    - :math:`\lambda` : the wavelength of light
    - :math:`\sigma` : width of Gaussian

    The intensity of a particular diffraction spike, :math:`I_0`, is given by,

    .. math::

        I_0 = \mathrm{sinc}^2\left(\frac{\theta w}{\lambda}\right)

    where,

    - :math:`w` : the width of the mesh wires
    - :math:`d` : spacing between two consecutive mesh wires

    The PSF for a given filter can then be calculated as,

    .. math::

        \mathrm{PSF} = \sum_{m=-\infty}^{+\infty}I_m(r,\theta)

    where, in practice, one can approximate the summation by simply summing
    over a sufficiently large number of diffraction orders. In this case, we
    sum from :math:`m=--100` to :math:`m=100`.

    Finally, the composite PSF of the entrance and focal plane filters is
    given by,

    .. math::

        \mathrm{PSF}_c = \left|\mathcal{F}\left\{
                            \mathcal{F}\{\mathrm{PSF}_f\}
                            \mathcal{F}\{\mathrm{PSF}_e\}
                          \right\}\right|

    where :math:`\mathcal{F}` denotes the Fourier transform,
    :math:`\mathrm{PSF}_f` is the PSF of the focal plane filter, and
    :math:`\mathrm{PSF}_e` is the PSF of the entrance filter. For a more
    detailed explanation of the PSF and the above calculation, see [1]_.

    Parameters
    ----------
    channel : `~astropy.units.Quantity`
        Wavelength of channel
    use_preflightcore : `bool`, optional
        If True, use the pre-flight values of the mesh width
    diffraction_orders : array-like, optional
        The diffraction orders to sum over. If None, the full
        range from -100 to +100 in steps of 1 will be used.
    use_gpu : `bool`, optional
        If True and `~cupy` is installed, do PSF deconvolution on the GPU
        with `~cupy`.

    Returns
    -------
    `~numpy.ndarray`
        The composite PSF of the entrance and focal plane filters.

    See Also
    --------
    filter_mesh_parameters
    deconvolve

    References
    ----------
    .. [1] `Grigis, P., Su, Y., Weber M., et al., 2012,
            AIA PSF Characterization and Deconvolution
            <https://sohoftp.nascom.nasa.gov/solarsoft/sdo/aia/idl/psf/DOC/psfreport.pdf>`__
    """
    if output_size is None:
        output_size = [1500, 1500]
    meshinfo = filter_mesh_parameters(channel, angle_arm = angle_arm, angles_focal_plane= angles_focal_plane, lpi = lpi)
    angles_entrance = meshinfo["angle_arm"]
    angles_focal_plane = meshinfo["angles_focal_plane"]
    if diffraction_orders is None:
        diffraction_orders = np.arange(-1500, 1500, 1)
    psf_entrance = _psf(meshinfo, angles_entrance, diffraction_orders, output_size, use_gpu=use_gpu)
    psf_focal_plane = _psf(
        meshinfo,
        angles_focal_plane,
        diffraction_orders,
        output_size,
        focal_plane=True,
        use_gpu=use_gpu,
    )
    # Composite PSF
    psf = abs(np.fft.fft2(np.fft.fft2(psf_focal_plane) * np.fft.fft2(psf_entrance)))
    # Center PSF in the middle of the image
    psf = np.roll(np.roll(psf, psf.shape[1] // 2, axis=1), psf.shape[0] // 2, axis=0)
    # Normalize by total number of pixels
    psf = psf / (psf.shape[0] * psf.shape[1])
    # If using cupy, cast back to a normal numpy array
    if HAS_CUPY and use_gpu:
        psf = cupy.asnumpy(psf)
    return psf


def _psf(meshinfo, angles, diffraction_orders, output_size, focal_plane=False, use_gpu=True):
    psf = np.zeros((output_size[0], output_size[1]), dtype=float)
    # If cupy is available, cast to a cupy array
    if HAS_CUPY and use_gpu:
        psf = cupy.array(psf)
    Nx, Ny = psf.shape
    width_x = meshinfo["width"].value
    width_y = meshinfo["width"].value
    # x and y position grids
    x = np.outer(np.ones(Ny), np.arange(Nx) + 0.5)
    y = np.outer(np.arange(Ny) + 0.5, np.ones(Nx))
    if HAS_CUPY and use_gpu:
        x = cupy.array(x)
        y = cupy.array(y)
    area_not_mesh = meshinfo["Area"]
    spacing = meshinfo["spacing_fp"] if focal_plane else meshinfo["spacing_e"]
    mesh_ratio = (meshinfo["mesh_pitch"] / meshinfo["mesh_width"]).decompose().value
    spacing_x = spacing * np.cos(angles)
    spacing_y = spacing * np.sin(angles)
    for order in diffraction_orders:
        if order == 0:
            continue
        intensity = np.sinc(order / mesh_ratio) ** 2  # I_0
        for dx, dy in zip(spacing_x.value, spacing_y.value):
            x_centered = x - (0.5 * Nx + dx * order + 0.5)
            y_centered = y - (0.5 * Ny + dy * order + 0.5)
            # NOTE: this step is the bottleneck and is VERY slow on a CPU
            psf += np.exp(-width_x * x_centered * x_centered - width_y * y_centered * y_centered) * intensity
    # Contribution from core
    psf_core = np.exp(-width_x * (x - 0.5 * Nx - 0.5) ** 2 - width_y * (y - 0.5 * Ny - 0.5) ** 2)
    psf_total = (1 - area_not_mesh) * psf / psf.sum() + area_not_mesh * psf_core / psf_core.sum()
    return psf_total
