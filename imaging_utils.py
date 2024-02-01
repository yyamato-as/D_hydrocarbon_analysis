import sys
sys.path.append("/home/yamato/Application/reduction_tools/")
import reduction_utils_mpi as utils
import casatasks
import casatools
import shutil
import os
import numpy as np
import analysis_utils as au
from astropy.convolution import convolve_fft
from scipy import optimize

tb = casatools.table()
ia = casatools.image()

def rad2arcsec(val):
    return val * 180 * 3600 / np.pi

def get_beam_info(imagename):
    ia.open(imagename)
    beam = ia.summary()["restoringbeam"]
    ia.close()

    bmaj = beam["major"]["value"]
    bmin = beam["minor"]["value"]
    bpa = beam["positionangle"]["value"]

    return bmaj, bmin, bpa

def get_xy_1D_coord(psf, imsize):
    peak_x, peak_y = np.unravel_index(np.argmax(psf), psf.shape)
    x = -(np.arange(imsize) - peak_x)
    y = np.arange(imsize) - peak_y
    return x, y

def get_xy_2D_coord(psf, imsize):
    x, y = get_xy_1D_coord(psf, imsize)
    xx, yy = np.meshgrid(x, y, indexing="ij") # indexing="ij" is needed to match the axis order to that of casa image
    return xx, yy

def rotate_coord(x, y, theta):
    rotangle = np.deg2rad(theta)
    x_rot = x * np.cos(rotangle) + y * np.sin(rotangle)
    y_rot = x * np.sin(rotangle) - y * np.cos(rotangle)
    return x_rot, y_rot

def Gaussian2D(x, y, sigma_x, sigma_y, PA):
    theta = 90. - PA # from negative x to positive y
    x_rot, y_rot = rotate_coord(x, y, theta)
    g = 1.*np.exp(-((x_rot/sigma_x)**2/2. + (y_rot/sigma_y)**2/2.))
    return g

def get_nulled_psf(psf_array, r, theta, ntheta=24):

    # azimuthal coordinate
    da = 2 * np.pi / ntheta
    azim = np.arange(-np.pi, np.pi, da)

    nulled_psf = psf_array.copy()

    for a in azim:
        # some overlapping azimuth: see https://github.com/MPoL-dev/MPoL/blob/main/src/mpol/gridding.py#L441
        wedge = (theta >= a - 0.3 * da) & (theta <= a + 1.3 * da)
        wedge_is_neg = wedge & (psf_array < 0)
        first_null = np.min(r[wedge_is_neg])
        nulled_psf[wedge & (r >= first_null)] = 0.0

    return nulled_psf

def calc_JvM_epsilon(xx, yy, target_beam_array, tapered_psf_array):
    r, theta = np.sqrt(xx**2 + yy**2), np.arctan2(yy, xx)
    epsilon = np.sum(target_beam_array) / np.sum(get_nulled_psf(tapered_psf_array, r, theta))
    print("JvM epsilon: {:.4f}".format(epsilon))
    return epsilon

def calc_uvtaper(vis, imagename, field="", specmode="mfs", weighting="superuniform", perchanweightdensity=True, cellsize="0.1arcsec", imsize=None, target_beam=(0.1, 0.1, 0.0), **tclean_kwargs):
    if not isinstance(target_beam, (tuple, list)):
        target_beam = (target_beam, target_beam, 0.0)
    target_bmaj, target_bmin, target_bpa = target_beam

    print("Starting to estimate the uvtaper which achieve the target beam of {:.4f} arcsec x {:.4f} arcsec (P.A. = {:.4f} deg)".format(*target_beam))

    # calculate image size
    if imsize is None:
        imsize = utils.get_imsize(vis=vis, field=field, cellsize=cellsize)
    
    print(f"Pixel size: {cellsize}")
    print(f"Image size: {imsize}")
    dpix = float(cellsize.replace("arcsec", ""))

    # calculate and fit psf
    msg = "Calculating and fitting PSF with the original weighting scheme (weighting = {:s}".format(weighting)
    if "briggs" in weighting:
        robust = tclean_kwargs.get("robust", 0.5)
        msg += f", robust = {robust}"
    msg += ")"
    print(msg)

    if specmode == "cube":
        tb.open(vis + "/SPECTRAL_WINDOW")
        nchan = tb.getcol("NUM_CHAN")
        tb.close()

        tclean_kwargs["nchan"] = 1
        tclean_kwargs["start"] = int(nchan/2)

    for ext in [".psf", ".pb", ".sumwt"]:
        if os.path.exists(imagename + ext):
            shutil.rmtree(imagename + ext)
    casatasks.tclean(
        vis=vis,
        imagename=imagename,
        specmode=specmode,
        weighting=weighting,
        perchanweightdensity=perchanweightdensity,
        cell=cellsize,
        imsize=imsize,
        calcres=False,
        calcpsf=True,
        **tclean_kwargs
    )
    print("Done.")

    # get the original psf array
    ia.open(imagename + ".psf")
    psf_orig = ia.getregion().squeeze()
    print(psf_orig.shape)
    ia.close()
    
    # if psf_orig.ndim > 2:
    #     psf_orig = psf_orig[]

    # fetch beam info
    bmaj, bmin, bpa = get_beam_info(imagename + ".psf")
    print(
        "Restoring beam shape: {:.4f} arcsec x {:.4f} arcsec (P.A. = {:.4f} deg)".format(
            bmaj, bmin, bpa
        )
    )

    # window out the central region for computational efficiency while keeping enough region for convolution to work well
    print("Windowing out the central region for computational efficiency")
    window_factor = 20 # factor to determine the window; multiplied to the target major axis
    start_idx = int(imsize*0.5) - int(target_bmaj * window_factor / dpix)
    end_idx = int(imsize*0.5) + int(target_bmaj * window_factor / dpix)
    psf_windowed = psf_orig[start_idx:end_idx, start_idx:end_idx]
    imsize_windowed = psf_windowed.shape[0]
    print(f"Image size of windowed psf: {imsize_windowed}")

    # minimizing the metric
    xx, yy = get_xy_2D_coord(psf_windowed, imsize_windowed)
    target_beam = Gaussian2D(xx, yy, sigma_x=au.FWHM_to_sigma(target_bmaj/dpix), sigma_y=au.FWHM_to_sigma(target_bmin/dpix), PA=target_bpa)
    psf_cutoff = tclean_kwargs.get("psf_cutoff", 0.35)

    def beam_chi2(params):
        tmaj, tmin, tpa = params

        kernel = Gaussian2D(xx, yy, sigma_x=au.FWHM_to_sigma(tmaj/dpix), sigma_y=au.FWHM_to_sigma(tmin/dpix), PA=tpa)
        tapered_beam = convolve_fft(psf_windowed, kernel=kernel)
        tapered_beam /= np.max(tapered_beam)

        metric = np.sum((target_beam[tapered_beam > psf_cutoff] - tapered_beam[tapered_beam > psf_cutoff])**2)

        return metric

    print("Calculating uvtaper parameter")
    x0 = (target_bmaj, target_bmin, target_bpa)
    res = optimize.minimize(beam_chi2, x0=x0, method="Nelder-Mead")
    uvtaper = [str(res.x[0]) + "arcsec", str(res.x[1]) + "arcsec", str(res.x[2]) + "deg"]
    print(f"Done. Best-fit uvtaper parameter: " + str(uvtaper))

    # check JvM epsilon
    print("Calculating the resulting beam shape after uvtaper")
    for ext in [".psf", ".pb", ".sumwt"]:
        shutil.rmtree(imagename + ext)
    casatasks.tclean(
        vis=vis,
        imagename=imagename,
        specmode=specmode,
        weighting=weighting,
        perchanweightdensity=perchanweightdensity,
        uvtaper=uvtaper,
        cell=cellsize,
        imsize=imsize,
        calcres=False,
        calcpsf=True,
        **tclean_kwargs
    )
    print("Done.")

    ia.open(imagename + ".psf")
    psf_tapered = ia.getregion().squeeze()
    ia.close()

     # fetch beam info
    bmaj, bmin, bpa = get_beam_info(imagename + ".psf")
    print(
        "Restoring beam shape after uvtaper: {:.4f} arcsec x {:.4f} arcsec (P.A. = {:.4f} deg)".format(
            bmaj, bmin, bpa
        )
    )

    xx, yy = get_xy_2D_coord(psf_tapered, imsize)
    target_beam = Gaussian2D(xx, yy, sigma_x=au.FWHM_to_sigma(target_bmaj/dpix), sigma_y=au.FWHM_to_sigma(target_bmin/dpix), PA=target_bpa)
    eps = calc_JvM_epsilon(xx, yy, target_beam_array=target_beam, tapered_psf_array=psf_tapered)
    if np.abs(eps - 1.0) > 0.1:
        print("Warning: dirty psf is deviated from the CLEAN beam by > 10%")

    return uvtaper


    
