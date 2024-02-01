# %%
import sys
sys.path.append("/home/yamato/Application/reduction_tools/")
sys.path.append("/home/yamato/Application/keplerian_mask/")
import reduction_utils_mpi as utils
from keplerian_mask import make_mask
import casatasks
import casatools
from linedictionary_v2 import line_dict
from diskdictionary import disk_dict
import shutil
import os
import analysis_utils as au
import imaging_utils as iutils
import numpy as np

tb = casatools.table()

source = "MWC_480"
msdir = "/works/yamato/MAPS_data/v1_product/measurement_set/"
vispath = "./data/measurement_set/"
impath = "./data/image/"
EBs = ["SB1", "SB2", "LB1", "LB2", "LB3", "LB4", "LB5", "LB6"]
ms_list = [msdir + f"{source}_{i}_spectral_line.ms" for i in EBs]
ms_wcont_list = [msdir + f"{source}_{i}_spectral_line_wcont.ms" for i in EBs]
lines = [f"C2H_3-2_hfs{i+1}" for i in range(5)]
velocity_range = 15  # maximum velocity w.r.t. vsys in km/s
velocity_width = 0.2

disk = au.DiskProperty(
    PA=disk_dict[source]["PA_gofish"],
    incl=disk_dict[source]["incl"],
    Mstar=disk_dict[source]["M_star"],
    vsys=disk_dict[source]["v_sys"],
    distance=disk_dict[source]["distance"],
)

# %%
# split out the respect spectral window and then concatenate
for line in lines:
    print("Processing " + line)
    spw = line_dict[source][line]["spw"]
    print("Transitions: ", line)
    print("SPW: ", spw)
    ms_list_split = []
    # split
    for ms in ms_list:
        print("Input: ", ms)
        outputvis = vispath + os.path.basename(ms).replace("_spectral_line.ms", f"_{line}.ms")
        print("Output: ", outputvis)
        os.system("rm -r " + outputvis)
        casatasks.split(vis=ms, outputvis=outputvis, spw=spw, datacolumn="data")
        ms_list_split.append(outputvis)

    # concat
    concatvis = vispath + f"{source}_{line}.ms"
    print("Concatenated MS: ", concatvis)
    os.system("rm -r " + concatvis)
    casatasks.concat(vis=ms_list_split, concatvis=concatvis)
    for vis in ms_list_split:
        shutil.rmtree(vis)

    # cvel
    print("cveled MS: ", concatvis + ".cvel")
    os.system("rm -r " + concatvis + ".cvel")
    casatasks.cvel2(
        vis=concatvis,
        outputvis=concatvis + ".cvel",
        mode="velocity",
        restfreq=line_dict[source][line]["restfreq"],
        outframe="LSRK",
        veltype="radio",
        width=f"{velocity_width}km/s",
        start=f"{disk.vsys - velocity_range}km/s",
        nchan=int(2 * velocity_range / velocity_width), 
    )

    # set the rest frequency
    tb.open(concatvis + ".cvel/SOURCE", nomodify=False)
    nu0 = float(line_dict[source][line]["restfreq"].replace("GHz", "")) * 1e9
    tb.putcol("REST_FREQUENCY", np.array([[nu0]]))
    tb.close()

# %%
# first make dirty image and then CLEAN
### basic parameters for tclean ###
cellsize = "0.02arcsec"
weighting = "superuniform"
scales = [0, 10, 30]
npixels = 6
target_bmaj = 0.3
deconvolver = "multiscale"
spw = "0"

### parameters for Keplerian mask ###
r_min = 0.25
r_max = 2.5
zr = 0.0
nbeams = 1.5

for line in lines:
    print("Processing " + line)
    # spw = line_dict[source][line]["spw"]
    # first calculate the taper parameter to acheive the desired beams
    vis = vispath + f"{source}_{line}.ms.cvel"
    imagename = impath + f"{source}_{line}_hybrid"
    imsize = utils.get_imsize(
        vis=vis,
        field=source,
        cellsize=cellsize,
    )

    for ext in [".alpha", ".image", ".mask", ".model", ".pb", ".psf", ".residual", ".sumwt"]:
        if os.path.exists(imagename + ext):
            shutil.rmtree(imagename + ext)
            
    uvtaper = iutils.calc_uvtaper(
        vis=vis,
        imagename=imagename,
        field=source,
        specmode="mfs",
        weighting=weighting,
        npixels=npixels,
        cellsize=cellsize,
        target_beam=(target_bmaj, target_bmaj * np.cos(np.deg2rad(disk.incl)), disk.PA),
    )

    for ext in [".alpha", ".image", ".mask", ".model", ".pb", ".psf", ".residual", ".sumwt"]:
        if os.path.exists(imagename + ext):
            shutil.rmtree(imagename + ext)

    casatasks.tclean(
        vis=vis,
        imagename=imagename,
        spw=spw,
        specmode="cube",
        outframe="LSRK",
        veltype="radio",
        imsize=imsize,
        cell=cellsize,
        scales=scales,
        weighting=weighting,
        npixels=npixels,
        uvtaper=uvtaper,
        restoringbeam="common",
        niter=0
    )

    # make mask
    rms = make_mask(
        imagename=imagename + ".image",
        inc=disk.incl,
        PA=disk.PA,
        mstar=disk.Mstar,
        dist=disk.distance,
        vlsr=disk.vsys * 1e3,
        nbeams=nbeams,
        r_min=r_min,
        r_max=r_max,
        zr=zr,
        restfreqs=line_dict[source][line]["maskfreqs"],
        export_FITS=True,
        overwrite=True,
    )
    # shutil.rmtree(imagename + ".mask.image")

    # CLEAN
    # rms = utils.calc_sensitivity(
    #     vis,
    #     cellsize=cellsize,
    #     imsize=imsize,
    #     weighting=weighting,
    #     npixels=npixels,
    #     specmode="cube",
    #     spw=[spw],
    #     chan=int(0.5 * (2 * velocity_range / velocity_width)),
    # )

    for ext in [".alpha", ".image", ".mask", ".model", ".pb", ".psf", ".residual", ".sumwt"]:
        if os.path.exists(imagename + ext):
            shutil.rmtree(imagename + ext)

    casatasks.tclean(
        vis=vis,
        imagename=imagename,
        spw=spw,
        specmode="cube",
        outframe="LSRK",
        veltype="radio",
        imsize=imsize,
        cell=cellsize,
        scales=scales,
        weighting=weighting,
        npixels=npixels,
        uvtaper=uvtaper,
        restoringbeam="common",
        niter=1000000,
        threshold=f"{3*rms}Jy",
        usemask="user",
        mask=imagename + ".mask.image"
    )

    # pbcor
    print("Correcting primary beam response")
    casatasks.impbcor(
        imagename=imagename + ".image",
        pbimage=imagename + ".pb",
        outfile=imagename + ".image.pbcor",
        overwrite=True,
    )

    # export fits
    print("Exporting to FITS")
    for ext in [".image", ".image.pbcor"]:
        casatasks.exportfits(
            imagename=imagename + ext,
            fitsimage=imagename + ext + ".fits",
            dropstokes=True,
            overwrite=True,
        )

# %%



