import os
import sys

sys.path.append("/home/yamato/Project/hydrocarbon_deuteration/")

execfile("/home/yamato/Project/hydrocarbon_deuteration/reduction_utils_mpi_v3.py")
# import reduction_utils_mpi as ru

delivered_data_path_SB = "/works/yamato/MAPS_data/2018.1.01055.L/science_goal.uid___A001_X133d_X19de/group.uid___A001_X133d_X19df/member.uid___A001_X133d_X19e2/calibrated/"
delivered_data_path_LB = "/works/yamato/MAPS_data/2018.1.01055.L/science_goal.uid___A001_X133d_X19de/group.uid___A001_X133d_X19df/member.uid___A001_X133d_X19e0/calibrated/"
processingdir = "/works/yamato/MAPS_data/processing/Band6/"
MSdir = "/works/yamato/MAPS_data/v1_product/measurement_set/"
imagedir = "/works/yamato/MAPS_data/v1_product/image/"
field = "MWC_480"
prefix = "MWC_480"
data_params_json_name = processingdir + prefix + "_data_params.json"

try:
    data_params = load_data_params(data_params_json_name)

except FileNotFoundError:
    data_params = {
        "SB1": {
            "vis": delivered_data_path_SB + "uid___A002_Xd42ec5_Xc98c.ms",
            "spws": "25,27,29,31,33,35",
            "field": field,
            "column": "corrected",
        },
        "SB2": {
            "vis": delivered_data_path_SB + "uid___A002_Xd42ec5_Xd87d.ms",
            "spws": "25,27,29,31,33,35",
            "field": field,
            "column": "corrected",
        },
        "LB1": {
            "vis": delivered_data_path_LB + "uid___A002_Xdfcc3f_X3b32.ms",
            "spws": "25,27,29,31,33,35",
            "field": field,
            "column": "corrected",
        },
        "LB2": {
            "vis": delivered_data_path_LB + "uid___A002_Xe014a2_X3caa.ms",
            "spws": "25,27,29,31,33,35",
            "field": field,
            "column": "corrected",
        },
        "LB3": {
            "vis": delivered_data_path_LB + "uid___A002_Xe014a2_X414a.ms",
            "spws": "25,27,29,31,33,35",
            "field": field,
            "column": "corrected",
        },
        "LB4": {
            "vis": delivered_data_path_LB + "uid___A002_Xe03886_X11879.ms",
            "spws": "25,27,29,31,33,35",
            "field": field,
            "column": "corrected",
        },
        "LB5": {
            "vis": delivered_data_path_LB + "uid___A002_Xe03886_X1781.ms",
            "spws": "25,27,29,31,33,35",
            "field": field,
            "column": "corrected",
        },
        "LB6": {
            "vis": delivered_data_path_LB + "uid___A002_Xe03886_Xc8ff.ms",
            "spws": "25,27,29,31,33,35",
            "field": field,
            "column": "corrected",
        },
    }


# split out the scientific data
for i in data_params.keys():
    outputvis = processingdir + prefix + f"_{i}.ms"
    os.system("rm -r " + outputvis)
    split(
        vis=data_params[i]["vis"],
        outputvis=outputvis,
        spw=data_params[i]["spws"],
        field=data_params[i]["field"],
        datacolumn=data_params[i]["column"],
    )
    data_params[i]["vis_split"] = outputvis
    data_params[i]["spws"] = ",".join(
        str(i) for i in range(len(data_params[i]["spws"].split(",")))
    )

save_data_params(data_params, data_params_json_name)

# contspws = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

# first make dirty image for spw 0 (originally 25) for visual inspection
spws = ["1", "2", "3", "4", "5"]
cellsize = "0.01arcsec"
scales = [0, 5, 15, 25]
# imsize = get_imsize(vis=[data_params[i]["vis_split"] for i in data_params.keys()], field=field, cellsize=cellsize)
su = casatools.synthesisutils()
imsize = su.getOptimumSize(
    1000
)  # to reduce the data size but corver approximately 10" for each direction
for spw in spws:
    print(f"Generating dirty image for spw {spw}")
    tclean_spectral_line_wrapper(
        vis=[data_params[i]["vis_split"] for i in data_params.keys()],
        imagename=processingdir + prefix + f"_spw{spw}.dirty",
        spw=spw,
        imsize=imsize,
        cellsize=cellsize,
        scales=scales,
        weighting="briggsbwtaper",
        perchanweightdensity=True,
        robust=0.5,
        niter=0,
        parallel=True,
    )

### flag the data based on the visual inspection on the dirty images produced above
flagchannels = "0:0~270;880~930;1660~1919,1:400~580,2:400~580,3:100~300;550~900,4:300~680,5:750~1200"
spw = "0~5"
width = [
    120,
    960,
    960,
    960,
    960,
    960,
]  # Note that the continuum window is averaged over only specified chans to avoid potential beam smearing effect; see https://casaguides.nrao.edu/index.php/Image_Continuum

### flagging
for i in data_params.keys():
    # # save the state of before flag
    flagmanager(
        vis=data_params[i]["vis_split"], mode="save", versionname="before_cont_flags"
    )
    # flag the line channels
    flagdata(
        vis=data_params[i]["vis_split"],
        mode="manual",
        spw=flagchannels,
        flagbackup=False,
        field=field,
    )
    # save the state of flag
    flagmanager(
        vis=data_params[i]["vis_split"],
        mode="save",
        versionname="cont_flags",
    )

################ check the flag state by plotms here ################

for i in data_params.keys():
    # average over the channels and split out the continuum visibility
    avg_vis = processingdir + prefix + f"_{i}_initcont.ms"
    os.system("rm -r " + avg_vis)
    split(
        vis=data_params[i]["vis_split"],
        outputvis=avg_vis,
        spw=spw,
        field=field,
        width=width,
        datacolumn="data",
    )
    # add data params entry
    data_params[i]["vis_split_avg"] = avg_vis
    # restore the original data
    flagmanager(
        vis=data_params[i]["vis_split"], mode="restore", versionname="before_cont_flags"
    )

save_data_params(data_params, data_params_json_name)

### image continuum for each EB
# mandatory parameters
scales = [0, 5, 15, 25]
cellsize = "0.01arcsec"
# imsize = get_imsize(
# 	vis=[data_params[i]["vis_split_avg"] for i in data_params.keys()],
# 	field=field,
# 	cellsize=cellsize,
# )
imsize = su.getOptimumSize(1000)
# continuum disk mask; based on the visual inspection into the delivered data
disk_mask = "ellipse[[{:s}, {:s}], [{:.1f}arcsec, {:.1f}arcsec], {:.1f}deg]".format(
    "04h58m46.274s", "+29d50m36.49s", 2.5, 2.1, 342.0
)

for i in data_params.keys():
    print(i)
    imagename = processingdir + prefix + "_" + i + "_initial"
    tclean_continuum_wrapper(
        vis=data_params[i]["vis_split_avg"],
        imagename=imagename,
        scales=scales,
        mask=disk_mask,
        imsize=imsize,
        cellsize=cellsize,
        nsigma=3,
        robust=0.5,
    )

### fit each image with Gaussian to estimate the peak
for i in data_params.keys():
    data_params[i]["phasecenter"] = fit_gaussian(
        imagename=processingdir + prefix + f"_{i}_initial.image.tt0",
        region=disk_mask,
    )

### There are slight shifts between the peaks of EBs
### To correct for this, we fix the common phase center and shift the peak of each EB to that
### common direction; this direction is come from the C2D dataset. This also corrects for possible proper motions
common_dir = "J2000 04h58m46.27938s +029d50m36.399246s"

for i in data_params.keys():
    data_params[i]["vis_split_avg_shift"] = (
        processingdir + prefix + f"_{i}_initcont_shift.ms"
    )
    fix_phasecenter(
        vis=data_params[i]["vis_split_avg"],
        outputvis=data_params[i]["vis_split_avg_shift"],
        phase_center=data_params[i]["phasecenter"],
        field=field,
        common_dir=common_dir,
    )

save_data_params(data_params, prefix + "_data_params.json")

# check if the center is indeed aligned
for i in data_params.keys():
    print(i)
    imagename = processingdir + prefix + "_" + i + "_initial_shifted"
    tclean_continuum_wrapper(
        vis=data_params[i]["vis_split_avg_shift"],
        imagename=imagename,
        scales=scales,
        mask=disk_mask,
        imsize=imsize,
        cellsize=cellsize,
        nsigma=3,
        robust=0.5,
    )

for i in data_params.keys():
    data_params[i]["phasecenter_aligned"] = fit_gaussian(
        imagename=processingdir + prefix + f"_{i}_initial_shifted.image.tt0",
        region=disk_mask,
    )
# confirm that the peak is well aligned

# 04h58m46.279384s +29d50m36.39923s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279384s +29d50m36.39923s
# PA of Gaussian component: 165.24 deg
# Inclination of Gaussian component: 35.32 deg
# Pixel coordinates of peak: x = 499.919 y = 500.041
# 04h58m46.279381s +29d50m36.39927s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279381s +29d50m36.39927s
# PA of Gaussian component: 153.58 deg
# Inclination of Gaussian component: 39.15 deg
# Pixel coordinates of peak: x = 499.923 y = 500.044
# 04h58m46.279387s +29d50m36.39937s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279387s +29d50m36.39937s
# PA of Gaussian component: 154.76 deg
# Inclination of Gaussian component: 32.57 deg
# Pixel coordinates of peak: x = 499.915 y = 500.054
# 04h58m46.279374s +29d50m36.39931s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279374s +29d50m36.39931s
# PA of Gaussian component: 144.32 deg
# Inclination of Gaussian component: 37.92 deg
# Pixel coordinates of peak: x = 499.933 y = 500.048
# 04h58m46.279374s +29d50m36.39928s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279374s +29d50m36.39928s
# PA of Gaussian component: 146.94 deg
# Inclination of Gaussian component: 35.42 deg
# Pixel coordinates of peak: x = 499.932 y = 500.045
# 04h58m46.279380s +29d50m36.39929s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279380s +29d50m36.39929s
# PA of Gaussian component: 150.26 deg
# Inclination of Gaussian component: 35.27 deg
# Pixel coordinates of peak: x = 499.924 y = 500.046
# 04h58m46.279379s +29d50m36.39927s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279379s +29d50m36.39927s
# PA of Gaussian component: 148.76 deg
# Inclination of Gaussian component: 32.81 deg
# Pixel coordinates of peak: x = 499.926 y = 500.044
# 04h58m46.279383s +29d50m36.39928s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279383s +29d50m36.39928s
# PA of Gaussian component: 149.57 deg
# Inclination of Gaussian component: 35.26 deg
# Pixel coordinates of peak: x = 499.921 y = 500.046

# looks ok

save_data_params(data_params, prefix + "_data_params.json")

### flux rescaling
# By inspecting the visibility profile, we check the flux consistency between EBs
export_vislist = []
for i in data_params.keys():
    visfilename = export_MS(data_params[i]["vis_split_avg_shift"])
    export_vislist.append(visfilename)

PA, incl = 148, 37  # based on MAPS diskdictionary
plot_deprojected(
    export_vislist,
    incl=incl,
    PA=PA,
    fluxscale=[1.0] * len(export_vislist),
    show_err=True,
    outfile="amp-vs-uvdistance-pre-scaling.pdf",
)

# all datasets looks really nicely aligned within ~10% offsets

### self-calibration
visID = "vis_split_avg_shift"
# mask for emission and noise
disk_mask = "ellipse[[{:s}, {:s}], [{:.1f}arcsec, {:.1f}arcsec], {:.1f}deg]".format(
    "04h58m46.27936s", "+29d50m36.398586s", 1.5, 1.3, 344.0
)
noise_mask = "annulus[[{:s}, {:s}],['{:.1f}arcsec', '{:.1f}arcsec']]".format(
    "04h58m46.27936s", "+29d50m36.398586s", 3.5, 5.0
)


### data params preparation
for i in data_params.keys():
    data_params[i]["reference_antennas"] = rank_reference_antennas(
        data_params[i][visID]
    )
save_data_params(data_params, prefix + "_data_params.json")

# solution intervals
for i in data_params.keys():
    print(i)
    if "SB" in i:
        solint = get_solints(data_params[i][visID], field)

# SB1
# Determined solints:  ['inf', '90.72s', '42.34s', '18.14s', 'int']
# SB2
# Determined solints:  ['inf', '90.72s', '42.34s', '18.14s', 'int']

# first make a test image using only SB data
scales = [0, 5, 15, 25]
cellsize = "0.01arcsec"
imsize = su.getOptimumSize(1000)
imagename = processingdir + prefix + "_SB-only_initial"
tclean_continuum_wrapper(
    vis=[data_params[i][visID] for i in data_params.keys() if "SB" in i],
    imagename=imagename,
    scales=scales,
    mask=disk_mask,
    imsize=imsize,
    cellsize=cellsize,
    nsigma=3,
    robust=0.5,
    nterms=1,
)
initial_SNR, _ = estimate_SNR(
    imagename=imagename + ".image.tt0", disk_mask=disk_mask, noise_mask=noise_mask
)
# nterms = 1
# /works/yamato/MAPS_data/processing/Band6/MWC_480_SB-only_initial.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.60 deg)
# Flux inside disk mask: 347.55 mJy
# Peak intensity of source: 95.22 mJy/beam
# rms: 2.52e-01 mJy/beam
# Peak SNR: 377.58

# nterms = 2
# /works/yamato/MAPS_data/processing/Band6/MWC_480_SB-only_initial.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.60 deg)
# Flux inside disk mask: 367.76 mJy
# Peak intensity of source: 99.77 mJy/beam
# rms: 2.41e-01 mJy/beam
# Peak SNR: 414.34

# below use nterms = 2

print(initial_SNR)
import numpy as np

nsigma_init = np.max(
    [initial_SNR / 15.0, 5.0]
)  # restricts initial nsigma to be at least 5
nsigma_per_solint = 10 ** np.linspace(np.log10(nsigma_init), np.log10(3.0), len(solint))
print(nsigma_per_solint)
# [25.17207839 14.79004925  8.69000778  5.10588125  3.        ]


for i in data_params.keys():
    bw = get_effective_bandwidth(
        data_params[i]["vis_split"], data_params[i]["field"], flagchannels
    )
    nu0 = get_mean_frequency(data_params[i][visID], data_params[i]["field"])
    print(bw, nu0)

# iteration 1
mode = "SB-only"
iteration = 1
solint = "inf"
combine = "spw"
calmode = "p"
nsigma = 25.17
run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    solint=solint,
    combine=combine,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
    showplot=False,
)
###### pre selfcal iteration 1 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_initial.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.60 deg)
# Flux inside disk mask: 376.74 mJy
# Peak intensity of source: 100.17 mJy/beam
# rms: 2.89e-01 mJy/beam
# Peak SNR: 346.88
###### post selfcal iteration 1 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_iter1_p.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.61 deg)
# Flux inside disk mask: 375.71 mJy
# Peak intensity of source: 105.14 mJy/beam
# rms: 1.00e-01 mJy/beam
# Peak SNR: 1048.63

# iteration 2
mode = "SB-only"
iteration = 2
solint = "90.72s"
combine = "spw"
calmode = "p"
nsigma = 14.79
run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    solint=solint,
    combine=combine,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
    showplot=False,
)
###### pre selfcal iteration 2 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_iter1_p.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.61 deg)
# Flux inside disk mask: 374.05 mJy
# Peak intensity of source: 104.87 mJy/beam
# rms: 8.34e-02 mJy/beam
# Peak SNR: 1258.14
###### post selfcal iteration 2 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_iter2_p.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.61 deg)
# Flux inside disk mask: 373.99 mJy
# Peak intensity of source: 105.87 mJy/beam
# rms: 7.91e-02 mJy/beam
# Peak SNR: 1337.75

# iteration 3
mode = "SB-only"
iteration = 3
solint = "42.34s"
combine = "spw"
calmode = "p"
nsigma = 8.690
run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    solint=solint,
    combine=combine,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
    showplot=False,
)
###### pre selfcal iteration 3 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_iter2_p.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.61 deg)
# Flux inside disk mask: 373.05 mJy
# Peak intensity of source: 105.79 mJy/beam
# rms: 7.20e-02 mJy/beam
# Peak SNR: 1468.52
###### post selfcal iteration 3 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_iter3_p.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.61 deg)
# Flux inside disk mask: 373.38 mJy
# Peak intensity of source: 106.85 mJy/beam
# rms: 7.08e-02 mJy/beam
# Peak SNR: 1508.85

# iteration 4
mode = "SB-only"
iteration = 4
solint = "18.14s"
combine = "spw"
calmode = "p"
nsigma = 5.106
run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    solint=solint,
    combine=combine,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
    showplot=False,
)
###### pre selfcal iteration 4 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_iter3_p.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.61 deg)
# Flux inside disk mask: 372.95 mJy
# Peak intensity of source: 106.82 mJy/beam
# rms: 6.66e-02 mJy/beam
# Peak SNR: 1603.47
###### post selfcal iteration 4 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_iter4_p.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.61 deg)
# Flux inside disk mask: 373.28 mJy
# Peak intensity of source: 107.72 mJy/beam
# rms: 6.62e-02 mJy/beam
# Peak SNR: 1628.45

# iteration 5
mode = "SB-only"
iteration = 5
solint = "int"
combine = "spw"
calmode = "p"
nsigma = 3.0
run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    solint=solint,
    combine=combine,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
    showplot=False,
)
###### pre selfcal iteration 5 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_iter4_p.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.61 deg)
# Flux inside disk mask: 373.11 mJy
# Peak intensity of source: 107.63 mJy/beam
# rms: 6.37e-02 mJy/beam
# Peak SNR: 1689.50
###### post selfcal iteration 5 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_iter5_p.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.61 deg)
# Flux inside disk mask: 373.87 mJy
# Peak intensity of source: 108.16 mJy/beam
# rms: 6.28e-02 mJy/beam
# Peak SNR: 1721.94

### phase + amplitude selfcal
# iteration 6
mode = "SB-only"
iteration = 6
solint = "inf"
combine = "spw"
calmode = "ap"
nsigma = 3.0
run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    solint=solint,
    combine=combine,
    solnorm_ap=False,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
    showplot=False,
)
###### pre selfcal iteration 6 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_iter5_p.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.61 deg)
# Flux inside disk mask: 373.87 mJy
# Peak intensity of source: 108.16 mJy/beam
# rms: 6.28e-02 mJy/beam
# Peak SNR: 1721.94
###### post selfcal iteration 6 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_iter6_ap.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.71 deg)
# Flux inside disk mask: 373.22 mJy
# Peak intensity of source: 108.09 mJy/beam
# rms: 6.35e-02 mJy/beam
# Peak SNR: 1701.76

# this step (iteration = 6) decreases the S/N and thus finalize at iteration = 5
iteration = 5
mode = "SB-only"
calmode = "p"
nsigma = 3.0
finalize_selfcal(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
)
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_SB-only_final.image.tt0
# Beam 0.339 arcsec x 0.235 arcsec (-6.61 deg)
# Flux inside disk mask: 373.87 mJy
# Peak intensity of source: 108.16 mJy/beam
# rms: 6.28e-02 mJy/beam
# Peak SNR: 1721.94

save_data_params(data_params, data_params_json_name)


### SB+LB self-calibration ###
for i in data_params.keys():
    print(i)
    if "LB" in i:
        solint = get_solints(data_params[i][visID], field)

# LB1
# Determined solints:  ['inf', '26.21s', '12.10s', '6.05s', 'int']
# LB2
# Determined solints:  ['inf', '26.21s', '12.10s', '6.05s', 'int']
# LB3
# Determined solints:  ['inf', '26.21s', '12.10s', '6.05s', 'int']
# LB4
# Determined solints:  ['inf', '26.21s', '12.10s', '6.05s', 'int']
# LB5
# Determined solints:  ['inf', '26.21s', '12.10s', '6.05s', 'int']
# LB6
# Determined solints:  ['inf', '26.21s', '12.10s', '6.05s', 'int']

scales = [0, 5, 15, 25]
cellsize = "0.01arcsec"
imsize = su.getOptimumSize(1000)
imagename = processingdir + prefix + "_combined_initial"
tclean_continuum_wrapper(
    vis=[data_params[i][visID] for i in data_params.keys()],
    imagename=imagename,
    scales=scales,
    mask=disk_mask,
    imsize=imsize,
    cellsize=cellsize,
    nsigma=3,
    robust=0.5,
    nterms=1,
)
initial_SNR, _ = estimate_SNR(
    imagename=imagename + ".image.tt0", disk_mask=disk_mask, noise_mask=noise_mask
)

# /works/yamato/MAPS_data/processing/Band6/MWC_480_combined_initial.image.tt0
# Beam 0.143 arcsec x 0.104 arcsec (-2.85 deg)
# Flux inside disk mask: 350.48 mJy
# Peak intensity of source: 29.47 mJy/beam
# rms: 5.10e-02 mJy/beam
# Peak SNR: 577.35

print(initial_SNR)
import numpy as np

nsigma_init = np.max(
    [initial_SNR / 15.0, 5.0]
)  # restricts initial nsigma to be at least 5
nsigma_per_solint = 10 ** np.linspace(np.log10(nsigma_init), np.log10(3.0), len(solint))
print(nsigma_per_solint)
# [38.49012511 20.33725367 10.74571428  5.67777622  3.        ]

# iteration 1
mode = "combined"
iteration = 1
solint = "inf"
combine = "spw"
calmode = "p"
nsigma = 38.49
run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    solint=solint,
    combine=combine,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
    showplot=False,
)
###### pre selfcal iteration 1 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_initial.image.tt0
# Beam 0.142 arcsec x 0.104 arcsec (-2.86 deg)
# Flux inside disk mask: 402.13 mJy
# Peak intensity of source: 32.12 mJy/beam
# rms: 4.91e-02 mJy/beam
# Peak SNR: 654.41
###### post selfcal iteration 1 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_iter1_p.image.tt0
# Beam 0.143 arcsec x 0.104 arcsec (-2.85 deg)
# Flux inside disk mask: 390.15 mJy
# Peak intensity of source: 34.27 mJy/beam
# rms: 3.40e-02 mJy/beam
# Peak SNR: 1008.63

# iteration 2
mode = "combined"
iteration = 2
solint = "26.21s"
combine = "spw"
calmode = "p"
nsigma = 20.34
run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    solint=solint,
    combine=combine,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
    showplot=False,
)
###### pre selfcal iteration 2 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_iter1_p.image.tt0
# Beam 0.143 arcsec x 0.104 arcsec (-2.85 deg)
# Flux inside disk mask: 382.71 mJy
# Peak intensity of source: 34.07 mJy/beam
# rms: 2.95e-02 mJy/beam
# Peak SNR: 1154.26
###### post selfcal iteration 2 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_iter2_p.image.tt0
# Beam 0.143 arcsec x 0.104 arcsec (-2.85 deg)
# Flux inside disk mask: 382.40 mJy
# Peak intensity of source: 34.38 mJy/beam
# rms: 2.92e-02 mJy/beam
# Peak SNR: 1176.11

# iteration 3
mode = "combined"
iteration = 3
solint = "12.10s"
combine = "spw"
calmode = "p"
nsigma = 10.75
run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    solint=solint,
    combine=combine,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
    showplot=False,
)
###### pre selfcal iteration 3 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_iter2_p.image.tt0
# Beam 0.143 arcsec x 0.104 arcsec (-2.85 deg)
# Flux inside disk mask: 378.47 mJy
# Peak intensity of source: 34.30 mJy/beam
# rms: 2.73e-02 mJy/beam
# Peak SNR: 1254.70
###### post selfcal iteration 3 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_iter3_p.image.tt0
# Beam 0.143 arcsec x 0.104 arcsec (-2.85 deg)
# Flux inside disk mask: 378.51 mJy
# Peak intensity of source: 34.56 mJy/beam
# rms: 2.71e-02 mJy/beam
# Peak SNR: 1273.05

"""
# iteration 4
mode = "combined"
iteration = 4
solint = "6.05s"
combine = "spw"
calmode = "p"
nsigma = 5.68
run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    solint=solint,
    combine=combine,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
    showplot=False,
)
###### pre selfcal iteration 4 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_iter3_p.image.tt0
# Beam 0.143 arcsec x 0.104 arcsec (-2.85 deg)
# Flux inside disk mask: 376.65 mJy
# Peak intensity of source: 34.54 mJy/beam
# rms: 2.66e-02 mJy/beam
# Peak SNR: 1300.59
###### post selfcal iteration 4 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_iter4_p.image.tt0
# Beam 0.143 arcsec x 0.104 arcsec (-2.85 deg)
# Flux inside disk mask: 376.78 mJy
# Peak intensity of source: 34.67 mJy/beam
# rms: 2.67e-02 mJy/beam
# Peak SNR: 1299.13

# this step decreases the S/N; instead try shorter solint and smaller nsigma
# iteration 4
mode = "combined"
iteration = 4
solint = "int"
combine = "spw"
calmode = "p"
nsigma = 3.0
run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    solint=solint,
    combine=combine,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
    showplot=False,
)
###### pre selfcal iteration 4 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_iter3_p.image.tt0
# Beam 0.143 arcsec x 0.104 arcsec (-2.85 deg)
# Flux inside disk mask: 375.68 mJy
# Peak intensity of source: 34.54 mJy/beam
# rms: 2.63e-02 mJy/beam
# Peak SNR: 1313.35
###### post selfcal iteration 4 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_iter4_p.image.tt0
# Beam 0.143 arcsec x 0.104 arcsec (-2.85 deg)
# Flux inside disk mask: 375.96 mJy
# Peak intensity of source: 34.67 mJy/beam
# rms: 2.69e-02 mJy/beam
# Peak SNR: 1287.09

# this also decreases the S/N; going to amplitude selfcal
# iteration 4
mode = "combined"
iteration = 4
solint = "inf"
combine = "spw"
calmode = "ap"
nsigma = 3.0
run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    solint=solint,
    combine=combine,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
    showplot=False,
)
###### pre selfcal iteration 4 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_iter3_p.image.tt0
# Beam 0.143 arcsec x 0.104 arcsec (-2.85 deg)
# Flux inside disk mask: 375.68 mJy
# Peak intensity of source: 34.54 mJy/beam
# rms: 2.63e-02 mJy/beam
# Peak SNR: 1313.35
###### post selfcal iteration 4 ######
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_iter4_ap.image.tt0
# Beam 0.141 arcsec x 0.103 arcsec (-3.60 deg)
# Flux inside disk mask: 375.79 mJy
# Peak intensity of source: 33.99 mJy/beam
# rms: 2.59e-02 mJy/beam
# Peak SNR: 1311.70
"""

# this step also decrease the S/N; finally stop at iteration = 3
iteration = 3
mode = "combined"
calmode = "p"
nsigma = 3.0
finalize_selfcal(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=mode,
    calmode=calmode,
    disk_mask=disk_mask,
    noise_mask=noise_mask,
    tclean_kwargs=dict(
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=2
    ),
)
# /works/yamato/MAPS_data/processing/Band6/MWC_480_selfcal_combined_final.image.tt0
# Beam 0.143 arcsec x 0.104 arcsec (-2.85 deg)
# Flux inside disk mask: 375.68 mJy
# Peak intensity of source: 34.54 mJy/beam
# rms: 2.63e-02 mJy/beam
# Peak SNR: 1313.35

save_data_params(data_params, data_params_json_name)


###### split off the final continuum data ######
for i in data_params.keys():
    print(i)
    outputvis = MSdir + prefix + "_" + i + "_continuum.ms"
    os.system("rm -r " + outputvis)
    split(
        vis=data_params[i][visID + "_selfcal"], outputvis=outputvis, datacolumn="data"
    )
    data_params[i]["vis_continuum"] = outputvis

save_data_params(data_params, data_params_json_name)


########### apply the solution to the line measurement sets ###########
# shift the center to the common_dir
for i in data_params.keys():
    print(i)
    data_params[i]["vis_split_shift"] = processingdir + prefix + "_" + i + "_shift.ms"
    os.system("rm -r " + data_params[i]["vis_split_shift"])
    fix_phasecenter(
        vis=data_params[i]["vis_split"],
        outputvis=data_params[i]["vis_split_shift"],
        phase_center=data_params[i]["phasecenter"],
        common_dir=common_dir,
        field=data_params[i]["field"],
    )

save_data_params(data_params, data_params_json_name)

# apply the solutions from continuum selfcal to the line MSs
for i in data_params.keys():
    print(i)
    ntables = len(data_params[i]["selfcal_table"])
    applycal(
        vis=data_params[i]["vis_split_shift"],
        spw="",
        gaintable=data_params[i]["selfcal_table"],
        spwmap=data_params[i]["selfcal_spwmap"],
        interp=["linearPD"] * ntables,
        calwt=True,
        applymode="calonly",
    )
    data_params[i]["vis_split_shift_selfcal"] = data_params[i][
        "vis_split_shift"
    ].replace(".ms", "_selfcal.ms")
    os.system("rm -r " + data_params[i]["vis_split_shift_selfcal"])
    split(
        vis=data_params[i]["vis_split_shift"],
        outputvis=data_params[i]["vis_split_shift_selfcal"],
        datacolumn="corrected",
    )

save_data_params(data_params, data_params_json_name)


########## continuum subtraction ##########
for i in data_params.keys():
    print(i, flagchannels)
    os.system("rm -r " + data_params[i]["vis_split_shift_selfcal"] + ".contsub")
    uvcontsub_old(
        vis=data_params[i]["vis_split_shift_selfcal"],
        fitspw=flagchannels,
        excludechans=True,
        combine="spw",
        fitorder=1,
    )
    data_params[i]["vis_split_shift_selfcal_contsub"] = (
        data_params[i]["vis_split_shift_selfcal"] + ".contsub"
    )

save_data_params(data_params, data_params_json_name)

# locate the final mesurement sets to the product directory
for i in data_params.keys():
    data_params[i]["vis_spectral_line_wcont"] = (
        MSdir + prefix + "_" + i + "_spectral_line_wcont.ms"
    )
    data_params[i]["vis_spectral_line"] = MSdir + prefix + "_" + i + "_spectral_line.ms"
    os.system(
        "cp -r "
        + data_params[i]["vis_split_shift_selfcal"]
        + " "
        + data_params[i]["vis_spectral_line_wcont"]
    )
    os.system(
        "cp -r "
        + data_params[i]["vis_split_shift_selfcal_contsub"]
        + " "
        + data_params[i]["vis_spectral_line"]
    )

save_data_params(data_params, data_params_json_name)


# # self-calibration parameters
# iteration = 1
# solint = "inf"
# combine = "spw,scan"
# spwmap = [0] * len(contspws)  # subparameter of combine
# mode = "p"

# print("Iteration: " + str(iteration))
# data_params = load_data_params(prefix + "_data_params.json")
# for i in data_params.keys():
# 	print(i)
# 	caltable = processingdir + prefix + "_" + i + "_selfcal_p" + str(iteration) + ".g"
# 	data_params[i]["solint"].append(solint)
# 	data_params[i]["combine"].append(combine)
# 	data_params[i]["spwmap"].append(spwmap)
# 	data_params[i]["selfcal_tables"].append(caltable)
# 	print("Solving gain...")
# 	gaincal(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		caltable=caltable,
# 		gaintype="T",
# 		refant=data_params[i]["reference_antennas"],
# 		calmode=mode,
# 		solint=solint,
# 		spw="0~9",
# 		minsnr=2.0,
# 		minblperant=4,
# 		combine=combine,
# 	)
# 	print("Plotting gain solutions...")
# 	plot_solution_SNR_hist(
# 		caltable,
# 		solint=solint,
# 		combine=combine,
# 		outfile=caltable.replace(".g", "_SNR_hist.pdf"),
# 	)
# 	plot_calibration_table(
# 		caltable,
# 		plotval="phase",
# 		outfile=caltable.replace(".g", "_gain.pdf"),
# 		show=False,
# 	)
# 	print("Applying solutions...")
# 	applycal(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		gaintable=[caltable],
# 		interp="linearPD",
# 		calwt=True,
# 		spw="0~9",
# 		spwmap=spwmap,
# 	)
# 	print("Split off corrected data...")
# 	outputvis = processingdir + prefix + "_" + i + "_selfcal_p" + str(iteration) + ".ms"
# 	os.system("rm -r " + outputvis + "*")
# 	split(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		outputvis=outputvis,
# 		datacolumn="corrected",
# 	)
# 	data_params[i][selected_vis_key + "_selfcal"] = outputvis
# 	print(
# 		"Iteration {:s} of self-calibration (mode = '{:s}') completed. Check the solved gain by inspecting the gain plot.".format(
# 			str(iteration), mode
# 		)
# 	)

# print("Image restoration to assess the improvement of SNR...")
# vis = [data_params[i][selected_vis_key + "_selfcal"] for i in data_params.keys()]
# imagename = processingdir + prefix + "_selfcal_p" + str(iteration)
# tclean_continuum_wrapper(
# 	vis=vis,
# 	imagename=imagename,
# 	scales=scales,
# 	mask=disk_mask,
# 	imsize=imsize,
# 	cellsize=cellsize,
# 	# threshold="0.7mJy",
# 	nsigma=3.0,
# 	robust=0.5,
# 	savemodel="modelcolumn",
# )

# estimate_SNR(imagename + ".image.tt0", disk_mask, noise_mask)
# # /works/yamato/C2D_disk/processing/MWC_480_selfcal_p1.image.tt0
# # Beam 0.307 arcsec x 0.204 arcsec (7.16 deg)
# # Flux inside disk mask: 522.91 mJy
# # Peak intensity of source: 107.83 mJy/beam
# # rms: 1.84e-01 mJy/beam
# # Peak SNR: 586.15

# # /works/yamato/C2D_disk/processing/MWC_480_selfcal_p1.image.tt0
# # Beam 0.307 arcsec x 0.204 arcsec (7.04 deg)
# # Flux inside disk mask: 531.84 mJy
# # Peak intensity of source: 121.80 mJy/beam
# # rms: 7.07e-02 mJy/beam
# # Peak SNR: 1722.32

# # if it's OK to apply this solution, save the data_params
# save_data_params(data_params, prefix + "_data_params.json")


# # self-calibration parameters
# iteration = 2
# solint = "inf"
# combine = "spw"
# spwmap = [0] * len(contspws)  # subparameter of combine
# mode = "p"

# print("Iteration: " + str(iteration))
# data_params = load_data_params(prefix + "_data_params.json")
# for i in data_params.keys():
# 	print(i)
# 	caltable = processingdir + prefix + "_" + i + "_selfcal_p" + str(iteration) + ".g"
# 	data_params[i]["solint"].append(solint)
# 	data_params[i]["combine"].append(combine)
# 	data_params[i]["spwmap"].append(spwmap)
# 	data_params[i]["selfcal_tables"].append(caltable)
# 	print("Solving gain...")
# 	gaincal(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		caltable=caltable,
# 		gaintype="T",
# 		refant=data_params[i]["reference_antennas"],
# 		calmode=mode,
# 		solint=solint,
# 		spw="0~9",
# 		minsnr=2.0,
# 		minblperant=4,
# 		combine=combine,
# 	)
# 	print("Plotting gain solutions...")
# 	plot_solution_SNR_hist(
# 		caltable,
# 		solint=solint,
# 		combine=combine,
# 		outfile=caltable.replace(".g", "_SNR_hist.pdf"),
# 	)
# 	plot_calibration_table(
# 		caltable,
# 		plotval="phase",
# 		outfile=caltable.replace(".g", "_gain.pdf"),
# 		show=False,
# 	)
# 	print("Applying solutions...")
# 	applycal(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		gaintable=[caltable],
# 		interp="linearPD",
# 		calwt=True,
# 		spw="0~9",
# 		spwmap=spwmap,
# 	)
# 	print("Split off corrected data...")
# 	outputvis = processingdir + prefix + "_" + i + "_selfcal_p" + str(iteration) + ".ms"
# 	os.system("rm -r " + outputvis + "*")
# 	split(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		outputvis=outputvis,
# 		datacolumn="corrected",
# 	)
# 	data_params[i][selected_vis_key + "_selfcal"] = outputvis
# 	print(
# 		"Iteration {:s} of self-calibration (mode = '{:s}') completed. Check the solved gain by inspecting the gain plot.".format(
# 			str(iteration), mode
# 		)
# 	)

# print("Image restoration to assess the improvement of SNR...")
# vis = [data_params[i][selected_vis_key + "_selfcal"] for i in data_params.keys()]
# imagename = processingdir + prefix + "_selfcal_p" + str(iteration)
# tclean_continuum_wrapper(
# 	vis=vis,
# 	imagename=imagename,
# 	scales=scales,
# 	mask=disk_mask,
# 	imsize=imsize,
# 	cellsize=cellsize,
# 	# threshold="0.7mJy",
# 	nsigma=3.0,
# 	robust=0.5,
# 	savemodel="modelcolumn",
# )

# estimate_SNR(imagename + ".image.tt0", disk_mask, noise_mask)
# # /works/yamato/C2D_disk/processing/MWC_480_selfcal_p2.image.tt0
# # Beam 0.307 arcsec x 0.204 arcsec (7.05 deg)
# # Flux inside disk mask: 533.10 mJy
# # Peak intensity of source: 124.51 mJy/beam
# # rms: 6.13e-02 mJy/beam
# # Peak SNR: 2029.66

# # /works/yamato/C2D_disk/processing/MWC_480_selfcal_p2.image.tt0
# # Beam 0.307 arcsec x 0.204 arcsec (7.17 deg)
# # Flux inside disk mask: 532.00 mJy
# # Peak intensity of source: 121.83 mJy/beam
# # rms: 7.01e-02 mJy/beam
# # Peak SNR: 1737.32

# # if it's OK to apply this solution, save the data_params
# save_data_params(data_params, prefix + "_data_params.json")


# # self-calibration parameters
# iteration = 3
# solint = "90.72s"
# combine = "spw"
# spwmap = [0] * len(contspws)  # subparameter of combine
# # spwmap = [""]
# mode = "p"

# print("Iteration: " + str(iteration))
# data_params = load_data_params(prefix + "_data_params.json")
# for i in data_params.keys():
# 	print(i)
# 	caltable = processingdir + prefix + "_" + i + "_selfcal_p" + str(iteration) + ".g"
# 	data_params[i]["solint"].append(solint)
# 	data_params[i]["combine"].append(combine)
# 	data_params[i]["spwmap"].append(spwmap)
# 	data_params[i]["selfcal_tables"].append(caltable)
# 	print("Solving gain...")
# 	gaincal(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		caltable=caltable,
# 		gaintype="T",
# 		refant=data_params[i]["reference_antennas"],
# 		calmode=mode,
# 		solint=solint,
# 		spw="0~9",
# 		minsnr=2.0,
# 		minblperant=4,
# 		combine=combine,
# 	)
# 	print("Plotting gain solutions...")
# 	plot_solution_SNR_hist(
# 		caltable,
# 		solint=solint,
# 		combine=combine,
# 		outfile=caltable.replace(".g", "_SNR_hist.pdf"),
# 	)
# 	plot_calibration_table(
# 		caltable,
# 		plotval="phase",
# 		outfile=caltable.replace(".g", "_gain.pdf"),
# 		show=False,
# 	)
# 	print("Applying solutions...")
# 	applycal(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		gaintable=[caltable],
# 		interp="linearPD",
# 		calwt=True,
# 		spw="0~9",
# 		spwmap=spwmap,
# 	)
# 	print("Split off corrected data...")
# 	outputvis = processingdir + prefix + "_" + i + "_selfcal_p" + str(iteration) + ".ms"
# 	os.system("rm -r " + outputvis + "*")
# 	split(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		outputvis=outputvis,
# 		datacolumn="corrected",
# 	)
# 	data_params[i][selected_vis_key + "_selfcal"] = outputvis
# 	print(
# 		"Iteration {:s} of self-calibration (mode = '{:s}') completed. Check the solved gain by inspecting the gain plot.".format(
# 			str(iteration), mode
# 		)
# 	)

# print("Image restoration to assess the improvement of SNR...")
# vis = [data_params[i][selected_vis_key + "_selfcal"] for i in data_params.keys()]
# imagename = processingdir + prefix + "_selfcal_p" + str(iteration)
# tclean_continuum_wrapper(
# 	vis=vis,
# 	imagename=imagename,
# 	scales=scales,
# 	mask=disk_mask,
# 	imsize=imsize,
# 	cellsize=cellsize,
# 	# threshold="0.7mJy",
# 	nsigma=3.0,
# 	robust=0.5,
# 	savemodel="modelcolumn",
# )

# estimate_SNR(imagename + ".image.tt0", disk_mask, noise_mask)
# # /works/yamato/C2D_disk/processing/MWC_480_selfcal_p3.image.tt0
# # Beam 0.313 arcsec x 0.212 arcsec (9.10 deg)
# # Flux inside disk mask: 535.08 mJy
# # Peak intensity of source: 131.99 mJy/beam
# # rms: 5.52e-02 mJy/beam
# # Peak SNR: 2389.28

# # /works/yamato/C2D_disk/processing/MWC_480_selfcal_p3.image.tt0
# # Beam 0.307 arcsec x 0.204 arcsec (7.17 deg)
# # Flux inside disk mask: 533.22 mJy
# # Peak intensity of source: 124.51 mJy/beam
# # rms: 6.20e-02 mJy/beam
# # Peak SNR: 2009.52

# # if it's OK to apply this solution, save the data_params
# save_data_params(data_params, prefix + "_data_params.json")


# # self-calibration parameters
# iteration = 4
# solint = "42.3s"
# combine = "spw"
# spwmap = [0] * len(contspws)  # subparameter of combine
# # spwmap = [""]
# mode = "p"

# print("Iteration: " + str(iteration))
# data_params = load_data_params(prefix + "_data_params.json")
# for i in data_params.keys():
# 	print(i)
# 	caltable = processingdir + prefix + "_" + i + "_selfcal_p" + str(iteration) + ".g"
# 	data_params[i]["solint"].append(solint)
# 	data_params[i]["combine"].append(combine)
# 	data_params[i]["spwmap"].append(spwmap)
# 	data_params[i]["selfcal_tables"].append(caltable)
# 	print("Solving gain...")
# 	gaincal(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		caltable=caltable,
# 		gaintype="T",
# 		refant=data_params[i]["reference_antennas"],
# 		calmode=mode,
# 		solint=solint,
# 		spw="0~9",
# 		minsnr=2.0,
# 		minblperant=4,
# 		combine=combine,
# 	)
# 	print("Plotting gain solutions...")
# 	plot_solution_SNR_hist(
# 		caltable,
# 		solint=solint,
# 		combine=combine,
# 		outfile=caltable.replace(".g", "_SNR_hist.pdf"),
# 	)
# 	plot_calibration_table(
# 		caltable,
# 		plotval="phase",
# 		outfile=caltable.replace(".g", "_gain.pdf"),
# 		show=False,
# 	)
# 	print("Applying solutions...")
# 	applycal(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		gaintable=[caltable],
# 		interp="linearPD",
# 		calwt=True,
# 		spw="0~9",
# 		spwmap=spwmap,
# 	)
# 	print("Split off corrected data...")
# 	outputvis = processingdir + prefix + "_" + i + "_selfcal_p" + str(iteration) + ".ms"
# 	os.system("rm -r " + outputvis + "*")
# 	split(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		outputvis=outputvis,
# 		datacolumn="corrected",
# 	)
# 	data_params[i][selected_vis_key + "_selfcal"] = outputvis
# 	print(
# 		"Iteration {:s} of self-calibration (mode = '{:s}') completed. Check the solved gain by inspecting the gain plot.".format(
# 			str(iteration), mode
# 		)
# 	)

# print("Image restoration to assess the improvement of SNR...")
# vis = [data_params[i][selected_vis_key + "_selfcal"] for i in data_params.keys()]
# imagename = processingdir + prefix + "_selfcal_p" + str(iteration)
# tclean_continuum_wrapper(
# 	vis=vis,
# 	imagename=imagename,
# 	scales=scales,
# 	mask=disk_mask,
# 	imsize=imsize,
# 	cellsize=cellsize,
# 	# threshold="0.7mJy",
# 	nsigma=3.0,
# 	robust=0.5,
# 	savemodel="modelcolumn",
# )

# estimate_SNR(imagename + ".image.tt0", disk_mask, noise_mask)
# # /works/yamato/C2D_disk/processing/MWC_480_selfcal_p4.image.tt0
# # Beam 0.313 arcsec x 0.213 arcsec (9.21 deg)
# # Flux inside disk mask: 537.20 mJy
# # Peak intensity of source: 134.49 mJy/beam
# # rms: 5.16e-02 mJy/beam
# # Peak SNR: 2604.54

# # /works/yamato/C2D_disk/processing/MWC_480_selfcal_p4.image.tt0
# # Beam 0.307 arcsec x 0.204 arcsec (7.18 deg)
# # Flux inside disk mask: 535.24 mJy
# # Peak intensity of source: 127.42 mJy/beam
# # rms: 5.47e-02 mJy/beam
# # Peak SNR: 2328.46

# # if it's OK to apply this solution, save the data_params
# save_data_params(data_params, prefix + "_data_params.json")


# # self-calibration parameters
# iteration = 5
# solint = "18.14s"
# combine = "spw"
# spwmap = [0] * len(contspws)  # subparameter of combine
# # spwmap = [""]
# mode = "p"

# print("Iteration: " + str(iteration))
# data_params = load_data_params(prefix + "_data_params.json")
# for i in data_params.keys():
# 	print(i)
# 	caltable = processingdir + prefix + "_" + i + "_selfcal_p" + str(iteration) + ".g"
# 	data_params[i]["solint"].append(solint)
# 	data_params[i]["combine"].append(combine)
# 	data_params[i]["spwmap"].append(spwmap)
# 	data_params[i]["selfcal_tables"].append(caltable)
# 	print("Solving gain...")
# 	gaincal(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		caltable=caltable,
# 		gaintype="T",
# 		refant=data_params[i]["reference_antennas"],
# 		calmode=mode,
# 		solint=solint,
# 		spw="0~9",
# 		minsnr=2.0,
# 		minblperant=4,
# 		combine=combine,
# 	)
# 	print("Plotting gain solutions...")
# 	plot_solution_SNR_hist(
# 		caltable,
# 		solint=solint,
# 		combine=combine,
# 		outfile=caltable.replace(".g", "_SNR_hist.pdf"),
# 	)
# 	plot_calibration_table(
# 		caltable,
# 		plotval="phase",
# 		outfile=caltable.replace(".g", "_gain.pdf"),
# 		show=False,
# 	)
# 	print("Applying solutions...")
# 	applycal(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		gaintable=[caltable],
# 		interp="linearPD",
# 		calwt=True,
# 		spw="0~9",
# 		spwmap=spwmap,
# 	)
# 	print("Split off corrected data...")
# 	outputvis = processingdir + prefix + "_" + i + "_selfcal_p" + str(iteration) + ".ms"
# 	os.system("rm -r " + outputvis + "*")
# 	split(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		outputvis=outputvis,
# 		datacolumn="corrected",
# 	)
# 	data_params[i][selected_vis_key + "_selfcal"] = outputvis
# 	print(
# 		"Iteration {:s} of self-calibration (mode = '{:s}') completed. Check the solved gain by inspecting the gain plot.".format(
# 			str(iteration), mode
# 		)
# 	)

# print("Image restoration to assess the improvement of SNR...")
# vis = [data_params[i][selected_vis_key + "_selfcal"] for i in data_params.keys()]
# imagename = processingdir + prefix + "_selfcal_p" + str(iteration)
# tclean_continuum_wrapper(
# 	vis=vis,
# 	imagename=imagename,
# 	scales=scales,
# 	mask=disk_mask,
# 	imsize=imsize,
# 	cellsize=cellsize,
# 	# threshold="0.7mJy",
# 	nsigma=3.0,
# 	robust=0.5,
# 	savemodel="modelcolumn",
# )

# estimate_SNR(imagename + ".image.tt0", disk_mask, noise_mask)
# # /works/yamato/C2D_disk/processing/MWC_480_selfcal_p5.image.tt0
# # Beam 0.313 arcsec x 0.213 arcsec (9.21 deg)
# # Flux inside disk mask: 539.26 mJy
# # Peak intensity of source: 135.29 mJy/beam
# # rms: 5.03e-02 mJy/beam
# # Peak SNR: 2691.80

# # if it's OK to apply this solution, save the data_params
# save_data_params(data_params, prefix + "_data_params.json")


# ### amplitude+phase self-calibration
# # self-calibration parameters
# iteration = 1
# solint = "inf"
# combine = "spw"
# spwmap = [0] * len(contspws)  # subparameter of combine
# # spwmap = [""]
# mode = "ap"

# print("Iteration: " + str(iteration))
# data_params = load_data_params(prefix + "_data_params.json")
# for i in data_params.keys():
# 	print(i)
# 	caltable = (
# 		processingdir + prefix + "_" + i + f"_selfcal_{mode}" + str(iteration) + ".g"
# 	)
# 	data_params[i]["solint"].append(solint)
# 	data_params[i]["combine"].append(combine)
# 	data_params[i]["spwmap"].append(spwmap)
# 	data_params[i]["selfcal_tables"].append(caltable)
# 	print("Solving gain...")
# 	gaincal(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		caltable=caltable,
# 		gaintype="T",
# 		refant=data_params[i]["reference_antennas"],
# 		calmode=mode,
# 		solint=solint,
# 		spw="0~9",
# 		minsnr=4.0,
# 		minblperant=4,
# 		combine=combine,
# 	)
# 	print("Plotting gain solutions...")
# 	plot_solution_SNR_hist(
# 		caltable,
# 		solint=solint,
# 		combine=combine,
# 		outfile=caltable.replace(".g", "_SNR_hist.pdf"),
# 	)
# 	plot_calibration_table(
# 		caltable,
# 		plotval="amplitude",
# 		outfile=caltable.replace(".g", "_gain.pdf"),
# 		show=False,
# 	)
# 	print("Applying solutions...")
# 	applycal(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		gaintable=[caltable],
# 		interp="linearPD",
# 		calwt=True,
# 		spw="0~9",
# 		spwmap=spwmap,
# 	)
# 	print("Split off corrected data...")
# 	outputvis = (
# 		processingdir + prefix + "_" + i + f"_selfcal_{mode}" + str(iteration) + ".ms"
# 	)
# 	os.system("rm -r " + outputvis + "*")
# 	split(
# 		vis=data_params[i][selected_vis_key + "_selfcal"],
# 		outputvis=outputvis,
# 		datacolumn="corrected",
# 	)
# 	data_params[i][selected_vis_key + "_selfcal"] = outputvis
# 	print(
# 		"Iteration {:s} of self-calibration (mode = '{:s}') completed. Check the solved gain by inspecting the gain plot.".format(
# 			str(iteration), mode
# 		)
# 	)

# print("Image restoration to assess the improvement of SNR...")
# vis = [data_params[i][selected_vis_key + "_selfcal"] for i in data_params.keys()]
# imagename = processingdir + prefix + f"_selfcal_{mode}" + str(iteration)
# tclean_continuum_wrapper(
# 	vis=vis,
# 	imagename=imagename,
# 	scales=scales,
# 	mask=disk_mask,
# 	imsize=imsize,
# 	cellsize=cellsize,
# 	# threshold="0.7mJy",
# 	nsigma=3.0,
# 	robust=0.5,
# 	savemodel="modelcolumn",
# )

# estimate_SNR(imagename + ".image.tt0", disk_mask, noise_mask)
# # /works/yamato/C2D_disk/processing/MWC_480_selfcal_ap1.image.tt0
# # Beam 0.315 arcsec x 0.213 arcsec (10.23 deg)
# # Flux inside disk mask: 539.02 mJy
# # Peak intensity of source: 135.66 mJy/beam
# # rms: 2.83e-02 mJy/beam
# # Peak SNR: 4790.02

# # if it's OK to apply this solution, save the data_params
# save_data_params(data_params, prefix + "_data_params.json")


# ### nearly the thermal noise, maybe it's OK
