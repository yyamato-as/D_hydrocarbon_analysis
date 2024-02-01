import os
import sys

sys.path.append("/home/yamato/Project/hydrocarbon_deuteration/")

execfile("/home/yamato/Project/hydrocarbon_deuteration/reduction_utils_mpi_v3.py")
# import reduction_utils_mpi as ru

# directory paths; don't forget trailing dashes!
# delivered_data_path_SB = "/works/yamato/MAPS_data/2018.1.01055.L/science_goal.uid___A001_X133d_X19de/group.uid___A001_X133d_X19df/member.uid___A001_X133d_X19e2/calibrated/"
delivered_data_path_LB = "/works/yamato/MAPS_data/MWC_480/Band3/pipeline-calibrated/"
processingdir = "/works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/"
MSdir = "/works/yamato/MAPS_data/MWC_480/Band3/self-calibrated/"
field = "MWC_480"
prefix = "MWC_480"
data_params_json_name = processingdir + prefix + "_data_params.json"

try:
    data_params = load_data_params(data_params_json_name)

except FileNotFoundError:
    data_params = {
        "LB1": {
            "vis": delivered_data_path_LB + "uid___A002_Xe0c2b4_X6b5.ms.split.cal",
            "spws": "25,27,29,31,33,35,37,39,41",
            "field": field,
            "column": "data",
        },
        "LB2": {
            "vis": delivered_data_path_LB + "uid___A002_Xe0cd4d_X5810.ms.split.cal",
            "spws": "25,27,29,31,33,35,37,39,41",
            "field": field,
            "column": "data",
        },
        "LB3": {
            "vis": delivered_data_path_LB + "uid___A002_Xe0cd4d_X5c6e.ms.split.cal",
            "spws": "25,27,29,31,33,35,37,39,41",
            "field": field,
            "column": "data",
        },
        "LB4": {
            "vis": delivered_data_path_LB + "uid___A002_Xe0cd4d_Xcd12.ms.split.cal",
            "spws": "25,27,29,31,33,35,37,39,41",
            "field": field,
            "column": "data",
        },
        "LB5": {
            "vis": delivered_data_path_LB + "uid___A002_Xe0cd4d_Xd1a5.ms.split.cal",
            "spws": "25,27,29,31,33,35,37,39,41",
            "field": field,
            "column": "data",
        },
        "LB6": {
            "vis": delivered_data_path_LB + "uid___A002_Xeee095_Xaa01.ms.split.cal",
            "spws": "25,27,29,31,33,35,37,39,41",
            "field": field,
            "column": "data",
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

# first make dirty image for all spws for visual inspection
# spws = ["0", "1", "2", "3", "4", "5", "6", "7", "8"]
spws = ["0"]
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
flagchannels = (
    "0:60~120;150~210;250~310;420~479,2:200~280,3:200~280,4:300~650,5:420~550,8:400~570"
)
spw = "0~8"
width = [
    480,
    480,
    480,
    480,
    960,
    960,
    60,
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

# LB1
# 04h58m46.279384s +29d50m36.39928s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279384s +29d50m36.39928s
# PA of Gaussian component: 149.37 deg
# Inclination of Gaussian component: 42.97 deg
# Pixel coordinates of peak: x = 499.920 y = 500.045
# LB2
# 04h58m46.279382s +29d50m36.39927s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279382s +29d50m36.39927s
# PA of Gaussian component: 149.39 deg
# Inclination of Gaussian component: 33.85 deg
# Pixel coordinates of peak: x = 499.922 y = 500.044
# LB3
# 04h58m46.279377s +29d50m36.39928s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279377s +29d50m36.39928s
# PA of Gaussian component: 140.67 deg
# Inclination of Gaussian component: 43.49 deg
# Pixel coordinates of peak: x = 499.928 y = 500.046
# LB4
# 04h58m46.279377s +29d50m36.39930s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279377s +29d50m36.39930s
# PA of Gaussian component: 51.33 deg
# Inclination of Gaussian component: 21.68 deg
# Pixel coordinates of peak: x = 499.928 y = 500.047
# LB5
# 04h58m46.279378s +29d50m36.39925s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279378s +29d50m36.39925s
# PA of Gaussian component: 128.78 deg
# Inclination of Gaussian component: 33.35 deg
# Pixel coordinates of peak: x = 499.928 y = 500.042
# LB6
# 04h58m46.279377s +29d50m36.39923s
# Peak of Gaussian component identified with imfit: J2000 04h58m46.279377s +29d50m36.39923s
# PA of Gaussian component: 150.02 deg
# Inclination of Gaussian component: 33.91 deg
# Pixel coordinates of peak: x = 499.929 y = 500.041

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

# LB5 looks really bad, with at most factor of ~5 difference from otehr datasets in its amplitude; LB4 also looks bad maybe due to the decorrelation
# consider if we can drop LB4 and LB5; these two are observed at bad weather (PWV > 2.0)
# first make a combined continuum image with all datasets and then compare with the case without LB4 and LB5

# image combined all datasets
imagename = processingdir + prefix + "_initcont_all"
tclean_continuum_wrapper(
    vis=[data_params[i]["vis_split_avg_shift"] for i in data_params.keys()],
    imagename=imagename,
    scales=scales,
    mask=disk_mask,
    imsize=imsize,
    cellsize=cellsize,
    nsigma=3,
    robust=0.5,
)

# without LB4 and LB5
imagename = processingdir + prefix + "_initcont_wo_LB4_LB5"
tclean_continuum_wrapper(
    vis=[
        data_params[i]["vis_split_avg_shift"]
        for i in data_params.keys()
        if i != "LB4" and i != "LB5"
    ],
    imagename=imagename,
    scales=scales,
    mask=disk_mask,
    imsize=imsize,
    cellsize=cellsize,
    nsigma=3,
    robust=0.5,
)

# without LB5
imagename = processingdir + prefix + "_initcont_wo_LB5"
tclean_continuum_wrapper(
    vis=[
        data_params[i]["vis_split_avg_shift"] for i in data_params.keys() if i != "LB5"
    ],
    imagename=imagename,
    scales=scales,
    mask=disk_mask,
    imsize=imsize,
    cellsize=cellsize,
    nsigma=3,
    robust=0.5,
)

# when without LB4 and LB5, resulting rms noise level slightly increases (1.9e-5 -> 2.1e-5), but without LB5, almost no difference is identified
# from here, we disregard LB5

# save original data_params and delete
save_data_params(data_params, prefix + "_data_params_original.json")

del data_params["LB5"]


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
    solint = get_solints(data_params[i][visID], field)

# LB1
# Determined solints:  ['inf', '26.21s', '12.10s', '6.05s', 'int']
# LB2
# Determined solints:  ['inf', '26.21s', '12.10s', '6.05s', 'int']
# LB3
# Determined solints:  ['inf', '26.21s', '12.10s', '6.05s', 'int']
# LB4
# Determined solints:  ['inf', '26.21s', '12.10s', '6.05s', 'int']
# LB6
# Determined solints:  ['inf', '26.21s', '12.10s', '6.05s', 'int']

# first make a test image
scales = [0, 5, 15, 25]
cellsize = "0.01arcsec"
imsize = su.getOptimumSize(1000)
imagename = processingdir + prefix + "_initial"
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
# nterms = 1
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_initial.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 29.73 mJy
# Peak intensity of source: 10.31 mJy/beam
# rms: 1.74e-02 mJy/beam
# Peak SNR: 591.75

# nterms = 2
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_initial.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 27.08 mJy
# Peak intensity of source: 9.45 mJy/beam
# rms: 1.94e-02 mJy/beam
# Peak SNR: 487.93

# below use nterms = 1

nterms = 1

print(initial_SNR)
import numpy as np

nsigma_init = np.max(
    [initial_SNR / 15.0, 5.0]
)  # restricts initial nsigma to be at least 5
nsigma_per_solint = 10 ** np.linspace(np.log10(nsigma_init), np.log10(3.0), len(solint))
print(nsigma_per_solint)
# [32.52866906 17.9258614   9.87856301  5.4438671   3.        ]


# iteration 1
mode = None
iteration = 1
solint = "inf"
combine = "spw"
calmode = "p"
nsigma = 32.53
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
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=nterms
    ),
    showplot=False,
)
###### pre selfcal iteration 1 ######
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_initial.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 31.22 mJy
# Peak intensity of source: 10.47 mJy/beam
# rms: 2.05e-02 mJy/beam
# Peak SNR: 509.51
###### post selfcal iteration 1 ######
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_iter1_p.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 31.09 mJy
# Peak intensity of source: 11.29 mJy/beam
# rms: 1.56e-02 mJy/beam
# Peak SNR: 723.54
save_data_params(data_params, prefix + "_data_params.json")


# iteration 2
mode = None
iteration = 2
solint = "26.21s"
combine = "spw"
calmode = "p"
nsigma = 17.93
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
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=nterms
    ),
    showplot=False,
)
###### pre selfcal iteration 2 ######
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_iter1_p.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 30.48 mJy
# Peak intensity of source: 11.24 mJy/beam
# rms: 1.39e-02 mJy/beam
# Peak SNR: 808.98
###### post selfcal iteration 2 ######
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_iter2_p.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 30.39 mJy
# Peak intensity of source: 11.28 mJy/beam
# rms: 1.39e-02 mJy/beam
# Peak SNR: 811.15
save_data_params(data_params, prefix + "_data_params.json")

# iteration 3
mode = None
iteration = 3
solint = "12.10s"
combine = "spw"
calmode = "p"
nsigma = 9.879
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
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=nterms
    ),
    showplot=False,
)
###### pre selfcal iteration 3 ######
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_iter2_p.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 30.11 mJy
# Peak intensity of source: 11.25 mJy/beam
# rms: 1.34e-02 mJy/beam
# Peak SNR: 838.46
###### post selfcal iteration 3 ######
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_iter3_p.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 30.22 mJy
# Peak intensity of source: 11.37 mJy/beam
# rms: 1.33e-02 mJy/beam
# Peak SNR: 851.80
save_data_params(data_params, prefix + "_data_params.json")

# iteration 4
mode = None
iteration = 4
solint = "6.05s"
combine = "spw"
calmode = "p"
nsigma = 5.444
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
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=nterms
    ),
    showplot=False,
)
###### pre selfcal iteration 4 ######
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_iter3_p.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 30.07 mJy
# Peak intensity of source: 11.33 mJy/beam
# rms: 1.31e-02 mJy/beam
# Peak SNR: 863.31
###### post selfcal iteration 4 ######
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_iter4_p.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 30.10 mJy
# Peak intensity of source: 11.34 mJy/beam
# rms: 1.30e-02 mJy/beam
# Peak SNR: 873.12
save_data_params(data_params, prefix + "_data_params.json")

# iteration 5
mode = None
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
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=nterms
    ),
    showplot=False,
)
###### pre selfcal iteration 5 ######
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_iter4_p.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 30.01 mJy
# Peak intensity of source: 11.31 mJy/beam
# rms: 1.29e-02 mJy/beam
# Peak SNR: 879.03
###### post selfcal iteration 5 ######
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_iter5_p.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 29.63 mJy
# Peak intensity of source: 10.47 mJy/beam
# rms: 1.60e-02 mJy/beam
# Peak SNR: 655.24

# S/N decreases, going to phase+amplitude selfcal

### phase + amplitude selfcal
# iteration 5
mode = None
iteration = 5
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
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=nterms
    ),
    showplot=False,
)
###### pre selfcal iteration 5 ######
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_iter4_p.image.tt0
# Beam 0.297 arcsec x 0.214 arcsec (-17.95 deg)
# Flux inside disk mask: 30.01 mJy
# Peak intensity of source: 11.31 mJy/beam
# rms: 1.29e-02 mJy/beam
# Peak SNR: 879.03
###### post selfcal iteration 5 ######
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_iter5_ap.image.tt0
# Beam 0.302 arcsec x 0.216 arcsec (-18.95 deg)
# Flux inside disk mask: 29.51 mJy
# Peak intensity of source: 11.49 mJy/beam
# rms: 1.24e-02 mJy/beam
# Peak SNR: 926.67
save_data_params(data_params, data_params_json_name)


# finalize at iteration = 5
iteration = 5
mode = None
calmode = "ap"
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
        scales=scales, mask=disk_mask, nsigma=nsigma, robust=0.5, nterms=nterms
    ),
)
# /works/yamato/MAPS_data/MWC_480/Band3/selfcal_processing/MWC_480_selfcal_final.image.tt0
# Beam 0.302 arcsec x 0.216 arcsec (-18.95 deg)
# Flux inside disk mask: 29.51 mJy
# Peak intensity of source: 11.49 mJy/beam
# rms: 1.24e-02 mJy/beam
# Peak SNR: 926.67

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
