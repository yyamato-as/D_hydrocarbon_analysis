import os
import sys

sys.path.append("/home/yamato/Project/C2D_disk/")

execfile("/home/yamato/Project/C2D_disk/reduction_utils_mpi_v3.py")
# import reduction_utils_mpi as ru

delivered_data_path = "/works/yamato/C2D_disk/2021.1.00535.S/science_goal.uid___A001_X15c6_X93/group.uid___A001_X15c6_X94/member.uid___A001_X15c6_X95/calibrated/working/"
processingdir = "/works/yamato/C2D_disk/processing/"
MSdir = "/works/yamato/C2D_disk/v1_product/measurement_set/"
imagedir = "/works/yamato/C2D_disk/v1_product/image/"
field = "MWC_480"
prefix = "MWC_480"
data_params_json_name = processingdir + prefix + "_data_params.json"

try:
	data_params = load_data_params(data_params_json_name)

except FileNotFoundError:
	data_params = {
		"EB1": {
			"vis": delivered_data_path + "uid___A002_Xfc8dc8_Xb1cf.ms",
			"spws": "23,25,27,29,31,33,35,37,39,41",
			"field": field,
			"column": "corrected",
		},
		"EB2": {
			"vis": delivered_data_path + "uid___A002_Xfca6fd_X2bf9.ms",
			"spws": "23,25,27,29,31,33,35,37,39,41",
			"field": field,
			"column": "corrected",
		},
		"EB3": {
			"vis": delivered_data_path + "uid___A002_Xfca6fd_X32c2.ms",
			"spws": "23,25,27,29,31,33,35,37,39,41",
			"field": field,
			"column": "corrected",
		},
		"EB4": {
			"vis": delivered_data_path + "uid___A002_Xfca6fd_X9c53.ms",
			"spws": "23,25,27,29,31,33,35,37,39,41",
			"field": field,
			"column": "corrected",
		},
		"EB5": {
			"vis": delivered_data_path + "uid___A002_Xfca6fd_Xb591.ms",
			"spws": "23,25,27,29,31,33,35,37,39,41",
			"field": field,
			"column": "corrected",
		},
	}

### Band 7 setup
# SpwID  #Chans   Frame   Ch0(MHz)        ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs
# 23     960      TOPO    288593.780      -122.070      117187.5   288535.2475        1  XX  YY -> CCD (N=4-3) *two lines
# 25     480      TOPO    286972.961      -122.070       58593.8   286943.7253        1  XX  YY -> c-HCCCD 9(0,9)-8(1,8)
# 27     480      TOPO    288183.380      -122.070       58593.8   288154.1440        1  XX  YY -> DCO+ (J=4-3)
# 29     480      TOPO    290662.940      -122.070       58593.8   290633.7045        2  XX  YY -> H2CO 4(0,4)-3(0,3)
# 31     480      TOPO    289248.664      -122.070       58593.8   289219.4284        2  XX  YY -> C34S (J=6-5)
# 33     480      TOPO    289979.041      -122.070       58593.8   289949.8056        2  XX  YY -> CH3OH 6(0,6)-5(0,5)
# 35     480      TOPO    289684.486      -122.070       58593.8   289655.2499        2  XX  YY -> DCN (J=4-3)
# 37    1920      TOPO    298894.271       976.562     1875000.0   299831.2831        3  XX  YY -> continuum (CH3OCH3, HC3N)
# 39     960      TOPO    301179.794       244.141      234375.0   301296.8593        4  XX  YY -> SO 3Sigma 7(7)-6(6)
# 41     960      TOPO    300730.331       244.141      234375.0   300847.3964        4  XX  YY -> H2CO 4(1,3)-3(1,2)


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

### split ms file spw
# SpwID #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs
# 0        960   TOPO  288593.780      -122.070    117187.5 288535.2475        1  XX  YY -> CCD (N=4-3) *two lines
# 1        480   TOPO  286972.961      -122.070     58593.8 286943.7253        1  XX  YY -> c-HCCCD 9(0,9)-8(1,8)
# 2        480   TOPO  288183.380      -122.070     58593.8 288154.1440        1  XX  YY -> DCO+ (J=4-3)
# 3        480   TOPO  290662.940      -122.070     58593.8 290633.7045        2  XX  YY -> H2CO 4(0,4)-3(0,3)
# 4        480   TOPO  289248.664      -122.070     58593.8 289219.4284        2  XX  YY -> C34S (J=6-5)
# 5        480   TOPO  289979.041      -122.070     58593.8 289949.8056        2  XX  YY -> CH3OH 6(0,6)-5(0,5)
# 6        480   TOPO  289684.486      -122.070     58593.8 289655.2499        2  XX  YY -> DCN (J=4-3)
# 7       1920   TOPO  298894.271       976.562   1875000.0 299831.2831        3  XX  YY -> continuum (CH3OCH3, HC3N)
# 8        960   TOPO  301179.794       244.141    234375.0 301296.8593        4  XX  YY -> SO 3Sigma 7(7)-6(6)
# 9        960   TOPO  300730.331       244.141    234375.0 300847.3964        4  XX  YY -> H2CO 4(1,3)-3(1,2)


# contspws = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

### flag the data based on the visual inspection on the pipeline calibrated data products
# Note that the continuum window has two lines and one absorption feature caused by a telluric line (see weblog spectra)
flagchannels = "0:150~380;650~800,1:160~320,2:160~320,3:160~320,4:160~320,5:160~320,6:160~320,7:950~1150;1250~1400;1780~1880,8:400~560,9:400~560"

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
	# average over the channels and split out the continuum visibility
	avg_vis = processingdir + prefix + f"_{i}_initcont.ms"
	os.system("rm -r " + avg_vis)
	split(
		vis=data_params[i]["vis_split"],
		outputvis=avg_vis,
		spw="0~9",
		field=field,
		width=[
			960,
			480,
			480,
			480,
			480,
			480,
			480,
			240,
			960,
			960,
		],  # Note that the continuum window is averaged over only 240 chans to avoid potential beam smearing effect; see https://casaguides.nrao.edu/index.php/Image_Continuum
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
scales = [0, 10, 30]
cellsize = "0.02arcsec"
imsize = get_imsize(
	vis=[data_params[i]["vis_split_avg"] for i in data_params.keys()],
	field=field,
	cellsize=cellsize,
)
# continuum disk mask; based on the visual inspection into the delivered data
disk_mask = "ellipse[[{:s}, {:s}], [{:.1f}arcsec, {:.1f}arcsec], {:.1f}deg]".format(
	"04h58m46.275s", "+29d50m36.43s", 2.5, 2.1, 342.0
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
		threshold="2mJy",
		robust=0.5,
	)

### fit each image with Gaussian to estimate the peak
for i in data_params.keys():
	data_params[i]["phasecenter"] = fit_gaussian(
		imagename=processingdir + prefix + f"_{i}_initial.image.tt0",
		region=disk_mask,
	)

# EB1
# 04h58m46.271792s +29d50m36.41795s
# #Peak of Gaussian component identified with imfit: ICRS 04h58m46.271792s +29d50m36.41795s
# 04:58:46.271792 +29:50:36.41795
# Separation: radian = 5.96181e-08, degrees = 0.000003 = 3.41586e-06, arcsec = 0.012297 = 0.0122971
# #Peak in J2000 coordinates: 04:58:46.27222, +029:50:36.406986
# #PA of Gaussian component: 149.86 deg
# #Inclination of Gaussian component: 35.07 deg
# #Pixel coordinates of peak: x = 510.010 y = 511.512

# EB2
# 04h58m46.272808s +29d50m36.40413s
# #Peak of Gaussian component identified with imfit: ICRS 04h58m46.272808s +29d50m36.40413s
# 04:58:46.272808 +29:50:36.40413
# Separation: radian = 5.97328e-08, degrees = 0.000003 = 3.42244e-06, arcsec = 0.012321 = 0.0123208
# #Peak in J2000 coordinates: 04:58:46.27324, +029:50:36.393166
# #PA of Gaussian component: 5.83 deg
# #Inclination of Gaussian component: 30.01 deg
# #Pixel coordinates of peak: x = 509.350 y = 510.821

# EB3
# 04h58m46.278955s +29d50m36.41021s
# #Peak of Gaussian component identified with imfit: ICRS 04h58m46.278955s +29d50m36.41021s
# 04:58:46.278955 +29:50:36.41021
# Separation: radian = 5.95326e-08, degrees = 0.000003 = 3.41097e-06, arcsec = 0.012279 = 0.0122795
# #Peak in J2000 coordinates: 04:58:46.27938, +029:50:36.399246
# #PA of Gaussian component: 151.53 deg
# #Inclination of Gaussian component: 40.03 deg
# #Pixel coordinates of peak: x = 505.351 y = 511.125

# EB4
# 04h58m46.273697s +29d50m36.39614s
# #Peak of Gaussian component identified with imfit: ICRS 04h58m46.273697s +29d50m36.39614s
# 04:58:46.273697 +29:50:36.39614
# Separation: radian = 5.97616e-08, degrees = 0.000003 = 3.42409e-06, arcsec = 0.012327 = 0.0123267
# #Peak in J2000 coordinates: 04:58:46.27413, +029:50:36.385176
# #PA of Gaussian component: 161.92 deg
# #Inclination of Gaussian component: 29.90 deg
# #Pixel coordinates of peak: x = 508.771 y = 510.422

# EB5
# 04h58m46.277417s +29d50m36.40478s
# #Peak of Gaussian component identified with imfit: ICRS 04h58m46.277417s +29d50m36.40478s
# 04:58:46.277417 +29:50:36.40478
# Separation: radian = 5.97616e-08, degrees = 0.000003 = 3.42409e-06, arcsec = 0.012327 = 0.0123267
# #Peak in J2000 coordinates: 04:58:46.27785, +029:50:36.393816
# #PA of Gaussian component: 145.76 deg
# #Inclination of Gaussian component: 42.19 deg
# #Pixel coordinates of peak: x = 506.351 y = 510.854


### There are slight shifts between the peaks of EBs
### To correct for this, we fix the common phase center and shift the peak of each EB to that
### common direction from the best-looking image (EB3)
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
		threshold="2mJy",
		robust=0.5,
	)

for i in data_params.keys():
	data_params[i]["phasecenter_aligned"] = fit_gaussian(
		imagename=processingdir + prefix + f"_{i}_initial_shifted.image.tt0",
		region=disk_mask,
	)
# confirm that the peak is well aligned

# EB1
# 04h58m46.279381s +29d50m36.39927s
# #Peak of Gaussian component identified with imfit: J2000 04h58m46.279381s +29d50m36.39927s
# #PA of Gaussian component: 149.68 deg
# #Inclination of Gaussian component: 35.11 deg
# #Pixel coordinates of peak: x = 511.962 y = 512.022

# EB2
# 04h58m46.279376s +29d50m36.39919s
# #Peak of Gaussian component identified with imfit: J2000 04h58m46.279376s +29d50m36.39919s
# #PA of Gaussian component: 5.56 deg
# #Inclination of Gaussian component: 29.95 deg
# #Pixel coordinates of peak: x = 511.965 y = 512.018

# EB3
# 04h58m46.279387s +29d50m36.39916s
# #Peak of Gaussian component identified with imfit: J2000 04h58m46.279387s +29d50m36.39916s
# #PA of Gaussian component: 151.33 deg
# #Inclination of Gaussian component: 40.13 deg
# #Pixel coordinates of peak: x = 511.958 y = 512.017

# EB4
# 04h58m46.279376s +29d50m36.39920s
# #Peak of Gaussian component identified with imfit: J2000 04h58m46.279376s +29d50m36.39920s
# #PA of Gaussian component: 161.79 deg
# #Inclination of Gaussian component: 29.86 deg
# #Pixel coordinates of peak: x = 511.965 y = 512.019

# EB5
# 04h58m46.279372s +29d50m36.39916s
# #Peak of Gaussian component identified with imfit: J2000 04h58m46.279372s +29d50m36.39916s
# #PA of Gaussian component: 145.70 deg
# #Inclination of Gaussian component: 42.23 deg
# #Pixel coordinates of peak: x = 511.967 y = 512.017

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

"""
# Visibilities are nicely aligned, but EB4 and EB5 are slightly offset from others; just check the scaling factor...
refdata = "EB3"

reference = [filename for filename in export_vislist if refdata in filename][0]
for i in data_params.keys():
	print(i)
	if i != refdata:
		comparison = [filename for filename in export_vislist if i in filename][0]
		data_params[i]["gencal_scale"] = estimate_flux_scale(
			reference=reference,
			incl=incl,
			PA=PA,
			comparison=comparison,
			outfile=f"scalefactor-vs-uvdistance-pre-scaling-{i}.png",
		)
	else:
		data_params[i]["gencal_scale"] = 1.0

# EB1
# The ratio of the fluxes of /works/yamato/C2D_disk/processing/MWC_480_EB1_initcont_shift.vis.npz to /works/yamato/C2D_disk/processing/MWC_480_EB3_initcont_shift.vis.npz is 1.12022
# The scaling factor for gencal is 1.058 for your comparison measurement
# The error on the weighted mean ratio is 1.903e-04, although it's likely that the weights in the measurement sets are too off by some constant factor

# EB2
# The ratio of the fluxes of /works/yamato/C2D_disk/processing/MWC_480_EB2_initcont_shift.vis.npz to /works/yamato/C2D_disk/processing/MWC_480_EB3_initcont_shift.vis.npz is 1.09443
# The scaling factor for gencal is 1.046 for your comparison measurement
# The error on the weighted mean ratio is 1.587e-04, although it's likely that the weights in the measurement sets are too off by some constant factor

# EB3

# EB4
# The ratio of the fluxes of /works/yamato/C2D_disk/processing/MWC_480_EB4_initcont_shift.vis.npz to /works/yamato/C2D_disk/processing/MWC_480_EB3_initcont_shift.vis.npz is 0.96998
# The scaling factor for gencal is 0.985 for your comparison measurement
# The error on the weighted mean ratio is 1.477e-04, although it's likely that the weights in the measurement sets are too off by some constant factor

# EB5
# The ratio of the fluxes of /works/yamato/C2D_disk/processing/MWC_480_EB5_initcont_shift.vis.npz to /works/yamato/C2D_disk/processing/MWC_480_EB3_initcont_shift.vis.npz is 0.92884
# The scaling factor for gencal is 0.964 for your comparison measurement
# The error on the weighted mean ratio is 1.484e-04, although it's likely that the weights in the measurement sets are too off by some constant factor


### The flux ratio as well as the scaling factor is within ~10%, which is reasonably within the flux calibration uncertainty of 10% (2sigma) in Band 7
### Therefore we do not apply any rescaling here
"""

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

# EB1
# Determined solints:  ['inf', '90.72s', '42.34s', '18.14s', 'int']
# EB2
# Determined solints:  ['inf', '90.72s', '42.34s', '18.14s', 'int']
# EB3
# Determined solints:  ['inf', '90.72s', '42.34s', '18.14s', 'int']
# EB4
# Determined solints:  ['inf', '90.72s', '42.34s', '18.14s', 'int']
# EB5
# Determined solints:  ['inf', '90.72s', '42.34s', '18.14s', 'int']
initial_SNR = 338.62
import numpy as np
nsigma_init = np.max([initial_SNR / 15.0, 5.0]) # restricts initial nsigma to be at least 5
nsigma_per_solint = 10 ** np.linspace(np.log10(nsigma_init), np.log10(3.0), len(solint))
print(nsigma_per_solint)
# [22.57466667 13.63001467  8.22945928  4.96874006  3.        ]


for i in data_params.keys():
	bw = get_effective_bandwidth(data_params[i]["vis_split"], data_params[i]["field"], flagchannels)
	nu0 = get_mean_frequency(data_params[i][visID], data_params[i]["field"])
	print(bw, nu0)

### The theoretical RMS
# bandwidth = 2.13 GHz
# nu0 = 298 GHz
# number of antenna = 41 (average over 5 executions)
# with excellent PWV (~0.3 mm in average!)
# int. time = 43.85 min/executions x 5 executions = 219.25 min
# ---> ~20.3 uJy by sensitivity calculator

# iteration 1
mode = None
iteration = 1
solint = "inf"
combine = "spw"
calmode = "p"
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
	tclean_kwargs=dict(scales=scales, mask=disk_mask, nsigma=22.57, robust=0.5, nterms=2),
	showplot=False,
)
###### pre selfcal iteration 1 ######
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_initial.image.tt0
#Beam 0.314 arcsec x 0.207 arcsec (7.18 deg)
#Flux inside disk mask: 513.67 mJy
#Peak intensity of source: 104.55 mJy/beam
#rms: 3.36e-01 mJy/beam
#Peak SNR: 310.99
###### post selfcal iteration 1 ######
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_iter1_p.image.tt0
#Beam 0.313 arcsec x 0.207 arcsec (7.22 deg)
#Flux inside disk mask: 533.59 mJy
#Peak intensity of source: 120.30 mJy/beam
#rms: 8.91e-02 mJy/beam
#Peak SNR: 1350.66

# ref.
# nterms=1: S/N ~ 1700
# nterms=3: S/N ~ 1200

# iteration 2
mode = None
iteration = 2
solint = "90.72s"
combine = "spw"
calmode = "p"
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
	tclean_kwargs=dict(scales=scales, mask=disk_mask, nsigma=13.63, robust=0.5, nterms=2),
	showplot=False,
)
###### pre selfcal iteration 2 ######
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_iter1_p.image.tt0
#Beam 0.313 arcsec x 0.207 arcsec (7.22 deg)
#Flux inside disk mask: 533.10 mJy
#Peak intensity of source: 120.49 mJy/beam
#rms: 7.61e-02 mJy/beam
#Peak SNR: 1583.79
###### post selfcal iteration 2 ######
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_iter2_p.image.tt0
#Beam 0.313 arcsec x 0.207 arcsec (7.22 deg)
#Flux inside disk mask: 534.04 mJy
#Peak intensity of source: 123.18 mJy/beam
#rms: 6.66e-02 mJy/beam
#Peak SNR: 1850.74

# iteration 3
mode = None
iteration = 3
solint = "42.3s"
combine = "spw"
calmode = "p"
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
	tclean_kwargs=dict(scales=scales, mask=disk_mask, nsigma=8.23, robust=0.5, nterms=2),
	showplot=False,
)
###### pre selfcal iteration 3 ######
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_iter2_p.image.tt0
#Beam 0.313 arcsec x 0.207 arcsec (7.22 deg)
#Flux inside disk mask: 533.86 mJy
#Peak intensity of source: 123.27 mJy/beam
#rms: 6.21e-02 mJy/beam
#Peak SNR: 1986.14
###### post selfcal iteration 3 ######
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_iter3_p.image.tt0
#Beam 0.313 arcsec x 0.207 arcsec (7.22 deg)
#Flux inside disk mask: 535.43 mJy
#Peak intensity of source: 126.16 mJy/beam
#rms: 5.43e-02 mJy/beam
#Peak SNR: 2322.06

# iteration 4
mode = None
iteration = 4
solint = "18.14s"
combine = "spw"
calmode = "p"
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
	tclean_kwargs=dict(scales=scales, mask=disk_mask, nsigma=4.969, robust=0.5, nterms=2),
	showplot=False,
)
###### pre selfcal iteration 4 ######
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_iter3_p.image.tt0
#Beam 0.313 arcsec x 0.207 arcsec (7.22 deg)
#Flux inside disk mask: 535.36 mJy
#Peak intensity of source: 126.23 mJy/beam
#rms: 5.26e-02 mJy/beam
#Peak SNR: 2399.74
###### post selfcal iteration 4 ######
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_iter4_p.image.tt0
#Beam 0.313 arcsec x 0.207 arcsec (7.22 deg)
#Flux inside disk mask: 536.75 mJy
#Peak intensity of source: 128.04 mJy/beam
#rms: 4.89e-02 mJy/beam
#Peak SNR: 2619.49

# iteration 5
mode = None
iteration = 5
solint = "int"
combine = "spw"
calmode = "p"
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
	tclean_kwargs=dict(scales=scales, mask=disk_mask, nsigma=3.0, robust=0.5, nterms=2),
	showplot=False,
)
###### pre selfcal iteration 5 ######
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_iter4_p.image.tt0
#Beam 0.313 arcsec x 0.207 arcsec (7.22 deg)
#Flux inside disk mask: 536.69 mJy
#Peak intensity of source: 128.11 mJy/beam
#rms: 4.82e-02 mJy/beam
#Peak SNR: 2655.95
###### post selfcal iteration 5 ######
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_iter5_p.image.tt0
#Beam 0.313 arcsec x 0.207 arcsec (7.22 deg)
#Flux inside disk mask: 538.32 mJy
#Peak intensity of source: 128.92 mJy/beam
#rms: 4.70e-02 mJy/beam
#Peak SNR: 2742.97

### phase + amplitude selfcal
# iteration 6
mode = None
iteration = 6
solint = "inf"
combine = "spw"
calmode = "ap"
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
	tclean_kwargs=dict(scales=scales, mask=disk_mask, nsigma=2.0, robust=0.5, nterms=2),
	showplot=False,
)
###### pre selfcal iteration 6 ######
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_iter5_p.image.tt0
#Beam 0.313 arcsec x 0.207 arcsec (7.22 deg)
#Flux inside disk mask: 538.28 mJy
#Peak intensity of source: 128.95 mJy/beam
#rms: 4.68e-02 mJy/beam
#Peak SNR: 2755.06
###### post selfcal iteration 6 ######
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_iter6_ap.image.tt0
#Beam 0.314 arcsec x 0.207 arcsec (8.17 deg)
#Flux inside disk mask: 538.91 mJy
#Peak intensity of source: 129.22 mJy/beam
#rms: 2.77e-02 mJy/beam
#Peak SNR: 4672.72


### nearly thermal noise level, stop here and finalize the selfcal
iteration = 6
mode = None
calmode = "ap"
finalize_selfcal(
	data_params,
	processingdir,
	prefix,
	iteration,
	visID,
	imsize,
	cellsize,
	mode=mode,
	calmode="ap",
	disk_mask=disk_mask,
	noise_mask=noise_mask,
	tclean_kwargs=dict(scales=scales, mask=disk_mask, nsigma=2.0, robust=0.5, nterms=2),
)
#/works/yamato/C2D_disk/processing/MWC_480_selfcal_final.image.tt0
#Beam 0.314 arcsec x 0.207 arcsec (8.17 deg)
#Flux inside disk mask: 538.91 mJy
#Peak intensity of source: 129.22 mJy/beam
#rms: 2.77e-02 mJy/beam
#Peak SNR: 4672.72

save_data_params(data_params, data_params_json_name)

###### split off the final continuum data ######
for i in data_params.keys():
	print(i)
	outputvis = MSdir + prefix + '_' + i + '_continuum.ms'
	os.system('rm -r '+ outputvis)
	split(vis=data_params[i][visID + '_selfcal'], outputvis=outputvis, datacolumn='data')
	data_params[i]['vis_continuum'] = outputvis

save_data_params(data_params, data_params_json_name)


########### apply the solution to the line measurement sets ###########
# shift the center to the common_dir
for i in data_params.keys():
	print(i)
	data_params[i]["vis_split_shift"] = (
		processingdir + prefix + "_" + i + "_shift.ms"
	)
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
	ntables = len(data_params[i]['selfcal_table'])
	applycal(vis=data_params[i]['vis_split_shift'], spw='', 
			 gaintable=data_params[i]['selfcal_table'],
			 spwmap=data_params[i]['selfcal_spwmap'],
			 interp=['linearPD']*ntables, calwt=True, applymode="calonly")
	data_params[i]["vis_split_shift_selfcal"] = data_params[i]['vis_split_shift'].replace(".ms", "_selfcal.ms")
	os.system("rm -r " + data_params[i]["vis_split_shift_selfcal"])
	split(vis=data_params[i]['vis_split_shift'],
		  outputvis=data_params[i]["vis_split_shift_selfcal"],
		  datacolumn='corrected')

save_data_params(data_params, data_params_json_name)


########## continuum subtraction ##########
for i in data_params.keys():
	print(i, flagchannels)
	os.system("rm -r " + data_params[i]["vis_split_shift_selfcal"] + ".contsub")
	uvcontsub(
		vis=data_params[i]["vis_split_shift_selfcal"],
		fitspw=flagchannels,
		excludechans=True,
		combine="spw",
		fitorder=1,
	)
	data_params[i]["vis_split_shift_selfcal_contsub"] = data_params[i]["vis_split_shift_selfcal"] + ".contsub"

save_data_params(data_params, data_params_json_name)

# locate the final mesurement sets to the product directory
for i in data_params.keys():
	data_params[i]["vis_spectral_line_wcont"] = MSdir + prefix + "_" + i + "_spectral_line_wcont.ms"
	data_params[i]["vis_spectral_line"] = MSdir + prefix + "_" + i + "_spectral_line.ms"
	os.system("cp -r " + data_params[i]["vis_split_shift_selfcal"] + " " + data_params[i]["vis_spectral_line_wcont"])
	os.system("cp -r " + data_params[i]["vis_split_shift_selfcal_contsub"] + " " + data_params[i]["vis_spectral_line"])

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
