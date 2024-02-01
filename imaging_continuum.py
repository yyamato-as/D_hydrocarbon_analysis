import os
import sys
sys.path.append("/home/yamato/Application/reduction_tools/")
from JvM_correction import apply_JvM_correction

execfile("/home/yamato/Application/reduction_tools/reduction_utils_mpi.py", globals())

prefix = "MWC_480"
processingdir = "/works/yamato/C2D_disk/processing/"
imagedir = "/works/yamato/C2D_disk/v1_product/image/"
field = "MWC_480"

data_params = load_data_params(processingdir + prefix + "_data_params.json")
vislist = [data_params[i]["vis_continuum"] for i in data_params.keys()]

### continuum imaging parameters
cellsize = "0.02arcsec"
imsize = get_imsize(
	vis=vislist,
	field=field,
	cellsize=cellsize,
)
scales = [0, 10, 30]
robust_list = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]
nterms = 2

########### continuum imaging ###########
for robust in robust_list:
	print(f"robust = {robust}")
	rms = calc_sensitivity(vislist, cellsize=cellsize, imsize=imsize, robust=robust)
	print(rms)
	imagename = imagedir + prefix + f"_continuum_robust{robust}"
	tclean_continuum_wrapper(
		vislist,
		imagename,
		imsize,
		cellsize,
		scales=scales,
		weighting="briggs",
		robust=robust,
		noisethreshold=5.0,
		sidelobethreshold=3.0,
		nterms=nterms,
		threshold=f"{rms*3}Jy",
		parallel=True,
	)
	# primary beam correction
	os.system("rm -r " + imagename + ".pbcor.tt0")
	impbcor(
		imagename=imagename + ".image.tt0",
		pbimage=imagename + ".pb.tt0",
		outfile=imagename + ".pbcor.tt0",
	)
	# export to fits
	for ext in [".image.tt0", ".pbcor.tt0"]:
		exportfits(imagename=imagename + ext, fitsimage=imagename + ext + ".fits", dropstokes=True, overwrite=True)
	# apply JvM correction
	apply_JvM_correction(imagename, mfs=True, pbcor=True)