import os
import sys
sys.path.append("/home/yamato/Project/C2D_disk/")
from JvM_correction import apply_JvM_correction

execfile("/home/yamato/Application/reduction_tools/reduction_utils_mpi.py", globals())

prefix = "MWC_480"
processingdir = "/works/yamato/C2D_disk/processing/"
imagedir = "/works/yamato/C2D_disk/v1_product/image/"
field = "MWC_480"

data_params = load_data_params(processingdir + prefix + "_data_params.json")
vislist = [data_params[i]["vis_spectral_line"] for i in data_params.keys()]
vislist_wcont = [data_params[i]["vis_spectral_line_wcont"] for i in data_params.keys()]

### continuum imaging parameters
cellsize = "0.02arcsec"
imsize = get_imsize(
	vis=vislist,
	field=field,
	cellsize=cellsize,
)
scales = [0, 10, 30]
robust_list = [0.5, 2.0]
spw_list = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]

########### continuum imaging ###########
for spw in spw_list:
	print(f"spw = {spw}")
	for robust in robust_list:
		print(f"robust = {robust}")
		rms = calc_sensitivity(vislist, cellsize=cellsize, imsize=imsize, robust=robust, specmode="cube", spw=[spw], chan=20)
		print(rms)
		imagename = imagedir + prefix + f"_spw{spw}_robust{robust}"
		tclean_spectral_line_wrapper(
			vislist,
			imagename,
			spw,
			imsize,
			cellsize,
			scales=scales,
			phasecenter="",
			weighting="briggsbwtaper",
			perchanweightdensity=True,
			robust=robust,
			mask="",
			noisethreshold=5.0,
			sidelobethreshold=3.0,
			threshold=f"{rms*3}Jy",
			parallel=True,
		)
		# primary beam correction
		os.system("rm -r " + imagename + ".pbcor")
		impbcor(
			imagename=imagename + ".image",
			pbimage=imagename + ".pb",
			outfile=imagename + ".pbcor",
		)
		# export to fits
		for ext in [".image", ".pbcor"]:
			exportfits(imagename=imagename + ext, fitsimage=imagename + ext + ".fits", dropstokes=True, overwrite=True)
		# apply JvM correction
		apply_JvM_correction(imagename, mfs=False, pbcor=True)





# ### line imaging
# vislist = {
#     "SB": [data_params[i]["vis_contsub"] for i in data_params.keys() if "SB" in i],
#     "LB": [data_params[i]["vis_contsub"] for i in data_params.keys() if "LB" in i],
#     "SBLB": [data_params[i]["vis_contsub"] for i in data_params.keys()],
# }

# image_list = {
#     "NH2D": dict(channelwidth=0.5, vrange=(-12, 20)),
#     "HC18O+": dict(channelwidth=0.5, vrange=(-12, 20)),
#     "H13CO+": dict(channelwidth=0.5, vrange=(-12, 20)),
#     "HOCO+_k0": dict(channelwidth=0.5, vrange=(-12, 20)),
#     "HOCO+_k1": dict(channelwidth=0.5, vrange=(-12, 20)),
#     "CH3OH_96.492GHz": dict(channelwidth=0.5, vrange=(-12, 20)),
#     "CH3OH_96.745GHz": dict(channelwidth=0.5, vrange=(-12, 20)),
#     "CH3OH_96.756GHz": dict(channelwidth=0.5, vrange=(-12, 20)),
#     "HCO": dict(channelwidth=0.5, vrange=(-12, 20)),
#     "CCH": dict(channelwidth=0.5, vrange=(-12, 20)),
# }

# # clean lines SB, LB, SBLB
# robust = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]
# for line in image_list.keys():
#     for i in vislist.keys():
#         if "LB" in i:
#             noisethreshold = 5.0
#             sidelobethreshold = 3.0
#         else:
#             noisethreshold = 4.25
#             sidelobethreshold = 2.0
#         for r in robust:
#             print(line, i, r)
#             imagename = (
#                 WD_path + "_".join([prefix, i, line, "robust{}".format(r)]) + "_dirty"
#             )
#             tclean_spectral_line_dirty(
#                 vis=vislist[i],
#                 imagename=imagename,
#                 spw=line_dict[line]["spw"],
#                 restfreq=line_dict[line]["restfreq"],
#                 weighting="briggs",
#                 robust=r,
#                 **image_list[line],
#             )
#             rms = imstat(imagename=imagename + ".image", region=noise_mask)["rms"][0]
#             imagename = WD_path + "_".join([prefix, i, line, "robust{}".format(r)])
#             tclean_spectral_line(
#                 vis=vislist[i],
#                 imagename=imagename,
#                 spw=line_dict[line]["spw"],
#                 restfreq=line_dict[line]["restfreq"],
#                 weighting="briggs",
#                 robust=r,
#                 noisethreshold=noisethreshold,
#                 sidelobethreshold=sidelobethreshold,
#                 threshold="{}Jy".format(3 * rms),
#                 **image_list[line],
#             )
#             # primary beam correction
#             os.system("rm -r " + imagename + ".pbcor")
#             impbcor(
#                 imagename=imagename + ".image",
#                 pbimage=imagename + ".pb",
#                 outfile=imagename + ".pbcor",
#             )
#             # export to fits
#             for ext in [".image", ".mask", ".pbcor"]:
#                 exportfits(
#                     imagename=imagename + ext,
#                     fitsimage=imagename + ext + ".fits",
#                     overwrite=True,
#                     dropdeg=True,
#                 )