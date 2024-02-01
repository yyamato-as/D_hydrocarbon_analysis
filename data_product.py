import argparse

parser = argparse.ArgumentParser(
    description="Script to generate value-added data products (VADP) of molecular line data toward three disks, AS 209, HD 163296, and MWC 480"
)

parser.add_argument(
    "--source",
    help="Source name for which the VADP will be generated",
    nargs="*",
)
parser.add_argument(
    "--species",
    help="Molecular species for which the VADP will be generated",
    nargs="*",
)
parser.add_argument(
    "--line",
    help="Lines for which the VADP will be generated. The molecular line keys in linedictionary.py can be used",
    nargs="*",
)
parser.add_argument("--product", help="Product types to be generated. ", nargs="*")


# velocity-integrated intensity maps
# for source in disk_dict
#     for line in lines:
#         print(line)
#         imagename = imagepath + f"{source}_{line}_hybrid.image.pbcor.fits"
#         image = FitsImage(imagename)
#         mask = FitsImage(imagename.replace(".image.pbcor", ".mask"))
#         image.moment_map(moment="0", mask=mask.data, save=True, savefilename=imagename)
