'''
This program reads in all tiff files from a given folder.
A rectangular region of interest is defined and extracted.
New tiff files are written.
New files are of type unsigned int 16 bit.
'''

# Python3.8
__version__ = "12 July 2024"
__author___ = "paul.klar@uni-bremen.de" # Paul B. Klar, Uni Bremen


import numpy as np
from pathlib import Path
from PIL import Image


# source path with tiff files
folder_from = Path(r"D:\ED_data_PlanetEarth\2024-02_Alke_Biotite\ed_rot_step_003_2_corr3\frames")

# destination for new tiff files
folder_to = Path(r"D:\ED_data_PlanetEarth\2024-02_Alke_Biotite\ed_rot_step_003_2_corr3\frames_cut")

files = folder_from.glob("*.tif")

for tif_file in folder_from.glob("*.tif"):
    print(tif_file.name)
    img = Image.open(tif_file)
    image_array = np.array(img)

    # define region of interest here <<<<<<<<<<           <<<<<<<<<<<<<<<<<<<           <<<<<<<
    image_array = image_array[144:903+1, 213:972+1]

    # make unsigned int 16-bit
    image_array[image_array > 2**16-1] = 2**16-1
    image_array[image_array < 0] = 0
    new_img = Image.fromarray(image_array.astype(np.uint16))
    new_img.save(folder_to / tif_file.name) #, software="python3_PIL")
    #break

input("... END OF PROGRAM ...")
