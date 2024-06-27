"""
Script to convert electron diffraction frame files to tif format.
Supported file types: CBF, IMG

This program requires Python and the package 'fabio'.
Install via pip: pip install fabaio

CALL PROGRAM:
If you run the program, you are asked to provide the path to the folder with the frames.
python Convert2tif.py

Alternatively, you can provide the path as an argument.
python Convert2tif.py C:/Path/to/frames

INPUT:      Path to folder with CBF or IMG files

This program does the following:
1. Scan folder for CBF or IMG files
2. Create subfolder "frames"
3. Convert any CBF / IMG file to 16-bit TIF and save it in the subfolder "frames"
"""

#__version__ = "13 June 2021" # Paul B. Klar, FZU Palatinus Group
__version__ = "27 June 2024"
__author___ = "paul.klar@uni-bremen.de" # Paul B. Klar, University of Bremen

from sys import argv as argument
from pathlib import Path
import fabio     # https://github.com/silx-kit/fabio  https://pythonhosted.org/fabio/api/modules.html

#import numpy as np
#import tifffile

folder = None

# check folder as provided argument
if len(argument) == 2:
    folder = Path(arguments[1])
    if not folder.exists():
        folder = None
        print('Folder does not exist:')
        print(folder)
    
while not folder:
    print("Please copy + paste here the path to the folder with frames in CBF or IMG format:")
    folder = Path(input())
    if not folder.exists():
        folder = None
        print("The folder does not exist.")
            
files_cbf = folder.glob('*.cbf')
files_img = folder.glob('*.img')
files = list(files_cbf) + list(files_img)

if files:
    print('Converting', len(files), 'files.')
    
    target_folder = folder / "frames"
    move_frames_folder = folder / "original"
    
    if not target_folder.exists():
        target_folder.mkdir()

    total = len(files)
    for i, file in enumerate(files):
        FrameFile = target_folder / (file.stem + ".tif")
        
        # open file
        Image = fabio.open(file)
        
        # convert and save
        Image.convert("tif").save(FrameFile)
        Image.close()
        
        print(FrameFile.name, " written.", i+1, "/", total)
    
print('The converted files are in the folder:')
print(target_folder)
print()  
input("... END OF PROGRAM ...")