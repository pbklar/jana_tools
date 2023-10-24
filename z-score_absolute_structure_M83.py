""" 27 July 2023
by Paul B. Klar, Department of Geosciences, University of Bremen

**Reference**:
Klar, P.B., Krysiak, Y., Xu, H. et al. Accurate structure models and absolute configuration determination
using dynamical effects in continuous-rotation 3D electron diffraction data. 
Nat. Chem. 15, 848–855 (2023). https://doi.org/10.1038/s41557-023-01186-1

This code makes use specifically of the subsection **Absolute structure determination** 
of the Methods part in the reference.  
"""

import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pathlib import Path
from sys import argv as arguments

def calcZP(df1, df2, selector=...):
    """ calculate Z score and related parameters.
    If z > 0, then the probability that the refinement corresponding to df1
    corresponds to the correct enantiomorph is p.
    
    Here, uncertainties of Io have no influence on the calculation of z.
    """    
    # N reflections in comparison
    if selector is ...:
        N = len(df1["DeltaI"])
    else:
        N = np.count_nonzero( selector )
    if N == 0:
        return float("NaN"), float("NaN"), 0, 0
    
    # k reflections are better fitted with refinement of df1
    k = np.count_nonzero( np.abs(df1["DeltaI"][selector]) < np.abs(df2["DeltaI"][selector]) )
    
    # z-score
    z = (k-N/2) / (math.sqrt(N)/2)
    
    # associated probability for correct assignment
    # p = Phi(z) = cumulative distribution function of standard normal distribution
    p = 0.5+0.5*math.erf(z/math.sqrt(2))
    return z, p, N, k


def calcWZP(df1, df2, selector=...):
    """
    Here, the uncertainties Isigma are used.
    Currently, it is not recommended to use this approach for limiting cases
    without further investigations of the methodology. (Work in progress.)
    """
    # N reflections in comparison
    if selector is ...:
        N = len(df1["DeltaI"])
    else:
        N = np.count_nonzero( selector )
    if N == 0: # bad selector
        return float("NaN"), float("NaN"), 0, 0
    DeltaIby2S = np.abs(df1["DeltaI"][selector] - df2["DeltaI"][selector])/(2*df1["Isigma"][selector])
    
    # w reflections that are likely to 
    w = np.sum([ 0.5-0.5*math.erf( DI2S/math.sqrt(2) ) for DI2S in DeltaIby2S ])
    
    # k reflections are better fitted with refinement of df1
    k = np.count_nonzero( np.abs(df1["DeltaI"][selector]) < np.abs(df2["DeltaI"][selector]) )
    
    # z-score corrected for w
    z = (k-N/2) / (math.sqrt(N-w)/2)
    
    # associated probability 
    p = 0.5+0.5*math.erf(z/math.sqrt(2))
    return z, p, N, k, w


def calcRall(df, selector=...):
    """ Conventional R factor
    selector: ... takes all data (no reflections filtered out)
    """
    if selector is not ... and np.count_nonzero( selector ) == 0: # bad selector
        print("Bad selector")
        return -1
    else: # R factor for selected subset
        return np.sum(np.abs(df["DeltaF"][selector]))/np.sum(np.abs(df["Fo"][selector]))

    
def read_m83(file):
    """ read Jana2006 or Jana2020 m83 file 
    returns a dictionary with
    h, k, l: integer Miller indices of reflection hkl
    Ic:      calculated intensity from refinement
    Io:      observed (measured) intensity from data reduction
    Isigma:  uncertainty of measured intensity Io from data reduction
    
    Fc, Fo, Fsigma are the corresponding values based on the square root of the intensity
    In the M83 file, these values are provided with limited precision. Therefore, they
    are calculated from Ic, Io, and Isigma.
    
    FrameID: ID of (virtual) frame used in the dynamical refinement
    BlockID: only relevant for multi-block refinements against multiple data sets
    """
    # read lines
    with open(file, "r") as fh:
        Lines = fh.readlines()
    if not Lines:
        return False
    
    # check if m83 file is from a dynamical refinement
    dynamical = False
    Line = Lines[2].split()
    if sum([ float(e) for e in Line[13:16]]) == 0: # check sum of A,B, which is only defined in kinematical refinement
        dynamical = True
        print("M83 file from dynamical refinement")
    elif Line[-1] == '0' and len(Lines[0]) in (150,165):
        dynamical = False
        print("M83 file from kinematical refinement")
    else:
        print("Unknown M83 format")
        print(Lines[0])
        print(Line)
        print("Line length (Python):", len(Lines[0]))
        
    # Prepare dictionary
    m83data = dict()
    for key in "h k l Ic Io Isigma Fc Fo Fsigma FrameID BlockID".split():
        m83data[key] = []
    
    # Read each line
    for Line in Lines:
        if Line.startswith("Block"):
            continue
        BlockID = 1
        if dynamical:
            FrameID = Line.split()[-1]
            if "%" in FrameID:
                FrameID, BlockID = FrameID.split("%")
        else:
            FrameID = -1
        m83data["BlockID"].append(int(BlockID))
        m83data["FrameID"].append(int(FrameID))

        # First three columns
        h,k,l = [ int(i) for i in Line[:12].split() ]
        m83data["h"].append(h)
        m83data["k"].append(k)
        m83data["l"].append(l)

        # Next three columns
        Ic = float(Line[14:27])
        Io = float(Line[29:42])
        Isigma = float(Line[45:57])
        m83data["Ic"].append(Ic) # Ic column is before Io column
        m83data["Io"].append(Io)  
        m83data["Isigma"].append(Isigma)

        # Calculate "structure factor amplitudes"
        if Io < 0.01*Isigma: # weak and negative intensities
            Fo = 0 if Io < 0 else math.sqrt(Io)
            Fsigma = 5*math.sqrt(Isigma) # this approach is how Jana treats weak reflections.
        else:
            Fo = math.sqrt(Io)
            Fsigma = Isigma/(2*Fo)
        m83data["Fc"].append(math.sqrt(Ic))
        m83data["Fo"].append(Fo)
        m83data["Fsigma"].append(Fsigma)
        
    return m83data
    
####################################################################################################################################

# USER INPUT / PROVIDED ARGUMENT

ask = True
if len(arguments) >  1:
    f1 = Path(arguments[1])
    
    if len(arguments) > 2:
        f2 = Path(arguments[2])
    else:
        f2 = f1.with_stem( f1.stem+"_INV" )
        
    if f1.exists() and f2.exists():
        ask = False
        
if ask:
    print("Please copy the path to the first M83 file")
    f1 = Path(input().strip('"'))
    print()
    print("Please copy the path to the second M83 file (leave empty if second file has label job_INV)")
    f2 = input().strip('"')
    print()
    if f2.strip() == '':
        f2 = f1.with_stem( f1.stem+"_INV" )
    else:
        f2 = Path(f2)

# READ FILES
   
# Read in data from M83 files
# Pandas data frame
print(f1.stem)
m83df = pd.DataFrame(data=read_m83(f1))
print(f2.stem)
m83dfINV = pd.DataFrame(data=read_m83(f2))

if np.all((m83df['h'] == m83dfINV['h']) & (m83df['Io'] == m83dfINV['Io'])):
    print("M83 files are compatible.")
else:
    print("M83 files are not compatible!")
    input("...")
print()    

# Define differences Io-Ic and Fo-Fc
# Iobs - Icalc
m83df["DeltaI"] = m83df["Io"] - m83df["Ic"]
m83dfINV["DeltaI"] = m83dfINV["Io"] - m83dfINV["Ic"]

# Fobs - Fcalc
m83df["DeltaF"] = m83df["Fo"] - m83df["Fc"]
m83dfINV["DeltaF"] = m83dfINV["Fo"] - m83dfINV["Fc"]


# CALCULATIONS AND OUTPUT

Blocks = set(m83df["BlockID"])
print("Block       N       k     N-k       z    p(#1)   Rall(#1)  Rall(#2)")

for Block in Blocks:
    Selector = (m83df["BlockID"]==Block)
    z, p, N, k = calcZP(m83df, m83dfINV, Selector)
    print(f"{Block:5d}{N:8d}{k:8d}{N-k:8d}{z:7.1f}σ{p*100:8.1f}% {calcRall(m83df,Selector):10.4f}{calcRall(m83dfINV,Selector):10.4f}")

# For multiblock data, that means if there are several data sets in the same refinement
if len(Blocks)>1:
    z, p, N, k = calcZP(m83df, m83dfINV)
    Rall = np.sum(np.abs(m83df["DeltaF"]))/np.sum(np.abs(m83df["Fo"]))
    RallINV = np.sum(np.abs(m83dfINV["DeltaF"]))/np.sum(np.abs(m83dfINV["Fo"]))
    print(f"comb.{N:8d}{k:8d}{N-k:8d}{z:7.1f}σ{p*100:8.1f}% {Rall:10.4f}{RallINV:10.4f}") # {p*100:10.6f}% 
    
input("... end of program ...")
