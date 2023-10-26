#!/usr/bin/env python
# coding: utf-8

"""
Remaining tasks:
- check weighting scheme / instability factor in M50 file
- multi-block refinement: block-based statistics
- optimise code for listing and merging sym. equiv. reflections
- use pyinputplus instead of input https://pyinputplus.readthedocs.io/en/latest/

- better description / comments
- installation / setup
"""

from math import sqrt
import numpy as np
import pandas as pd
from pathlib import Path

from sys import argv as arguments

__VERSION__ = "26 Oct 2023"

# [HKL, Ic, Io, Isigma, Fc, Fo, Fsigma, Fweight, wDF, Zone, Block]

# M83 format
# kinematical: M80factor = sqrt(1/scale)
#   h   k   l   Ic             Io             Isigma      obs   Twin  w(Fo-Fc) s*sqrt(Io) s*sqrt(Ic) 1/weight M80factor   sqrt(A^2+B^2) A         B 
#   0   1   0   0.259254E+02   0.289628E+02   0.466400E+00 o    1     4.197     21.00     19.87      0.27    0.65690E-01      5.09     -5.09      0.00
# dynamical:
#   h   k   l   Ic             Io             Isigma      obs   Twin  w(Fo-Fc) s*sqrt(Io) s*sqrt(Ic) 1/weight M80factor       0         0         0       Zone%Block
#   2   2  -7   0.112630E+03   0.299800E+03   0.203200E+03 <    1     1.142     17.31     10.61      5.87    0.20610E-01      0.00      0.00      0.00     1

# Read M50: get symmetry operations, calculate symmetry operations in reciprocal space
# Read M83: calculate wRall, blockwise wRall, symm.avg. wRall

def xyz2matrix(xyz):
    '''get matrix representation of symmetry operators'''
    if type(xyz) is str:
        e3 = xyz.split()
    elif len(xyz) == 3:
        e3 = xyz
    else:
        return False
    matrix = np.zeros([3,3],dtype=int)
    s = "xyz"
    for i in (0,1,2):
        for j in (0,1,2):
            if "-"+s[j] in e3[i]: # "-x" -> -1
                matrix[i,j] = -1
            elif s[j] in e3[i]:
                matrix[i,j] = 1
    return matrix

def read_m50(file):
    " Read M50 file, symmetry, UC "
    with open(file, "r") as fh:
        Lines = fh.readlines()
    if not Lines:
        return False
        
    m50 = dict()
    m50["File"] = file
    m50["Symmetry"] = []
    m50["SymmetryRotation"] = []
    m50["SymmetryTranslation"] = []
    m50["SymmetryCentering"] = []
    m50["SymmetryReciprocal"] = []
    m50["Refinement"] = ""    
    for i, line in enumerate(Lines):
        if line.startswith("title"):
            m50["Title"] = line[6:-1]
        elif line.startswith("cell"):
            m50["UC"] = tuple([ float(e) for e in line.split()[1:] ])
        elif line.startswith("esdcell"):
            m50["UCesd"] = tuple([ float(e) for e in line.split()[1:] ])            
        elif line.startswith("spgroup"):
            m50["SpaceGroup"], m50["SpaceGroupID"] = line.split()[1:3]
            m50["SpaceGroupID"] = int(m50["SpaceGroupID"])
        elif line.startswith("symmetry"):
            SymXYZ = line[8:]
            m50["Symmetry"].append( SymXYZ ) # tuple of strings, use eval + replace # e for e in line[8:].split())
            m50["SymmetryRotation"].append( xyz2matrix(SymXYZ) )
            m50["SymmetryReciprocal"].append( np.linalg.inv(m50["SymmetryRotation"][-1].transpose()).astype("int") )
        elif line.startswith("unitsnumb"):
            m50["Z"] = int(line.split()[1])
        elif line.startswith("atlist"):
            m50["Formula"] = line[7:-1]
        elif line.startswith("refine"):
            m50["Refinement"] = "".join(Lines[i+1:]).replace("\n &", " ").split("\n")
            if m50["Refinement"][-1] == "":
                m50["Refinement"].pop() # remove last element
            if m50["Refinement"][-1] == "end refine":
                m50["Refinement"].pop() # remove last element
            if "snlmx" in Lines[i+1]:
                m50["gmax"] = 2*float(Lines[i+1].split("snlmx")[1].split()[0])
        elif line.startswith("end refine"):
            break
    else:
        return False
    
    # CENTERING VECTORS
    CC = m50["SpaceGroup"][0] # centering
    m50["SymmetryCentering"] = [(0,0,0)]
    if CC != "P":
        add = 0.5
        a = (0, add, add)
        b = (add,0,add)
        c = (add,add,0)
        if CC == "F":
            m50["SymmetryCentering"] += [a,b,c]
        elif CC == "C":
            m50["SymmetryCentering"].append(c)            
        elif CC == "I":
            m50["SymmetryCentering"].append( (add,add,add) )
        elif CC == "A":
            m50["SymmetryCentering"].append(a)
        elif CC == "B":
            m50["SymmetryCentering"].append(b)
    #m50["Symmetry"] += [ tuple([x[i] + c[i] for i in (0,1,2)]) for c in CCvectors for x in m50["Symmetry"] ] # expand symmetry 
    return m50


def read_m83(file):
    """
    hkl h k l Icalc Iobs Isigma Fcalc Fobs Fsigma w(Fo-Fc), 1/weight~=sigma^2 Zone%Block
    """
    #m83df = pd.DataFrame(m83, columns="HKL, Ic, Io, Isigma, Fc, Fo, Fsigma, Fweight, wDF, Zone, Block".split(", "))
    with open(file, "r") as fh:
        Lines = fh.readlines()
    if not Lines:
        return False
    m83 = list()
    #m83["HKL"] = list()
    #m83["Zone"] = list()
    # Line length: Jana2006 KIN=150 DYN=156
    #              Jana2020 KIN=165 DYN=171
    BlockID = 1
    Dynamical = False
    if len(Lines[0]) in (150,165): # kinematical refinement
        Zone = 0
        print("Kinematical refinement")
    elif len(Lines[0]) in (154, 155, 156,157, 169, 171,172,173,174): # dynamical refinement
        print("Dynamical refinement")
        Dynamical = True
    else:
        print("Unknown M83 format.", len(Lines[0]))
        
    # Adapt code so that it first determines the instabilty factor.
    # Then use instability factor to get proper Fo Fc w wDF
        
    for Line in Lines:
        BlockID = 1
        if Dynamical:
            Zone = Line.split()[-1]
            if "%" in Zone:
                Zone, Block = Zone.split("%")
                BlockID = int(Block)
        
        h,k,l = [ int(i) for i in Line[:12].split() ]
        HKL = (h,k,l)
                
        Ic = float(Line[14:27]) # Ic column is before Io column
        Io = float(Line[29:42]) # MINUS sign not read until 02.11.2020
        Isigma = float(Line[45:57])
        
        wDF = float(Line[64:74])
        #Fo = float(Line[76:84])
        #Fc = float(Line[86:94])        
        Fweight = 1/float(Line[96:104])
        
        # kinematical case: read scale from M40 and include scale in calculation of Fo
        if Io < 0.01*Isigma: # weak and negative intensities
            Fo = 0 if Io < 0 else sqrt(Io)
            Fsigma = 5*sqrt(Isigma)
        else:
            Fo = sqrt(Io)
            Fsigma = Isigma/(2*Fo)
        Fc = np.sqrt(Ic)
        #Fweight = 1/sqrt(Fsigma**2 + (0.02*Fo)**2) 
        #wDF = Fweight*(Fo-Fc)
        m83.append( (HKL, Ic, Io, Isigma, Fc, Fo, Fsigma, Fweight, wDF, int(Zone), BlockID) )
    return m83



# Check provided arguments for Jana files
ask = True
if len(arguments) >  1: # and ".m83" in arguments[1]:
    f1 = Path(arguments[1].strip('"').strip("'"))
    f_m50 = f1.with_suffix(".m50")
    f_m83 = f1.with_suffix(".m83")
    
    print(f_m50)
    print(f_m83)
        
    if f_m50.exists() and f_m83.exists():
        ask = False
        
if ask:
    print("Please copy the path to the M83 file")
    f_m83 = Path(input().strip('"').strip("'")).with_suffix(".m83")
    f_m50 = f_m83.with_suffix(".m50")
    
    if not f_m50.exists():
        print(f_m50.stem, ' not found!')
        input('... end of program ...')
        



print("""
            _____             _      _                _ 
           |  __ \           (_)    | |              | |
 __      __| |__) | ___  ___  _   __| | _   _   __ _ | |
 \ \ /\ / /|  _  / / _ \/ __|| | / _` || | | | / _` || |
  \ V  V / | | \ \|  __/\__ \| || (_| || |_| || (_| || |
   \_/\_/  |_|  \_\\___| |___/|_| \__,_| \__,_| \__,_||_|""") #http://patorjk.com/software/taag/#p=display&h=1&f=Big&t=DynKin   
print("\n\tF I L E   C H E C K")
#print("Program called with argument:", f_jana)

m50 = read_m50(f_m50) # reciprocal space symmetry operators now in m50["SymmetryReciprocal"]
m83 = read_m83(f_m83) # [HKL, Ic, Io, Isigma, Fc, Fo, Fsigma, Fweight, wDF, Zone, Block]

# read_m83 determines Fweight = 1/np.sqrt(Fsigma**2 + (0.01*Fo)**2)

if not m50 or not m83:
    print("Error reading files...")
    input("... end or program")
    input("... end or program")
    input("... end or program")

m83df = pd.DataFrame(m83, columns="HKL, Ic, Io, Isigma, Fc, Fo, Fsigma, Fweight, wDF, Zone, Block".split(", "))

Blocks = set(m83df["Block"])

IobsLimit = 3 # Fobs > 2*Iobs

Fobs = []
Fcalc = []
Weights = []
Sigmas = []
UsedIDs = []

Iobs = []
Icalc = []
Isigmas = []

# Data containers, key = hkl, value 
FO = dict()
FC = dict()
FS = dict()
#FS2 = dict() # 1/sqrt(weigtht)
FW = dict()

IO = dict() # Iobs
IC = dict() # Icalc
IS = dict() # Isigma
#IW = dict() # Iweight

H = dict() # hkl

# hkl Fo Fc Fs Fw

# Prepare dict to find symmetrically equivalent reflections
# Prepare data containers
Equiv = dict()
for hkl in m83df["HKL"]:
    if hkl in Equiv:
        continue
    # prepare entry in dictionaries for hkl-group
    FO[hkl] = [] # Fobs
    FC[hkl] = [] # Fcalc
    FS[hkl] = [] # Fsigma
    FW[hkl] = [] # Fweight
    IO[hkl] = []
    IC[hkl] = []
    IS[hkl] = []
    H[hkl] = [] # hkl list
    for symop in m50["SymmetryReciprocal"]:
        sHKL = tuple([ int(e) for e in np.matmul(symop, hkl) ]) # includes identity operation
        if sHKL not in Equiv:
            Equiv[sHKL] = hkl
# Equiv: Symmetrically equivalent reflections (as key) have the same value hkl (as value).
#UniqueHKL = set(Equiv.values()) # unique in point group

# PRINT GROUPS
#for key, value in Equiv.items():
#    print(key, value)
    
#HKL, Ic, Io, Isigma, Fc, Fo, Fsigma, Fweight, wDF, Zone, Block
# Loop again. Group data together
                # 0    1     2     3       4         5     6      7
for r in m83df[["HKL","Fo","Fc","Fsigma","Fweight", "Ic", "Io", "Isigma"]].to_numpy():
    #print(f"Fc {r[2]:6.2f} Fo {r[1]:6.2f}  sigma {r[3]:4.1f}     ", ("{:3d}"*3).format(*r[0]), ">>>", ("{:3d}"*3).format(*Equiv[r[0]]))
    FO[Equiv[r[0]]].append(r[1]) # Equiv[r[0]] is the hkl index of the "mother" reflection of the same reflection group
    FC[Equiv[r[0]]].append(r[2])
    FS[Equiv[r[0]]].append(r[3]) # r[3] is Fsigma
    FW[Equiv[r[0]]].append(r[4]) # r[4] is Fweight
    
    IO[Equiv[r[0]]].append(r[6])
    IC[Equiv[r[0]]].append(r[5])
    IS[Equiv[r[0]]].append(r[7])
    
    H[Equiv[r[0]]].append(r[0])
    #FS2[Equiv[r[0]]].append(1/np.sqrt(r[3]))
    
NoDuplicates = sum([ len(set(g)) for g in H.values() ])

""" Correct: working with intensities """
# INTENSITIES: Process reflection groups --> arithmetic mean
IOM = np.array([ sum(g)/len(g) for g in IO.values() ]) # g = group, IO = <Iobs>
ICM = np.array([ sum(g)/len(g) for g in IC.values() ]) # IC = <Icalc>
ISM = np.array([ sqrt( sum([s*s for s in g]) )/len(g) for g in IS.values() ]) # IS = sigma(I) = sqrt(sigma1**2 + sigma2**2 + ... +sigma_i**2)/n

FCM = np.sqrt(ICM)
FSM = np.array(ISM)

IOMt = np.array(IOM) # temp...
IOMt[IOMt < 0 ] = 0
FOM = np.sqrt(IOMt)
FSM[IOM>=0.01*ISM] = ISM[IOM>=0.01*ISM]/(2*np.sqrt(IOM[IOM>=0.01*ISM])) # only correct for "strong" reflections
FSM[IOM<0.01*ISM] = 5*np.sqrt(ISM[IOM<0.01*ISM]) # JANA treatmeant of weak reflections

InstabilityFactor = 0.01

FWM = 1/np.sqrt(FSM*FSM + InstabilityFactor**2*FOM*FOM)

""" JANA treatmeant of weak reflections
if Io < 0.01*Isigma:
        Fsigma = 5*sqrt(Isigma)*Scale
        Fo = 0 if Io < 0 else sqrt(Io)*Scale
    else:
        Fsigma = Isigma/(2*sqrt(Io))*Scale
        Fo = sqrt(Io)*Scale
    Fc = sqrt(Ic)*Scale
    w = 1/sqrt(Fsigma**2 + (0.01*Fo)**2)
"""    

# Calculate mR-factors
MObsIDs = IOM > 3*ISM
MRobs = sum( abs(FOM[MObsIDs] - FCM[MObsIDs]) )/sum( abs(FOM[MObsIDs]) )
MRall = sum( abs(FOM - FCM) )/sum( abs(FOM) )
MwRall = sqrt( sum( (FWM*(FOM-FCM))**2 ) / sum( (FWM*FOM)**2 ) )

# Dynamical R-factors as calculated in JANA
ObsIDs = m83df["Fo"] > 6*m83df["Fsigma"]
Robs = sum(abs(m83df["Fo"][ObsIDs] - m83df["Fc"][ObsIDs]))/sum(abs(m83df["Fo"][ObsIDs]))
#wRobs = np.sqrt( np.sum( m83df["wDF"][ObsIDs]**2 ) / np.sum( (m83df["Fweight"][ObsIDs]*m83df["Fo"][ObsIDs])**2 ) )
Rall = sum( abs(m83df["Fo"] - m83df["Fc"]) )/sum( abs(m83df["Fo"]) )
wRall = np.sqrt( np.sum( m83df["wDF"]**2 ) / np.sum( (m83df["Fweight"]*m83df["Fo"])**2 ) )

print("Number of duplicates:", len(m83df["Fo"])-NoDuplicates)
print("All reflections without duplicates:", NoDuplicates)
print("'P1':  ", np.count_nonzero(ObsIDs),"/", len(m83df["Fo"]), " observed/all  reflections.")

print()
print(f"Merged reflections that are symmetrically equivalent in the space group type {m50['SpaceGroup']}:")
print("Merge: ", np.count_nonzero(MObsIDs),"/", len(FOM), " observed/all  reflections.")
#print(f"MRobs  {100*MRobs:6.2f}  Robs  {100*Robs:6.2f}")
#print(f"MRall  {100*MRall:6.2f}  Rall  {100*Rall:6.2f}")
#print(f"MwRall {100*MwRall:6.2f}  wRall {100*wRall:6.2f}")

print()

print(f"         obs       all  weighted")
print(f" R  {100*Robs:8.2f}  {100*Rall:8.2f}  {100*wRall:8.2f}")
print(f"MR  {100*MRobs:8.2f}  {100*MRall:8.2f}  {100*MwRall:8.2f}")

input("... end of program ... ")
