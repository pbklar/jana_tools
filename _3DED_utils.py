# -*- coding: utf-8 -*-
'''
Written by Paul Klar, University of Bremen
'''

__version__ = "13 August 2024"

from pathlib import Path
import numpy as np
import pandas as pd

ToRad = np.pi/180

# Jana
# ### J A N A
# 

def read_m40(file):
    ''' Read M40 file
    
    return list of tuples: Label Occ X Y Z ADPs
    
    Does NOT read in s.u.s
    
    Does not work with modulated structures
    '''
    with open(file, 'r') as fh:
        Lines = fh.readlines()
        
    if not Lines:
        return False
    NAtoms = int(Lines[0][:5])
    Atoms = []
    for i, Line in enumerate(Lines):
        if Line[0].isupper(): # atom line
            Label = Line[0:9].rstrip()
            Occ = float(Line[18:27])
            X = float(Line[27:36])
            Y = float(Line[36:45])
            Z = float(Line[45:54])
            ADPs = [ float(Lines[i+1][j*9:(j+1)*9]) for j in range(6) ]
            if sum(ADPs[1:]) == 0: 
                ADPs = [ ADPs[0] ]                      # LABEL occ x y z Uiso
            Atoms.append( (Label,Occ,X,Y,Z,*ADPs) ) # LABEL occ x y z u11 u22 u33 u12 u13 u23        
        elif len(Atoms) == NAtoms or "s.u. block" in Line or Line.startswith("--------------"):
            # one method should be enough...
            break
    return Atoms
    
def read_m42(file, read_frames=True, frames_as_pandas=True):
    ''' updated 13 August 2024
    read (almost) all parameters from M42 file.
    
    By default, a dictionary with the block settings and a pandas dataframe with frame parameters is returned.
    
    read_frames: True = read complete M42 file
                 False = only read the configuration block at the very top of the M42 file
                 
    frames_as_pandas: True = return a pandas data frame and a dictionary
                      False = return one dictionary with all settings
    '''
    
    with open(file, 'r') as fh:
        Lines = fh.readlines()    

    end_config_line = Lines.index('end\n')
    N_blocks = (end_config_line-2) // 7

    m42 = dict()
    m42["File"]    = file
    m42["commands"] = Lines[0].strip()
    m42["dynamical"] = True if Lines[1][8] == '1' else False       # calcdyn 1
    m42["scale_to_Fcalc"] = True if Lines[1][27] == '1' else False # scalefc 1

    m42["iedt"]    = []
    m42["tiltcorr"] = []
    m42["nzones"]  = []
    m42["intsteps"]= []
    m42["gmax_BW"]    = []
    m42["Sgmax_BW"]   = []
    m42["Sgmax"]   = []
    m42["RSgmax"]  = []
    m42["DSgmin"]  = []
    m42["OM"]  = []

    if read_frames:
        m42["BlockID"] = []
        m42["FrameID"] = []
        m42["UVW"]   = []
        m42["alpha"] = []
        m42["beta"]  = []
        m42["phi"]   = []   # rotation semiangle or precession angle
        m42["EDphi"]   = [] # optimised frame orientation correction direction
        m42["EDtheta"] = [] # optimised frame orientation correction tilt
        m42["use_flag"]    = []
        m42["thickness"]   = []
        m42["scale"]       = []

    # refinement data block settings
    BlockID = 1
    DSgmin = -1
    for i, line in enumerate(Lines[:end_config_line]):
        if line.startswith('Refblock'):
            BlockID = int( line.split('Block')[1] )
        elif line.startswith('threads'):
            data = line.split()
            tiltcorr = True if data[3]=='1' else False
            iedt = True if data[5]=='1' else False
        elif line.startswith('nzones'):
            data = line.split()
            nzones = int(data[1])
            intsteps = int(data[3])
        elif line.startswith('ormat'):
            OM = np.array( [ l.split() for l in Lines[i+1:i+4] ] ).astype(float)
        elif line.startswith('omega'):
            data = line.split()
            offset = 0
            if 'sca' in line:
                offset = 2
            gmax_BW = float(data[3])
            Sgmax_BW = float(data[5+offset])
            Sgmax = float(data[7+offset])
            RSgmax = float(data[9+offset])
            if 'dsgmin' in line:
                DSgmin = float(data[11+offset])

            m42["iedt"].append(iedt)
            m42["tiltcorr"].append(tiltcorr)
            m42["nzones"].append(nzones)
            m42["intsteps"].append(intsteps)
            m42["gmax_BW"].append(gmax_BW)
            m42["Sgmax_BW"].append(Sgmax_BW)
            m42["Sgmax"].append(Sgmax)
            m42["RSgmax"].append(RSgmax)
            m42["DSgmin"].append(DSgmin)
            m42["OM"].append(OM)

    if not read_frames:
        return m42

    # frame parameters
    before = 0
    BlockID = 1
    for i, line in enumerate(Lines): 
        if line.startswith('# Zone'):
            frameID = int(line[6:].strip().split()[0])

            UVW = [ float(Lines[i+1][k*9:k*9+9]) for k in (0,1,2)]
            alpha = float(Lines[i+1][27:36])
            beta = float(Lines[i+1][36:45])
            phi = float(Lines[i+1][45:54])
            UseFlag = True if Lines[i+1][60] == "T" else False

            scale = float(Lines[i+2][:9])
            thickness = float(Lines[i+2][9:18])
            EDphi = float(Lines[i+3][:9])
            EDtheta = float(Lines[i+3][9:18])

            if frameID < before:
                BlockID += 1

            m42["BlockID"].append( BlockID )
            m42["FrameID"].append( frameID )
            m42["UVW"].append(UVW)
            m42["alpha"].append(alpha)
            m42["beta"].append(beta)
            m42["phi"].append(phi)
            m42["use_flag"].append(UseFlag)
            m42["thickness"].append(thickness)
            m42["scale"].append(scale)

            m42["EDphi"].append(EDphi)
            m42["EDtheta"].append(EDtheta)

            before = frameID
        elif line.startswith('---------------'): # s.u. block
            break
            
    if frames_as_pandas:
        keys = 'use_flag BlockID FrameID alpha phi thickness scale beta EDphi EDtheta UVW'.split()
        df = pd.DataFrame(data=m42, columns=keys)
        for key in keys:
            del(m42[key])
        return df, m42
    else:
        return m42
        
       
        
def read_m50(file):
    " Read M50 file, symmetry, UC "
    with open(file, 'r') as fh:
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
    m50["gmax"] = 9
    m50["Refinement"] = ""
    m50["skipbad"] = False
    m50["skipflag42"] = False
    for i, line in enumerate(Lines):
        if line.startswith("title"):
            m50["Title"] = line[6:-1]
        elif line.startswith("cell"):
            m50["UC"] = tuple([ float(e) for e in line.split()[1:] ])
            m50["V"] = np.sqrt( np.linalg.det( metric_tensor(*m50["UC"]) ) )
        elif line.startswith("esdcell"):
            m50["UCesd"] = tuple([ float(e) for e in line.split()[1:] ])            
        elif line.startswith("spgroup"):
            m50["SpaceGroup"], m50["SpaceGroupID"] = line.split()[1:3]
            m50["SpaceGroupID"] = int(m50["SpaceGroupID"])
        elif line.startswith("symmetry"):
            SymXYZ = line[8:].strip()
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
            m50["skipbad"] = "skipbad 1" in Lines[i+1]
            m50["skipflag42"] = "  skipflag 42" in m50["Refinement"]
            if "snlmx" in Lines[i+1]:
                m50["gmax"] = 2*float(Lines[i+1].split("snlmx")[1].split()[0])
        elif line.startswith("end refine"):
            break
    else:
        return False
    
    # CENTERING VECTORS
    # R centering missing
    
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
        elif CC == "A":
            m50["SymmetryCentering"].append(a)
        elif CC == "B":
            m50["SymmetryCentering"].append(b)
        elif CC == "I":
            m50["SymmetryCentering"].append( (add,add,add) )
    #m50["Symmetry"] += [ tuple([x[i] + c[i] for i in (0,1,2)]) for c in CCvectors for x in m50["Symmetry"] ] # expand symmetry 
    return m50
  
    
def read_m83(file, return_pandas=True):
    """
    hkl h k l Icalc Iobs Isigma Fcalc Fobs Fsigma w(Fo-Fc), 1/weight~=sigma^2 Zone%Block
    
    Example M83 line:
    h   k   l   Ic             Io             Isigma      obs   Twin  w(Fo-Fc) s*sqrt(Io) s*sqrt(Ic) 1/weight M80factor       0         0         0       Zone%Block
    1  -4  -5   0.464230E+03   0.402300E+03   0.168300E+03 <    1    -0.354     20.06     21.55      4.20    0.30474E-01      0.00      0.00      0.00     1
   
    For kinematical refinement, proper use of scale factors not implemented yet.
    """    
    with open(file, 'r') as fh:
        Lines = fh.readlines()

    if not Lines:
        return False

    m83 = list()
    # Line length: Jana2006 KIN=150 DYN=156
    #              Jana2020 KIN=165 DYN=171
    BlockID = 1
    Dynamical = False
    if len(Lines[0]) in (150,165): # kinematical refinement
        Zone = 0
        print("Kinematical refinement!")
    elif len(Lines[0]) in (13,154, 155, 156,157, 169, 171,172,173,174): # dynamical refinement
        #print("Dynamical refinement")
        Dynamical = True
    else:
        print("Unknown M83 format.", len(Lines[0]))

    # Adapt code so that it first determines the instabilty factor.
    # Then use instability factor to get proper Fo Fc w wDF
    Zone = 1
    for Line in Lines:
        BlockID = 1
        if Line[0] == 'B': # BlockN begin / end
            continue
        elif Dynamical:
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

        # kinematical case: read scale from M40 and include in calculation of Fo
        if Io < 0.01*Isigma: # weak and negative intensities
            Fo = 0 if Io < 0 else Io**0.5
            Fsigma = 5*Isigma**0.5
        else:
            Fo = Io ** 0.5
            Fsigma = Isigma/(2*Fo)
        Fc = Ic ** 0.5
        #Fweight = 1/sqrt(Fsigma**2 + (0.02*Fo)**2) 
        #wDF = Fweight*(Fo-Fc)
        m83.append( (HKL, Ic, Io, Isigma, Fc, Fo, Fsigma, Fweight, wDF, int(Zone), BlockID) )

    if return_pandas:
        return pd.DataFrame(m83, columns="HKL Ic Io Isigma Fc Fo Fsigma Fweight wDF Zone Block".split())
    else:
        # list of tuples
        return m83


def merge_reflections_post_refinement(m83df, symmetry_operations_RS):    
    '''
    merge symmetrically equivalent reflections
    return pandas data frame with merged reflections, respective intensities and amplitudes
    
    symmetry_operations_RS is a list of 3x3 matrices with reciprocal space symmetry operations,
    i.e., (R^T)^(-1), where R is the rotational part of a space group operation
    
    symmetry_operations_RS is determined during read_m50()
    
    Example:
    f_m50 = Path(file_m50)
    f_m83 = f_m50.with_suffix(".m83")
    m50 = read_m50(f_m50)
    m83df = read_m83(f_m83)
    m83_merged = merge_reflections_post_refinement(m83df, m50["SymmetryReciprocal"])
    '''
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
        for symop in symmetry_operations_RS:
            sHKL = tuple([ int(e) for e in np.matmul(symop, hkl) ]) # includes identity operation
            if sHKL not in Equiv:
                Equiv[sHKL] = hkl
    # Equiv: Symmetrically equivalent reflections (as key) have the same value hkl (as value).
    #UniqueHKL = set(Equiv.values())

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

    #NoDuplicates = sum([ len(set(g)) for g in H.values() ])

    """ Correct: working with intensities """
    # INTENSITIES: Process reflection groups --> arithmetic mean
    IOM = np.array([ sum(g)/len(g) for g in IO.values() ]) # g = group, IO = <Iobs>
    ICM = np.array([ sum(g)/len(g) for g in IC.values() ]) # IC = <Icalc>
    ISM = np.array([ sum([s*s for s in g])**0.5/len(g) for g in IS.values() ]) # IS = sigma(I) = sqrt(sigma1**2 + sigma2**2 + ... +sigma_i**2)/n
    # 07.10.: len(g) moved outside the sqrt().

    FCM = np.sqrt(ICM)
    FSM = np.array(ISM)

    IOMt = np.array(IOM) # temp...
    IOMt[IOMt < 0 ] = 0
    FOM = np.sqrt(IOMt)
    FSM[IOM>=0.01*ISM] = ISM[IOM>=0.01*ISM]/(2*np.sqrt(IOM[IOM>=0.01*ISM])) # only correct for "strong" reflections
    FSM[IOM<0.01*ISM] = 5*np.sqrt(ISM[IOM<0.01*ISM]) # JANA treatmeant of weak reflections

    FWM = 1/np.sqrt(FSM*FSM + 0.0001*FOM*FOM)

    merged_reflections = pd.DataFrame( [H.keys(), ICM, IOM, ISM, FCM, FOM, FSM, FWM, IOM > 3*ISM] ).transpose()
    merged_reflections.columns = "HKL Ic Io Isigma Fc Fo Fsigma Fweight obs_flag".split()
    return merged_reflections        
    


def calc_R_factors(reflections, observed_sigma_factor=3, instability_factor=0.01): # R_in_percent=True # R_in_percent:           R factors are returned in % 
    """ calculate R factors based on structure factor amplitudes 
    reflections must be an object like a dictionary or a pandas DataFrame with columns
    
    Io:   measured intensity
    Ic:   calculated intensity
    Isigma: standard uncertainty of measured intensity
    Fc:   calculated amplitude
    Fo    observed amplitude
    Fsigma: standard uncertainty of observed amplitude
    
    Optional:
    Fweight: weighting factor for Fo
    
    WARNING:
    For kinematical refinement, an output might be calculated without error message.
    It is very likely, that the output will be wrong! (EXTI, scale factors)
    
    Parameters:
    observed_sigma_factor:  Reflections with Iobs > observed_sigma_factor * Isigma are considered 'observed'
    instability_factor:     if provided, weighting scheme is used with instability factor
                            if is None and Fweight provided, Fweight is used
                            otherwise, wR factors are not calculated
                            
    Weightigh scheme:
    see Jana98 manual, p. 223, equations E57 and E59
    """
       
    if 'Io' and 'Isigma' in reflections:
        observed_flag = ( reflections['Io'] > observed_sigma_factor*reflections['Isigma'] )
        Nobs = np.count_nonzero(observed_flag)
        
        Robs = ( np.sum( np.abs(reflections['Fo'][observed_flag] - reflections['Fc'][observed_flag]) ) / 
                 np.sum( np.abs(reflections['Fo'][observed_flag]) ) )
    else:
        #    observed_flag = ( reflections['Io'] > observed_sigma_factor*reflections['Isigma'] )
        Robs = np.nan
        Nobs = np.nan
    
    Rall = ( np.sum( np.abs(reflections['Fo'] - reflections['Fc']) ) / 
             np.sum( np.abs(reflections['Fo']) ) )
    
    # weighted R factors
    Fweight = False
        
    if instability_factor:
        Fweight = 1 / (   reflections['Fsigma']**2 + (instability_factor*reflections['Fo'])**2   ) ** 0.5
        Iweight = 1 / ( ( reflections['Fsigma']**2 + (instability_factor*reflections['Fo'])**2 ) * np.abs(4*reflections['Io']) ) ** 0.5        
    elif 'Fweight' in reflections:
        Fweight = reflections['Fweight']
        Iweight = ( reflections['Fweight'] ** 2 / (4*reflections['Io']) ) ** 0.5
        
    if np.any(Fweight):
        # replace infinity
        Iweight[ Iweight == np.inf ] = 0
        
        wRobs =  ( np.sum( ( Fweight[observed_flag]*(reflections['Fo'][observed_flag]-reflections['Fc'][observed_flag]) )**2 ) / 
                   np.sum( ( Fweight[observed_flag]* reflections['Fo'][observed_flag] )**2 )   )**0.5
        
        wRall =  ( np.sum( ( Fweight*(reflections['Fo']-reflections['Fc']) )**2 ) / 
                   np.sum( ( Fweight* reflections['Fo'] )**2 )   )**0.5
                   
        wR2obs = ( np.sum( ( Iweight[observed_flag]*(reflections['Io'][observed_flag]-reflections['Ic'][observed_flag]) )**2 ) / 
                   np.sum( ( Iweight[observed_flag]* reflections['Io'][observed_flag] )**2 )   )**0.5
        
        wR2all = ( np.sum( ( Iweight*(reflections['Io']-reflections['Ic']) )**2 ) / 
                   np.sum( ( Iweight* reflections['Io'] )**2 )   )**0.5
    else:
        wRobs = np.nan
        wRall = np.nan
        wR2all = np.nan
        wR2obs = np.nan
    
    return {'Robs': Robs,
            'Rall': Rall,
            'wRall': wRall,
            'wR2all': wR2all,
            'wRobs': wRobs,
            'wR2obs': wR2obs,
            'Nobs': Nobs,
            'Nall': len(reflections)}

    #GoF = ( np.sum( ( reflections['Fweight']**2 * (reflections['Io']-reflections['Ic']) )**2 ) /
    #        (M_reflections - N_parameters) ) ** 0.5
    
    # JANA98 manual: p. 223, E57, E59
    # Fweight = 1 / ( reflections['Fsigma']**2 + (instability_factor*reflections['Fo'])**2 )
    # Iweight = Fweight / (4*reflections['Fo']**2)   
    

def generate_equivalent_reflections(reflections, symmetry_operations_RS):
    '''
    1. apply symmetry operations to list of reflections
    2. general full cuboid of reflections
    
    Useful to determine the completeness.
    
    Systematic absences NOT considered. Separate routine needed, but not written yet (see notes in Jupyter notebook)
    '''
    
    # make sure that h_exp has the required format
    h_exp = np.array( [ list(h) for h in list(reflections) ] )
    
    # Theoretical set of reflections
    h_max = np.max( np.abs( h_exp ) ) + 5 # is this bullet proof? better work with gmax?
    h_range = range(-h_max,h_max+1)
    h_all = np.array( [ (h,k,l) for h in h_range for k in h_range for l in h_range ] )

    # Generate all symmetrically equivalent reflections from input list of reflections
    h_exp_point_group_symmetry = set()

    for hkl in reflections:
        # apply symmetry
        for symop in symmetry_operations_RS:
            sHKL = tuple([ int(e) for e in np.matmul(symop, hkl) ]) # includes identity operation
            h_exp_point_group_symmetry.add( sHKL )
            
    [ list(h) for h in h_exp_point_group_symmetry ]
            
    return [ list(h) for h in h_exp_point_group_symmetry ], [ list(h) for h in h_all ]

# Crystal geometry
def metric_tensor(a,b,c,alpha=90,beta=90,gamma=90): # 
    '''
    return the metric tensor based on 6 unit cell parameters.
    
    Example:
    GG = metric_tensor(*UC)
    GGr = np.linalg.inv(GG)
    d2 = HKL.dot(GGr).dot(HKL) # dstar^2
    d = 1 / np.sqrt( d2 )
    '''
    alpha = alpha*ToRad
    beta = beta*ToRad
    gamma = gamma*ToRad
    return np.array([[ a*a,               a*b*np.cos(gamma), a*c*np.cos(beta)  ],
                     [ a*b*np.cos(gamma), b*b,               b*c*np.cos(alpha) ],
                     [ a*c*np.cos(beta),  b*c*np.cos(alpha), c*c               ]])

def general_orientation_matrix(a,b,c,alpha=90,beta=90,gamma=90):
    '''
    generates one orientation matrix with
    a parallel to x,
    b perpendicular to z, and
    c completes the right-handed basis
    '''    
    V = np.linalg.det( metric_tensor(a,b,c,alpha,beta,gamma) )**0.5
    
    alpha = alpha*ToRad
    beta = beta*ToRad
    gamma = gamma*ToRad
    return np.linalg.inv( np.array([
        [ a,               0,                  0 ],
        [ b*np.cos(gamma), b*np.sin(gamma),    0 ],
        [ c*np.cos(beta), -c*(np.cos(beta)*np.cos(gamma) -np.cos(alpha))/np.sin(gamma), V/(a*b*np.sin(gamma)) ] ]) )
        
        
def read_OM(OMstring):
    " convert a string of 9 numbers to a 3x3 rotation matrix "
    OM = [ float(e) for e in OMstring.split() ]
    if len(OM) != 9:
        return False
    #return ([ tuple([OM[i+j] for j in (0,1,2)]) for i in (0,3,6) ])        
    return np.array( OM ).reshape( (3,3) )


# I/O read write tools

def split_value_su(string):
    if "(" in string:
        value, su = string.split("(")
        value = value.strip()
        TenTo = 0
        if "." in value:
            TenTo = len(value) - value.index(".") - 1
        su = float(su.replace(")","")) * 10**-TenTo
    else:
        su = 0
        value = string
    return (float(value), su)

# file tols

def argument(required="m50"):
    " RETURN argument that fits desired extension"
    if extension.endswith("m*"):
        extension = ".m50"
    if extension[0] != ".":
        extension = "."+extension
    for arg in sys.argv: # check all arguments
        # direct match: desired extension in list of arguments
        if arg.endswith(extension):
            p = Path(arg)
            if p.exists():
                return p
        # indirect match: desired extension is a JANA file, but different JANA file provided
        elif arg[-3:] in "m40 m42 m50 m70 m80 m81 m83 m85 m90 m95 ref".split():
            p = Path(arg)
            p = p.with_suffix(extension)
            if p.exists():
                return p    
    return input_filepath()
    
def input_filepath():
    print("Please enter the path. You can use 'Copy path' from Windows Explorer.")
    file = input()
    file = Path(file.strip("\""))
    if file.exists():
        return file
    else:
        input_filepath()  
        
def xyz2matrix(xyz):
    '''
    convert symmetry operation in 'x,y,z' form into 3x3 matrix
    '''
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


# Jana M42 utils

def m42_profile(alpha0,beta0=0.0,phi0=1.0,scale=1000,thickness=400,EDphi=0,EDtheta=0, step=0.01):
    '''
    prepare M42 file to get simulated rocking curve profile, intstep = 1
    
    Example:
    Get rocking curves for alpha in the range 30.0-1.4 = 28.6 deg to 30.0+1.4 = 31.4 deg,
    beta 0.2, thickness set to 1000 Ang in 0.01 deg steps:
    m42_profile(alpha0=30.0, beta0=0.2, phi0=1.4, thickness=1000, step=0.01)     
    
    output given on screen ...
    '''
    alpha = np.linspace(alpha0-phi0, alpha0+phi0, int(2*phi0/step))
    for i, a in enumerate(alpha):
        print(f"# Zone {i+1:d}")
        print(Jana9(1),Jana9(1),Jana9(1),Jana9(a),Jana9(beta0), Jana9(0),"      T   1",sep="") # u v w alpha beta phi
        print(Jana9(scale),Jana9(thickness),Jana9(0),Jana9(0),"                        1000",sep="")
        print(Jana9(EDphi),Jana9(EDtheta),"                                          00",sep="")


def m42_precession_profile(alpha0,beta0=0.0,phi0=1.0,scale=1000,thickness=400,EDphi=0,EDtheta=0, step=1.0): # 3.6 = 100 steps, 1.0 = 360 steps, 0.36 = 1000
    ''' similar to m42_profile
    Here, phi0 is interpreted as precession angle. The orientation is set by adapting EDphi and EDtheta.
    Input EDphi and EDtheta define the frame orientation without beam precession.
    Output EDphi and EDtheta describe the frame orientation + precession movement.
    '''
    EDphi_range = np.arange(0+EDphi, 360+EDphi, step)
    for i, p in enumerate(EDphi):
        print(f"# Zone {i+1:d}")
        print(Jana9(1),Jana9(1),Jana9(1),Jana9(alpha0),Jana9(beta0), Jana9(0),"      T   1",sep="") # u v w alpha beta phi
        print(Jana9(scale),Jana9(thickness),Jana9(0),Jana9(0),"                        1000",sep="")
        print(Jana9(p),Jana9(phi0+EDtheta),"                                          00",sep="")
#m42_precession_profile(alpha0=-49.891, beta0=-0.065, phi0=0.92, scale=502.2951, thickness=507.0736, step=1.0)        


# chemistry

ELEMENTS = """                                                                             H He 
Li Be                                                                          B  C  N  O  F Ne
Na Mg                                                                         Al Si  P  S Cl Ar
 K Ca Sc                                           Ti  V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
Rb Sr  Y                                           Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te  I Xe
Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta  W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
Fr Ra Ac Th Pa  U Np Pu"""
ELEMENTS = { i+1:E for i,E in enumerate(ELEMENTS.split()) }




# REQUIRES OPTIMISATION / CAREFUL CHECK:
def calc_sg(wavevector, xyz):
    """ Calculate excitation error Sg, which is the distance beteen the
    ideal Ewald sphere defined by the wavevector and the coordinates xyz.
    Sg is defined in Equation (7) in Palatinus et al (2015), Acta Cryst A71, 235-244, "Structure refinement using precession electron diffraction tomography and dynamical diffraction: theory and implementation"
    doi:10.1107/S2053273315001266
    """ 
    K = wavevector
    g = xyz
    KplusG = [ K[i]+g[i] for i in range(3) ]
    Sg = ( length(K)**2 - length(KplusG)**2 ) / (2*length(K))
    return Sg


def calc_rsg(wavevector, xyz, phi, geometry="continuous", absolute=False):
    """ RSg <= 1: reflection passes through exact Bragg condition during frame recording
        RSg > 1: exact Bragg condition not reached during precession movement
    RSg for precession is defined in Section 2.4(d) in Palatinus et al (2015), Acta Cryst A71, 235-244, "Structure refinement using precession electron diffraction tomography and dynamical diffraction: theory and implementation"
    RSg for continuous rotation method is similar, but |g| must be replaced by the distance to the projected rotation axis of the goniometer.
    """
    if phi == 0 or sum(xyz) == 0 or (geometry=="continuous" and sum(xyz[1:3])==0):
        return float("NaN")
    semiangle = phi*ToRad # Precession: precession semiangle. Continuous rotation: rotation semiangle
    K = wavevector
    g = xyz
    Sg = calc_sg(K,g)
    #print(Sg, g, semiangle)
    if geometry == "precession":
        g = xyz
    elif geometry == "continuous":
        g = (0, xyz[1], xyz[2])
    else:
        return False
    RSg = Sg / ( sqrt(sum([e**2 for e in g]))*semiangle )
    if absolute:
        return abs(RSg)
    else:
        return RSg
    
def calc_dsg(wavevector, xyz, phi, geometry="continuous"):
    """ DeltaSg = DSG
    DSg is the distance of a reflection from the closest limiting Ewald sphere.
    Calculation: 
        a) Calculate Sg of reflection AFTER rotating the crystal by an additional +- phi. The smaller value is the DSg.
    or  b) Calculate Sg of reflection. DSg = g*phi - abs(Sg)
    """ # Let's call it "Delta SG" = DSg
    semiangle = phi*ToRad
    K = wavevector
    g = xyz
    Sg = calc_sg(K, g)
    if phi == 0:
        return -abs(Sg)
    if geometry == "precession":
        g = xyz
    elif geometry == "continuous":
        g = (0, xyz[1], xyz[2]) # project to yz plane, g is now the shortest line between xyz and the rotation axis.
    else:
        return False
    return sqrt(sum([e**2 for e in xyz]))*semiangle - abs(Sg) # DSg
    
    
def calc_distance_to_x_axis(xyz):
    """ Distance to rotation axis. z || electron beam, x || goniometer tilt axis (alpha rotation), y || second tilt axis (beta rotation)
    The d2axe is the distance to the "shadow" (projection) of the rotation axis on the diffraction pattern defined in the Cartesian coordinate system by the axis (0,y,z).
    The unit is in reciprocal Angstrom.
    """
    return length(xyz[1:3])

def calc_azimuth(xyz):
    """ returns the azimuth in degrees
    x    y    z    -180 ... -90 ... 0 ... 90 ... 180 azimuth angle
    +    +                             X                          
    -    +                                    X
    +    -                       X
    -    -               X    
    
    Azimuth independent of z
    atan2(azimuth) = sgn(y)*d2x/x 
    """
    #print(asin( calc_distance_to_x_axis(xyz) / length(xyz) )*180/np.pi )
    sgn = abs(xyz[1])/xyz[1] if xyz[1] != 0 else 1 # sgn(y) either 1 or -1
    return atan2(sgn*calc_distance_to_x_axis(xyz), xyz[0])*180/np.pi       
    
    
    
    
    
def read_dyntmp(file):
    dyntmp = []
    Lines = readlines(file)
    for line in Lines:       
        h,k,l = int(line[0:4]), int(line[4:8]), int(line[8:12])
        Io = float(line[12:28])
        Isigma = float(line[28:44])
        ZoneID = int(line[84:89])
        g = float(line[94:104]) 
        Sg = float(line[104:114])
        RSg = float(line[114:124]) # |RSg|
        if Isigma == 0:
            Isigma = float("NaN")
        dyntmp.append( ((h,k,l), Io, Isigma, int(ZoneID), g, Sg, RSg) )
    return dyntmp    
    
    


def read_cif_pets(file):
    " Read cif_pets file, RC width, mosaicity, UC "
    cif_pets = dict()
    Lines = readlines(file, 2000)
    if not Lines:
        return False
    
    for line in Lines:
        if line.startswith("_cell_length_a"): a = float(line.split()[1])
        elif line.startswith("_cell_length_b"): b = float(line.split()[1])
        elif line.startswith("_cell_length_c"): c = float(line.split()[1])
        elif line.startswith("_cell_angle_alpha"): alpha = float(line.split()[1])
        elif line.startswith("_cell_angle_beta"): beta = float(line.split()[1])
        elif line.startswith("_cell_angle_gamma"): gamma = float(line.split()[1])
        elif line.startswith("_cell_volume"): cif_pets["V"] = float(line.split()[1])
        elif line.startswith("dstarmax"): cif_pets["gmax"] = float(line.split()[-1])
        elif line.startswith("RC width"): cif_pets["RC"] = float(line.split()[-1])
        elif line.startswith("mosaicity"): cif_pets["mosaicity"] = float(line.split()[-1])
    cif_pets["UC"] = (a,b,c,alpha,beta,gamma)
    g = 0.2
    d2x = 0.2
    FWHM = cif_pets["RC"] + np.sqrt(cif_pets["RC"]**2 + 8*np.log(2)*(cif_pets["mosaicity"]*np.pi/180*g)**2)
    Lifetime = FWHM*6/(2*np.sqrt(2*np.log(2))) # in reciprocal angstrom
    LifetimeDeg = np.arctan(Lifetime/d2x)*180/np.pi
    
    cif_pets["PhiEstimate"] = LifetimeDeg/2
    
    return cif_pets