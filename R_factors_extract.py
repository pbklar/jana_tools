# -*- coding: utf-8 -*-
'''
Written by Paul Klar, University of Bremen

Search for *.ref files from JANA refinement jobs and list R factors.

Useful tool to compare R factors of e.g. several refinement steps and
get quickly the R factor of the latest refinement in a folder.

Example:
- 
'''

from pathlib import Path
from sys import exit, argv as arguments

#__version__ = "17 Dec 2019"
#__version__ = "26 Jul 2021"
__version__ = "20 August 2024"

JANA_ref_files = []

print("- - - - - - - - - - ")
print(" J A N A  R factors:")
print("- - - - - - - - - - ")
print("This program reads in JANA refinement *.ref files and gives an overview on selected parameters.")
print()

for folder in arguments:
    folder = Path(folder)
    if folder.is_dir():
        ref_files = list(folder.glob('*.ref'))
        N_files = len(ref_files)
        if N_files > 0:
            JANA_ref_files += ref_files
            print(N_files, 'JANA refinements found in', folder)
            
if not JANA_ref_files:
    print("Current working directory:")
    folder = Path(".")
    print(folder.resolve())
    print("User input required!")
    print("- If you just hit ENTER (empty input), all subfolders of the working directory are searched for JANA refinement files.")
    print("- If you input the path to a folder with JANA refinement file(s), only refinements of that folder are read.")
    print('PATH:')
    choice = input()
    
    if choice.strip() == '':
        # Recursive search for ref files
        JANA_ref_files = list(folder.glob('**/*.ref'))
    else:
        folder = Path(choice)
        if folder.is_dir():
            JANA_ref_files = list(folder.glob('*.ref'))
    
    print(len(JANA_ref_files), 'JANA refinements found.')

# PREPARE
empty = []
this_folder = False

# this variable will be printed at the end
output = ''

# this variable collects information from each relevant folder
output_tmp = ''

# READ
for file in JANA_ref_files:
    multiblock = False
    Nobs = 0
    
    if file.parent != this_folder:
        if output_tmp.count("\n") > 2:
            # if last folder contained useful ref files, use the generated output
            output += output_tmp + '\n'
            
        this_folder = file.parent
        print('       ', '.' * len(str(file.parent)))
        print('FOLDER:', file.parent)
        print('       ', '.' * len(str(file.parent)), '\n')
        output_tmp = ' > > ' + str(file.parent) + '\\ < < \n'
        output_tmp += "File                              Parameters    Nobs    Nall  GOFobs  GOFall    Robs    Rall   wRall\n"
    
    print('READING:', file.name)
    with open(file, "r") as fh:
        Lines = fh.readlines() 
    
    for i, line in enumerate(Lines):
        if line.startswith("Last screen information window"):
            counter = 2
            j = 0
            while counter:
                LL = Lines[i+j]
                if LL[0] == "-":
                    counter -= 1
                    print(LL, end='')
                elif LL[0] == '|':
                    print(LL, end='')
                    if LL.startswith('|R factors'):
                        # SINGLE BLOCK
                        N_data = LL
                    elif LL.startswith('|Overall R factors'):
                        # MULTIBLOCK
                        multiblock = True
                        N_data = LL                        
                    elif LL.startswith('|GOF(obs)='):
                        GoF_data = LL
                    elif LL.startswith('|R(obs)='):
                        R_data = LL
                elif LL.startswith('***'): # avoid infinite loop
                    break
                j += 1

            # Extract Nobs, Nall, N_params
            start = N_data.index('[')+1
            end = N_data.index(']')
            data = N_data[start:end]
            Nall, data = data.split('=')
            Nobs, data = data.split('+')
            _, Nparameters = data.split('/')
            
            # Goodness of fit
            _, GoFobs, _, GoFall, _ = GoF_data.split()
                        
            # R factors
            _, Robs, _, wRobs, _, Rall, _, wRall, _ = R_data.split()
            
            output_tmp += f"{str(file.stem)[-40:]:<36}{Nparameters:>8}{Nobs:>8}{Nall:>8}{GoFobs:>8}{GoFall:>8}{Robs:>8}{Rall:>8}{wRall:>8}\n" 
            
            if not multiblock:
                break
                
        elif line.startswith('* R-factors overview *'):
            # only needed if multiple blocks are present in refinement
            BlockID = 0
            for line in Lines[i+2:]:
                data = line.split()
                if line.startswith('Block'):
                    BlockID = line.split('->')[0]
                elif line.startswith('#') and ":" in line:
                    BlockID = "Block" + line[1:].split(":")[0]
                elif len(data) == 17:
                    # Extract parameters for individual data blocks
                    cycle, Robs, wRobs, Rall, wRall, Nobs, Nall, Nparameters, ratio, Nskip, _, damp, GoFobs, GoFall, change_over_su_avg, change_over_su_max, _ = data
                elif len(data) == 9:
                    cycle, Robs, wRobs, Rall, wRall, Nobs, Nall, Nparameters, ratio = data
                    output_tmp += f"|{BlockID:->35}{' ':>8}{Nobs:>8}{Nall:>8}{' ':>8}{' ':>8}{Robs:>8}{Rall:>8}{wRall:>8}\n" 
                elif line.startswith('***'):
                    break 
            output_tmp += '\n'
    if not Nobs:
        empty.append(file)
        
    print()

for file in empty:
    print('INCOMPLETE REFINEMENT OUTPUT: ', file)
print()    
    
print('O V E R V I E W :')
print('"""""""""""""""""')
print(output)
if output_tmp.count("\n") > 2:
    print(output_tmp)
    
print()
input(' ... end of program ...')