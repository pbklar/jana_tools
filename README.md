# jana_tools
Tools related to 3D ED, dynamical refinement and especially reading/analysing files written by Jana2020 (or Jana2006)

### R_factors_extract.py
> Provides a covenient overview over all JANA refinement outputs (refinements against single crystal.
> Input: Folder with **complete** REF file(s) written after the last refinement cycle.

**Example 1**: python R_factors_extract.py
> Input the path to a folder with JANA refinement files. If no input is provided, the current working directory and all subdirectories are searched.

**Example 2**: python R_factors_extract.py C:/refinements
> Looks for refinement files in the given folder (but not in subdirectories).

**Example output**:
```
- - - - - - - - - -
 J A N A  R factors:
- - - - - - - - - -
This program reads in JANA refinement *.ref files and gives an overview on selected parameters.

Current working directory:
D:\3DED\rocks\
User input required!
- If you just hit ENTER (empty input), all subfolders of the working directory are searched for JANA refinement files.
- If you input the path to a folder with JANA refinement file(s), only refinements of that folder are read.
PATH:
D:\refinements

O V E R V I E W :
"""""""""""""""""

 > > D:\refinements\single_block_refinement\ < <
File                              Parameters    Nobs    Nall  GOFobs  GOFall    Robs    Rall   wRall
DYN2_PbSb                                 74     988    1794    3.03    2.30    9.65   11.79   10.58

 > > D:\refinements\amazing_compound\multiblock\refinement\ < <
File                              Parameters    Nobs    Nall  GOFobs  GOFall    Robs    Rall   wRall
FINAL_ADP_cations                        184    2361    6299    3.23    2.03   10.56   17.04   11.59
|-----------------------------Block1             796    2194                   12.24   19.38   13.10
|-----------------------------Block2             661    2187                    9.63   18.61   10.93
|-----------------------------Block3             904    1918                    9.73   14.12   10.79

FINAL_IDP                                179    2361    6299    3.28    2.07   10.84   17.47   11.83
|-----------------------------Block1             796    2194                   12.29   19.51   13.13
|-----------------------------Block2             661    2187                    9.93   19.09   11.22
|-----------------------------Block3             904    1918                   10.19   14.78   11.18
```

### dyn+cif.py
> Update one CIF file with parameters relevant for 3D ED and especially dynamical refinement.
> Many parameters, including thickness and merged R factors, are are then written into the field _refine_special_details.
> In the header of dyn+cif.py, individual parameters (author information, microscope information) can be defined.
> Run "CIF make" in JANA2020 before running dyn+cif.py
> The main output is the updated CIF file. A new CIF file is written without overwriting the initial CIF file.
> **REQUIRES _3DED_utils.py**

**Example 1**: python dyn+cif.py C:/path_to_refinement/jana_job.m83
> Reads jana_job.m42, jana_job.m50, jana_job.m83, and jana_job.cif.
> The file jana_job_updated.cif is written.

**Example 2**: python dyn+cif.py
> The programs asks to copy the path to a JANA jobfile. If the CIF file does not have the standard name like the JANA job, you are asked to provide the path to the CIF file as well.

**Example output**:
```
F I L E S:
m42 (frames, ED settings)      D:\test\1-OO.m42
m50 (general settings)         D:\test\1-OO.m50
m83 (reflection list)          D:\test\1-OO.m83
cif (CIF from Jana CIF make)   D:\test\1-OO.cif
new, updated cif file          D:\test\1-OO_updated.cif


;
Calculated intensities based on dynamical theory of electron diffraction.
Number of individual data blocks in refinement: 1
Block     Thickness     Nobs    Nall    Robs    Rall   wRall
    1        786.81      417     513  0.1405  0.1570  0.1603
Thickness given in Angstrom.
A crystal with a thickness distribution of a 'wedge' is assumed.
The refined thickness corresponds to the maximum thickness of the shape.
The effective crystal thickness was assumed to be independent of the goniometer angle.

Refinement statistics relevant for the non-linear least-squares minimisation of wR(all):
Number of reflections in refinement (obs/all):   417 / 513
Number of reflections present more than once:    66
Number of reflections unique in point group 1:   447
Robs:   0.1405
Rall:   0.1570
wRall:  0.1603

Post-refinement analysis of symmetrically-equivalent reflections:
Number of unique reflections (obs/all): 150 / 190
MRobs:  0.1279
MRall:  0.1524
MwRall: 0.1426

REMOVED LOOP: _geom_hbond_atom_site_label_D
REMOVED LOOP: _geom_torsion_atom_site_label_1
REMOVED LOOP: _geom_angle_atom_site_label_1
REMOVED LOOP: _geom_bond_atom_site_label_1
REMOVED LOOP: _jana_atom_site_ADP_F_label
REMOVED LOOP: _jana_atom_site_ADP_E_label
REMOVED LOOP: _jana_atom_site_ADP_D_label
REMOVED LOOP: _jana_atom_site_ADP_C_label
REMOVED LOOP: _atom_site_aniso_label
REMOVED LOOP: _restr_equal_torsion_atom_site_label_1
REMOVED LOOP: _restr_equal_angle_atom_site_label_1
REMOVED LOOP: _restr_equal_distance_atom_site_label_1
REMOVED LOOP: _restr_torsion_atom_site_label_1
REMOVED LOOP: _restr_angle_atom_site_label_1
REMOVED LOOP: _restr_distance_atom_site_label_1
REMOVED LOOP: _diffrn_standard_refln_index_h
REMOVED LOOP: _exptl_crystal_face_index_h
REMOVED LOOP: _twin_individual_id
REMOVED: _diffrn_standards_interval_count         ?
REMOVED: _diffrn_standards_interval_time          ?
REMOVED: _diffrn_standards_decay_%                ?
REMOVED: _diffrn_radiation_monochromator          ?
REMOVED: _diffrn_source_power                     ?
REMOVED: _diffrn_source_current                   ?
Removed: 16 empty lines
CIF file written: D:\test\1-OO_updated.cif      
```

### R_factors_merge_postref.py
> Calculate merged R-factors for dynamical refinement against 3D ED data based on m50 and m83 files.
> For details see pages S2 and S3 in [Klar et al. 2023, Nat. Chem. 15, 848, SI](https://static-content.springer.com/esm/art%3A10.1038%2Fs41557-023-01186-1/MediaObjects/41557_2023_1186_MOESM1_ESM.pdf).

**Example 1**: python R_factors_merge_postref.py C:/path_to_refinement/jana_job.m50
> Reads jana_job.m50 and jana_job.m83. Prints merged R factors on screen.

**Example 2**: python R_factors_merge_postref.py
> The programs asks to copy the path to a JANA jobfile. Then prints merged R factors on screen.

**Example output**:
```
        F I L E   C H E C K
Dynamical refinement
Number of duplicates: 240
All reflections without duplicates: 1554
'P1':   988 / 1794  observed/all  reflections.
Merge:  176 / 218  observed/all  reflections.
MRobs    6.91  Robs    9.65
MRall    7.39  Rall   11.79
MwRall   7.72  wRall  10.58
```
