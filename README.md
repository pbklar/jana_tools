# Jana Tools
A collection of tools to read, update, write, and call files related to the software programs [Jana2020](http://jana.fzu.cz/) and [PETS2](http://pets.fzu.cz/). The tools are currently developed by Paul Klar (University of Bremen). Note that Paul Klar is **not** a developer of [Jana2020](http://jana.fzu.cz/) and [PETS2](http://pets.fzu.cz/). 

# Contact
Paul Klar, University of Bremen

# Install and update package via pip
Install via pip:
```
pip install git+https://github.com/pbklar/jana_tools.git
```

Update via pip:
```
pip install --upgrade git+https://github.com/pbklar/jana_tools.git
```

# Install and update package using git
Install using git:
```
cd D:/path/to/desired_destiation
git clone https://github.com/pbklar/jana_tools.git
cd jana_tools
pip install -e .
```

Update:
```
cd D:/path/to/desired_destiation
git clone https://github.com/pbklar/jana_tools.git
```


#Preliminary usage guide
In some better future, this section will move to readthedocs using Sphinx. Until then:

After installation a script called jana-tools.exe (on Windows) is available. If the PATH variables are properly set, this tool can be called from any folder.

The first argument is the tool you want to use. Currently, there is only one tool labelled *runtime* implemented. This tool takes one additional argument for the file path to a refinement \*.ref file written by Jana2020 or a folder with several Jana2020 jobs.

Example call of Jana Tools:
```
jana-tools.exe runtime -p D:\path\to\refinement.ref
```

Example output:
```
Refinement:
***********
 # of parameters:             65
 # of reflections:          5932
 # of cycles N:                2 <-- with derivatives and parameter shifts
 # of R-factor cycles:         1 <-- no derivatives, no parameter shifts
CPU time:                0:06:37 <-- JANA + DYNGO
Non-CPU time:            0:00:05 <-- Read/Write files

Refinement timestamps:
*******************
Start Refinement:      2025-04-21 08:17:49
Start R-factor cycle:  2025-04-21 08:24:01
End:                   2025-04-21 08:24:31

Refinement runtime:
*******************
 Runtime N cycles:       0:06:12 <-- without R-factor cycle
 R-factor cycle:    +    0:00:30 <-- no derivatives, no parameter shifts
                      ----------
 Total runtime:     =    0:06:42

 Runtime per cycle:      0:03:06 <-- with derivatives and parameter shifts
```

# List of tools

## Implemented
- runtime: Jana2020 refinement runtime analysis

## Planned
- Jana2020
    - batch modify Jana refinement files
    - command-line based call of JANA
    - absolute structure hypothesis testing
    - detailed R factor reporting

- PETS2
    - jana job generator from PETS2 output + CIF
    - virtual frames identifier + visualiser
    - Aperpixel guide
    - autopets 
    - frame correlation analysis
    - frame or image file format converter
    - mask generator
    
- Crystallography
    - mini toolbox
    - scattering form factors
    - electron wavelength
    - rocking curve visualiser
    
    
<!--
## Envisioned
- ShelX <-> Jana restraints/constraints conversion
- structure model comparison tool
- ESP tools (xplor)
-->    