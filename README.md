# orca-mb
A Python 3 script for (hassle-free) plotting of Mößbauer spectra from [ORCA](https://orcaforum.kofo.mpg.de) 
output files. The script calculates the isomeric shift δ from ρ(0) according to the equation δ = ɑ*(ρ(0) - C) + β, 
ΔE<sub>Q</sub> is taken from the ORCA output. 
The parameters ɑ, C and β can be changed in the script or from the command line. The Mößbauer spectrum is calculated as a summation
of Lorentz functions. The script can plot and also export single or multiple spectra. 

## External modules
 `numpy` 
 `matplotlib`
 
## Quick start
 Start the script with:
```console
python3 orca-mb.py filename
```
it will save the plot as PNG:
`filename-mb.png`
