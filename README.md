# plot-mb
A Python 3 script for (hassle-free) plotting of Mößbauer (MB) spectra from parameter files (text files) or from [ORCA](https://orcaforum.kofo.mpg.de) 
output files.  

In case of ORCA output files, the script calculates the isomeric shift δ from ρ(0) according to the equation δ = ɑ*(ρ(0) - C) + β, 
ΔE<sub>Q</sub> is taken from the ORCA output. The parameters ɑ, C and β can be changed in the script or from the command line.  

The script can also read δ, ΔE<sub>Q</sub>, line width (fwhm, optional) and the ratio (optional) of components from a text file and plot the Mößbauer (MB) spectra.  

The Mößbauer spectrum is calculated as a summation of Lorentz functions. The script can plot and also export single or multiple spectra. 

## External modules
 `numpy` 
 `matplotlib`
 
## Quick start
 Start the script with:
```console
python3 plot-mb.py filename
```
it will save the plot as PNG:
`filename-mb.png`

## Command-line options
- `filename` , required: filename
- `-s` , optional: shows the `matplotlib` window
- `-n` , optional: do not save the spectrum
- `-e` , optional: export the line spectrum or spectra in a csv-like fashion; filename of the export is Atomname + ".dat" for each atom and "all.dat" for the whole spectrum
- `-w` `N` , optional: line width (fwhm) of the Lorentzian (default is `N = 0.2`) - ignored if 'fwhm' is specified in a text file with MB parameters
- `-xmin`  `N` , optional: start spectrum at N mm/s (automatic scaling if not specified)
- `-xmax`  `N` , optional: end spectrum at N mm/s (automatic scaling if not specified)
- `-a` `N`, optional: ɑ for δ = ɑ*(ρ(0) - C) + β (if not specified, ɑ is taken from the script) - ignored in the case of MB parameters from a text file
- `-b` `N`, optional: β for δ = ɑ*(ρ(0) - C) + β (if not specified, β is taken from the script) - ignored in the case of MB parameters from a text file
- `-c` `N`, optional: C for δ = ɑ*(ρ(0) - C) + β (if not specified, C is taken from the script) - ignored in the case of MB parameters from a text file
- `-sh` `N`, optional: shift  δ by +/- N mm/s

## Script options
There are numerous ways to configure the spectrum in the script:
Check `# plot config section - configure here` in the script. 

The delimiter for the line spectrum export can be changed by changing the value of `export_delim =`.

## Code options
Colors, line thickness, line styles and 
more can be changed in the code directly.

## Remarks
The PNG file will be replaced everytime you start the script with the same output file. 
If you want to keep the file, you have to rename it. 

The script can only handle singlets and doublets.

The key for the calculation of the isomeric shift ist the equation δ = ɑ*(ρ(0) - C) + β. You should adjust ɑ, β and C according to functional and basis to achieve reasonable results.

For MB parameters from a text file, have a look at [`mb-param.txt`](mb-param.txt). 

## Examples
### Example 1:
```console
python3 plot-mb.py test.out -w 0.6 -s
```
Open the ORCA output file `test.out`, select a line width of 0.6 (`-w 0.6`), show the spectrum (`-s`). A PNG of the spectrum is saved as `test-mb.png`.
Output:
```console
=============================================================
δ = ɑ(ρ(0) - C) + β (+ shift)
ɑ = -0.366 C = 11810 β = 2.852 shift = 0
=============================================================
0Fe: δ = 0.47 mm/s ΔEQ = 3.45 mm/s ρ(0)=11816.51263 a.u.⁻³
1Fe: δ = 0.27 mm/s ΔEQ = 1.96 mm/s ρ(0)=11817.06310 a.u.⁻³
=============================================================
```
![Example 1](/examples/example1a.png)

### Example 2:
```console
python3 plot-mb.py test.out -s -e
```
Open the ORCA output file `test.out`, show the spectrum (`-s`), export the data (`-e`). The files `Fe0.dat`, `Fe1.dat`, `Fe4.dat` and `all.dat` containing the data of the doublets and the entire spectrum are saved. 
A PNG of the spectrum is saved as `test-mb.png`.
Output:
```console
=============================================================
δ = ɑ(ρ(0) - C) + β (+ shift)
ɑ = -0.366 C = 11810 β = 2.852 shift = 0
=============================================================
0Fe: δ = 0.61 mm/s ΔEQ = 2.37 mm/s ρ(0)=11816.13420 a.u.⁻³
1Fe: δ = 0.88 mm/s ΔEQ = -2.89 mm/s ρ(0)=11815.37543 a.u.⁻³
4Fe: δ = 1.74 mm/s ΔEQ = -2.90 mm/s ρ(0)=11813.03611 a.u.⁻³
=============================================================
```
![Example 2](/examples/example2a.png)

### Example 3:
```console
python3 plot-mb.py test.out -xmin -12 -xmax 12 -sh 2 -s
```
Open `test.out`, start the spectrum at -12 mm/s and end at 12 mm/s (`-xmin -12 -xmax 12`), shift δ by 2 mm/s (`-sh 2`), show the spectrum (`-s`). A PNG of the spectrum is saved as `test-mb.png`.
Output:
```console
==========================================================
δ = ɑ(ρ(0) - C) + β (+ shift)
ɑ=-0.366 C=11810 β=2.852 shift=2.0
==========================================================
Fe0: δ=3.78 mm/s ΔEQ=-0.00 mm/s ρ(0)=11812.93015 a.u.⁻³
==========================================================
```
![Example 3](/examples/example3.png)

### Example 4:
```console
python3 plot-mb.py test.out -s
```
Open `test.out` and show the spectrum (`-s`). A PNG of the spectrum is saved as `test-mb.png`.
Output:

![Example 4](/examples/show-use2.gif)
