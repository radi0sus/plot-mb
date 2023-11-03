# -*- coding: utf-8 -*-
'''

plot-mb

'''

import sys                                      #sys files processing
import os                                       #os file processing
import re                                       #regular expressions
import numpy as np                              #summation and other math
import matplotlib.pyplot as plt                 #plots
from matplotlib.ticker import AutoMinorLocator  #for minor ticks
import argparse                                 #argument parser

#for windows console
sys.stdout.reconfigure(encoding = 'utf-8')  

'''
Only for ORCA output files and 57Fe:

  M. Römelt, S. Ye, F. Nesse, Inorg. Chem., 2009, 48, 784–785. 
  DOI: 10.1021/ic801535v
 
  delta (I.S.) = alpha*(rho - C) + beta
 
 =========================================
 functional basis   alpha   beta    C   
 =========================================
 BP86       CP(PPP) -0.425  7.916   11810
            TZVP    -0.340  1.034   11580 
            TZVPa   -0.362  4.957   13800
 B3LYP      CP(PPP) -0.366  2.852   11810 
            TZVP    -0.298  1.118   11580 
            TZVPa   -0.307  4.045   13770 
 TPSS       CP(PPP) -0.421  5.154   11810
            TZVP    -0.336  1.327   11580 
            TZVPa   -0.365  1.385   13800 
 TPSSh      CP(PPP) -0.376  4.130   11810 
            TZVP    -0.321  1.466   11580
            TZVPa   -0.322  1.830   13780 
 B2PLYP     CP(PPP) -0.336  2.642   11810 
            TZVP    -0.261  1.483   11580 
            TZVPa   -0.311  2.256   13790 
 =========================================
 a) ZORA recontracted basis sets    
 =========================================
 
 change parameters here or use -a -b -c options
 (this is B3LYP & CP(PPP))
'''
##########################################
alpha    = -0.366
beta     =  2.852
C        =  11810
##########################################
#
'''
 The script can also read δ and ΔEQ directly from a text file.
 The parameters must be entered in the following order:
 
 compound/atom_name_1 δ_1 ΔEQ_1 fwhm_1 ratio_1
 compound/atom_name_2 δ_2 ΔEQ_2 fwhm_2 ratio_2
 compound/atom_name_3 δ_3 ΔEQ_3 fwhm_3 ratio_3
 ...
 Lines starting with '#' are ignored by the script.
 
 'fwhm' and 'ratio' are optional. If there is no value for 'fwhm', 
 but there is a value for 'ratio', 'ratio' is considered to be 'fwhm', 
 because the third parameter in the line is assumed to be 'fwhm'.
 
 If 'fwhm' is not specified it is set to 0.2 and can be changed
 with the -w option. If 'fwhm' is not specified, then 'ratio' 
 cannot be specified (see remark above).
 
 If 'ratio' is not specified it is set to 1.
 Defining a ratio different from 1 or 100% is useful if there is 
 a main component and an impurity which is also MB active
 or a mixture of two or more compounds with MB active nuclei.

 'ratio' can be defined as 1:2:1:2 for example or 0.15:0.50:0.15:0.20 or 75:25...

 The -a, -b, -C options are ignored in case of MB parameters from a text file.
 The -w option is ignored if 'fwhm' is specified in the text file.
'''
#
##########################################
is_shift  =  0                              #shift delta (I.S.) by +x mm/s or -x mm/s
##########################################
scale_mul =  1                              #y range for export; for scale_mul = '10' ymin is 0, ymax is 10 
export_delim = " "                          #delimiter for data export

# plot config section - configure here
colors = ['red','blue','green','grey']      #colors for plots in case of multiple MB-active atoms
show_area = True                            #show area (filled) plots in case of multiple MB-active atoms if True
alpha_line = 0.3                            #value for the alpha channel for line plots in case of multiple MB-active atoms (0-1)
alpha_area = 0.1                            #value for the alpha channel for area plots in case of multiple MB-active atoms (0-1)
linewidth = 0.2                             #line width for broadening (FWHM)
draw_v_line = False                         #draw a vertical line at delta (I.S.) if True
show_legend = True                          #show legend if True
show_minor = False                          #show minor ticks on x-axis if True
spectrum_title = r'Mößbauer spectrum'       #title
spectrum_title_weight = "light"             #weight of the title font: 'normal' | 'bold' | 'heavy' | 'light' | 'ultrabold' | 'ultralight'
y_label = "relative transmission"           #label of the y-axis 
x_label = r'velocity /mm$\cdot$s$^{-1}$'    #label of the x-axis
figure_dpi = 300                            #DPI of the PNG
# end plot config section

#global lists
deltaeqlist = list()                        #list with Delta-EQ
rho0list = list()                           #list with rho(0) from orca.out
nucnamelist = list()                        #list of atom / compound names
ishiftlist = list()                         #list of isomeric shifts (deltas)
lorentz_sum = list()                        #list for the sum of single lorentzians = the convoluted spectrum
xminmaxlist = list()                        #list for max and min delta (I.S.) + Delta-EQ/2 for auto range in x
ratiolist = list()                          #list with ratios of several MB active nuclei; set to 1 in case of orca.out
wlist = list()                              #list with individual w (fwhm) or from -w option

#gauss and pseudo-voigt profiles ... in case
#but must be renamed to 'lorentz' or more code
#has to be changed
#
#def gauss(a,m,x,w,ratio):
#    # calculation of the Gaussian line shape
#    # a = amplitude (max y, intensity)
#    # x = position
#    # m = maximum/median (stick position in x, wave number)
#    # w = line width, FWHM
#    # ratio = extra parameter for ratio of several MB active species
#    # ratio is always 1 in case of data calculated with ORCA
#    return a*np.exp(-(np.log(2)*((m-x)/w)**2))*abs(ratio)

#def pvoigt(a,x,m,w,ratio):
#    # calculation of Pseudo-Voigt: pvoigt = h*lorentz + (1-h)*gauss; 0 < h < 1
#    h = 0.5
#    return h*(a/(1+((m-x)/w)**2))+(1-h)*(a*np.exp(-(np.log(2)*((m-x)/w)**2)))*abs(ratio)

def lorentz(a,m,x,w,ratio):
    # calculation of the Lorentzian line shape
    # a = amplitude (max y, intensity)
    # x = position
    # m = maximum/median (I.S. (delta) in x, mm/s)
    # w = line width, FWHM
    # ratio = extra parameter for ratio of several MB active species
    # ratio is always 1 in case of data calculated with ORCA
    #return abs(a/(1+((m-x)/w)**2))*abs(ratio)
    return abs((a/np.pi)*((w/2)/((m-x)**2+(w/2)**2)))*abs(ratio)
    
def rho_to_ishift(rho):
    #rho(0) to I.S. (delta)
    return args.alpha*(rho - args.C) + args.beta + args.shift

#considered unnecessary    
#def switchme(label):
#    #switches label from orca.out, 0Fe --> Fe0
#    newlabel = re.split('(\d+)',label)
#    return newlabel[2]+newlabel[1]

def export_csv(filename,scale_mul,index):
    #export spectra in csv-like fashion
    #get data from plot (window)
    plotdata = ax.lines[index]
    xdata = plotdata.get_xdata()
    #for transmission style 1 (or scale_mul) - data
    ydata = scale_mul-plotdata.get_ydata()
    try:
        with open(filename + ".dat","w") as output_file:
            for elements in range(len(xdata)):
                output_file.write(str(xdata[elements]) + export_delim + str(ydata[elements]) +'\n')
    #write error -> exit here
    except IOError:
        print("Write error. Exit.")
        sys.exit(1)

#parse arguments
parser = argparse.ArgumentParser(prog = 'orca-mb', description = 'Easily plot Mößbauer spectra')

#filename is required
parser.add_argument("filename", help = "file with Mößbauer parameters or ORCA output file")

#show the matplotlib window
parser.add_argument('-s','--show',
    default = 0, action = 'store_true',
    help = 'show the plot window')

#do not save the png file of the spectra
parser.add_argument('-n','--nosave',
    default = 1, action = 'store_false',
    help = 'do not save the spectra as PNG')

#change line width for line broadening, default is 0.2
#override default value from the top of the script
#ignored if 'fwhm' is specified in a text file with MB parameters
parser.add_argument('-w','--linewidth',
    type = float,
    default = linewidth,
    help = 'line width (fwhm) for broadening - 0.2 if not specified')

#turn off automatic scaling, start spectrum at xmin
parser.add_argument('-xmin','--xmin',
    type = float,
    default = 0,
    help = 'start spectrum - automatic scaling if not specified')

#turn off automatic scaling, end spectrum at xmax
parser.add_argument('-xmax','--xmax',
    type = float,
    default = 0,
    help = 'end spectrum - automatic scaling if not specified')

#alpha for delta (I.S.) = alpha*(rho - C) + beta
#override default values from the top of the script
#ignored in the case of MB parameters from a text file
parser.add_argument('-a','--alpha',
    type = float,
    default = alpha,
    help = 'ɑ for δ (I.S.) = ɑ*(ρ(0) - C) + β')

#beta for delta (I.S.) = alpha*(rho - C) + beta
#override default values from the top of the script
#ignored in the case of MB parameters from a text file
parser.add_argument('-b','--beta',
    type = float,
    default = beta,
    help = 'β for δ (I.S.) = ɑ*(ρ(0) - C) + β')

#C for delta (I.S.) = alpha*(rho - C) + beta
#override default values from the top of the script
#ignored in the case of MB parameters from a text file
parser.add_argument('-c','--C',
    type = float,
    default = C,
    help = 'C for δ (I.S.) = ɑ*(ρ(0) - C) + β')

#shift delta (I.S.) by +x or -x mm/s
parser.add_argument('-sh','--shift',
    type = float,
    default = is_shift,
    help = 'shift δ (I.S.) by +x or -x mm/s')

#export data for the line spectrum in a csv-like fashion
parser.add_argument('-e','--export',
    default = 0, action = 'store_true',
    help = 'export data')

#pare arguments
args = parser.parse_args()

#just tell the script that w is line width (fwhm)
w = args.linewidth

#open an orca.out file
#check existence
try:
    with open(args.filename, "r") as input_file:
        for line in input_file:
            #start exctract atom names, Delta-EQ and rho(0)
            #in 3 lists
            if "Nucleus:" in line:
                nucnamelist.append(line.strip().split()[1])
            if "Delta-EQ" in line:
                deltaeqlist.append(float(line.strip().split()[5]))
            if "RHO(0)" in line:
                rho0list.append(float(line.strip().split()[1]))
            #add w (fwhm) to wlist 
            #add 1 to ratiolist, no ratio or ratio always 1 in calculated data
            #from orca.out
            wlist.append(args.linewidth)
            ratiolist.append(1)
#file not found -> exit here
except IOError:
    print(f"'{args.filename}'" + " not found")
    sys.exit(1)

#check if MB parameters are from orca.out; existence of deltaeqlist
if not deltaeqlist:
    print("This file does not contain MB parameters from an ORCA calculation.")
    print("Trying to read MB parameters directly.")
    #if not, assume text file with MB parameters
    #reset the above wlist
    #reset the above ratiolist
    wlist.clear()
    ratiolist.clear()
    try:
        with open(args.filename, "r") as input_file:
            for line in input_file:
                if not line.startswith("#"):
                    nucnamelist.append(line.strip().split()[0])
                    ishiftlist.append(float(line.strip().split()[1])+args.shift)
                    deltaeqlist.append(float(line.strip().split()[2]))
                    try:
                        wlist.append(float(line.strip().split()[3]))
                    except IndexError:
                        wlist.append(args.linewidth)
                    try:
                        ratiolist.append(float(line.strip().split()[4]))
                    except IndexError:
                        ratiolist.append(1)
#file not found -> exit here
    except IOError:
        print(f"'{args.filename}'" + " not found")
        sys.exit(1)
#strings where should be numbers -> exit here
    except ValueError:
        print("Warning! Numerical value expected. Exit.")
        sys.exit(1)
#no values found -> exit here
    except IndexError:
        print("Warning! Value missing. Exit.")
        sys.exit(1)
    
#no MB parameters in either orca.out or text file with MB parameters -> exit here
if not deltaeqlist:
    print("Warning! No Mößbauer parameters found in " + f"'{args.filename}'. Exit.")
    sys.exit(1)

#values from rho(0) list to list with isomeric shifts
#function delta (I.S.) = alpha*(rho - C) + beta is applied
#only for values from orca.out
if not ishiftlist:
    ishiftlist = [rho_to_ishift(rho) for rho in rho0list]

#autoscale in case of xmin and xmax args are not given (or zero)
if args.xmin == 0 and args.xmax == 0:
    #calculate all delta (I.S.) +/- Delta-EQ/2 to obtain the max and min values
    #for x autoscale
    for index, ishift in enumerate(ishiftlist):
        xminmaxlist.append(ishift-abs(deltaeqlist[index]/2))
        xminmaxlist.append(ishift+abs(deltaeqlist[index]/2))

    #obtain the max and min values for auto scale
    #add or subtract 1, also take w or mean w from wlist into account
    x_min = min(xminmaxlist) - np.mean(wlist) * 5 - 1
    x_max = max(xminmaxlist) + np.mean(wlist) * 5 + 1
else:
    #dont use autoscale
    x_min = args.xmin 
    x_max = args.xmax 

#plotrange in x from x_min and x_max
#change resolution with the last parameter if necessary (default = 0.01) 
plt_range_x = np.arange(x_min, x_max, 0.01)

#print parameters for delta (I.S.) calculation
if rho0list:
    #for MB parameters from orca.out
    print('=============================================================')
    print('δ = ɑ(ρ(0) - C) + β (+ shift)')
    print('ɑ = '+str(args.alpha), 'C = '+str(args.C), 'β = '+str(args.beta), 'shift = '+str(args.shift))
print('=============================================================')
for index, ishift in enumerate(ishiftlist):
    #generate summation of single lorentz functions
    #every lorentzian of a doublet
    #consists of two single lorentzian functions 
    #delta (I.S.) + Delta-EQ/2 and delta (I.S.) - Delta-EQ/2
    #it is here to save an extra iteration
    #ratio of different compounds is considered, ratio is always 1 for orca.out
    lorentz_sum.append(lorentz(1, plt_range_x, ishift - abs(deltaeqlist[index]/2), 
    wlist[index], ratiolist[index]))
    lorentz_sum.append(lorentz(1, plt_range_x, ishift + abs(deltaeqlist[index]/2), 
    wlist[index], ratiolist[index]))
    ####
    #print the results of the calculation: delta (I.S.) = alpha*(rho - C) + beta
    #and Delta-EQ and rho(0) from orca.out
    if rho0list:
        #for MB parameters from orca.out
        print(nucnamelist[index]+':', 
            'δ = {:.2f}'.format(ishift), 
            'mm/s',
            'ΔEQ = {:.2f}'.format(deltaeqlist[index]), 
            'mm/s',
            'ρ(0)={:.5f}'.format(rho0list[index]),
            'a.u.⁻³')
    else:
        #for MB parameters from text file
        print(nucnamelist[index]+':', 
        'δ ={: .2f}'.format(ishift), 
        'mm/s',
        'ΔEQ ={: .2f}'.format(deltaeqlist[index]), 
        'mm/s',
        'fwhm ={: .2f}'.format(wlist[index]),
        'ratio ={: .2f}'.format(ratiolist[index]))
print('=============================================================')

#prepare plot
fig, ax = plt.subplots()

#get colors from config section on top
ax.set_prop_cycle(color = colors)

#summation of all lorentz functions values in y
plt_range_lorentz_sum_y = np.sum(lorentz_sum, axis = 0)

#plot single lorentz function for every frequency freq

for index, ishift in enumerate(ishiftlist):
    #every lorentz function of a doublet
    #consists of two single lorentzian functions 
    #delta (I.S.) + Delta-EQ/2 and delta (I.S.) - Delta-EQ/2
    #ratio of different compounds is considered, ratio is always 1 for orca.out
    #generate an empty list for every doublet
    lorentz_doublet = []
    #generate the doublet from two singlets
    lorentz_doublet.append(lorentz(1/plt_range_lorentz_sum_y.max() * scale_mul, 
                           plt_range_x, ishift - abs(deltaeqlist[index]/2), 
                           wlist[index], ratiolist[index]))
    lorentz_doublet.append(lorentz(1/plt_range_lorentz_sum_y.max() * scale_mul, 
                            plt_range_x, ishift + abs(deltaeqlist[index]/2), 
                            wlist[index], ratiolist[index]))
                            
    #if only one atom is present: no filled plot, black line color for plot 
    if len(ishiftlist) == 1:
        ax.plot(plt_range_x,np.sum(lorentz_doublet, axis = 0),
                label = nucnamelist[index] + 
                r': $\delta$ = '+'{:.2f}'.format(ishift) + r' mm$\cdot$s$^{-1}$' +
                r', $ΔE_Q =$' +'{:.2f}'.format(deltaeqlist[index]) + 
                ' mm$\cdot$s$^{-1}$', color = 'black') 
        if draw_v_line:
            #draw vertical line at I.S. if True
            ax.axvline(x = ishift, linestyle = '--', color = 'black')
        if args.export:
            #save data on request
            export_csv(nucnamelist[index],scale_mul,index)
    else:
        #plots for more than one atom 
        line, = ax.plot(plt_range_x,np.sum(lorentz_doublet, axis = 0),
                label = nucnamelist[index] +  
                r': $\delta$ = '+'{: .2f}'.format(ishift) + r' mm$\cdot$s$^{-1}$' + 
                r', $ΔE_Q =$' +'{: .2f}'.format(deltaeqlist[index]) + 
                ' mm$\cdot$s$^{-1}$',alpha = alpha_line) 
        if draw_v_line:
            #draw vertical line at I.S. if True
            ax.axvline(x = ishift, linestyle = '--',color = line.get_color(), alpha = 0.8)
        if args.export:
        #save data on request
            export_csv(nucnamelist[index],scale_mul,index)
            
    if len(ishiftlist) > 1 and show_area:
        #fill plot area on request, but only in case of more than one atom
        ax.fill_between(plt_range_x,np.sum(lorentz_doublet, axis = 0), alpha = alpha_area) 



#increase figure size N times
N = 1.5
params = plt.gcf()
plSize = params.get_size_inches()
params.set_size_inches((plSize[0]*N, plSize[1]*N))

#plot the overall (convoluted) spectrum
ax.plot(plt_range_x,plt_range_lorentz_sum_y/plt_range_lorentz_sum_y.max()*scale_mul,color = "black",linewidth = 1)

#set y and x limits; y + 10% for legend etc.
plt.ylim(plt_range_lorentz_sum_y.max()/plt_range_lorentz_sum_y.max()*scale_mul+scale_mul*0.1,0)
plt.xlim(x_min,x_max)

#save data for the overall (convoluted) spectrum on request, but only in case of more than one atom
#in case of one atom, this is already achieved in the above section
if args.export and len(ishiftlist) > 1:
    export_csv('all',scale_mul,len(ishiftlist))

#show legend
if show_legend:
    ax.legend(fancybox = True, shadow = True, prop = {'size': 9})

ax.set_title(spectrum_title,fontweight = spectrum_title_weight)   #title from config section
ax.get_yaxis().set_ticks([])                                      #remove ticks from y axis
if show_minor:
    ax.xaxis.set_minor_locator(AutoMinorLocator())                #if True, show minor ticks on x-axis
ax.set_ylabel(y_label)                                            #set ylabel from config section
ax.set_xlabel(x_label)                                            #set xlabel from config section

#use tight layout
plt.tight_layout()

#(do not) save the PNG
if args.nosave:
    filename, file_extension = os.path.splitext(args.filename)
    plt.savefig(f"{filename}-mb.png", dpi = figure_dpi)
    
#(do not) show the plot
if args.show:
    plt.show()
