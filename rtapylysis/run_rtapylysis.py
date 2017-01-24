# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 15:41:05 2017

@author: Sonal Dey
"""
############################### START: USER PARAMETERS #########################################
import os
import func_rtapylysis
############################### START: USER PARAMETERS #########################################
projectdir=os.path.join(os.getcwd(), "example")

# enter multiple filenames seperated by comma and within quotes
fnames=[ 'example.mat' ] 

t1=75       # START: x-coordinate of the RTA map
t2=850      # END: x-coordinate of the RTA map
d1=1.95     # START: y-coordinate of the RTA map     
d2=3.65     # END: y-coordinate of the RTA map

# 3 x-values at which XRD plot is obtained
txrd=[ 100, 600, 800 ] 
# Label text for the phases for 3 x-values above
txrdlabel=[ "phase-a", "phase-b", "phase-c" ] 
# for 3 phases above, dhkl file names
txrdphase=[ "phaseA.txt", "phaseB.txt", "phaseC.txt" ] 
# Label text for the phases for 3 x-values above during indexing
txrdphaselabel=[ "a", "b", "c" ] 

mphvalin=0.05   # for detection of phase transition temp. in rate-peak-area-plot
mpdvalin=10     # for detection of phase transition temp. in rate-peak-area-plot

# below t is for True, 'f' is for false
save_='t'                         # for saving the final plot in the current working directory
smooth_XRD_='t'                   # for smoothing the 3 XRD plots
image_smooth_bkg_='t'             # option to smooth the RTA after background subtraction
smooth_peak_area_and_grad_='t'    # option to smooth the peak area and gradient plots

color_=32 # chose the color scheme for the RTA plots 32, 33, 2 works best
map_='auto' # auto or man for chosing the color limits for RTA:Fill range
min_=2.2  # minimum limit of the RTA color map (only works if map_ is set to man)    
max_=3.3  # maximum limit of the RTA color map (only works if map_ is set to man)

# Range of d-spacings along which peak area is calculated after background subtraction
drange=[ 1.955, 2.005, 2.050, 2.127, 2.160, 2.245, 2.290, 2.490, 2.750, 2.815, 2.865, 2.975 ]
################################ END: USER PARAMETERS ##########################################
for filename in fnames:
    args=[projectdir, filename, t1, t2, d1, d2 ] + \
         txrd + txrdlabel + txrdphase + txrdphaselabel + [mphvalin, mpdvalin] + \
         [save_, smooth_XRD_, image_smooth_bkg_, smooth_peak_area_and_grad_] + \
         [color_, map_, min_, max_] + \
         drange    
    func_rtapylysis.rtapylysis(args)
############################# BELOW: Usage: rtapylysis(*args) ##################################    
'''
Usage: rtapylysis(*args)
-----
*args or function parameters must appear in the following order, separated by comma:
*args or function parameters must appear in the following order, separated by comma:
-----
        1st: projectdir: the project directory
        2nd: filename*: the name of the matfile (*.mat) stored in above projectdir.
            *Below is a list of tuples as returned by scipy.io.whosmat("example.mat"), 
            one for each array (or other object) in the file. Each tuple contains 
            the name, shape and data type of the array.
            [('Temp', (272, 1), 'double'),
            ('time', (272, 1), 'double'),
            ('tth', (1, 638), 'double'),
            ('ekev', (1, 1), 'double'),
            ('q', (1, 638), 'double'),
            ('intensity', (272, 638), 'double')]
            
            * The input filename needs to be of same format as above example where the above structure format is maintained 
            * (272 and 638 are specific to above example.mat; they can be changed).
            * Atleast the 'Temp', 'q', 'data' numpy arrays needs to be present in filename. 
            * Log10 of the 'data' numpy array is taken; so 'data' should not have values < 0.
            
        3-rd- 4-th:  Temperature minimum and maximum.
        5-th- 6-th:  d-spacing minimum and maximum (in Ang).
        6-th- 8-th:  3 x-values or Temperature values at which XRD plot is obtained.
        9-th-11-th:  3 Label text for the phases for 3 x-values above.
        12-th-14-th: 3 dhkl_file names for 3 phases above, all contained inside dhkl_dir. 
        15-th-17-th: 3 Label text for the phases for 3 x-values above during indexing.
        18-th-19-th: for detection of phase transition temperatures;
                     18-th is mph: detect peaks that are greater than minimum peak height;
                     19-th is mpd: detect peaks that are at least separated by minimum peak distance (in number of data).
                     mph=0.05 and mpd=10 seem to work well for the example file; change as necessary.

        * For 20-th - 23-rd: 
                't', 'T', 'True', 'true' for setting save option as true;
                'f', 'F', 'False', 'false' for setting save option as false.
        *20-th: option for saving the final plot in the current working directory.                    
        *21-st: option for smoothing the 3 XRD plots.
        *22-nd: option to smooth the RTA after background subtraction. 
        *23-rd: option to smooth the peak area and gradient plots.
        
        24-th: chose the colormap for the RTA plots by any integers from 0 -35;
                32 ('cubehelix') and 33 ('viridis') work best
                For details of the colormap values, 
                In a python command line or in spyder, enter "from func_rtapylysis import *"
                followed by, col_dict and hit enter.
                32 and 33 represents colormaps that are closet to grayscale;
                the default 'jet' may not be used; read here: 
                https://jakevdp.github.io/blog/2014/10/16/how-bad-is-your-colormap
        
        25-th: 'auto' or 'man' for choosing the color limits for RTA: Full range sub-plot
        26th, 27th: minimum, maximum limits of the RTA color map** 
                    **Give values regardless of 25-th set to 'man' or 'auto';
                    will only be used in the code if 25-th is set to 'man' though;
                    for 25-th set to 'auto', minimum and maximum are obtained from the
                    minimum and maximum values of the data  
                    
        28th onwards: Must appear in pairs of increasing d-spacing values (in Angstrom), 
                    e.g, 1.955, 2.005, 2.050, 2.127 and so on.
    
    -----           
    Returns: No return object but displays a 14x18 figure using seaborn with the following sub-plots
             (Saves the figure to present working directory if save_option is True):
                sub-plot 1 (*RTA: Full range); sub-plot 2 (*RTA: Partial);
                sub-plot 3 (XRD: phase 1);    sub-plot 4 (Peak area plots for given d-spacing ranges);
                sub-plot 5 (XRD: phase 2);    sub-plot 6 (Peak area gradients for sub-plot 4 with transition temperatures marked from peak positions);
                sub-plot 7 (XRD: phase 3); 
                
                * RTA: Rapid thermal annealing map. For RTA,
                * The x-axis is the given temperature range or 't' range;
                * The y-axis is the whole or partial d-spacing range as calculated from whole or partial range of 'q';
                * The z-axis contains the log10 of 'data' and/or sub-range of 'data'. 
    -----
'''
