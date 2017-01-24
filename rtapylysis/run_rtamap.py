# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 15:41:05 2017

@author: Sonal Dey
"""
############################### START: USER PARAMETERS #########################################
import os
from func_rtapylysis import rtamap
############################### START: USER PARAMETERS #########################################
projectdir=os.path.join(os.getcwd(), "example")

# enter multiple filenames seperated by comma and within quotes
fnames=[ 'example.mat' ] 

xmin=75       # START: x-coordinate of the RTA map
xmax=850      # END: x-coordinate of the RTA map


# below t is for True, 'f' is for false
save_='t'     # for saving the plot in the current working directory

color_=32 # chose the color scheme for the RTA plots 32, 33, 2 works best
map_='auto' # auto or man for chosing the color limits for RTA:Fill range
min_=2.2  # minimum limit of the RTA color map (only works if map_ is set to man)    
max_=3.3  # maximum limit of the RTA color map (only works if map_ is set to man)

################################ END: USER PARAMETERS ##########################################
for filename in fnames:
    rtamap(projectdir, filename, xmin, xmax, save_, color_, map_, min_, max_)
############################### BELOW: Usage: rtapmap(*args) ###################################
'''
Usage: rtamap(projectdir, filename, xmin=0, xmax=10, save_str, color_choice, color_map_control, man_col_min, man_col_max)
-----
Function parameters:
-----
        projectdir: the project directory
        filename*: the name of the matfile (*.mat) stored in above projectdir.
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
            
        xmin: Temperature minimum.
        xmax: Temperature maximum.
        *save_str: option for saving the RTA map in the present working directory.
        *save_str: 
            't', 'T', 'True', 'true' for setting save option as true;
            'f', 'F', 'False', 'false' for setting save option as false.

        color_choice: chose the colormap for the RTA plots by any integers from 0 -35;
                32 ('cubehelix') and 33 ('viridis') work best
                For details of the colormap values, 
                In a python command line or in spyder, enter "from func_rtapylysis import *"
                followed by, col_dict and hit enter.
                32 and 33 represents colormaps that are closet to grayscale;
                the default 'jet' may not be used; read here: 
                https://jakevdp.github.io/blog/2014/10/16/how-bad-is-your-colormap
        
        color_map_control: 'auto' or 'man' for choosing the color limits for RTA: Full range sub-plot
        
        man_col_min, man_col_max: minimum, maximum limits of the RTA color map** 
                    **Give values regardless of 25-th set to 'man' or 'auto';
                    will only be used in the code if 25-th is set to 'man' though;
                    for 25-th set to 'auto', minimum and maximum are obtained from the
                    minimum and maximum values of the data             
-----
Returns: no return object but displays a 8x6 figure using seaborn with the
                rapid thermal annealing map (rtamap) in 'cubehelix' colormap.
                The x-axis is the given temperature range or 't' range;
                The y-axis is the whole d-spacing range of as calculated from range of 'q';
                The z-axis contains the log10 of 'data'.
'''