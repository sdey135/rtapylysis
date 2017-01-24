# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 16:04:41 2017

@author: Sonal Dey
"""
import os
#import sys
import numpy as np
import scipy.io
import scipy.linalg as LA
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pylab import pi, log10, shape
import seaborn as sns
from detect_peaks import detect_peaks

#==============================================================================
## some useful variables: dhkl_dir and col_dict
#==============================================================================
# the dhkl files are called for indexing purposes from dhkl_dir
dhkl_dir = os.path.join(os.getcwd(), 'dhkl_dir')
####
#colormap to be used: http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
#Also see for colormap: http://matplotlib.org/users/colormaps.html
#Interesting read on colormap and unwanted luminance spikes in default 'jet': 
#https://jakevdp.github.io/blog/2014/10/16/how-bad-is-your-colormap/
# #32: 'cubehelix' and #33: 'viridis' are the top color schemes to be used
col_dict={ 0:'rainbow', 1:'binary', 2:'gist_yarg', 3:'gray', 4:'bone', 5:'gist_earth', 6:'brg',\
           7:'gist_stern', 8:'copper', 9:'jet', 10:'rainbow_r', 11:'coolwarm_r', 12:'hot',\
           13:'spectral', 14:'nipy_spectral', 15:'Paired', 16:'coolwarm', 17:'PiYG_r',\
           18:'seismic', 19:'Greys', 20:'afmhot', 21:'RdBu', 22:'YlGn', 23:'RdGy', 24:'cool',\
           25:'autumn', 26:'YlOrRd', 27:'YlOrBr', 28:'YlGnBu', 29:'YlGn',30:'gist_rainbow',\
           31:'hsv', 32: 'cubehelix', 33: 'viridis', 34: 'prism' , 35: 'flag'  }
#
#==============================================================================
## some useful functions
#==============================================================================
def rtamap(projectdir, filename, xmin, xmax, save_str, color_choice, color_map_control, man_col_min, man_col_max):
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
    ###################################################################################################
    fig = plt.figure(facecolor='w', figsize=[8.0,6.0])
    sns.set_style("ticks")    # using seaborn (http://web.stanford.edu/~mwaskom/software/seaborn/)
    sns.set_context("talk", font_scale=1.8, rc={"lines.linewidth": 2.0})
    ####################################################################################################
    filepath = os.path.join(projectdir, filename)
    matfile = scipy.io.loadmat(filepath)
    t = matfile['Temp'][:,0] #converting Temperature to 1D array    
    q = matfile['q'][0,:] #converting q vector to 1D array
    data = matfile['intensity'] # 2D array; 
    #time = matfile['time'][:,0] #converting time to 1D array
    #tth = matfile['tth'][0,:] #converting 2-theta to 1D array
    #ekev = matfile['ekev'][0,0] # getting the x-ray energy in keV 

    # find the indices near user given values for abscissa and ordinate 
    xmin, xmax = float(xmin), float(xmax)
    id_x_min = np.argmin(np.abs(t - xmin))
    id_x_max = np.argmin(np.abs(t - xmax))

    # select parts of the matfile variables for plotting
    # x-range or Temperature range here
    x = t[id_x_min:id_x_max]
    # for complete RTA plot in specified x-range
    y = q
    z = data[id_x_min:id_x_max, : ]
    
    # Full RTA for the chosen x-range, no BKG subtraction
    if color_map_control == 'auto':        
        int_min = np.log10(np.min(np.min(z.T)))
        int_max = np.log10(np.max(z.T))
        print "\nRTA: Full range: cmap: auto min and max -> %.2f and %.2f" % (int_min, int_max)
    else:
        int_min = float(man_col_min) 
        int_max = float(man_col_max) 
        print "\nRTA: Full range: cmap: man min and max -> %.2f and %.2f" % (int_min, int_max)    

        
    ax = plt.subplot(1,1,1)
    ax.pcolormesh(x,y,log10(z.T), cmap=col_dict[int(color_choice)], vmin=int_min, vmax=int_max)
    ax.axis('tight');
    
    # put inverse levels (i.e., d=2*pi/q instead of q, for the y-axis
    get_y_ticks = np.array(ax.get_yticks().tolist())
    new_y_ticks = np.round(2*pi/get_y_ticks, 2)
    ax.set_yticklabels(new_y_ticks);
    ax.set_xlabel(r'Temperature [$^{\circ}$C]')
    ax.set_ylabel(r'd [$\AA$]')
    #ax.set_title('RTA map', y=1.02)
    
    # visualizing and saving the plot 
    out_plot_title = filename.split('.')[0] + '_rtamap'
    ax.set_title(out_plot_title, y=1.03)
    fig.tight_layout()

    save_option = bool_sd(save_str)
    if save_option == True:
        # output image -- save it as png or pdf
        output_image_name = out_plot_title + '.png'
        plt.savefig(output_image_name,bbox_inches='tight')
        print "RTA map saved in present working directory as: %s\n" % (output_image_name)
    
    #plt.close(fig)
    #plt.show() 
    return
#
#==============================================================================
def rtapylysis(*args):
    '''
    Usage: rtapylysis(*args)
    -----
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
    ###################################################################################################
    fig = plt.figure(facecolor='w', figsize=[14.0,18.0])
    sns.set_style("ticks")    # using seaborn (http://web.stanford.edu/~mwaskom/software/seaborn/)
    sns.set_context("talk", font_scale=1.8, rc={"lines.linewidth": 2.0})
    ####################################################################################################
    # input START ---
    inv = args[0]
    invk=0;  projectdirin, filename = inv[invk:invk+2]
    invk=2;  xmin, xmax, dmin, dmax = map(float, inv[invk:invk+4])
    invk=6;  t_in = map(float, inv[invk:invk+3]) # temperature for XRD plots
    invk=9;  t_in_label = inv[invk:invk+3] # label for XRD plots
    invk=12; phase = inv[invk:invk+3] # dhkl file names 
    invk=15; phasetxt = inv[invk:invk+3] # dhkl label text
    invk=18; mphvalin, mpdvalin = map(float, inv[invk:invk+2]) # for detection of phase transition temp.
    invk=20; save_option  = bool_sd(inv[invk]) # 'True' or 'False' 
    invk=21; smooth_XRD = bool_sd(inv[invk]) # 'True' or 'False'; for the 3 XRD plots
    invk=22; SG_smooth_and_bkg_subtr_in_image= bool_sd(inv[invk]) # 'True' or 'False'
    invk=23; smoothing_option = bool_sd(inv[invk]) # 'True' or 'False'; for peak area & rate-area plots
    invk=24; color_choice = int(inv[invk]) # chose the colormap for col_dict in myfunctions
    invk=25; color_map_control = inv[invk] # 'man' or 'auto'
    invk=26; man_col_min, man_col_max = map(float, inv[invk:invk+2])
    invk=28;
    # create dminArea and dmaxArea ranges -- multiple ranges
    dminArea = {}
    dmaxArea = {}
    sub_inv = inv[invk:len(inv)]
    _len = len(sub_inv)
    if _len%2 == 0: # for even
        #print _len
        for idd in range( 0, int(_len/2) ):
            #print i
            dminArea[idd] = float(sub_inv[ 2*idd ])
            dmaxArea[idd] = float(sub_inv[ 2*idd +1 ])
    # input END ---
    ####################################################################################################
        projectdir = projectdirin # get the project directory path as a system variable
        filepath = os.path.join(projectdir, filename)
    
        matfile = scipy.io.loadmat(filepath)
        t = matfile['Temp'][:,0] #converting Temperature to 1D array    
        q = matfile['q'][0,:] #converting q vector to 1D array
        data = matfile['intensity'] # 2D array; 
        #time = matfile['time'][:,0] #converting time to 1D array
        #tth = matfile['tth'][0,:] #converting 2-theta to 1D array
        #ekev = matfile['ekev'][0,0] # getting the x-ray energy in keV
    
        # find the indices near user given values for abscissa and ordinate 
        id_x_min = np.argmin(np.abs(t - xmin))
        id_x_max = np.argmin(np.abs(t - xmax))
        #val_near_target_x_min = t[id_x_min]
        #val_near_target_x_max = t[id_x_max]
        
        target_y_min = 2*pi/float(dmax)
        target_y_max = 2*pi/float(dmin)
        id_y_min = np.argmin(np.abs(q - target_y_min))
        id_y_max = np.argmin(np.abs(q - target_y_max))
        #val_near_target_y_min = q[id_y_min]
        #val_near_target_y_max = q[id_y_max]
        
        # for peak area plots
        target_y_min2 = {}
        target_y_max2 = {}
        id_y_min2 = {}
        id_y_max2 = {}
        val_near_target_y_min2 = {}
        val_near_target_y_max2 = {}
        for i in range( 0, len(dminArea) ): # using one of dminArea or dmaxArea for this loop
            #print i
            target_y_min2[i] = 2*pi/dmaxArea[i]
            target_y_max2[i] = 2*pi/dminArea[i]
            id_y_min2[i] = np.argmin(np.abs(q - target_y_min2[i]))
            id_y_max2[i] = np.argmin(np.abs(q - target_y_max2[i]))
            val_near_target_y_min2[i] = q[id_y_min2[i]]
            val_near_target_y_max2[i] = q[id_y_max2[i]]
        
        # select parts of the matfile variables for plotting
        # x-range or Temperature range here
        x = t[id_x_min:id_x_max]
        
        # for complete RTA plot in specified x-range
        y = q
        z = data[id_x_min:id_x_max, : ]
        
        # for bkg subtracted RTA plots
        y_bkg = q[id_y_min:id_y_max]
        z_bkg = data[id_x_min:id_x_max, id_y_min:id_y_max]
        
        
        # for calculating the peak areas & plotting them in a single plot
        y_bkg2 = {}
        z_bkg2 = {}
        for i in range( 0, len(id_y_min2) ): # using one of id_y_min2 or id_y_max2 for this loop
            y_bkg2[i] = q[id_y_min2[i]:id_y_max2[i]]
            z_bkg2[i] = data[id_x_min:id_x_max, id_y_min2[i]:id_y_max2[i]]
            
        # names of the dhkl files 
        dhkl =  [
                    {'f': phase[0], 'color': 'r', 'id': phasetxt[0]},\
                    {'f': phase[1], 'color': 'g', 'id': phasetxt[1]},\
                    {'f': phase[2], 'color': 'b', 'id': phasetxt[2]},\
                ]        
    ####################################################################################################
        # create subplots
        ax = {}
        for i in range(1,7+1):
            #print i
            ax[i] = plt.subplot(4,2,i)
    
        #'''
        # Full RTA for the chosen x-range, no BKG subtraction
        if color_map_control == 'auto':        
            int_min = np.log10(np.min(np.min(z.T)))
            int_max = np.log10(np.max(z.T))
            print "\nRTA: Full range: cmap: auto min and max -> %.2f and %.2f" % (int_min, int_max)
        else:
            int_min = float(man_col_min) 
            int_max = float(man_col_max) 
            print "\nRTA: Full range: cmap: man min and max -> %.2f and %.2f" % (int_min, int_max)
            
        i=1
        ax[i].pcolormesh(x,y,log10(z.T), cmap=col_dict[color_choice], vmin=int_min, vmax=int_max)
        ax[i].axis('tight');
        
        
        # put inverse levels (i.e., d=2*pi/q instead of q, for the y-axis
        get_y_ticks = np.array(ax[i].get_yticks().tolist())
        new_y_ticks = np.round(2*pi/get_y_ticks, 2)
        ax[i].set_yticklabels(new_y_ticks);
        ax[i].set_xlabel(r'Temperature [$^{\circ}$C]')
        ax[i].set_ylabel(r'd [$\AA$]')
        ax[i].set_title('RTA: Full range', y=1.02)
        #''' 
    ####################################################################################################
        #'''
        # XRD Plots from 3 regions to be put in three subplots
        id_xrd_plot = [3, 5, 7]  # for chosing the subplots ---
        colors = cm.gnuplot(np.linspace(0, 0.5, len(t_in)))
        for k, j, c in zip(id_xrd_plot, range(0, len(t_in)), colors):
            idval = np.argmin( np.abs( t - t_in[j] ) )
            y_plot = data[idval,:]
            
            if smooth_XRD:
                y_plot = savitzky_golay(y_plot, 5, 1)
                print 'Smoothing done for the XRD plot in subplot ', k
                
            #ax[k].plot( q, log10(y_plot), color=c, label = '%.f $^{\circ}$C; \n%s' % ( t[idval], t_in_label[j] ) )
            ax[k].plot( q, y_plot, color=c, label = '%.f $^{\circ}$C; \n%s' % ( t[idval], t_in_label[j] ) )
            ax[k].legend(loc='best', fontsize = 20, shadow=True)
            
            ax[k].axis('tight');
            get_x_ticks = np.array(ax[k].get_xticks().tolist())
            new_x_ticks = np.round(2*pi/get_x_ticks,1)
            ax[k].set_xticklabels(new_x_ticks);
            ax[k].set_yscale('log')
            #ax[k].ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
            ax[k].set_xlabel(r'd [$\AA$]')
            ax[k].set_ylabel(r'Norm. I')
    
            #==============================================================================
            # annotate with phase name and hkl information
            if k == id_xrd_plot[0]:
                _a = 0
            if k == id_xrd_plot[1]:
                _a = 1
            if k == id_xrd_plot[2]:
                _a = 2       
            qlist = get_qhkl( dhkl[_a]['f'], dmin, dmax-0.5 )
            draw_vert_qlines( ax[k], qlist, c )
            trans = ax[k].get_xaxis_transform()
            for j in range( 0, len(qlist) ):
                ax[k].annotate( dhkl[_a]['id'] + '($' + qlist[j]['label']  + '$)', xy=( float(qlist[j]['q']), 0.9 ), xycoords=trans, rotation=90, size = 15, color=c )
            #==============================================================================
            
        #'''
        ####################################################################################################
        #'''
        # BKG subtracted partial RTA
        i=2    
        # for SG smooting and background subtraction in the image itself, this section is used
        #z_bkg22 = z_bkg.copy() # this is used in case the plot is shown without bkg subtraction but peak area is calculated after bkg.... this is a temprary solution....
        if SG_smooth_and_bkg_subtr_in_image:
            for idx in range( 0, shape(z_bkg.T)[1] ):
                #print idx
                z2 = z_bkg.T[:, idx]
                z2_sg = savitzky_golay(z2, 5, 1) 
                base = baseline(z2_sg, 2)[0]
                z2 = z2_sg - base
                z_bkg.T[:, idx] = z2 + abs(z2.min()) + 1. #10.
                #print 'Smoothing & Bkg. subtr. done for 2nd RTA plot'
    
        color_map_control2 = 'auto'
        if color_map_control2 == 'auto':        
            int_min2 = np.mean( [np.min(z_bkg.T[np.where(z_bkg.T > 0)]), np.median(z_bkg.T > 0)] )
            int_max2 = (np.max(z_bkg.T))
            print "\nRTA: partial: cmap: auto min and max -> %.2f and %.2f" % (int_min2, int_max2)
        else:
            int_min2 = float(man_col_min)
            int_max2 = float(man_col_max)
            print "\ncRTA: partial: cmap: man min and max -> %.2f and %.2f" % (int_min2, int_max2)
    
        
        ax[i].pcolormesh(x,y_bkg, z_bkg.T, cmap=col_dict[color_choice], vmin=int_min2, vmax=int_max2)    
        #ax[i].pcolormesh(x,y_bkg, log10(z_bkg22.T), cmap=col_dict[color_choice], vmin=int_min, vmax=int_max)     
        ax[i].axis('tight');
        
        # put inverse levels (i.e., d=2*pi/q instead of q, for the y-axis
        get_y_ticks = np.array(ax[i].get_yticks().tolist())
        new_y_ticks = np.round(2*pi/get_y_ticks, 2)
        ax[i].set_yticklabels(new_y_ticks);
        ax[i].set_xlabel(r'Temperature [$^{\circ}$C]')
        ax[i].set_ylabel(r'd [$\AA$]')
        ax[i].set_title(r'RTA: Partial', y=1.02)
        #'''
        ####################################################################################################
        #'''
        # Peak area plots below the BKG subtracted partial RTA
        i=4;
        peak_area = {}
        peak_name = {}
        for j in range( 0, len(z_bkg2) ): 
            peak_area[j] = []
            peak_name[j] = 'P_' + str(j)
            for idx in range( 0, shape(z_bkg2[j].T)[1] ):
                peak_area[j].append(trapz(y = z_bkg2[j].T[:, idx], x = y_bkg2[j]))
            
            peak_area[j] = np.array(peak_area[j])
            
            #smoothing
            if smoothing_option:
                peak_area[j] = savitzky_golay(peak_area[j], 5, 1)
                print 'Smoothing done for Peak Area: ', peak_name[j]
        
        plot_colors = cm.spectral(np.linspace(0., 0.5, len(peak_area)))
        for k in range( 0, len(peak_area) ):
            ax[i].plot(x, peak_area[k], color=plot_colors[k])
    
        trans = ax[i].get_xaxis_transform()
        annot_y_gap = 1.0/(len(peak_area) - 1); # dividing total y-height
        for k in range( 0, len(peak_area) ):
            ax[i].annotate( peak_name[k], xy=(xmax, (len(peak_area) - k - 1)*annot_y_gap), xycoords=trans, annotation_clip=False, size = 20, color=plot_colors[k] )        
        
        #ax[i].set_yscale('log')
        ax[i].axis('tight');
        #ax[i].set_xlabel(r'Temperature [$^{\circ}$C]')
        ax[i].set_ylabel(r'P.A. [a.u.]')
        ####################################################################################################
        #'''
        # Derivative of Peak area plots below the Peak Area plots
        i=6;
        peak_area_grad_mag = {}
        for j in range( 0, len(peak_area) ):
            peak_area_grad_mag[j] = np.abs( np.gradient(peak_area[j], np.gradient(x)) )    
        
            #smoothing
            smoothing_option_PAG = smoothing_option
            if smoothing_option_PAG: # not using this at this point
                peak_area_grad_mag[j] = savitzky_golay(peak_area_grad_mag[j], 5, 1)
                print 'Smoothing done for Peak Area Gradient: ', peak_name[j]    
    
        for k in range( 0, len(peak_area_grad_mag) ):
            ax[i].plot(x, peak_area_grad_mag[k], color=plot_colors[k])
    
        trans = ax[i].get_xaxis_transform()
        annot_y_gap = 1.0/(len(peak_area_grad_mag) - 1); # dividing total y-height
        for k in range( 0, len(peak_area_grad_mag) ):
            ax[i].annotate( peak_name[k], xy=(xmax, (len(peak_area_grad_mag) - k - 1)*annot_y_gap), xycoords=trans, annotation_clip=False, size = 20, color=plot_colors[k] )        
        
        #ax[i].set_yscale('log')
        ax[i].axis('tight');
        ax[i].set_xlabel(r'Temperature [$^{\circ}$C]')
        ax[i].set_ylabel(r'$|\Delta$(P.A.)/$\Delta$T| [a.u./$^{\circ}$C]')
    
        #'''
        ####################################################################################################
        # detect transition temeprature by max in peak area changes
        #annotate text, drawing lines etc.
        # put vline on the subplots 
        #The parameters mphval, mpdval follow the convention of the Matlab function findpeaks.m
        # http://nbviewer.jupyter.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb
        mphval = mphvalin #0.05 
        mpdval = mpdvalin #10
        peak_area_max = peak_area_grad_mag[1];
        idn = detect_peaks( peak_area_max, mph=mphval, mpd=mpdval, show=False ) 
    
        #idn_id = range(len(idn))
        #print 'Peaks: idn = ', idn
        #print 'idn_id = ', idn_id
        print 'Transition points along x = ', np.around(x[idn],1)
        
        idv = [2,4,6,8]
        for k in idv:
            if k == 2:
                c = 'w'
            else:
                c = 'k'
            if k in (2,4,6):
            #if k == 4: 
                #print k
                for idx in range(0, len(idn) ):
                    #print idx
                    ax[k].axvline(x=float( x[idn[idx]] ), color=c, ls='--' ) #,linewidth=1.2)
    
        # put hline on the subplots to show range used for peak-area 
        k=2
        # for area plots
        for i in range( 0, len(val_near_target_y_min2) ): 
            ax[k].axhline(y=val_near_target_y_min2[i], color='w', ls='-', linewidth=1.0)
            ax[k].axhline(y=val_near_target_y_max2[i], color='w', ls='-', linewidth=1.0)  
            annot_y_postn2 = 0.5 * ( val_near_target_y_min2[i] + val_near_target_y_max2[i] )
            ax[k].annotate( peak_name[i], xy=(xmax, annot_y_postn2), xycoords='data', annotation_clip=False, size = 20, color=plot_colors[i] )
        
        # annotate text
        k=6
        trans = ax[k].get_xaxis_transform()
        for idx in range(0, len(idn) ):
            annotate_text = '%.0f$^{\circ}$C' % ( float( x[idn[idx]] ) )
            ax[k].annotate(annotate_text, xy=( float( x[idn[idx]] ), 1.17 ), xycoords=trans, rotation=45, size = 20, color='r')
    
    
#==============================================================================
#         # annotate textlabel for Si peak    
#         dSi = 1.64
#         annotate_text = 'Si'
#         k=1
#         trans = ax[k].get_yaxis_transform()
#         ax[k].annotate(annotate_text, xy=(0.02, 2*pi/dSi), xycoords=trans, rotation=0, size = 25, color='k')
#         
#         k=3
#         trans = ax[k].get_xaxis_transform()
#         ax[k].annotate(annotate_text, xy=(2*pi/dSi, 0.55 ), xycoords=trans, rotation=0, size = 25, color='k')
#==============================================================================
     
        # visualizing and saving the plot 
        out_plot_title = filename.split('.')[0] + '_rtapylysis'
        fig.tight_layout()
        fig.suptitle('%s' % (out_plot_title), y = 1.02, fontsize=30)
    
        if save_option == True:
            # output image -- save it as png or pdf
            output_image_name = out_plot_title + '.png'
            plt.savefig(output_image_name,bbox_inches='tight')
            print "RTA analysis saved in present working directory as: %s\n" % (output_image_name)
        
        #plt.close(fig)
        #plt.show()
    
    else:
        print 'Error in d-range; Give multiple d-spacing ranges in pairs and try again\n'
    return
#
#==============================================================================
def bool_sd(_str): # convert string to bool; e.g.,when a bool value is passed as sys.argv
    '''
    Usage: bool_sd(_string) returns:
        True  (type Bool) if _string is one of 'True', 'true', 'T', 't'
        False (type Bool) if _string is one of 'False', 'false', 'F', 'f'
        The string itself (type string) if _string does not belong to any of the above options. 
    '''
    if _str in ('True', 'true', 'T', 't'):
        _str = True
    if _str in ('False', 'false', 'F', 'f'):
        _str = False
    return _str 
# 
#==============================================================================
def get_qhkl( dhkl_file, xmin=0, xmax=10, dhkl_dir1=dhkl_dir ):
    '''
    Usage: get_qhkl( dhkl_file, xmin=0, xmax=10, dhkl_dir1=dhkl_dir )
        Function parameters:
        dhkl_file: 
            a space delimited text file with at least 4 columns with values of
            h, k, l, and d-spacing(in Angstrom; DECREASING with line numbers in the file). 
            The file must contain a space delimited header, say h k l d_in_Angstrom.
            This file may be obtained from ICSD or generated using Vesta or made by hand. 
            The (hkl) values and d-spacing values correspond to the particular 
            crystallographic phase under consideration.
        xmin: Lowest d-spacing value (in Angstrom) extracted from dhkl_file; defualt 0.
        xmax: Highest d-spacing value (in Angstrom) extracted from dhkl_file; defualt 10.
            The file is scanned for (hkl) values between xmin <= d-spacing <= xmax.
        dhkl_dir1: This is a directory. By default, it is set to dhkl_dir where 
            the dhkl_file is stored for indexing purposes. 
            This may be set to any other directory where dhkl_file is stored.
            It is best to keep all the dhkl_file files in a single directory and call as necessary.
            
        Returns: a list qlist where each element is a dictionary consisting of:
            q-value (type string; units Ang^-1) corresponding to the d-spacing (unit Ang)
            hkl string corresponding to the q-value. These values are used for indexing purposes.
            E.g., qlist[0]['q'] = q-value in Ang^-1
                  qlist[0]['label'] = hkl values corresponding to qlist[0]['q']
    '''
    # open the indexing file with h, k, l, d values in a space delimited text file containing one header
    indexfile_with_path = os.path.join(dhkl_dir1, dhkl_file)

    ###begin: load the indexfiles one by one
    filename = open(indexfile_with_path, 'r')
    filename.readline() # reading the 1st line as the header
    #initialize lists for h k l and d-spacing values
    h_lst = []
    k_lst = []
    l_lst = []
    h_lst_label = []
    k_lst_label = []
    l_lst_label = []            
    d_lst = []
    for line in filename:
        h_lst.append(line.split()[0])
        k_lst.append(line.split()[1])
        l_lst.append(line.split()[2])
        h_lst_label.append(overline_txt(line.split()[0]))
        k_lst_label.append(overline_txt(line.split()[1]))
        l_lst_label.append(overline_txt(line.split()[2]))
        d_lst.append(float(line.split()[3]))
    filename.close()    # close the open file before proceeding further
    ###individual indexfile loading completed    

    # convert the lists into numpy array 
    d_array = np.flipud(np.array(d_lst)) # d-arry was flipped to make it an ascending array instead of a descending one
    #h_array = np.array(h_lst)
    #k_array = np.array(k_lst)
    #l_array = np.array(l_lst)
    #h_array_label = np.array(h_lst_label)
    #k_array_label = np.array(k_lst_label)
    #l_array_label = np.array(l_lst_label)            
    #hkl_array = np.flipud(np.array([a + b + c for a, b, c in zip(h_lst, k_lst, l_lst)])) # flipped
    hkl_array_label = np.flipud(np.array([a + b + c for a, b, c in zip(h_lst_label, k_lst_label, l_lst_label)])) # flipped            
    # conversion to numpy array completed

    x_array = d_array
    #identify the minm and maxm indices needed for plotting the vertical lines
    id_x_min = np.argmin(np.abs(x_array - xmin))    
    id_x_max = np.argmin(np.abs(x_array - xmax))
    #print "Before: MIN(id, val)= (%d, %f), MAX(id, val)=(%d, %f)" %(id_x_min, x_array.item(id_x_min), id_x_max, x_array.item(id_x_max))
    #if x_array.item(id_x_min) >= xmax:
        #id_x_min = id_x_min + 1
    #if x_array.item(id_x_max) <= xmin:
        #id_x_max = id_x_max - 1
    if x_array.item(id_x_min) <= xmin:
        id_x_min = id_x_min + 1
    if x_array.item(id_x_max) >= xmax:
        id_x_max = id_x_max - 1
    #print "After: MIN(id, val)= (%d, %f), MAX(id, val)=(%d, %f)" %(id_x_min, x_array.item(id_x_min), id_x_max, x_array.item(id_x_max))

    dictval = {}    
    qlist = []
    for j in range( id_x_min, id_x_max + 1 ):
        dictval['q'] = 2*pi/x_array[j] # converting d to q
        dictval['label'] = hkl_array_label[j]
        qlist.append(dictval.copy())

    return qlist
#
#==============================================================================
def overline_txt(txt):
    '''
    This function is used for labelling the negative hkl indices.
    Usage: overline_txt(string) returns string or \overline{string} depending if the 
               numerical value of the string is positive or negative.
               
               E.g., overline_txt(-2.94) will return 2 with a bar on top; 
              This will be displayed in figure labels, figure captions etc.
              
    Caution: 
        Ideally, only integer (positive, 0, or negative) should be used with this function.
        The function converts the text into int, i.e, 
        using any of 2.984, 2.0, 2 will yield '\overline{2}'
    '''
    from math import fabs    
    val = int(txt)
    if val >= 0:
        txt = '%d' % (val)
    else:
        val = fabs(val)
        txt = '\overline{%d}' % (val)

    return txt
#
#==============================================================================
def draw_vert_qlines( pltobject, qlist, _color='k' ):
    '''
    Usage: draw_vert_qlines( pltobject, qlist, _color='k' )
    pltobject: a matplotlib object where vertical lines of _color will be drawn at
    x-positions defined in q_list and (hkl) indices defined in qlist.
    Here, qlist is a list which is the output of get_qhkl function such that 
            each element in qlist is a dictionary consisting of:
            q-value (type string; units Ang^-1) and hkl string corresponding to the q-value.
    The line color may be changed by changing _color to other values, e.g.,
    _color='r' for red, _color + 'g' for green etc.; any matplotlib color would work.
    
    The function returns the same pltobject with the vertical lines drawn at the q-values.
    '''
    for i in range(0, len(qlist) ):
        pltobject.axvline(x=float(qlist[i]['q']), color=_color, ls='--',linewidth=1.2)
    return pltobject
#
#==============================================================================
def baseline(y, deg=3, max_it=100, tol=1e-3):
    '''
    Notes: this function is copied from PeakUtils.
    http://pythonhosted.org/PeakUtils/    
    (PeakUtils targets Python 3.4 and depends on numpy, scipy, and optionally on matplotlib.) 
    *This module is written in python 2.7 where the baseline function could still be used.*
    
    Computes the baseline of a given data.

        Iteratively performs a polynomial fitting in the data to detect its
        baseline. At every iteration, the fitting weights on the regions with
        peaks is reduced to identify the baseline only.

        Parameters
        ----------
        y : ndarray
            Data to detect the baseline.
        deg : int
            Degree of the polynomial that will estimate the data baseline. A low
            degree may fail to detect all the baseline present, while a high
            degree may make the data too oscillatory, especially at the edges.
        max_it : int
            Maximum number of iterations to perform.
        tol : float
            Tolerance to use when comparing the difference between the current
            fit coefficient and the ones from the last iteration. The iteration
            procedure will stop when the difference between them is lower than
            *tol*.

        Returns
        -------
        ndarray
            Array with the baseline amplitude for every original point in *y*
    '''
    import math    
    order = deg+1
    coeffs = np.ones(order)

    # try to avoid numerical issues
    cond = math.pow(y.max(), 1./order)
    x = np.linspace(0., cond, y.size)
    base = y.copy()

    vander = np.vander(x, order)
    vander_pinv = LA.pinv2(vander)

    for _ in range(max_it):
        coeffs_new = np.dot(vander_pinv, y)

        if LA.norm(coeffs_new-coeffs) / LA.norm(coeffs) < tol:
            break

        coeffs = coeffs_new
        base = np.dot(vander, coeffs)
        y = np.minimum(y, base)

    return base, coeffs
#
#==============================================================================
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    '''
    Usage: savitzky_golay(y, window_size, order, deriv=0, rate=1)
    ----------
    http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
    ----------    
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688    
    
    '''
    
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
#
#==============================================================================


