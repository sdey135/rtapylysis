# rtapylysis
#### Updated to _version 0.3_: 2017-01-24
> [**_rtapylysis_**](https://github.com/sdey135/rtapylysis.git) is a data-analysis tool written in *Python 2.7.12* for fast analysis of big datasets collected by **synchrotron rapid thermal annealing (_RTA_) X-ray diffraction (_XRD_) experiments**. 
> This README file is intended for someone who may use the functions inside [**_rtapylysis_**](https://github.com/sdey135/rtapylysis.git) for analyzing similar kinds of synchrotron RTA data. [**_rtapylysis_**](https://github.com/sdey135/rtapylysis.git) has been tested on Linux, Mac, and Windows OS.

## Prerequisites
1. It is strongly recommended that the [Anaconda python](https://www.continuum.io/downloads/) is installed which would include all but one of the required python packages needed to use [**_rtapylysis_**](https://github.com/sdey135/rtapylysis.git).  
2. The python package [seaborn](http://seaborn.pydata.org) needs to be installed:

  * If [Anaconda python](https://www.continuum.io/downloads/) is installed, then 

  ```
$ conda install seaborn
```

  * If [pip](https://pip.pypa.io/en/stable/) is installed:

  ```
$ pip install --user seaborn
```

3. For a bare-bone version of _Python 2.7.12_, the following packages must be installed before running the script: `numpy`, `scipy`, `matplotlib`, `seaborn`, in addiiton to already installed python packages `os`, and `sys`.

  ```
$ pip install --user package_name # where, package_name: numpy, scipy, matplotlib, seaborn 
```

## Installation
1. It may be downloaded directly from the [**_github repository of rtapylysis_**](https://github.com/sdey135/rtapylysis.git) to any user directory called **_TARGET_**. No installation required.

2. Alternatively, [**_rtapylysis_**](https://github.com/sdey135/rtapylysis.git) may be installed using `pip`:

  ```
$ pip install --user rtapylysis
```

  * For a **Linux OS**, the above command would install [**_rtapylysis_**](https://github.com/sdey135/rtapylysis.git) in **_TARGET_** = `$HOME/.local/lib/pythonx.y/site-packages/` where _pythonx.y_ = _python2.7_ for _for _python 2.7.12_ (Use `python --version` for _python version_). 

  * For a **Mac OS**, see this [stackoverflow link](http://stackoverflow.com/questions/7143077/how-can-i-install-packages-in-my-home-folder-with-pip) for value of **_TARGET_**.

## Quickstart
* Change present working directory to **_TARGET_** by using the following command from terminal `cd $TARGET` (replace $TARGET with **_TARGET_** path). Then use the following commands to run the example scripts:

  * Obtain the RTA map for example.mat file contained inside **_$TARGET/example/_**
```
$ python run_rtamap.py 
```

  * Obtain the RTA analysis for example.mat file contained inside **_$TARGET/example/_**
```
$ python run_rtapylysis.py 
```
  Note: Any previous files with same name as **_example_rtamap.png_** and **_example_rtapylysis.png_** inside the present working directory will be overwritten without any prompt upon execution of the above python scripts.

## Description
[**_rtapylysis_**](https://github.com/sdey135/rtapylysis.git) `tree`: 
```
rtapylysis
├── detect_peaks.py
├── dhkl_dir
│   ├── phaseA.txt
│   ├── phaseB.txt
│   └── phaseC.txt
├── example
│   ├── example.mat
│   ├── example_rtamap.png
│   ├── example_rtapylysis.png
│   └── example_structure.txt
├── func_rtapylysis.py
├── run_rtamap.py
└── run_rtapylysis.py
```

`func_rtapylysis.py` contains the following functions: 
  1. rtamap
  2. rtapylysis
  3. bool_sd
  4. get_qhkl
  5. overline_txt
  6. draw_vert_qlines
  7. [baseline](https://bitbucket.org/lucashnegri/peakutils/src/cf22985c2cc1b4ea32cbd201e368593caf598d71/peakutils/baseline.py?at=master&fileviewer=file-view-default) (copied from [PeakUtils](http://pythonhosted.org/PeakUtils/))
  8. [savitzky_golay](http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay) (copied from [scipy](https://www.scipy.org/))

`detect_peaks.py`, written by [**Marcos Duarte**](https://github.com/demotu) contains the following functions:  
  1. detect_peaks
  2. _plot

Sub-directory _example_ contains a matfile **_example.mat_**, a text file **_example_structure.txt_**, and two image files **_example_rtamap.png_** and **_example_rtapylysis.png_** showing examples of typical data-analysis that could be achieved using `run_rtamap.py' and `run_rtapylysis.py' respectively on the RTA data present in **_example.mat_**. The file **_example_structure.txt_** contains information about the structure of the input datafile **_example.mat_**. The new matfiles from the user must be constructed exactly after **_example.mat_** following the **_example_structure.txt_** for the data-analysis to work smoothly. It is best to keep the new matfiles under the **_example_** sub-directory in which case so that `run_rtamap.py` and `run_rtapylysis.py` could be run with minimal changes.  

Sub-directory _dhkl_dir_ contains the three text files needed for phase identification, viz., **_phaseA.txt, phaseB.txt, and phaseC.txt_**. Presumably, the data in new matfiles will represent different phases and the user would also need to construct new indexing files and save them inside the **_dhkl_dir_** as text files of similar structure as one of the _dhkl_files_, say, after **_phaseA.txt_** which is a space delimited text file with at least 4 columns with values of h, k, l, and d-spacing(in Angstrom; DECREASING with line numbers in the file). The file must contain a space delimited header, say h k l d_in_Angstrom. This file may be obtained from [ICSD](https://www.fiz-karlsruhe.de/icsd.html) or generated using [Vesta](http://jp-minerals.org/vesta) or prepared manually. The (hkl) values and d-spacing values are specific to the particular crystallographic phase under consideration.

The functions `rtamap` and `rtapylysis` are the two functions the user would use for RTA XRD data analysis, rest of the functions are supporting the use of these two functions. 

In addition, two user-friendly python script files are included: `run_rtamap.py`, `run_rtapylysis.py`. These script-files are meant to be edited by the user for data-analysis. The files contain instructions which are self-explanatory. The script may be written on a single file or a series of files collected under same experimental conditions. 

The input data (e.g., example.mat) consists of a series of _Intensity vs. q-vector_ data collected at various temperatures while the sample is being heated up rapidly.   
The `rtamap` function takes this data as its input and creates a plot with Temperature as x-axis, d-spacing (calculated from q-vector) as y-axis, and Intensity values as the z-axis for initial data visualization (see `run_rtamap.py` for details).  

The `rtapylysis` function takes the same input data (e.g., example.mat) as its input and creates the first sub-plot with Temperature as x-axis, d-spacing (calculated from q-vector) as y-axis, and Intensity values as the z-axis. This 2D map is named “RTA: Full Range” during analysis.  

The user could select the sub-regions of temperature and d-spacing values to get the second sub-plot for which _on-the-fly polynomial background subtraction_ can be applied to all scans. This second sub-plot is named “RTA: Partial” in the final figure. The on-the-fly background subtraction step is particularly useful to isolate meaningful signals from the noise while studying thin and ultra-thin films as thin as 2-3 nanometers.  

The `rtapylysis` function also outputs 3 “Intensity vs. d-spacing” sub-plots at 3 temperature values chosen by the user, i.e., representing three different phases encountered with temperature. Assuming the user is able to supply three text files with information about the anticipated XRD peak locations for that particular phase (see `run_rtapylysis.py` for details), `rtapylysis` function will also annotate the 3 XRD sub-plots with the user-provided indexing information.  

The `rtapylysis` function integrates the area under the background-subtracted-peaks from d-spacing ranges given by the user and plots ALL of the integrated-peak-areas into a separate sub-plot with a peak identifier on the side.  
Finally, the `rtapylysis` function determines the transition temperatures from the gradient of the individual integrated-peak-areas and annotates the gradient-of-integrated-peak-areas sub-plot with these values.  

The `rtapylysis` function outputs a single publication quality figure combining all the seven sub-plots.  

The example output of `run_rtamap.py` and `run_rtapylysis.py` are included below.

## Example plots
  1. ### RTA map for initial data visualization  
  
    ![example_rtamap](https://cloud.githubusercontent.com/assets/20307497/22233758/8e445f52-e1c0-11e6-9e47-bec9bccada86.png)

  2. ### Analysis of the RTA map 
   
    ![example_rtapylysis](https://cloud.githubusercontent.com/assets/20307497/22233772/bd09c084-e1c0-11e6-8961-fc42a37656d0.png)

## Troubleshoot
  1. In case of issues, please make sure that [**_rtapylysis_**](https://github.com/sdey135/rtapylysis.git) is residing in the `PYTHONPATH` of your system.

## References
1. Dey, Sonal, Kai-Hung Yu, Steven Consiglio, Kandabara Tapily, Takahiro Hakamata, Cory S. Wajda, Gert J. Leusink, et al. “Atomic Layer Deposited Ultrathin Metal Nitride Barrier Layers for Ruthenium Interconnect Applications.” Journal of Vacuum Science & Technology A: Vacuum, Surfaces, and Films 35, no. 3 (May 2017): 03E109. https://doi.org/10.1116/1.4979709.
2. Savitzky, Abraham., and M. J. E. Golay. “Smoothing and Differentiation of Data by Simplified Least Squares Procedures.” Analytical Chemistry 36, no. 8 (July 1, 1964): 1627–39. https://doi.org/10.1021/ac60214a047.
3. Press, William H., Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. Numerical Recipes 3rd Edition: The Art of Scientific Computing. 3 edition. Cambridge, UK ; New York: Cambridge University Press, 2007.
4. [_detect_peaks_](https://github.com/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb), written by [**Marcos Duarte**](https://github.com/demotu), has been included with [**_rtapylysis_**](https://github.com/sdey135/rtapylysis.git).
5. [PeakUtils](http://pythonhosted.org/PeakUtils/).
6. [scipy](https://www.scipy.org/).

