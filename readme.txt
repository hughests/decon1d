The program decon1d consists of two files, decon1d.py and decon_set.py. Settings that can be changed to allow analysis of your specific
spectrum are found in decon_set.py, nothing should need to be changed in the file decon1d.py. For details about this software see Deconvolution of Complex 1D NMR Spectra Using Objective Model Selection
Hughes TS, Wilson HD, de Vera IMS, Kojetin DJ (2015) Deconvolution of Complex 1D NMR Spectra Using Objective Model Selection. PLoS ONE 10(8): e0134474. doi: 10.1371/journal.pone.0134474

In addition to python and it's standard library and commonly included modules (e.g. random, math, csv) the modules numpy, pylab, matplotlib, seaborn, scipy.signal are required along with the module lmfit. 
All of these modules (except lmfit) are included in the Enthought python distribution (https://www.enthought.com), which is likely the easiest method to get 
decon1d to work. The Enthought free version contains all these modules. The download and installation instructions for lmfit can be found at (http://lmfit.github.io/lmfit-py/).

After installing the required software create a folder and place decon1d.py, decon_set.py and your data in the folder. Adjust the variables
in decon_set.py as required for your data. Next open a terminal window and navigate to the folder that contains your data and type python (or your alias for python)
decon1d.py. A single pdf document should be created within the same folder that displays the deconvolution models.


