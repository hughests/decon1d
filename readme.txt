I have changed several things about the program, including now condensing both the settings and compuational script into one file. decon2.0.py allows specfic named files to be deconvoluted and requires you to enter the current path to the directory you are working in, whereas decon2.0b.py (and 3.0+) differs only in that it will attempt to deconvolute all files in the working directory that match the specifiers within the '' on line this line: "for filename in glob.glob(os.path.join(path, '*')):"

For details about this software see Deconvolution of Complex 1D NMR Spectra Using Objective Model Selection
Hughes TS, Wilson HD, de Vera IMS, Kojetin DJ (2015) Deconvolution of Complex 1D NMR Spectra Using Objective Model Selection. PLoS ONE 10(8): e0134474. doi: 10.1371/journal.pone.0134474

In addition to python and it's standard library and commonly included modules (e.g. random, math, csv) the modules numpy, pylab, matplotlib, seaborn, scipy.signal are required along with the module lmfit. 
All of these modules (except lmfit) are included in the Enthought python distribution (https://www.enthought.com), which is likely the easiest method to get 
decon1d to work. The Enthought free version contains all these modules. The download and installation instructions for lmfit can be found at (http://lmfit.github.io/lmfit-py/).

After installing the required software create a folder and place decon1d.py and your data in the folder. Adjust the variables
on the top portion of decon1d.py as required for your data. Next open a terminal window and navigate to the folder that contains your data and type python (or your alias for python) and decon1d.py. A single pdf document should be created within the same folder that displays the deconvolution models.

S1 Data.txt contains example simulated data from the above mentioned paper.
