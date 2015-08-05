import math as m
###Fitting variables that are routinely changed for your specific data set
filenames=['S1 Data.txt']### This is where you should put the name of the file(s) you want to analyze
savemetime=[m.pi/4,m.pi/50] # max limit for phase (automatically sets min as -max) of fit peaks if less than or equal to pi/50 then no BIC penalty. Will run all the phases listed here.
peaks_to_fit=7 ## this is the number of peaks that are initially attempted to fit
freqsig=376.4731756 ##need to change this to whatever the precession frequency is of the atom you are analyzing in Mhz (e.g for fluorine on a '400' NMR it is ~376.5 for proton it would be ~400)
leftright=False #this applies to fitting if this equals False then all points are used in fit  if True fits only data between the limits specified in the left and right variables, remember this is NMR so x axis is reversed. 
left=20.0 #desired left limit of analysis
right=-20.0 #desired right limit of analysis
whichcol=1#if you have a multi-column text file (with separate data sets in each column) set this to the column you want to analyze. NOTE!! 0 is the first column
mc0123=True  #put True here to run montecarlo
num_of_mc_peaks=5 #number of monte carlo fits to run
maxFWHM=2000 ###this sets the maximum width of the fitted peaks in Hz
minFWHM=4 ##this sets the minimum width of the fitted peaks in Hz

########Graph output settings that are routinely changed for your specific data set
showgraphs=False ##if this =True then matplotlib interactive graphs will pop up, but if not only pdfs will be saved.
spectrum=True##if true this adds a color bar that shows color of various ppms.
mancentcol=False## This can be changed to a floating point number to manually set the center of the color spectrum. This is to manually set center of peak coloring, this should ideally be centered at left or right end of peaks
manrangecol=False # This can be changed to a floating point number to manually set the width of the peak coloring spectrum should be about half the width of the peak cluster in ppm
setyh=False #set to False for auto y scale
setyl=-0.003 #sets lower limit of y axis WHEN setyh IS SET TO FALSE.
tickfont=30 #set the font size of the tick labels (e.g. values of ppm and intensity)
labelfont=30 #set the font size of the x-y axis labels (ppm and intensity)
legendfont=12 #set the font size of the legends
xviewset=False ##determines whether left and right are controlled for the view graphed
leftview=20.0#desired left limit of graph (for view only)
rightview=-20.0#desired right limit of the graph (for view only)
showsmooth=True ##set this to True to show smoothed raw data, False shows raw data
colorblind=False## set to true to produce graphs better for colorblind people
boty=False #controls whether end axis label is shown
topy=True#controls whether end axis label is shown
leftx=False#controls whether end axis label is shown
rightx=False#controls whether end axis label is shown
fac=15 #number to multiply rms by to determine where to put the spectrum rainbow if the rainbow is too high or two low on the graph change this number.
plotresid=True ###whether residuals are plotted on graphs
fares='n'	
widthmin=2##FOR DISPLAY purposes ONLY in Hz this controls information display on the first graph (the one that shows ind peaks) for the fit peaks if peaks are less then area of these will not be taken into account
widthmax=2000##FOR DISPLAY purposes ONLY in Hz this controls information display  on the first graph (the one that shows ind peaks) for the fit peaks if peaks are greater then area of these will not be taken into account
linewidthgraph=4 ##set the linewidth of some of the graphs
inputpeaks=False #Whether there are input peaks of simulated data you want to graph
inputpeakloc='./run6_param_file.text'
show_title=True #this turns on or off the title on the graphs

###########################below here are advanced commands that in general don't need to be changed
showfitsphased=False##this determines whether peaks with 0 phase are shown on the second graph where just residual and data are usually found.
correct_baseline=True##this adjusts the baseline (zero order correction, that is just adjusts all values up or down by a constant)
maxsplit=10 #this sets the maximum number of split peaks that are attempted to fit over ten gets pretty slow 16 will last about a half hour on a average computer.
units=1 #(changes fitting equation depending on if x axis is in rad or Hz/ppm) 
split_peak=True ##set this to false if you don't want to run the split peak part of the process
Plimit=15 # Two BICs that differ less than this in value are considered equally good. If you put this value very high it will fit less peaks as adding more peaks won't improve the BIC more than a high value.
