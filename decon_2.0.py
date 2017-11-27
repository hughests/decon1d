import numpy as n
import seaborn as sns
import sys
import random
import pylab as p
import math as m
from matplotlib.backends.backend_pdf import PdfPages
import scipy.signal as sig
import lmfit as lm
print 'lm version',lm.__version__
import csv
from matplotlib import mlab as ml
import matplotlib.ticker as mticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
#from matplotlib import font_manager as fontm
########here are some things that may need to be adjusted for your specific application.
filenames=['WT DMSO Mock_20_2x']
#savemetime=[0,1,2,3,4,5,6,7,8,9,10,11]
savemetime=[False,True]

for jikn in filenames:
	for thgiu in savemetime:
		showgraphs=False ##if this =True then matplotlib interactive graphs will pop up, but if not only pdfs will be saved.
		spectrum=True##if true this adds a color bar that shows color of various ppms.
		mancentcol=-84.5##False #84.75#85 # MUST be a floating point number! 67.0 This is to manually set center of peak coloring, this should ideally be centered at left or right end of peaks
		manrangecol=1.0 #3.0 # MUST be a floating point number! This is to manually set the width of the peak coloring spectrum should be about half the width of the peak cluster in ppm
		peaks_to_fit=8 ## this is the number of peaks that are initially attempted to fit
		sirc=-20
		fn=jikn+".txt"### This is where you should put the name of the file you want to analyze str(raw_input('Enter name of text file to analyze '))
		well_file=""+fn ##this is the path to the file you want to analyze
		print 'analyzing this file',well_file
		freqsig=658.8021882 ##need to change this to whatever the precession frequency is of the atom you are analyzing (e.g for fluorine on a '400' NMR it is 376.5 for proton in would be ~400)
		Plimit=15 # Two BICs that differ less than this in value are considered equally good. If you put this value very high it will fit less peaks as adding more peaks won't improve the BIC more than a high value.
		linewidthgraph=4 ##set the linewidth of some of the graphs
		inputpeaks=False #Whether there are input peaks of simulated data you want to graph
		inputpeakloc='./run6_param_file.text'
		show_title=True #this turns on or off the title on the graphs
		leftright=True #this applies to fitting if this = False then all points are used in fit  if True fits only data between the limits specified in the left and right variables, remember this is NMR so x axis is reversed. 
		left=-82.6 #desired left limit of analysis
		right=-84.6 #desired right limit of analysis
		
		setyh=False #set to False for auto y scale
		setyl=-0.003 #sets lower limit of y axis when setyh is a float (number) set setyh to False for auto y scale.
		tickfont=30 #set the font size of the tick labels (e.g. values of ppm and intensity)
		labelfont=30 #set the font size of the x-y axis labels (ppm and intensity)
		legendfont=12 #set the font size of the legends
		xviewset=True ##determines whether left and right are controlled for the view graphed
		leftview=-82.6 #desired left limit of graph (for view only)
		rightview=-84.6 #desired right limit of the graph (for view only)
		showsmooth=True ##set this to True to show smoothed 
		
		colorblind=False
		maxph=m.pi/50 #max limit for phase (automatically sets min as -max) of fit peaks if less than or equal to pi/50 then no BIC penalty
		boty=False #controls whether end axis label is shown
		topy=True#controls whether end axis label is shown
		leftx=False#controls whether end axis label is shown
		rightx=False#controls whether end axis label is shown
		num_of_mc_peaks=1 #number of monte carlo fits to run
		whichcol=0#if you have a multi-column text file (with separate data sets in each column) set this to the column you want to analyze. NOTE!! 0 is the first column
		fac=6 #number to multiply rms by to determine where to put the spectrum rainbow
		exhaust=False##set this to True if you want the most exhaustive method run for deleting unecessary peaks, however it often takes a lot longer to run. If set to False then a quicker close to as good method runs. 
		plotresid=True ###whether residuals are plotted on graphs
		fares='n'
		if plotresid==True:
			fares='r'
			
		savgol=True
		smoothwindow=1 #smoothing window in Hz
		widthmin=2##FOR DISPLAY purposes ONLY in Hz this controls information display on the first graph (the one that shows ind peaks) for the fit peaks if peaks are less then area of these will not be taken into account
		widthmax=1000##FOR DISPLAY purposes ONLY in Hz this controls information display  on the first graph (the one that shows ind peaks) for the fit peaks if peaks are greater then area of these will not be taken into account
		split_peak=thgiu
		fit_with_smooth=True
		###########################below here are advanced commands that in general don't need to be changed#########################################################################################################################
		min_or_not=0 ##set to 0 to exhaustively search all possible good peaks 
		mc0123=True  #put True here to run montecarlo
		mcperturb=1000 ##default should be 300 a random number up to mcperturb divided by 10 is multiplied by the standard error of each parameter and added to the fitted parameter and then the fit it tried again
		A_or_B='BIC'
		mod=1  #this changes the penalty for each new peak, 1 does nothing, >1 increases penalty, <1 decreases penalty
		maxFWHM=2000 ###this sets the maximum width of the fitted peaks
		maxFWHM=maxFWHM/(freqsig*2)
		minFWHM=4 ##this sets the minimum width of the fitted peaks.
		minFWHM=minFWHM/(freqsig*2)
		showfitsphased=False##this determines whether peaks with 0 phase are shown on the second graph where just residual and data are usually found.
		correct_baseline=True##this adjusts the baseline (zero order correction, that is just adjusts all values up or down by a constant)
		maxsplit=10 #this sets the maximum number of split peaks that are attempted to fit over ten gets pretty slow 16 will last about a half hour on a average computer.
		units=1 #(changes fitting equation depending on if x axis is in rad or Hz/ppm) 
		threedfig=False
		
		#########################################################
		OR=False #if this is set to true then the three settings below are manual (determined by the settings below) instead of dependent on the signal to noise value.
		fit_with_smoothOR=False## if OR is set to True then this determines whether the center and height are set using raw residual error (fit-raw) if this is anything but 'true', if true then the fit-smoothed data is used.
		split_peakOR=True ##if OR is set to True then this will determine whether peak splitting is attempted: once the initial fit is made each fitted peak is attempted to be split into two peaks and then refit and compared against the previous fit.
		lowsignoiseOR=True ##if OR is set to True, and lowsingnoiseOR is set to True then length of consecutive positive runs are used exclusively as the method to place new peaks: as opposed to alternating between runs and picking the highest residual to place the peak. When the highest value is used in low signal to noise data the program fits noise spikes quite a bit.
		###############################################################
		varyphi=True #whether phase is allowed to be fit
		minph=-maxph
		hfactor=1.1 ##this is the limit for the height of the fits (this number times the highest signal value)
		###########################################Start of program####################################################  
		class autovivification(dict): #autovivification taken from http://stackoverflow.com/questions/651794/whats-the-best-way-to-initialize-a-dict-of-dicts-in-python
			"""Implementation of perl's autovivification feature."""
			def __getitem__(self, item):
				try:
					return dict.__getitem__(self, item)
				except KeyError:
					value = self[item] = type(self)()
					return value
		
		########################################################################################################
		####### Here we extract the data from the text file and make the Raw1 array containing the Y data.######
		########################################################################################################
		def find_nearest(array,value):
			idx = (n.abs(array-value)).argmin()
			return idx
		Raw1 =[]
		with open(well_file, 'rb') as f:
			reader = csv.reader(f)
			counter7=0
		#    counter5=0
			for row in reader:
				counter7=counter7+1
				if counter7==4:
					startend=row[0].split()
					s=float(startend[3])
					e=float(startend[7])
				if counter7==6:
					points=row[0].split()     		
					dp=int(points[3])
				if counter7>10: ##get rid of header
					hipster=row[0].split()
					Raw1.append(float(hipster[whichcol]))
		Raw2=n.array(Raw1)
		print 's',s,'e',e ##s is the value on the left of the textfile line number four (the high value) and e is the value of the right (low value) because nmr x axis are weird
		specwidth=n.abs(s-e)
		
		print 'specwidth',specwidth
		print 'number of data points',dp
		step=specwidth/(dp-1) ##ppm per step this is the step size for the ppm array
		stepsperhz=1/(step*freqsig)
		print 'stepsperhz',stepsperhz
		ppm2=n.array([])
		current=s
		for i in range(dp):
			ppm2=n.append(ppm2,current)
			current-=step
		Onewindow_length=m.ceil(stepsperhz*smoothwindow)
		if Onewindow_length <=2:
			Onewindow_length=3
		if Onewindow_length % 2 == 0:
			Onewindow_length+=1
		print 'Onewindow_length',Onewindow_length ## this is the SMOOTHING WINDOW LENGTH for Raw data for display purposes
		Twowindow_length=Onewindow_length ## this is the Smoothing widow length for the residual for fit purposes
		
		limit=len(Raw2)/30 ## this is how far from the left on x axis the noise will be calculated to for estimating noise standard deviation
		Rawsmooth2=sig.savgol_filter(Raw2, Onewindow_length, 2, mode='mirror')
		print 'Rawsmooth2',Rawsmooth2
		print 'Rawsmooth2 2:25',Rawsmooth2[2:25]
		noisey=Rawsmooth2[3:limit]
		#print 'noisy',noisey
		noiseyright=Rawsmooth2[::-1]
		noiseyright=noiseyright[3:limit]
		baseco=n.concatenate((noisey,noiseyright))	
		mean=n.mean(baseco)
		print 'mean',mean
		
		if correct_baseline==True:
			Raw2=Raw2-mean
		#build the ppm array
		if leftright==True:
			leftindex=find_nearest(ppm2,left)
			rightindex=find_nearest(ppm2,right)
			print 'leftindex',leftindex
			print 'rightindex',rightindex
			Raw1=Raw2[leftindex:rightindex]
			ppm=ppm2[leftindex:rightindex]
			print 'lenRaw1',len(Raw1)
			print 'lenppm',len(ppm)
		else:
			Raw1=Raw2
			ppm=ppm2
		minc=ppm[-1]
		print 'minc',minc
		maxc=ppm[0]
		print 'maxc',maxc
		specwidth=n.abs(ppm[-1]-ppm[0])
		print 'specwidth2',specwidth
		print 'dp',dp
		print 'step',step
		middle=n.amax(ppm)-(specwidth/2)
		print 'middle', middle
		centcol=middle
		rangecol=specwidth/10
		print 'center',centcol,'range',rangecol
		if isinstance(mancentcol, float):
			centcol=mancentcol
		if isinstance(manrangecol, float):
			rangecol=manrangecol
		resolution=step*freqsig
		print 'resolution',resolution
		print 'minFWHM',minFWHM
		ppmlen=len(ppm)
		print 'ppmlen',ppmlen
		print 'limit', limit
		limit_in_ppm=round((s-(limit*step)),2)
		print 'limit_in_ppm',limit_in_ppm
		maxindex1=Raw1.argmax() #This finds the point of maximum height in the data and uses this as first guess for center of first peak
		print 'maxindex', maxindex1
		c1=ppm[maxindex1] #center of first peak=point of max height
		print 'c1 start',c1
		phpos=n.abs(n.amax(Raw1)) #guess of height of first peak is height of highest point
		phneg=n.abs(n.amin(Raw1))
		listph=[phpos,phneg]
		print 'phpos,phneg',listph
		ph1=max(listph)#guess of height of first peak is height of most extreme y point
		print 'ph1 start',ph1
		danh=int(rangecol/step)  ###rangecol is in ppm step is ppm per step sp this is range in steps
		color_all=sns.color_palette("husl", danh)
		if colorblind==True:
			color_all=sns.cubehelix_palette(danh)
		noisey2=Raw2[0:limit] ##Raw2 is the end of the full spectrum no matter if only a portion is fitted
		rms = n.sqrt(n.mean(noisey2**2))
		lowsignoise=False
		
		sigtonoise=(ph1/rms)
		mcperturb=int(mcperturb/sigtonoise)
		print 'signtonoise',sigtonoise
		if sigtonoise<65:
			lowsignoise=True
		if OR==True:
			lowsignoise=lowsignoiseOR
			split_peak=split_peakOR
			fit_with_smooth=fit_with_smoothOR
		print 'rms',rms
		integratearray=n.sum(Raw1*step)
		pw=integratearray/(ph1*n.pi) #guess of width of first peak is based on total area of signal and highest peak
		print 'pw',pw
		endvalue=(ph1*0.6)*((pw/2)**2/((pw/2)**2+(n.abs(s-e)/2)**2)) ##this gives the y value of a single peak of height ph1 with pw width at the edge of the spectrum
		print 'endvalue',endvalue
		noisestd=n.std(noisey2) ## here is the noise standard deviation
		if correct_baseline==True:
			Raw1=Raw1+endvalue ## add back what should be added for an estimated peak values at edge of spectra
		Rawsmooth=sig.savgol_filter(Raw1, Onewindow_length, 2, mode='mirror')
		hcw_names=['height1','center1','width1','phase1','height2','center2','width2','phase2','height3','center3','width3','phase3','height4','center4','width4','phase4','height5','center5','width5','phase5','height6','center6','width6','phase6','height7','center7','width7','phase7','height8','center8','width8','phase8','height9','center9','width9','phase9','height10','center10','width10','phase10','height11','center11','width11','phase11','height12','center12','width12','phase12','height13','center13','width13','phase13','height14','center14','width14','height15','center15','width15','height16','center16','width16','height17','center17','width17','height18','center18','width18','height19','center19','width19','height20','center20','width20','height21','center21','width21','height22','center22','width22','height23','center23','width23','height24','center24','width24']
		################################################################################################
		#################here are the functions#########################################################
		################################################################################################
		def avoidbounds(temp): ##this makes sure the parameters that are input are not = to or greater than the bounds because if they are equal weirdness ensues
			mrbubbs=0
			for item in range(len(temp)/4):
				if temp[mrbubbs]<=0.0:
					temp[mrbubbs]=phmax/1000
				if temp[mrbubbs]>=phmax:
					temp[mrbubbs]=phmax-(phmax/1000)
				if temp[mrbubbs+1]<=minc:
					temp[mrbubbs+1]=minc+(specwidth/10000)
				if temp[mrbubbs+1]>=maxc:
					temp[mrbubbs+1]=maxc-(specwidth/10000)
				if temp[mrbubbs+2]<=minFWHM:
					temp[mrbubbs+2]=minFWHM+(minFWHM/1000)
				if temp[mrbubbs+2]>=maxFWHM:
					temp[mrbubbs+2]=maxFWHM-(minFWHM/1000)
				if m.ceil(temp[mrbubbs+3]*1000)<=m.ceil(minph*1000):
					temp[mrbubbs+3]=minph+(n.abs(minph)/10)
				if m.ceil(temp[mrbubbs+3]*1000)>=m.ceil(maxph*1000):
					temp[mrbubbs+3]=maxph-(n.abs(maxph)/10)
				mrbubbs+=4
			return temp
		def lorentzDK(params, ppm, Raw1,hcw_names):##this function works with the curve_fit function to turn curve_fit guesses for peak height center widths into an array of height values that is the model based on the HCW values.
					dlen=len(params)
					dstep=dlen/4
					c=0
					d2=0
					for i in range(0,dstep):
						h1=params[hcw_names[c]].value
						c1=params[hcw_names[c+1]].value
						w1=params[hcw_names[c+2]].value
						phi=params[hcw_names[c+3]].value
						d2=d2+(h1*(m.cos(phi)*(w1**2/(w1**2+(ppm-c1)**2))-m.sin(phi)*((-w1*(ppm-c1))/(w1**2+(ppm-c1)**2)))) ##w1 is HWHM
						c=c+4	
					
					return (Raw1-d2)
		
		def fit(params, *args): #this function contains curve_fit which minimizes the difference between the data and the summaation of the lorentzian peaks that are specified by the Parameter_init dictionary of lists.
			try: ## there are three tries given before giving up on fitting something.
				result=lm.minimize(lorentzDK, params, args=(ppm,Raw1,hcw_names)) 
				return result 
			except RuntimeError:
				try:
					clen=len(params)
					print 'length of failed parameters'
					print 'params',params
					params[hcw_names[clen-3]].value=params[hcw_names[clen-3]].value+(random.uniform(-0.2,0.2)) ##this takes the last center try and changes it a small random amount to be used in the next attempt
					params[hcw_names[clen-2]].value=params[hcw_names[clen-2]]/(random.uniform(5,20)) ##this takes the last width and decreases the width by dividing by 5 to 20 for next try
					result=lm.minimize(lorentzDK, params, args=(ppm,Raw1,hcw_names)) 
					print '####################runtimeerror'
					return result
					
				except RuntimeError:
					try:
						params[hcw_names[clen-2]].value=pw*(random.uniform(2,5)) ## this is the third try and increases the width beyond the original attempt
						result=lm.minimize(lorentzDK, params, args=(ppm,Raw1,hcw_names)) 
						print '##############################runtimerror'
						return result
					except RuntimeError:
						print 'did not converge'
						pass
					except:
						print 'did not fit for some random reason try3 '
						why=' did not fit b/c? try3'
						pass
					
				except:
					print 'did not fit for some random reason try2 '
					why=' did not fit b/c? try2'
					pass
					
			except:
				print 'did not fit for some random reason try1 '
				why=' did not fit b/c?'
		
		def extract_ind_peaks(something2): ##this function makes a dictionary of arrays of the y values of individual fitted peaks (the composite peaks)
			dlen=len(something2)
			ind_peaks={}
			dstep=(dlen)/4
			c=0
			for i in range(0,dstep):
				ind_peaks[i]=something2[c]*(m.cos(something2[c+3])*(something2[c+2]**2/(something2[c+2]**2+(ppm-something2[c+1])**2))-m.sin(something2[c+3])*(-something2[c+2]*(ppm-something2[c+1]))/(something2[c+2]**2+(ppm-something2[c+1])**2))
				c=c+4
			return ind_peaks
			
		def extract_ind_peaks_b(something2): ##this function makes a list of arrays of the y values of individual fitted peaks (the composite peaks)
			dlen=len(something2)
			ind_peaks=[]
			dstep=(dlen)/4
			c=0
			for i in range(0,dstep):
				ind_peaks.append(something2[c]*(m.cos(something2[c+3])*(something2[c+2]**2/(something2[c+2]**2+(ppm-something2[c+1])**2))-m.sin(something2[c+3])*(-something2[c+2]*(ppm-something2[c+1]))/(something2[c+2]**2+(ppm-something2[c+1])**2)))
				c=c+4
			return ind_peaks
			
		def lorentzDKp(ppm, dataArray): ##this is basically the same as lorentzDK it just can handle a dictionary of arrays input instead of a dictionary of lists which lorentzDK uses
					dlen=len(dataArray)
					dstep=dlen/4
					#print 'dstep p',dstep
					c=0
					d2=0
					for i in range(0,dstep):
						h1=dataArray[c]
						c1=dataArray[c+1]
						w1=dataArray[c+2]
						phi=dataArray[c+3]
						d2=d2+(h1*(m.cos(phi)*(w1**2/(w1**2+(ppm-c1)**2))-m.sin(phi)*((-w1*(ppm-c1))/(w1**2+(ppm-c1)**2)))) ##w1 is HWHM
						c=c+4
					return d2
					
		def lorentzDKz(dataArray): ##this takes and array of the heights, centers and widths (hcwarray) and seperates out into a list of three dictionarys that contain either heights centers or widths.
					wi=[]
					ce=[]
					he=[]
					pha=[]
					dlen=len(dataArray)
					dstep=(dlen)/4
					#print 'dstep',dstep
					c=0
					for i in range(0,dstep):
						h1=dataArray[c]
						c1=dataArray[c+1]
						w1=dataArray[c+2]
						p1=dataArray[c+3]
						c=c+4
						wi.append(w1)
						ce.append(c1)
						he.append(h1)
						pha.append(p1)
					hcw=[he,ce,wi,pha]
					return hcw
		
		def remove_peaks(dataArray): ## this removes each individual peak from a dictionary of arrays 
					dlen=len(dataArray)
					dstep=dlen/4
					c=dlen-1
					sub_a_peak={}
					for i in range(0,dstep):
						sub_a_peak[i]=n.delete(dataArray,[c-3,c-2,c-1,c])
						c=c-4
					return sub_a_peak
					
					
		def add_peaks(dataArray): ## this adds a peak to an hcw list 
					dlen=len(dataArray)
					dataArray=n.asarray(dataArray)
					dstep=dlen/4
					c=0
					add_a_peak=[]
					for i in range(0,dstep):
						# if dataArray[c+2]<0.01: ##don't split really thin peaks
# 							add_a_peak[i]=dataArray
# 							c=c+4
# 							continue
						#else:
						#print 'dataArray[c]',dataArray[c]
						split_h=0.7*dataArray[c]
						split_c_1=(dataArray[c+1])+((dataArray[c+2])/4)
						split_c_2=(dataArray[c+1])-((dataArray[c+2])/4)
						split_w=0.7*dataArray[c+2]
						split_p=0
						deleteadd_peak=dataArray
						#print 'deleteadd_peak',deleteadd_peak
						deleteadd_peak=n.append(deleteadd_peak,[split_h,split_c_1,split_w,split_p,split_h,split_c_2,split_w,split_p])
						#print 'deleteadd_peak',deleteadd_peak
						deleteadd_peak=n.delete(deleteadd_peak,[c,c+1,c+2,c+3])
						#print 'deleteadd_peak',deleteadd_peak
						add_a_peak.append(deleteadd_peak)
						c=c+4
					return add_a_peak
		def split_widest_peak(dataArray): ## this adds a peak to an hcw list 
					dlen=len(dataArray)
					dstep=dlen/4
					c=0
					b=0
					#print'dataArray',dataArray
					widths=n.array([])
					add_a_peak={}
					for i in range(0,dstep):
						widths=n.append(widths,dataArray[c+3])
						c+=4
					#print'widths',widths
					maxi=n.argmax(widths)
					b=2+(maxi*4)
					#print 'index of maximum width',b
						
					if dataArray[b]<0.01: ##don't split really thin peaks
						add_a_peak[0]=dataArray
					else:
						#print'gothere'
						split_h=0.7*dataArray[b-2]
						split_c_1=(dataArray[b-1])+((dataArray[b])/4)
						split_c_2=(dataArray[b-1])-((dataArray[b])/4)
						split_w=0.7*dataArray[b]
						split_p=0
						print dataArray
						add_a_peak[0]=n.delete(dataArray,[b-2,b-1,b,b+1])
						print add_a_peak[0]
						add_a_peak[0]=n.insert(add_a_peak[0],[b-2,b-2,b-2,b-2,b-2,b-2,b-2,b-2],[split_h,split_c_1,split_w,split_p,split_h,split_c_2,split_w,split_p])
						print add_a_peak[0]
					return add_a_peak
			
		def extract_val(params):
			hcw=[]
			standard_err=[]
			c=0
			lenpa=len(params)
			for i in range (lenpa/4):
				hcw.append(params[hcw_names[c]].value)
				hcw.append(params[hcw_names[c+1]].value)
				hcw.append(params[hcw_names[c+2]].value)
				hcw.append(params[hcw_names[c+3]].value)
				#print 'height',params[hcw_names[c]].value
				#print 'phase',params[hcw_names[c+3]].value
				
				
				standard_err.append(params[hcw_names[c]].stderr)
				standard_err.append(params[hcw_names[c+1]].stderr)
				standard_err.append(params[hcw_names[c+2]].stderr)
				standard_err.append(params[hcw_names[c+3]].stderr)
				c+=4
			wrapped=[hcw,standard_err]
			return wrapped
		#####################################################################################################################	
		###below is the initial fitting where peaks are sequentially fit until the desired amount of peaks is reached########
		###and then later in the next section peaks not significantly contributing to the fit per BIC are deleted######################################
		#####################################################################################################################
		poptn={}
		pcovn={}
		numberb=0
		twolnL={}
		ind_peaks=[]
		resi_err={}
		resi_err_smooth={}
		resi_err_smooth_abs={}
		AICc={}
		BIC={}
		k=4
		phmax=ph1*hfactor
		stder=[]
		STDER=[]
		#### creating the paramter file for use in the fitting routine##########
		print "c1",c1
		print 'minc',minc
		print 'maxc',maxc
		params=lm.Parameters()
		params.add('height1', value=ph1, min=0.0, max=phmax)
		params.add('width1', value=pw, min=minFWHM, max=maxFWHM)
		params.add('center1', value=c1, min=minc, max=maxc)
		params.add('phase1', value=0, min=minph, max=maxph, vary=varyphi)
		#print 'params', params
		evenodd=1
		fittedpeaks=0
		for number in range(peaks_to_fit):
			currentnum=k
			print 'Attempting to fit %.f peaks'%(number+1)
			result=fit(params) ##this is the fitting routine see def
			pen=currentnum+1-(currentnum/4)
			#print 'pen',pen
			if maxph > 0.063:
				kphi=maxph/m.pi
				pen=pen+((currentnum/4)*kphi)
				#print 'pen+kphi',pen
			fitted=extract_val(params)  ##see def
			
				
			poptn[number]=fitted[0]
			#print 'fitted[0]',fitted[0]
			stder.append(fitted[1])
			print 'stder', stder[number]
			if result.success==False:
				print 'NOTE up: fit only up to %.f peaks'%number
				break
			if result.success==True:
				fittedpeaks+=1
			try:
				resi_err[number]=result.residual
				resi_err_smooth[number]=sig.savgol_filter(resi_err[number], Twowindow_length, 2, mode='mirror')
				resi_err_smooth[number]=resi_err_smooth[number]
				resi_err_smooth_abs[number]=n.absolute(resi_err_smooth[number])
				maxindexES1_smooth=resi_err_smooth[number].argmax()
				maxindexES1_smooth_abs=resi_err_smooth_abs[number].argmax()
				maxindexES1=resi_err[number].argmax()
				if fit_with_smooth==True:
					resi_err1=resi_err_smooth[number]
					c2=ppm[maxindexES1_smooth]
				else:
					resi_err1=resi_err[number]
					c2=ppm[maxindexES1]
				y21 = n.ma.masked_less(resi_err1,0) ### this returns an array of all the numbers greater than or equal to 0 with dashes in place of numbers that are less than zero.
				y22=n.ma.flatnotmasked_contiguous(y21) #### this returns a list of slices of the form [slice(0,7,none),slice(10,15,none)] this means 0 to 6 indices are greater than zero and 10-14 are also. 
				difference=n.array([])
				
				for item in y22:
					difference=n.append(difference,(item.indices(len(resi_err1))[1]-item.indices(len(resi_err1))[0])) ###adds the lengths of all the positive runs to the array called difference. this takes the right slice value and subtracts that left slice value for example 7-0=7 so this run is 7 long (actually 6 I think but they are all off by one and it is the relative values that matters so OK)
				longestrunstart_ind=y22[n.argmax(difference)].indices(len(resi_err1))[0] ###finds index of start of longest run in residual
				longestrunend_ind=y22[n.argmax(difference)].indices(len(resi_err1))[1] ###finds index of end of longest run in residual
				longestrunvalues=resi_err1[longestrunstart_ind:longestrunend_ind] ### array? of longest residual run
				runlength=len(longestrunvalues)
				max_run_index=n.argmax(longestrunvalues)+longestrunstart_ind+20 ###this is the index of the largest (highest) residual in the longest run
				max_run_value=n.amax(longestrunvalues) ###this is the value of the largest (highest) residual in the longest run
				if evenodd % 2 == 0:
					c2=ppm[max_run_index]####this is the new initial center value, which is the ppm of the maximum value within the longest run of + residuals.
				if lowsignoise==True:
					c2=ppm[max_run_index]####this is the new initial center value, which is the ppm of the maximum value within the longest run of + residuals.
				ph2=max_run_value  ####this is the new initial peak height value, which is the value of the largest (highest) residual in the longest run
				twolnL[number]=-ppmlen*m.log(result.chisqr/ppmlen) ###part of BIC/AICc value calculation despite it's name result.chisqr is not chisquare but is instead the sum of the squared residuals.
				AICc[number]=(2*pen-twolnL[number])+((2*pen*(k+1))/(ppmlen-pen-1))
				BIC[number]=(mod*pen*m.log(ppmlen))-twolnL[number]
				p_1_width=(runlength*step)/4 ###this sets the initial value of the width of the peak at 1/4 the width of the longest run.
				numberb+=4
				somelist=[ph2,c2,p_1_width,0]
				somelist=avoidbounds(somelist)
				#print 'c2',c2
				params.add(hcw_names[numberb], value=somelist[0], min=0.0, max=phmax)
				params.add(hcw_names[numberb+1], value=somelist[1], min=minc, max=maxc)
				params.add(hcw_names[numberb+2], value=somelist[2], min=minFWHM, max=maxFWHM)
				params.add(hcw_names[numberb+3], value=somelist[3], min=minph, max=maxph, vary=varyphi)
				evenodd+=1
				k+=4
			except:
				print 'NOTE: fit only up to %.f peaks'%number
				break
		
		
		BestBICindex=min(BIC, key=BIC.get)
		print 'BIC',BIC
		print 'BestBICindex',BestBICindex
		best_fitindex=BestBICindex
		BestAICcindex=min(AICc, key=AICc.get)
		BestBICvalue=BIC[BestBICindex]
		BestAICcvalue=AICc[BestAICcindex]
		BICarray=n.asarray(BIC.values())
		diffBICarray=n.diff(BICarray)
		print 'diffBICarray',diffBICarray
		the_index=n.where(diffBICarray<2)
		
		if the_index[0].size:
			use_this_fit_index=n.max(the_index[0])+1
			print 'use this fit',use_this_fit_index
		if min_or_not==0 and the_index[0].size:
			best_fitindex=use_this_fit_index
			print 'best_fitindex1',best_fitindex
		#print 'the index', the_index
		statsame=([])
		c34=0
		for item in BICarray:
			if item <=BestBICvalue+Plimit:
				statsame=n.append(statsame,c34)
			c34+=1
		print 'statsame',statsame			
		AICcarray=n.asarray(AICc.values())
		STDER=stder[best_fitindex]
		print 'firstSTDER',STDER
		bestHCWs=poptn[best_fitindex]
		#print 'bestHCWsfirst',bestHCWs
		number_of_peaks_start=best_fitindex+1
		bestpeak=lorentzDKp(ppm,bestHCWs)
		NP=n.arange(1,(numberb/4)+1)
		bestC=BIC[best_fitindex]
		bestAICc=AICc[best_fitindex]
		if A_or_B=='AICc':
			bestC=AICc[best_fitindex]	
		number_of_peaks_in_best_fit=best_fitindex+1
		print 'starting with %.f peaks'%number_of_peaks_in_best_fit
		hcw=lorentzDKz(bestHCWs)
		#############################################################################################################################################################
		##this next part gets rid of peaks that do not significantly add to the quality of fit as judged by BIC with the BIC limit set at start above same as limit for above code.
		##############################################################################################################################################################
		overall=autovivification()
		poptv={}
		stanerr={}
		C_compare=-1
		one_less_peak={}
		bestHCWsdict={}
		splitBIC=0
		set=0
		k=5
		best_refit_index3=[]
		
		
		###This next bit of code takes the fit with the lowest BIC value from the original fit and deletes peaks that don't significantly add to the fit.###This fitting routine is conservative as peaks must add more value than Plimit in order to be kept.
		if number_of_peaks_in_best_fit>0:
			print'now refitting and deleting peaks that do not signficantly contribute to fit. Process=1)delete peak 2)if new BIC is less than old BIC then peak needs to be deleted'
			Crefits=n.ones(40)*1E10
			overall[set]['stanerror'][0]=STDER
			overall[set]['peaks_refit'][0]=bestpeak
			overall[set]['hcwarray_refit'][0]=bestHCWs
			overall[set]['BIC_refit'][0]=bestC
			overall[set]['AICc_refit'][0]=bestAICc
			best_refit_index3.append(0)
			one_less_peak[0]=bestHCWs
			k=len(one_less_peak[0])
			#print'lenonelesspeak',k
			Crefits[(k/4)-1]=bestC
			#print'Crefits before',Crefits
			while C_compare<Plimit: 
				C_com=1
				set=set+1
				params=lm.Parameters()
				if set>=1:
					one_less_peak=remove_peaks(one_less_peak[best_refit_index3[set-1]])
				k=len(one_less_peak[0])
				carty=len(one_less_peak)-1
				cart=-1
				pen=k+1-(k/4)
				#print 'pen',pen
				if maxph>0.063:
					kphi=maxph/m.pi
					pen=pen+((k/4)*kphi)
					#print 'pen+kphi',pen	
				while C_com>-2 and cart<carty:
					cart+=1
					c=0
					#print 'onelesspeak[cart]before',one_less_peak[cart]
					one_less_peak[cart]=avoidbounds(one_less_peak[cart])
					#print 'onelesspeak[cart]after',one_less_peak[cart]
					for i in range(k/4):
						params.add(hcw_names[c], value=one_less_peak[cart][c], min=0.0, max=phmax)
						params.add(hcw_names[c+1], value=one_less_peak[cart][c+1], min=minc, max=maxc)
						params.add(hcw_names[c+2], value=one_less_peak[cart][c+2], min=minFWHM, max=maxFWHM)
						params.add(hcw_names[c+3], value=one_less_peak[cart][c+3], min=minph, max=maxph, vary=varyphi)
						c+=4
					print 'before fit %.i'%cart
					result2=fit(params)
					print 'after fit  %.i'%cart
					fitted_temp=extract_val(params)
					poptv[cart]=fitted_temp[0]
					stanerr[cart]=fitted_temp[1]
					print 'cart',cart
					print 'standard error',fitted_temp[1]
					if poptv[cart]==None:
						print 'somecart went wrong on refit # %.f'%cart
						break
					try:
						poptz=poptv[cart]
						overall[set]['stanerror'][cart]=stanerr[cart]
						overall[set]['peaks_refit'][cart]=lorentzDKp(ppm,poptz)
						overall[set]['hcwarray_refit'][cart]=poptz
						twolnL_refit=-ppmlen*m.log(result2.chisqr/ppmlen) ###part of BIC/AICc value calculation despite it's name result.chisqr is not chisquare but is instead the sum of the squared residuals.
						overall[set]['BIC_refit'][cart]=(mod*pen*m.log(ppmlen))-twolnL_refit
						overall[set]['AICc_refit'][cart]=(2*pen)-twolnL_refit+((2*pen*(pen+1))/(ppmlen-pen-1))
						if exhaust == False:
							compare=overall[set]['BIC_refit'][cart]-bestC
							print 'Compare BIC inside',compare
							if compare<-5:
								C_com=compare						
					except:
						print 'something went wrong lower on refit # %.f'%cart
						break
				overall[set]['BIC_refit_array']=n.asarray(overall[set]['BIC_refit'].values())
				overall[set]['AICc_refit_array']=n.asarray(overall[set]['AICc_refit'].values())
				if A_or_B=='BIC':
					best_refit_index3.append(n.argmin(overall[set]['BIC_refit_array']))##here are the indexes (cart) of the lowest BIC value fit for each set
				if A_or_B=='AICc':
					best_refit_index3.append(n.argmin(overall[set]['AICc_refit_array']))
				if A_or_B=='BIC':
					Crefits[(k/4)-1]=n.amin(overall[set]['BIC_refit_array'])	####the (k/4)-1 in the Crefits dictionary? key is so that the BIC value gets associated with 0 for 1 peak and 8 for 9 peaks etc.
					C_compare=n.amin(overall[set]['BIC_refit_array'])-bestC ####if this is less than zero it means it is a nominal improvement if it is less than Plimit then it is a "significant" improvemtn
					if set>0:
						print 'C_compare',C_compare
						print 'best standard error %i'%(k/4)
						print overall[set]['stanerror'][best_refit_index3[set]]
						if C_compare<Plimit:
							print '%.f peaks deleted'%set
						if C_compare<0:
							bestC=n.amin(overall[set]['BIC_refit_array'])
							if k==0:
								sys.exit('No significant peaks found')
					
				if A_or_B=='AICc':
					Crefits[(k/4)-1]=n.amin(overall[set]['AICc_refit_array'])
					C_compare=n.amin(overall[set]['AICc_refit_array'])-bestC
					if set>0:
						print 'C_compare',C_compare
						if C_compare<0:
							bestC=n.amin(overall[set]['BIC_refit_array'])
							print '%.f peaks deleted'%set
							if k==0:
								sys.exit('No significant peaks found')
					
			
			else:
				if set>0:
					yuio=n.int_([])
					refitmin_val=n.amin(Crefits)
					refitmin_index=n.argmin(Crefits)#this is the index of the best refit by BIC
					compare=Crefits-refitmin_val #take the array of the refit BIC values and subract the best BIC value
					print 'comparison after delete',compare
					for i in range (fittedpeaks):  ## pick out the indices of the fits that are statistically the same as best fit (less than Plimit away)
						if compare[i] <=Plimit: ###right here we are getting all the indices within Plimit of each other.
							yuio=n.append(yuio,int(i))
					print 'fits that are less than Plimit different from best fit',yuio
					split_index=refitmin_index#n.amin(yuio) ###picking out the split fit index within Plimit that has the least number of peaks
					bestset=best_fitindex-split_index ###best_fitindex is the index of the starting number of peaks with the best BIC from the initial fit before deleting any
					splitset=best_fitindex-split_index
					bestfitindexarray=n.ones_like(yuio)
					bestfitindexarray=bestfitindexarray*best_fitindex
					bestsets=bestfitindexarray-yuio
					#print 'bestsets',bestsets #these are the sets that are within Plimit of the best refit (which includes the best initial fit)
					split_refit_index=best_refit_index3[splitset]
					splithcw=overall[splitset]['hcwarray_refit'][split_refit_index]
					splitBIC=overall[splitset]['BIC_refit'][split_refit_index]
					best_refit_index=best_refit_index3[bestset]
					bestpeak=overall[bestset]['peaks_refit'][best_refit_index]	
					best_AICc_refit=overall[bestset]['AICc_refit'][best_refit_index]
					bestC=overall[bestset]['BIC_refit'][best_refit_index]
					gh=split_index
					bestHCWs=overall[bestset]['hcwarray_refit'][best_refit_index]
					STDER=overall[bestset]['stanerror'][best_refit_index]
					print 'bestCfinalyuio',bestC
					number_of_peaks_in_best_fit=gh+1
					best_fitindex=refitmin_index
					hcw=lorentzDKz(bestHCWs)##this is the shortened version
					if A_or_B=='AICc':
						bestC=best_AICc_refit
					peaksdelete=bestset	
					print 'Best fit is with %.f unecessary peak(s) deleted'%peaksdelete
					print 'The Best fit has %.f peaks'%number_of_peaks_in_best_fit
				elif set==1:
					print 'Did not delete any peaks, best unsplit fit is with %.f peaks'%number_of_peaks_in_best_fit		
		mc_hcw_list=[] ##good
		mc_BIC_list=[] ##good
		split_hcw_list=[] ##good
		split_STDER=[] ##good
		split_BIC_list=[] ##good
		#print 'best_refit_index3',best_refit_index3
		#print 'yuio', yuio
		statsameset=[]
		BICdel=[]
		HCWdel=[]
		serr=[]
		for item in yuio:
			statsameset.append(best_fitindex-item)
		statsameset=bestsets
		#print 'statsamesets', statsameset
		for item in statsameset:
			BICdel.append(overall[item]['BIC_refit'][best_refit_index3[item]])
			HCWdel.append(overall[item]['hcwarray_refit'][best_refit_index3[item]])	
			serr.append(overall[item]['stanerror'][best_refit_index3[item]])	
		print 'HCWdel',HCWdel
		#print 'BICdel',BICdel
		print 'serr',serr
		########################################################################################################################################################
		#####this definition (split_peaks), takes the best fit from above and trys to split each individual peak and refit and get a better fit by BIC##########
		########################################################################################################################################################
		def split_peaks(input_HCW,BICin):
			split=False
			next=True
			BICsplit=[]
			BIC_compare_split=-1
			BIC_compare_split2=-100
			one_more_peak={}
			one_more_peak=add_peaks(input_HCW)
			#print 'onemorepeak',one_more_peak
			set=-1
			k=0
			print 'now splitting each peak and if the new BIC+Plimit< unsplit BIC then the peak is kept and try to split more'
			while next==True and k<((maxsplit)*4):
				set=set+1
				next=False
				if set>=1:
					#one_more_peak={}
					one_more_peak=add_peaks(next_peak)
					#print 'one morepeak before shuffle',one_more_peak
					random.shuffle(one_more_peak) ###this is necessary and helps speed things up because without it you sequentially work through all of the peaks each time you add a new one whereas with it you can hit here and there peaks that need to be split at the start
					#print 'one morepeak after shuffle',one_more_peak
					
				tempBIClist=[]
				currentnum=len(one_more_peak[0])
				k=len(one_more_peak[0])
				pen=currentnum+1-(currentnum/4)
				#print 'pen',pen
				if maxph>0.063:
					kphi=maxph/m.pi
					pen=pen+((currentnum/4)*kphi)
					#print 'pen+kphi',pen
				hcwpsplit=[]
				counteri=-1
				for thing in one_more_peak:
					params=lm.Parameters()
					c=0
					#print 'onemorepeak thing before',thing
					thing=avoidbounds(thing)
					#print 'onemorepeak thing after',thing
					for i in range(k/4):
						#print 'c',c
						params.add(hcw_names[c], value=thing[c], min=0.0, max=phmax)
						params.add(hcw_names[c+1], value=thing[c+1], min=minc, max=maxc)
						params.add(hcw_names[c+2], value=thing[c+2], min=minFWHM, max=maxFWHM)
						params.add(hcw_names[c+3], value=thing[c+3], min=minph, max=maxph, vary=varyphi)
						c+=4
					result2=fit(params)	
					twolnL=-ppmlen*m.log(result2.chisqr/ppmlen) ###part of BIC/AICc value calculation despite it's name result.chisqr is not chisquare but is instead the sum of the squared residuals.
					BIC_split=(mod*(pen)*m.log(ppmlen))-twolnL ###the +1 on the k is emperical.
					BIC_compare_split=BIC_split-BICin
					print 'BIC_compare_split',BIC_compare_split
					print 'number of peaks=%.f'%(k/4)
					fitted_temp3=extract_val(params)
					hcwpnow=fitted_temp3[0]
					hcwpsplit.append(hcwpnow)
					tempBIClist.append(BIC_split)
					counteri+=1
					if BIC_compare_split<sirc:
						split=True
						fitted_temp=extract_val(params)
						if fitted_temp[0]==None:
							print 'something went wrong on refit split'
						else:
							split_hcw_list.append(fitted_temp[0])
							split_STDER.append(fitted_temp[1])
							split_BIC_list.append(BIC_split)	
						tempc=BICin-Plimit
						BICin=BIC_split
						next=True
						next_peak=fitted_temp[0]
						print 'Winner!!!!!!!!!!!number of peaks=%.f'%(k/4)	
						break
				BIC_compare_split2=min(tempBIClist)-BICin
				best_refit_index=tempBIClist.index(min(tempBIClist))
				next_peak=hcwpsplit[best_refit_index]
					
			return [split,BICin]
		
		def montecarlo(inputsplit):
			mcbic_best=[]
			mchcw_best=[]
			counter54=0
			successmcpeaks=0
			trigger=False
			delete_peak=False
			#print 'inputsplit',inputsplit
			currentnum=len(inputsplit)
			pen=currentnum+1-(currentnum/4)
			#print 'pen',pen
			if maxph>0.063:
				kphi=maxph/m.pi
				pen=pen+((currentnum/4)*kphi)
				#print 'pen+kphi',pen
			if STDER[0]==0:
				stder56=[]
				for be in range (16):
					stder56=stder56+stder[0]+stder[0]
				print 'standard error used instead for mc because normal is all zeros',stder56
			for x in range(num_of_mc_peaks):
				print 'mc peaks tried',x
				paramsmc=lm.Parameters()
				counter3=0
				temptry=[]
				randnum=random.randrange(0, mcperturb)
				#print randnum
				#print 'inputsplit',inputsplit
				if STDER[0]==0:
					absth=0
					for i in inputsplit:
						mu=(stder56[absth]/50)*randnum
						temptry.append(random.gauss(i, mu))
						absth=+1
						#print 'egggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggstderrorzero'
				else:	
					for i in inputsplit:
						mu=STDER[counter3]*randnum
						temptry.append(random.gauss(i, mu))
						counter3+=1
				c=0
				#print 'temptry',temptry
				
				#print 'minc',minc
				#print 'maxc',maxc
				
				avoidbounds(temptry)
				#print 'temptryafter',temptry
					
					
					
				pn=len(temptry)
				for i in range(pn/4):
						#print 'c',c
						paramsmc.add(hcw_names[c], value=temptry[c], min=0.0, max=phmax)
						paramsmc.add(hcw_names[c+1], value=temptry[c+1], min=minc, max=maxc)
						paramsmc.add(hcw_names[c+2], value=temptry[c+2], min=minFWHM, max=maxFWHM)
						paramsmc.add(hcw_names[c+3], value=temptry[c+3], min=minph, max=maxph, vary=varyphi)
						c+=4
				#print 'params before mc', paramsmc
				result3=fit(paramsmc)
				mctwolnL=-ppmlen*m.log(result3.chisqr/ppmlen) ###part of BIC/AICc value calculation despite it's name result.chisqr is not chisquare but is instead the sum of the squared residuals.
				mcBIC=(mod*(pen)*m.log(ppmlen))-mctwolnL
				comp2=mcBIC-bestC
				print 'mcBIC-bestC',comp2
				if comp2<Plimit:
					trigger=True
					counter54+=1
					print 'Number of succesful monte carlo peak sets=%.f'%counter54
					successmcpeaks=counter54
					#print 'mcBIC=%.f'%mcBIC
					fitted_temp2=extract_val(paramsmc)
					popte=fitted_temp2[0]
					#print 'fitted hcw',popte
					mc_hcw_list.append(popte)
					mc_BIC_list.append(mcBIC)
					hcwp.append(popte)
					print 'added to hcwp'
					bic.append(mcBIC)
					if comp2<2:
						mchcw_best.append(popte)
						mcbic_best.append(mcBIC)
							
			return [trigger,delete_peak,successmcpeaks,mchcw_best,mcbic_best]
		def separate(lister):
			all=[]
			if any(isinstance(el, list) for el in lister):
				list_o_lists=[]
				for row in lister:
					splitpeakh=[]
					splitpeakc=[]
					splitpeakw=[]
					splitpeakp=[]
					c=0
					len24=((len(row))/4)#this is the number of peaks in the set
					for i in range(len24):
						splitpeakh.append(row[c])
						splitpeakc.append(row[c+1])
						splitpeakw.append(row[c+2]*freqsig*2)
						splitpeakp.append(row[c+3])
						c+=4
					list_o_lists.extend([splitpeakh,splitpeakc,splitpeakw,splitpeakp])	
				return list_o_lists
			else:
				#print 'not list of lists'
				splitpeakh=[]
				splitpeakc=[]
				splitpeakw=[]
				splitpeakp=[]
				c=0
				len24=((len(lister))/4)#this is the number of peaks in the set
				for i in range(len24):
					splitpeakh.append(lister[c])
					splitpeakc.append(lister[c+1])
					splitpeakw.append(lister[c+2]*freqsig*2)
					splitpeakp.append(lister[c+3])
					c+=4
				return [splitpeakh,splitpeakc,splitpeakw,splitpeakp]
		
		def stringize(lister):
			allstr=[]
			if any(isinstance(el, list) for el in lister):
				for row in lister:
					str=[]
					for thing in row:
						str.append("%.3f" %thing)
					allstr.append(str)
				return allstr
			else:
				str=[]
				for thing in lister:
					str.append("%.3f" %thing) ###!!!!!!!!!!!!?????This controls number of digits
				return str
		def area(h_list_of_lists,w_list_of_lists):
			counter57=-1
			all_areas=[]
			all_pc_areas=[]
			all_pc_areas_str=[]
			FWHM=[]
			if any(isinstance(el, list) for el in h_list_of_lists):
				for row in h_list_of_lists:
					counter57+=1
					peakset_area=[]
					peakset_pc_area=[]
					fwhm_list=[]
					for i in range(len(row)):
						peakset_area.append(row[i]*w_list_of_lists[counter57][i]*n.pi)
					total_a=sum(peakset_area)
					for thing in peakset_area:
						item2=(thing/total_a)*100
						peakset_pc_area.append(item2)
					for i in range(len(row)):
						fwhm_list.append((peakset_pc_area[i]/100)*w_list_of_lists[counter57][i])
					FWHM.append(sum(fwhm_list))
						
					
					all_areas.append(peakset_area)
					all_pc_areas.append(peakset_pc_area)
					all_pc_areas_str.append(stringize(peakset_pc_area))
				return [all_areas,all_pc_areas,all_pc_areas_str,FWHM]
			else:
				counter57+=1
				peakset_area=[]
				peakset_pc_area=[]
				fwhm_list=[]
				for i in range(len(h_list_of_lists)):
					peakset_area.append(h_list_of_lists[i]*w_list_of_lists[i]*n.pi)
				total_a=sum(peakset_area)
				for thing in peakset_area:
					item2=(thing/total_a)*100
					peakset_pc_area.append(item2)
				for i in range(len(w_list_of_lists)):
					fwhm_list.append((peakset_pc_area[i]/100)*w_list_of_lists[i])
				FWHM=sum(fwhm_list)
				pc_areas_str=(stringize(peakset_pc_area))
				return [peakset_area,peakset_pc_area,pc_areas_str,FWHM]
		######################################################################################################################################################
		####This code uses the split peaks definition from above to try to improve the fit by splitting the peaks in the best fit from above##################		
		######################################################################################################################################################
		usesplit=False
		if split_peak==True:
			if number_of_peaks_in_best_fit==1:
				splithcw=bestHCWs
				splitBIC=BestBICvalue
			if set==0:
				splithcw=bestHCWs
				splitBIC=BestBICvalue
			print 'splithcw',splithcw
			print 'splitBIC',splitBIC
			worked=split_peaks(splithcw,splitBIC)
			print 'worked[0]',worked[0]
			if worked[0]==True:###the stuff above is legit but the stuff in this line and below may not be necessary
				split_minBIC=min(split_BIC_list)
				split_minBIC_index=split_BIC_list.index(split_minBIC)
				split_BIC_array=n.asarray(split_BIC_list)
				dif_BIC=split_BIC_array-split_minBIC
				counter45=-1
				for item in dif_BIC:
					counter45+=1
					if item<Plimit:
						BICdel.append(split_BIC_list[counter45])
						HCWdel.append(split_hcw_list[counter45])
						serr.append(split_STDER[counter45])
				minsplitcomp=split_minBIC-bestC
				if minsplitcomp<-0.01:
					bestC=split_minBIC
					usesplit=True
				print 'bestC',bestC
				#print 'usesplit',usesplit
		hcwp=[]
		bic=[]
		sterr=[]
		bestC=min(BICdel)
		clarkind=0
		for thing in BICdel:
			clarkwat=thing-bestC
			if clarkwat<Plimit:
				hcwp.append(HCWdel[clarkind])
				bic.append(thing)
				sterr.append(serr[clarkind])
			clarkind+=1	
		indexminc=bic.index(min(bic))	
####################################################################
####################### mc routine##################################
####################################################################
		improved=False
		if mc0123==True:
			#print 'bestHCWs',bestHCWs
			#print 'STDER',STDER
			STDER=sterr[indexminc]
			print 'split_STDER usedin mc',STDER
			inputsplit=hcwp[indexminc]
			#print 'inputsplit',inputsplit
			monte_carlo=montecarlo(inputsplit)
			print 'monte carlo ran!!!!!!!'
			#print 'monte_carlo',monte_carlo
			
				
			if monte_carlo[0]==True:
				#print 'mc_BIC_list',mc_BIC_list
				mc_peaks=[]
				#print 'mc_hcw_list just after fit', mc_hcw_list
				for item in mc_hcw_list:
					mc_peaks.append(lorentzDKp(ppm,item))
				minBIC=min(mc_BIC_list)
				maxBIC=max(mc_BIC_list)
				print 'minBIC',minBIC
				print 'maxBIC',maxBIC
				bestmchcwp=separate(monte_carlo[3])
				bestmcbic=separate(monte_carlo[4])
				print len(bestmchcwp)
				#print 'bestmchcwp',bestmchcwp
				# if monte_carlo[2]>0:
# 					jhnm=len(bestmchcwp[0])#number of peaks in fit
# 					
# 					bestmchcwph=[]
# 					bestmchcwpc=[]
# 					bestmchcwpw=[]
# 					bestmchcwpp=[]
# 					c=0
# 					for i in range(len(bestmchcwp)/4):
# 						for item in bestmchcwp[c]:
# 							bestmchcwph.append(item)
# 						for item in bestmchcwp[c+1]:
# 							bestmchcwpc.append(item)
# 						for item in bestmchcwp[c+2]:
# 							bestmchcwpw.append(item)
# 						for item in bestmchcwp[c+3]:
# 							bestmchcwpp.append(item)
# 						c+=4
# 					#print 'bestmchcwph',bestmchcwph
# 					#print 'bestmchcwpc',bestmchcwpc
# 					totalfits=len(bestmchcwpc)
# 					if monte_carlo[2]>1:
# 						bestmchcwpc, bestmchcwph, bestmchcwpw, bestmchcwpp = zip(*sorted(zip(bestmchcwpc,bestmchcwph,bestmchcwpw,bestmchcwpp)))#sort the data by the centers
# 					noepp=monte_carlo[2] #number of entries per peak
# 					print 'noepp',noepp
# 					mch=[]
# 					mcc=[]
# 					mcw=[]
# 					mcp=[]
# 					cghi=0
# 					for i in range(jhnm):
# 						start=noepp*cghi
# 						end=noepp*(cghi+1)
# 						mch.append(bestmchcwph[start:end])
# 						mcc.append(bestmchcwpc[start:end])
# 						mcw.append(bestmchcwpw[start:end])
# 						mcp.append(bestmchcwpp[start:end])
# 						cghi+=1
# 					#print 'mch',mch
# 					#print 'mcc',mcc
# 					#print 'mcw',mcw
# 					#print 'mcp',mcp
# 					mch=n.asarray(mch)
# 					mcc=n.asarray(mcc)
# 					mcw=n.asarray(mcw)
# 					mcp=n.asarray(mcp)
# 					mcwmed=[]
# 					for i in range(len(mcw)):
# 						mcwmed.append(n.median(mcw[i]))
# 					#print 'mcwmed',mcwmed
# 					# for j in range(len(mch)/noepp):
# 	# 					for i in range(len(mch[j])):
# 	# 						soyu=mcw[j][i]
# 	# 						if soyu>2*mcwmed:
# 	# 							mch=n.delete(mch,j,[i])
# 	# 							mcc=n.delete(mcc,j,[i])
# 	# 							mcw=n.delete(mcw,j,[i])
# 	# 							mcp=n.delete(mcp,j,[i])
# 	# 						elif soyu<0.5*mcwmed:
# 	# 							mch=n.delete(mch,j,[i])
# 	# 							mcc=n.delete(mcc,j,[i])
# 	# 							mcw=n.delete(mcw,j,[i])
# 	# 							mcp=n.delete(mcp,j,[i])
# 					mchmed=[]
# 					mccmed=[]
# 					mcwmed=[]
# 					mcpmed=[]
# 					mchstd=[]
# 					mccstd=[]
# 					mcwstd=[]
# 					mcpstd=[]
# 					huty=len(mch)
# 					for i in range(huty):
# 						mchmed.append(n.median(mch[i]))
# 						mchstd.append(n.std(mch[i]))
# 					for i in range(huty):
# 						mccmed.append(n.median(mcc[i]))
# 						mccstd.append(n.std(mcc[i]))
# 					for i in range(huty):
# 						mcwmed.append(n.median(mcw[i]))
# 						mcwstd.append(n.std(mcw[i]))
# 					for i in range(huty):
# 						mcpmed.append(n.median(mcp[i]))
# 						mcpstd.append(n.std(mcp[i]))
					#print 'mchmed',mchmed
					#print 'mchstd',mchstd
					#print 'mch',mch
					
							
						
					
					# print 'bestmchcwph',bestmchcwph
	# 				print 'bestmchcwpc',bestmchcwpc
	# 				print 'bestmchcwpw',bestmchcwpw
	# 				print 'bestmchcwpp',bestmchcwpp
			
				
				mcdif1=[]
				cal_mc_BIC_list=[]
				for item in mc_BIC_list:
					mcdif=item-bestC
					mcdif1.append(item-bestC)
					if mcdif<-2:
						vluw=-mcdif+0.5
						cal_mc_BIC_list.append(vluw)
						print 'sig improved BIC by mc',mcdif
						improved=True
					elif -2<=mcdif<2:
						cal_mc_BIC_list.append(2.5)
					else:
						cal_mc_BIC_list.append((2/(mcdif))+1.0)
				#print 'cal_mc_BIC_list',cal_mc_BIC_list
				#print 'mc_hcw_list',mc_hcw_list
				monte_carlo_indpeaks=[]
				for thing in mc_hcw_list:
					monte_carlo_indpeaks.append(extract_ind_peaks_b(thing))
				#print 'monte_carlo_indpeaks',monte_carlo_indpeaks
				mcpeakh=[]
				mcpeakc=[]
				mcpeakw=[]
				mcpeakp=[]
				for row in mc_hcw_list:
					c=0
					len24=((len(row))/4)#this is the number of peaks in the set
					for i in range(len24):
						mcpeakh.append(row[c])
						mcpeakc.append(row[c+1])
						mcpeakw.append(row[c+2]*376.5*2)
						mcpeakp.append(row[c+3])
						c+=4
				# couytu=-1
# 				couytu1=-1
# 				minmcdif1=min(mcdif1)
# 				for thing in mcdif1:
# 					couytu+=1
# 					rzc=minmcdif1-thing
# 					if rzc>-2 and couytu1==-1:
# 						hcwp.append(mc_hcw_list[couytu])
# 						bic.append(mc_BIC_list[couytu])	
# 						couytu1+=1
		indexminc=bic.index(min(bic))
		bestC=bic[indexminc]
		#print 'hcwp',hcwp
		if len(bic)>1:
			bic, hcwp= zip(*sorted(zip(bic,hcwp)))#sort the data by the centers
		#print 'hcwp',hcwp
		river=-1
		kjh1=-1
		kjh2=-1
		kjh3=-1
		kjh4=-1
		kjh5=-1
		hcwfilt=[]
		bicfilt=[]
		for item in bic:
			river+=1
			BICcompared=bic[river]-bestC
			if BICcompared <3 and kjh1==-1:
				hcwfilt.append(hcwp[river])
				bicfilt.append(bic[river])
				kjh1+=1
			if 3<BICcompared <6 and kjh2==-1:
				hcwfilt.append(hcwp[river])
				bicfilt.append(bic[river])
				kjh2+=1
			if 6<BICcompared <9 and kjh3==-1:
				hcwfilt.append(hcwp[river])
				bicfilt.append(bic[river])
				kjh3+=1
			if 9<BICcompared <12 and kjh4==-1:
				hcwfilt.append(hcwp[river])
				bicfilt.append(bic[river])
				kjh4+=1
			if 12<BICcompared <15 and kjh5==-1:
				hcwfilt.append(hcwp[river])
				bicfilt.append(bic[river])
				kjh5+=1
		#print 'hcwfilt before graphing',hcwfilt
		##################################################################################################################################
		### code below here is for making graphs of the data. #######################################################################
		##################################################################################################################################
		sns.set_style("white")
		sns.set_style("ticks")
		lenfilt=len(hcwfilt[0])/4
		lenfilts='%.i'%lenfilt
		phase='%.3fph'%maxph
		minorLocator1   = MultipleLocator(1.0)
		whichcols="%i"%whichcol
		howdo='nosplit'
		if split_peak==True:
			howdo='split'    
		pdf=PdfPages(fn+whichcols+'_'+phase+fares+lenfilts+howdo+'.pdf')
		###Figure 1 this figure is split and shows the BIC/AICc values vs number of peaks on bottom and the best unsplit fit on top#########
		p.figure( figsize=(15, 9))
		#p.figtext(0.1,0.94,BICarray,size='small')
		ax=p.subplot(311)
		ax.invert_xaxis()
		ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
		if leftright==True:
			ax.set_autoscalex_on(False)
			ax.set_xlim([left,right])
		else:
			ax.set_autoscalex_on(False)
			ax.set_xlim([s,e])
		splitbest='No'
		numpeaksbest=len(hcwp[indexminc])/4
		if usesplit==True:
			splitbest='Yes'
		p.title((fn,'Splitbest?',splitbest,'peak(s)',numpeaksbest,'BIC=%.2f'%bicfilt[0],'Width of data used in fit displayed'))
		curves1=('Raw data','Best fit','Residual')
		p.plot(ppm,Raw1,linewidth=1,color='k')
		sumpeak=lorentzDKp(ppm,hcwfilt[0])
		p.plot(ppm,sumpeak,linewidth=1,color='green')
		if plotresid==True:
			residualmain=Raw1-sumpeak
			#residualmain_smooth=sig.savgol_filter(residualmain, Twowindow_length, 2, mode='mirror')
			p.plot(ppm,residualmain,linewidth=1.5,color='b')
		statsame_ip=extract_ind_peaks(hcwfilt[0])
		for peak in statsame_ip:
			p.plot(ppm,statsame_ip[peak],linewidth=0.5)
		leg=p.legend(curves1,loc=1)
		for t in leg.get_texts():
			t.set_fontsize(legendfont)
		p.setp(ax.get_xticklabels(), fontsize=tickfont)
		p.setp(ax.get_yticklabels(), fontsize=tickfont)
		
		 
		sns.despine()		
		p.locator_params(nbins=4)
		
		#p.setp(ax.get_xticklabels()[0], visible=leftx) ##controls left x axis label 
		p.setp(ax.get_xticklabels()[-1], visible=rightx) ##controls right x axis label
		
		p.xlabel('ppm', fontsize=labelfont)
		p.ylabel('Intensity', fontsize=labelfont)
		
		p.subplot(313)
		p.title((fn,'Initial BIC(AICc) values','Starting number of peaks',number_of_peaks_start,'S/N=%.2f'%sigtonoise))
		curves1=('AICc','BIC',)
		#print 'BICarray',BICarray
		p.plot(NP,AICcarray,linewidth=1)
		p.plot(NP,BICarray, linewidth=1)
		leg=p.legend(curves1,loc=1)
		for t in leg.get_texts():
			t.set_fontsize(legendfont)
		sns.despine()
		p.xlabel('Peaks in fit')
		p.ylabel('Value')
		p.savefig(pdf, format='pdf')
		
		ppmshort=ppm[40:-40:5] ##########This adjusts the number of colored dots in the spectrum
		
		# def colorme(di):
# 			one=rangecol/6.0
# 			two=2*one
# 			three=3*one
# 			four=4*one
# 			five=5*one
# 			six=6*one
# 			if di<one:
# 				color_red=1
# 			if di>=one:
# 				color_red=1-((di-one)/one)
# 			if di>=two:
# 				color_red=0
# 			if di>=four:
# 				color_red=(di-four)/one
# 			if di>=five:
# 				color_red=1
# 			if di>=six:
# 				color_red=1
# 			if di<two:
# 				color_green=0
# 			if di>=two:
# 				color_green=(di-two)/one
# 			if di>=three:
# 				color_green=1
# 			if di>=five:
# 				color_green=1-((di-five)/one)
# 			if di>=six:
# 				color_green=0
# 			if di<one:
# 				color_blue=(di/one)
# 			if di>=one:
# 				color_blue=1
# 			if di>=three:
# 				color_blue=1-((di-three)/one)
# 			if di>=four:
# 				color_blue=0
# 			color_all=(color_red,color_green,color_blue)
# 			return color_all
			
			
		
		
		##############graphs all stat similar models######################################################################################################
		#need two listst [hcwp] and [bic] that have the infor for graphing.
		def graphthis(hcwp,bic):
			#print 'graph ran!!!!!!!!!!!!!!!'
			#annacar=-1
			#print 'hcwfilt in graphthis',hcwp
			river=0
			#kjh=-1
			HCWdellen=len(hcwp)
			print 'HCWdellen %.i'%HCWdellen
			#while annacar<6 and annacar<(HCWdellen-1):
			for item in range(HCWdellen):
				#annacar+=1
				#print 'item',item
				
				sorted13=separate(hcwp[item])
				zadia11=0
				print 'hcwp[item]',hcwp[item]
				#print 'sorted13',sorted13
				
				statsameh=sorted13[zadia11]
				statsamec=sorted13[zadia11+1]
				numpeak=len(statsamec)
				statsamew=sorted13[zadia11+2]
				futurec=statsamec
				print 'futurec',futurec
				print 'statsamew',statsamew
				statsamep=sorted13[zadia11+3]
				i=-1
				statsamehfilt=[]
				statsamecfilt=[]
				statsamewfilt=[]
				statsamepfilt=[]
				
				for anghlk in statsamew:
					i+=1
					if widthmin<anghlk<widthmax:
						statsamehfilt.append(statsameh[i])
						statsamecfilt.append(statsamec[i])
						statsamewfilt.append(statsamew[i])
						statsamepfilt.append(statsamep[i])
						print 'appendi=',i
					
				print 'statsamewfilt',statsamewfilt
					
				numgoodpeaks=len(statsamecfilt)
				
				statsamec_str=stringize(statsamecfilt)
				print 'statsame center %i'%item,statsamec
				print 'statsame width %i'%item, statsamew
				print 'statsame height %i'%item,statsameh
				print 'statsame phase %i'%item,statsamep
				statsamew_str=stringize(statsamewfilt)
				statsamep_str=stringize(statsamepfilt)
				BICcompared=bic[river]-bestC
				p.figure( figsize=(15, 9))
				ax=p.subplot(111)
				ax.xaxis.set_minor_locator(minorLocator1)
				ax.invert_xaxis()
				p.tick_params(which='both', width=1)
				p.tick_params(which='major', length=5)
				p.tick_params(which='minor', length=2, color='k')
				p.figtext(0.01,0.98,statsamec_str,size='small')
				p.figtext(0.01,0.96,statsamew_str,size='small')
				p.figtext(0.5,0.98,statsamep_str,size='small')
				eviratio=1/(m.exp(0.5*BICcompared))
				#p.figtext(0.5,0.96,,size='small')
				p.title((fn,'BIC=%.2f'%bic[item],'%.i peaks'%numpeak,'Rel. BIC= %.2f'%BICcompared,'Rel. likelihood of model=%.3f'%eviratio))##relative likelihood is from page 74 of Burnham and Anderson "model selection and multimodel inference" 2002 second edition Springer
				if showsmooth==True:
					p.plot(ppm,Rawsmooth,linewidth=linewidthgraph,color='K')
				else: 
					p.plot(ppm,Raw1,linewidth=linewidthgraph,color='K')
				sumpeak=lorentzDKp(ppm,hcwp[item])
				p.plot(ppm,sumpeak,linewidth=linewidthgraph,color='green')
				if plotresid==True:
					residualmain=Raw1-sumpeak
					if showsmooth==True:
						residualmain_smooth=sig.savgol_filter(residualmain, Twowindow_length, 2, mode='mirror')
						p.plot(ppm,residualmain_smooth,linewidth=1,color='0.8')
					else:
						p.plot(ppm,residualmain,linewidth=1,color='0.8')	
				statsame_ip=extract_ind_peaks(hcwp[item])
				statsameharray=n.asarray(statsamehfilt)
				statsamewarray=n.asarray(statsamewfilt)
				
				statsamewarray=statsamewarray/(freqsig*2)
				Areass=statsamewarray*statsameharray*n.pi
				sumAreass=n.sum(Areass)
				fAreass=(Areass/sumAreass)
				pcAreass=fAreass*100
				#print 'Best unsplit pcAreass',pcAreass
				pcAreassstr=[]
				for thing in pcAreass:
					pcAreassstr.append("%.1f" %thing)
				p.figtext(0.01,0.94,pcAreassstr,size='small')
				#print 'statsame_ip',statsame_ip
				print len(statsame_ip)
				#print 'futurec',futurec
				for nuth in range(len(statsame_ip)):
						#print nuth
						currentc=futurec[nuth]
						#print currentc
						di=int(n.abs(centcol-currentc)/step)
						if di>(danh-2):
							di=-1
						#print statsame_ip[nuth]
						p.plot(ppm,statsame_ip[nuth],linewidth=linewidthgraph,color=color_all[di])
				clarkh=setyl*(2.0/3.0)
				#print 'clarkh',clarkh
				if setyh==False:
					clarkh=-fac*rms
								
				if spectrum==True:			
					#print 'colorall',color_all
					for hgi in ppmshort:
						di=int(n.abs(centcol-hgi)/step)
						if di>(danh-2):
							di=-1
						#print 'di',di
						p.scatter(hgi,clarkh, linewidth=1, color=color_all[di], marker='o')			
				for t in leg.get_texts():
					t.set_fontsize(legendfont)
				if isinstance(setyh, float):
					ax.set_autoscaley_on(False)
					ax.set_ylim([setyl,setyh])
				if xviewset==True:
					ax.set_autoscalex_on(False)
					ax.set_xlim([leftview,rightview])
				p.setp(ax.get_xticklabels(), fontsize=tickfont)
				p.setp(ax.get_yticklabels(), fontsize=tickfont)
				ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
				p.setp(ax.get_yticklabels()[0], visible=boty) ##controls bottom y axis label   
				p.setp(ax.get_yticklabels()[-1], visible=topy) ##controls top y axis label
				p.setp(ax.get_xticklabels()[0], visible=leftx) ##controls left x axis label 
				p.setp(ax.get_xticklabels()[-1], visible=rightx) ##controls right x axis label
				
				river+=1
				
				# print 'start',start
				# print 'end',end
				# 
				sns.despine()
				p.xlabel('ppm', fontsize=labelfont)
				p.ylabel('Intensity', fontsize=labelfont)
				p.savefig(pdf, format='pdf')
		t564=graphthis(hcwfilt,bicfilt)
		
		def graphresidual(hcwp,bic):
			# annacar=-1
 			river=0
 			HCWdellen=len(hcwp)
# 			while annacar<6 and annacar<(HCWdellen-1):
			for item in range(HCWdellen):
				#annacar+=1
				#print 'item',item
				
				sorted13=separate(hcwp[item])
				zadia11=0
				futurec=sorted13[zadia11+1]
				#print 'hcwp[item]',hcwp[item]
				#print 'sorted13',sorted13
				statsameh=sorted13[zadia11]
				statsamec=sorted13[zadia11+1]
				statsamew=sorted13[zadia11+2]
				statsamep=sorted13[zadia11+3]
				numpeak=len(statsamec)
				
				statsamec_str=stringize(statsamec)
				# print 'statsame center %i'%item,statsamec
# 					print 'statsame width %i'%item, statsamew
# 					print 'statsame height %i'%item,statsameh
# 					print 'statsame phase %i'%item,statsamep
				statsamew_str=stringize(statsamew)
				statsamep_str=stringize(statsamep)
				BICcompared=bic[river]-bestC
				p.figure( figsize=(15, 9))
				ax=p.subplot(111)
				ax.xaxis.set_minor_locator(minorLocator1)
				ax.invert_xaxis()
				p.tick_params(which='both', width=1)
				p.tick_params(which='major', length=5)
				p.tick_params(which='minor', length=2, color='k')
				p.figtext(0.01,0.98,statsamec_str,size='small')
				p.figtext(0.01,0.96,statsamew_str,size='small')
				
				p.figtext(0.5,0.98,statsamep_str,size='small')
				
				p.title((fn,'data for all peaks shown','BIC=%.2f'%bic[item],'%.i peaks'%numpeak,'BIC-bestBIC= %.2f'%BICcompared,))
				if showsmooth==True:
					p.plot(ppm,Rawsmooth,linewidth=linewidthgraph,color='K')
				else:
					p.plot(ppm,Raw1,linewidth=linewidthgraph,color='K')
				sumpeak=lorentzDKp(ppm,hcwp[item])
				p.plot(ppm,sumpeak,linewidth=linewidthgraph,color='green')
				#if plotresid==True:
				residualmain=Raw1-sumpeak
				if showsmooth==True:
					residualmain_smooth=sig.savgol_filter(residualmain, Twowindow_length, 2, mode='mirror')
					p.plot(ppm,residualmain_smooth,linewidth=1,color='#FFA500')
				else:
					p.plot(ppm,residualmain,linewidth=1,color='#FFA500')
				hcwpphased=hcwp[item]
				if showfitsphased==True:
					hcwpphased=hcwp[item]
					cnh=0
					for i in range(len(hcwpphased)/4):
						hcwpphased[cnh+3]=0
						cnh+=4
				
				statsame_ip=extract_ind_peaks(hcwpphased)
				statsameharray=n.asarray(statsameh)
				statsamewarray=n.asarray(statsamew)
				statsamewarray=statsamewarray/(freqsig*2)
				Areass=statsamewarray*statsameharray*n.pi
				sumAreass=n.sum(Areass)
				sumAreasstr='%.2E'%sumAreass
				fAreass=(Areass/sumAreass)
				pcAreass=fAreass*100
				#print 'Best unsplit pcAreass',pcAreass
				pcAreassstr=[]
				Areassstr=[]
				for thing in pcAreass:
					pcAreassstr.append("%.1f" %thing)
				for thing in Areass:
					Areassstr.append('%.2E'%thing)
				p.figtext(0.01,0.94,pcAreassstr,size='small')
				p.figtext(0.5,0.96,Areassstr,size='small')
				p.figtext(0.5,0.94,sumAreasstr,size='small')
				if showfitsphased==True:
					c2=0
					for nuth in range(len(statsame_ip)):
						#print nuth
						currentc=statsamec[nuth]
						#print currentc
						di=int(n.abs(centcol-currentc)/step)
						if di>(danh-2):
							di=-1
						#print statsame_ip[nuth]
						p.plot(ppm,statsame_ip[nuth],linewidth=linewidthgraph,color=color_all[di])
				clarkh=setyl*(2.0/3.0)
				print 'clarkh',clarkh
				if setyh==False:
 					clarkh=-fac*rms	
				for t in leg.get_texts():
					t.set_fontsize(legendfont)
				if isinstance(setyh, float):
					ax.set_autoscaley_on(False)
					ax.set_ylim([setyl,setyh])
				if xviewset==True:
					ax.set_autoscalex_on(False)
					ax.set_xlim([leftview,rightview])
				p.setp(ax.get_xticklabels(), fontsize=tickfont)
				p.setp(ax.get_yticklabels(), fontsize=tickfont)
				ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
				p.setp(ax.get_yticklabels()[0], visible=boty) ##controls bottom y axis label   
				p.setp(ax.get_yticklabels()[-1], visible=topy) ##controls top y axis label
				p.setp(ax.get_xticklabels()[0], visible=leftx) ##controls left x axis label 
				p.setp(ax.get_xticklabels()[-1], visible=rightx) ##controls right x axis label
				
				river+=1
				
				# print 'start',start
				# print 'end',end
				# 
				sns.despine()
				p.xlabel('ppm', fontsize=labelfont)
				p.ylabel('Intensity', fontsize=labelfont)
				p.savefig(pdf, format='pdf')
		residualplots=graphresidual(hcwfilt,bicfilt)
	
				
		
		feed=[]
		Zadia=[]
		feedf=[]
		#########################################################################################################################
		######################### Input peaks figure ############################################################################
		#########################################################################################################################
		if inputpeaks==True:
			inputindpeaks=[]
			with open(inputpeakloc, 'rb') as f:
				
				reader = csv.reader(f)
				for row in reader:
					Zadia=row[0].split()
					for thing in Zadia:
						feed.append(float(thing))
					feedf=[feed]
			
			integratearray=n.sum(Raw1*step)
			z=0
			#print 'feedf',feedf
			dstep2=(len(feedf[0]))/4
			for i in range(0,dstep2):
				feedf[0][z+2]=feedf[0][z+2]/(freqsig/2)
				feedf[0][z]=feedf[0][z]*integratearray/(feedf[0][z+2]*n.pi)
				z+=4
			feedpeak=[]
			for thing in feedf:
				inputindpeaks.append(extract_ind_peaks_b(thing))
			#print 'inputindpeaks',inputindpeaks
			len456=len(inputindpeaks[0][0])
			#print 'len inputindpeaks',len456
			#print 'inputindpeaks[0][maxindex1]',inputindpeaks[0][0][maxindex1]
			
			#print 'integratearray',integratearray
			
			for item in feedf:
				feedpeak.append(lorentzDKp(ppm,item))
			
			feedpeakmax=n.max(feedpeak[0])
			bestpeakmax=n.max(bestpeak)
			scale=bestpeakmax/feedpeakmax
			feedpeak[0]=feedpeak[0]*scale
			p.figure(figsize=(15, 9))
			ax=p.subplot(111)
			ax.xaxis.set_minor_locator(minorLocator1)
			ax.invert_xaxis()
			p.tick_params(which='both', width=1)
			p.tick_params(which='major', length=5)
			p.tick_params(which='minor', length=2, color='k')
			if show_title==True:
				p.title((inputpeakloc,'Input peaks','S/N=%.2f'%sigtonoise,))
			best_together=separate(feedf)
			feedh=[]
			feedc=[]
			feedw=[]
			feedp=[]
			c23=0
			for i in range(len(best_together)/4):
				feedh.append(best_together[c23])
				feedc.append(best_together[c23+1])
				feedw.append(best_together[c23+2])
				feedp.append(best_together[c23+3])
				c23+=4
			
			
			
			feedc_str=stringize(feedc)
			#print 'input center',feedc
			feedw_str=stringize(feedw)
			#print 'input FWHM',feedw
			feed_area=area(feedh,feedw)
			#print 'input pc area',feed_area[1]
			feed_area_str=stringize(feed_area[1])
			feedp_str=stringize(feedp)
			p.figtext(0.1,0.98,feedc_str,size='small')
			p.figtext(0.1,0.96,feedw_str,size='small')
			p.figtext(0.1,0.94,feed_area_str,size='small')
			p.figtext(0.5,0.98,feedp_str,size='small')
			p.figtext(0.5,0.96,scale,size='small')
			p.plot(ppm,Rawsmooth,linewidth=linewidthgraph,label='Smoothed Data',color='K')
			c2=0
			counter67=-1
			for item in inputindpeaks:
				counter67+=1
				nopey=1
				p.plot(ppm,feedpeak[c2],linewidth=linewidthgraph,color="green")
				counter68=-1
				for i in range(len(item)):
					counter68+=1
					currentc=float(feedc_str[counter67][counter68])
					di=int(n.abs(centcol-currentc)/step)
					if di>(danh-2):
							di=-1
					item[i]=item[i]*scale
					if nopey==1:
						p.plot(ppm,item[i],linewidth=linewidthgraph,color=color_all[di],label='%.f peaks'%len(item))
						nopey+=1
					else:
						p.plot(ppm,item[i],linewidth=linewidthgraph,color=color_all[di])
				c2+=1
			clarkh=setyl*(2.0/3.0)
			
			if setyh==False:
				clarkh=-fac*rms
			if spectrum==True:				
				for hgi in ppmshort:
					di=int(n.abs(centcol-hgi)/step)
					if di>(danh-2):
						di=-1
					p.scatter(hgi,clarkh, linewidth=1, color=color_all[di], marker='o')	
			p.setp(ax.get_xticklabels(), fontsize=tickfont)
			p.setp(ax.get_yticklabels(), fontsize=tickfont)
			ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
			p.setp(ax.get_yticklabels()[0], visible=boty) ##controls bottom y axis label   
			p.setp(ax.get_yticklabels()[-1], visible=topy) ##controls top y axis label
			p.setp(ax.get_xticklabels()[0], visible=leftx) ##controls left x axis label   
			p.setp(ax.get_xticklabels()[-1], visible=rightx) ##controls right x axis label
			sns.despine()
			p.xlabel('ppm', fontsize=labelfont)
			p.ylabel('Intensity', fontsize=labelfont)
			if isinstance(setyh, float):
				ax.set_autoscaley_on(False)
				ax.set_ylim([setyl,setyh])
			if xviewset==True:
				ax.set_autoscalex_on(False)
				ax.set_xlim([leftview,rightview])
			p.savefig(pdf, format='pdf')
		
		print 'done'
		
		print 'c1end',c1
		print 'ph1 end',ph1
		
		
		#################################  montecarlo ########################################################
		
		from mpl_toolkits.mplot3d import Axes3D
		includeinput=False
		if mc0123==True:	
			if monte_carlo[0]==True:	
				p.figure( figsize=(15, 9))
				sns.despine()
				pcsucess=(monte_carlo[2]*100)/num_of_mc_peaks
				print pcsucess
				p.title(('MinBIC=%.2f'%minBIC,'MaxBIC=%.2f'%maxBIC,'Number tried=%i'%num_of_mc_peaks,'pc success=%i'%pcsucess,))
				if improved==True:
					p.title(('Better mc!','MinBIC=%.2f'%minBIC,'MaxBIC=%.2f'%maxBIC,'MC tries=%i'%num_of_mc_peaks,'pc success=%0i'%pcsucess,))
				
				ax=p.subplot(111)
				ax.invert_xaxis()
				p.plot(ppm,Rawsmooth,linewidth=1,color='0.3')
				if leftright==True:
					ax.set_autoscalex_on(False)
					ax.set_xlim([left,right])
				else:
					ax.set_autoscalex_on(False)
					ax.set_xlim([s,e])
				c2=0
				betterpeakHCW1=[]
				betterpeakhcw1=[]
				betterbic1=[]
				#print 'mc_hcw_list',mc_hcw_list
				for item in mc_hcw_list:
					mc_ip=[]
					summcpeak=lorentzDKp(ppm,item)
					
					mc_ip=extract_ind_peaks_b(item)
					print 'mcip',mc_ip
					if mcdif1[c2]<-2:
						betterpeakHCW1.append(mc_ip)
						betterpeakhcw1.append(item)
						betterbic1.append(mcdif1[c2])
					
					if mcdif1[c2]>=2:
						colory='c'
					elif 2>mcdif1[c2]>=-2:
						colory='0.1'
					elif -2>mcdif1[c2]>=-5:
						colory='g'
						#print 'mc found a better peak',colory
						print 'mcdif1',mcdif1[c2]
						print item
						print '2.5to4.5better',item
						
					elif -5>mcdif1[c2]>=-10:
						colory='#FFA500'
						print 'mc found a better peak',colory
						print 'mcdif1',mcdif1[c2]
						print '4.5to8.5better',item 
						#p.figtext(0.1,0.96,item,size='small')
					elif mcdif1[c2]<-10:
						colory='r'
						print 'mc found a better peak',colory
						print 'mcdif1',mcdif1[c2]
						print 'more than 8.5 better',item
						#p.figtext(0.1,0.94,item,size='small')
						
					p.plot(ppm,summcpeak,linewidth=0.3,color=colory)
					for peak in mc_ip:
						print 'peak',peak
						p.plot(ppm,peak,linewidth=0.3,color=colory)
							
					c2+=1
					#print 'c1',c1
					
				
				if includeinput==True and inputpeaks==True:
					c2=0
					counter67=-1
					
					for item in inputindpeaks:
						counter67+=1
						nopey=1
						p.plot(ppm,feedpeak[c2],linewidth=linewidthgraph,color="green")
						counter68=-1
						for i in range(len(item)):
							counter68+=1
							currentc=float(feedc_str[counter67][counter68])
							di=int(n.abs(centcol-currentc)/step)
							if di>(danh-2):
									di=-1
							item[i]=item[i]*scale
							if nopey==1:
								p.plot(ppm,item[i],linewidth=linewidthgraph,color=color_all[di],label='%.f peaks'%len(item))
								nopey+=1
							else:
								p.plot(ppm,item[i],linewidth=linewidthgraph,color=color_all[di])
						c2+=1
				p.xlabel('ppm',fontsize=labelfont)
				p.ylabel('Intensity',fontsize=labelfont)
				p.setp(ax.get_xticklabels(), fontsize=tickfont)
				p.setp(ax.get_yticklabels(), fontsize=tickfont)
				ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
				p.setp(ax.get_yticklabels()[0], visible=boty) ##controls bottom y axis label   
				p.setp(ax.get_yticklabels()[-1], visible=topy) ##controls top y axis label
				p.setp(ax.get_xticklabels()[0], visible=leftx) ##controls left x axis label   
				p.setp(ax.get_xticklabels()[-1], visible=rightx) ##controls right x axis label	
				
				
				p.savefig(pdf, format='pdf')
				if showgraphs==True:
					p.show()
				doit=False
				if betterpeakHCW1:
					doit=True
					print 'doit=True##########'
				lenbic1ed=len(betterbic1)
				if doit==True:
				
					p.figure( figsize=(15, 9))
					ax=p.subplot(111)
					p.plot(ppm,Rawsmooth,linewidth=1,color='0.3')
					if lenbic1ed > 1:
						betterbic1, betterpeakhcw1, betterpeakHCW1= zip(*sorted(zip(betterbic1,betterpeakhcw1,betterpeakHCW1)))#sort the data by the centers
					counteruy=-1
					hytg=-1
					while counteruy<3 and counteruy<(len(betterpeakHCW1)-1):
						for item in betterpeakHCW1:
							counteruy+=1
							comp734=betterbic1[0]-betterbic1[counteruy]
							if comp734>-2:
								hytg+=1
								if counteruy==0:
									for i in item:
										print 'i',i
										print 'len(i)',len(i)
										print 'len(ppm)',len(ppm)
										p.plot(ppm,i,linewidth=2,color='r')##BIC differend is between -2 and -4 (mc model may be better than best fitted evidence ratio between 2.7 and 7.4)
								else:
									for i in item:
											print 'i',i
											print 'len(i)',len(i)
											print 'len(ppm)',len(ppm)
											p.plot(ppm,i,linewidth=1,color='0.5')##BIC differend is between -2 and -4 (mc model may be better than best fitted evidence ratio between 2.7 and 7.4)
					p.title(('All the montecarlo fits (%.i -grey) that are stat same as best mc fit (red)'%hytg))		
					p.xlabel('-ppm',fontsize=labelfont)
					p.ylabel('Intensity',fontsize=labelfont)
					
					ax.invert_xaxis()
					ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
					p.setp(ax.get_xticklabels(), fontsize=tickfont)
					p.setp(ax.get_yticklabels(), fontsize=tickfont)
					p.setp(ax.get_yticklabels()[0], visible=boty) ##controls bottom y axis label   
					p.setp(ax.get_yticklabels()[-1], visible=topy) ##controls top y axis label
					p.setp(ax.get_xticklabels()[0], visible=leftx) ##controls left x axis label   
					p.setp(ax.get_xticklabels()[-1], visible=rightx) ##controls right x axis label
					if leftright==True:
						ax.set_autoscalex_on(False)
						ax.set_xlim([left,right])
					else:
						ax.set_autoscalex_on(False)
						ax.set_xlim([s,e])
					
					if showgraphs==True:
						p.show()
					p.savefig(pdf, format='pdf')
					
					
				pdf.close()
				if threedfig==True:		
					pdf=PdfPages('3DMC'+fn+whichcols+'_'+phase+fares+'.pdf')
					p.figure( figsize=(15, 9))
					
					ax=p.subplot(111,projection='3d')
					fx124=mcpeakw
					
					fy126=mcpeakc
					fz125=mcpeakh
					apple=n.array([])
		# 			print 'fx',fx124
		# 			print 'fz',fz125
		# 			print 'fy',fy126
		# 			print 'cal_mc_BIC_list',cal_mc_BIC_list
					for i in cal_mc_BIC_list:
						apple=n.append(apple,30*i)
					# print 'apple',apple
		# 			print 'feedh',feedh
		# 			print 'feedw',feedw
		# 			print 'feedc',feedc
					big=n.amax(apple)
				
					ax.scatter(fx124, fy126, fz125,s=apple,c='0.5',marker='o')
					if inputpeaks==True:
						ax.scatter(feedw[0], feedc[0], feedh[0],s=big,c='r',marker='^')
				
					ax.set_xlabel('FWHM (Hz)',fontsize=labelfont*0.66)
					ax.set_ylabel('Center (ppm)',fontsize=labelfont*0.66)
					ax.set_zlabel('Height (Rel. units)',fontsize=labelfont*0.66)
					ax.ticklabel_format(axis='z', style='sci', scilimits=(-2,2))
					p.setp(ax.get_xticklabels(), fontsize=tickfont*0.66)
					p.setp(ax.get_yticklabels(), fontsize=tickfont*0.66)
					p.setp(ax.get_zticklabels(), fontsize=tickfont*0.66)
					p.savefig(pdf, format='pdf')
					if showgraphs==True:
						p.show()
					
					pdf.close()
			else:
				if showgraphs==True:
					p.show()
				pdf.close()
		else:
			if showgraphs==True:
				p.show()
			pdf.close()
		p.close('all')
		###########you can export data to a csv file using the code below###################################		
		# import csv
# 		data_output=(ppm,Raw1)#####this is what data you want to export#############
# 		writer = csv.writer(open(fn+"  "+".csv", "wb"))
# 		for item in data_output:
# 			writer.writerow(item)