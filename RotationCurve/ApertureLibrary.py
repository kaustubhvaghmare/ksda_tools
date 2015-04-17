'''
This program contains tools for aperture extraction.
It is based on the original Aperture Extraction code written by Rajin.

The following are the primary modifications done to the program.
* Use spatial profile itself as a weighing function.
* Don't actually create aperture summed FITS file but merely return apertures

Version 1.0
'''

import numpy as np
import pylab as plt
import pyfits as pft
import os,sys,glob,string,time,subprocess
import math as mth
import scipy as sc
from scipy.integrate import simps
from scipy import interpolate
from astropy.table import Table, Column


#------------------------>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----------------------------------------------
#			 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#------------------------>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----------------------------------------------
# path to where code is running 
pathrun = os.getcwd() # '/home/kaustubh/SALT/Experiments/Rajin_Apex/'
# name of fits file
#filename = 'NGC1553_2df.fits'
#file_err =  'NGC1553_2df.fits'
#setting pylab to interactive mode
plt.ion()
#------------------------>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----------------------------------------------
#			 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#------------------------>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----------------------------------------------


#------------------------------------------------------------------------------------------------------------------------------------
#function to get the coordinates of the windows through GUI
def wincoord(winno):
	cor1 = 1 ; cor2 = 0
	while cor1 > cor2:
		print("Please click at the start (left-side) of the window %d" % (winno))
		cor1 = plt.ginput(1)[0][0]
		print(cor1)
		print("Please click at the end (right-side) of the window %d" % (winno))
		cor2 = plt.ginput(1)[0][0]
		print(cor2)
		if cor1 > cor2:
			#print('pixel coordinates of the start of the window is greater than the end! So lets restarts input for window ' % (winno)) 
			# Check if coordinates were entered in Arabic fashion! Swap if needed.
			cor1, cor2 = cor2, cor1
	return round(cor1,0),round(cor2,0)

# Like wincoord but gets coordinates by accepting width.
def cuicoord(winno, centroid):
	while True:
		width = int(raw_input("Enter width: "))
		centroid = int(centroid)
		cor1 = centroid - width/2
		cor2 = centroid + width/2
		if cor1 > 0:
			return cor1, cor2
	
#------------------------------------------------------------------------------------------------------------------------------------
#function to get the input value of a user (and deal in the cases of typing mistakes)
def input_val(script,lim1,lim2):
        while True:
           try:
               value = int(raw_input('%s : ' % (script)))
           except ValueError: # just catch the exceptions you know!
               print 'That\'s not a number!'
           else:
               if lim1 <= value < lim2: 
                   break
               else:
                   print 'Out of range. Try again'
        return value

#------------------------------------------------------------------------------------------------------------------------------------

#function written to ask the questions to the user!
def input_str(script):
        while True:
                userInput = str(raw_input('%s' % (script)))
                if len(userInput) == 1:
                        if userInput in string.letters:
                                if userInput.lower() == 'y' or userInput.lower() == 'n':                
                                        break
                                print 'Please enter only "y" or "n"'
                        else:
                                print 'Please enter only letters!'
                elif len(userInput) == 0:
                        print 'Please enter at least 1 character!'
                elif len(userInput) >1 and userInput.isalpha():
                        print 'Please enter only 1 character!'
                else:
                        print 'Please enter only letters and no numbers'
        return userInput.lower()

#------------------------------------------------------------------------------------------------------------------------------------

def zscale(image, contrast=1.0):
	'''Implementation of the IRAF zscale algorithm to find vmin and vmax parameters for the dynamic range of a display. It finds the 		image values near the median image value without the time consuming process of computing a full image histogram.
	'''
	from scipy import optimize
	# Get ordered list of points
	I=np.sort(image.flatten())

	# Get number of points
	npoints=len(I)
	# Find the midpoint (median)
	midpoint=(npoints-1)/2

	# Fit a linear function
	# I(i) = intercept + slope * (i - midpoint)
	fitfunc = lambda p, x: p[0]*x+p[1]
	errfunc = lambda p, x, y: fitfunc(p, x) - y

	# Initial guess for the parameters
	p0 = [(I[-1]-I[0])/npoints,I[midpoint]] 
	# Fit
	i=np.arange(len(I))
	p1, success = optimize.leastsq(errfunc, p0[:], args=(i, I))

	if success in [1,2,3,4]:
	       slope=p1[0]
	       z1=I[midpoint]+(slope/contrast)*(1-midpoint)
	       z2=I[midpoint]+(slope/contrast)*(npoints-midpoint)
	else:
	       z1=np.min(image)
	       z2=np.max(image)

	return z1, z2
#------------------------------------------------------------------------------------------------------------------------------------
#function to calculate area of a region with simpson's rule
def area(x,y):
	return simps(y=y,x=x)

#------------------------------------------------------------------------------------------------------------------------------------
#
def weigh_func(ans,xnew):
	size = len(xnew)
	if ans == 'y':
		x,y = np.loadtxt('selfdeffunc.txt')
	else:
		x = np.arange(-1,1.05,0.05)
		y = np.zeros(len(x),float)
		for j in range(0,len(x)):
			if x[j] < 0: y[j] = -x[j]+1
			if x[j]== 0: y[j] = 1
			if x[j] > 0: y[j] = x[j]+1
	#interpolate to resize x,y array to proper length
	#xinterp = np.arange(0,1,1./size)
	f = interpolate.interp1d(x, y)
	#yinterp = f(xinterp)
	ynew = f(xnew)
#	ynew = 	yinterp*len(xnew)
#	return xinterp,yinterp
	trap_area=[]; diff_x = x[1]-x[0]
	for j in range(0,size-1):
		trap_area.append(0.5*(ynew[j]+ynew[j+1])*diff_x)
	trap_area = np.array(trap_area)
	total_area = np.sum(trap_area)
	frac_area = trap_area/total_area
	return frac_area

#------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------

def GetCentroid( spatial_points, intensities):
	return np.sum(spatial_points * intensities)/np.sum(intensities)

#change directory to where image file is
os.chdir(pathrun)
def get_aperture(filename, file_err=None):

	if file_err is None:
		file_err = filename
	#getting fits file info
	data = pft.getdata(filename)
	header = pft.getheader(filename)
	dataerr = pft.getdata(file_err)
	errheader = pft.getheader(file_err)
	#Summing in the column direction to get sum of flux in the cross section
	xsection = np.sum(data,axis=1)
	temp = np.arange(len(xsection))
	np.savetxt("selfdeffunc.txt", np.vstack( (temp, xsection) ) )
	
	centroid = GetCentroid( range(header["NAXIS2"]), np.sum(data, axis=1) )
	
	#plotting the x-section for window selection
	plt.plot(xsection)
	plt.xlabel('spatial pixel direction')
	plt.ylabel('Flux')
	plt.title('Cross Section Plot; Central Row = %d' % int(centroid))
	windowno = input_val('Enter the number of windows(targets) you want to have aperture extraction',0,10)
	
	cui_mode = input_str('Enter window boundaries through keyboard? (y/n): ')
	coord = np.zeros((windowno,2),int)
	#window selection 
	for i in range(0,windowno):
		if cui_mode == "n":
			coord[i,:] = wincoord(i+1)
		else:
			coord[i,:] = cuicoord(i+1, centroid)
	
	plt.close()
	
	#plotting the selected windows
	plt.plot(xsection)
	plt.xlabel('spatial pixel direction')
	plt.ylabel('Flux')
	plt.title('Plot selected windows')
	xsec_max = xsection.max()
	for i in range(0,windowno):
		low = coord[i,0] ; up = coord[i,1]
		plt.broken_barh([ (low, up-low) ] , (0, xsec_max), facecolors='red',alpha=0.2,label='win_1')
	plt.close()
	
	#add a flagging module here
	#---
	
	#---
		
	#calculating area of the windows in the x-section -> as the x-section is a sum of the fluxes along the row direction, 
	#area below the x-section curve represent fluxes in the row direction as well
	cum_area = np.zeros((windowno,len(data)),float)
	ins_area = np.zeros((windowno,len(data)),float)
	max_area = []
	for i in range(0,windowno):
		startx = coord[i,0] ; endx = coord[i,1]
		for j in range(startx+1,endx):
			datax = np.arange(startx,j)
			datay = xsection[startx:j]
			cum_area[i,j] = area(datax,datay)
			datax1 = np.arange(j-1,j)
			datay1 = xsection[j-1:j]
			ins_area[i,j] = area(datax1,datay1)
		max_area.append(cum_area[i,j])
	
	
	#need to use weighing function to select apertures 
	
	aperture_coord = np.zeros((windowno,40),float)
	for i in range(0,windowno):
		#plotting the selected windows
		plt.plot(xsection)
		plt.xlabel('spatial pixel direction')
		plt.ylabel('Flux')
		plt.title('Plot selected windows')
		low = coord[i,0] ; up = coord[i,1]
		xsec_max = xsection[low:up].max()
		plt.broken_barh([ (low, up-low) ] , (0, xsec_max), facecolors='red',alpha=0.1,label='win_1')
		num_apert = int(raw_input('Please enter the number of apertures for window (max=40) '+str(i)+': '))
		plt.close()
		x_in = np.linspace(-1,1,num_apert+1)
		y_in = weigh_func('n',x_in)
		aperture_areas = max_area[i]*y_in
		calc_cum_areas = np.cumsum(aperture_areas)
		for j in range(0,num_apert+1):
			if j == 0:
				dum = np.where(cum_area[i] > 0)[0][0]
				aperture_coord[i,j] = dum
				plt.plot([dum,dum,dum],[0,xsec_max/2.0,xsec_max])
			if j > 0 and j < num_apert:
				dum = np.where(cum_area[i] > calc_cum_areas[j-1])[0][0]
				aperture_coord[i,j] = dum
				plt.plot([dum,dum,dum],[0,xsec_max/2.0,xsec_max])
			if j == num_apert:
				dum = np.where(cum_area[i] > 0)[0][-1]
				aperture_coord[i,j] = dum
				plt.plot([dum,dum,dum],[0,xsec_max/2.0,xsec_max])
		plt.plot(xsection)
	
	return aperture_coord	
#	for i in range(0,windowno):
#		for j in range(0,num_apert):
#			x1 = int(aperture_coord[i,j])
#			x2 = int(aperture_coord[i,j+1])
#			dataout = np.sum(data[x1:x2,:],axis=0)
#			dataerrout = (np.sum(dataerr[x1:x2,:]*dataerr[x1:x2,:],axis=0))**0.5
#			header['APERTURE'] = 'Win'+str(i)+'ap'+str(j)
#			header['AP-VAL'] = str(x1)+':'+str(x2)
#			header['AP-MID'] = str((x1+x2)/2.)
#			errheader['APERTURE'] = 'Win'+str(i)+'ap'+str(j)
#			errheader['AP-VAL'] = str(x1)+':'+str(x2)
#			errheader['AP-MID'] = str((x1+x2)/2.)
#	
#			object_name = header['OBJECT']
#	                split_name = object_name.split()
#	                new_name = ''
#	                for ii in range(0,len(split_name)):
#	                        if ii == 0: new_name = new_name+split_name[ii]
#	                        else: new_name = new_name+'_'+split_name[ii]
#		
#			pft.writeto(new_name+'win'+str(i)+'ap'+str(j)+'.fits',data = dataout, header = header)
#			pft.writeto(new_name+'err_win'+str(i)+'ap'+str(j)+'.fits',data = dataerrout, header = errheader)
#



