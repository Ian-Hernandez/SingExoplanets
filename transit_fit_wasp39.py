# import packages to use
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.io import fits
import glob 
import batman
import lmfit
import corner


# directory='/Users/annaburkholder/exp_det_scripts/visit23_defringed/'
directory='/home/ian/Desktop/WebbData/visit23_defringed/'
#directory='/visit23_defringed/' #Change directory to proper location
number_of_images=43

#Load images into a list
list=glob.glob(directory+"*.fits")
#print, first image in list
print(list[0])
####Exmaple load first fits image
hdul=fits.open(list[1])
#Get MJD mid time of exposure from Header, which has start and end MJD times
mjd_start=hdul[0].header['EXPSTART']
mjd_end=hdul[0].header['EXPEND']
mjd=(mjd_end+mjd_start)/2.
#maybe no period?
print(mjd)

#load fits file image into an array called 'data'
data = hdul[0].data
data.shape      #size of image
data.dtype.name #type of image
print(np.sum(data)) #total counts in image

#close fits after loading in data needed
hdul.close()

#plot loaded image
#plt.imshow(data, cmap='hot')
#plt.colorbar()
#show()


#load all fits images
#Arrays created for MJD time, and the white light curve total_counts
index_of_images=np.arange(number_of_images) #
mjd=np.zeros((number_of_images))
total_counts=np.zeros((number_of_images))
for i in index_of_images:
    img=list[i]
    print(img)
    hdul=fits.open(img)
    mjd_start=hdul[0].header['EXPSTART']
    mjd_end=hdul[0].header['EXPEND']
    mjd_image=(mjd_end+mjd_start)/2.
    mjd[i]=mjd_image
    print(mjd[i])
    data = hdul[0].data
    print(np.sum(data)) #total counts in image
    total_counts[i]=np.sum(data[60-6:60+6,0:1024]) #total counts in 12 pix wide aperature around pixel 60 image
	
	#Want to find max of each column to get a 12 pix range, more accurate than just picking pixel 60
	#use these pixels to create a line of best fit (linear, second order, etc)
	#phase marks times where hubble is looking at the star, b/c it's on the wrong side of earth half the time
	#can plot total_counts vs pixels to see spectra

total_error=np.sqrt(total_counts)


#TRANSIT model batman package https://astro.uchicago.edu/~kreidberg/batman/
#Setup inital parameters (can get from http://exoplanets.org/detail/WASP-39_b)
params = batman.TransitParams()        #object to store transit parameters
params.t0 = 55342.46880                 #mjd time of inferior conjunction  56368.454950
params.per = 4.0552590                #orbital period
params.rp = 0.14491                    #planet radius (in units of stellar radii)
params.a = 11.70                      #semi-major axis (in units of stellar radii)
params.inc = 87.83                    #orbital inclination (in degrees)
params.ecc = 0.                        #eccentricity
params.w = 90.                         #longitude of periastron (in degrees)
params.limb_dark = "nonlinear"         #limb darkening model
params.u = [0.6642, -0.4701, 1.0466,-0.5097]      #limb darkening coefficients
#>>> dir(params)    #View attributes
#['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'a', 'ecc', 'fp', 'inc', 'limb_dark', 'per', 'rp', 't0', 't_secondary', 'u', 'w']
#print(params.rp)
#0.1453


#Batman Transit model
t_fine = np.linspace(np.min(mjd), np.max(mjd), 1000) #times at which to calculate light curve
m = batman.TransitModel(params, t_fine)     #initializes model
transit_fine = m.light_curve(params)        #Calculates Fine-Grid Tranist model

m = batman.TransitModel(params, mjd)     #initializes model
transit_mjd = m.light_curve(params)        #Calculates Fine-Grid Tranist model


#Calculate phase of HST telescope, 96.36 minutes, used for detrending
#phase_hst array values should be between -0.5 and 0.5
phase_hst=(mjd-np.min(mjd))/(6.691666E-2) #phase in days relative to first exposure
phase_hst_fix=np.fix(phase_hst)      #rounds to nearest zero
phase_hst=phase_hst-phase_hst_fix    #all between 0 to 1.0
phase_hst[np.where(phase_hst > 0.5)]=phase_hst[np.where(phase_hst > 0.5)]-1.0 #all between -0.5 and 0.5
#matplotlib.pyplot.scatter(mjd,phase_hst)  #View


#Calculate Orbital Phase
phase=(mjd-params.t0)/(params.per)  #phase in days relative to T0 ephemeris
phase=phase-np.fix(phase[number_of_images-1]) # Have current phase occur at value 0.0

# create a set of Parameters
p = lmfit.Parameters()    #object to store L-M fit Parameters
p.add('rprs', value=params.rp)   # HST phase^4 parameter
p.add('t0', value=params.t0)   # HST phase^4 parameter
p.add('f0', value=10736338.8, min=0 , max=total_counts[0]*1.2)    #Baseline Flux
p.add('p1', value=0.0)   # Orbital phase parameter
p.add('p2', value=0.0)   # HST phase^1 parameter
p.add('p3', value=0.00)   # HST phase^4 parameter
p.add('p4', value=0.0)   # HST phase^3 parameter
p.add('p5', value=0.0)   # HST phase^4 parameter

x=mjd
y=total_counts
err=total_error

def residual(p):
    params.rp = p['rprs'].value                # Set Batman rprs to new fit rprs
    params.t0 = p['t0'].value                # Set Batman rprs to new fit rprs
    m = batman.TransitModel(params, mjd)       #initializes model
    transit_mjd = m.light_curve(params)        #Calculates Fine-Grid Tranist model
    model=transit_mjd*p['f0'] * (p['p1']*phase + p['p2']*phase_hst + p['p3']*phase_hst**2. + p['p4']*phase_hst**3. + p['p5']*phase_hst**4. + 1.0)              #Simple transit model is baseline flux X transit model
    chi2now=np.sum((y-model)**2/err**2)
    #print("current chi^2=",chi2now)
    return (y-model)/err

def model(p):
    params.rp = p['rprs'].value                # Set Batman rprs to new fit rprs
    params.t0 = p['t0'].value                # Set Batman rprs to new fit rprs
    m = batman.TransitModel(params, mjd)       #initializes model
    transit_mjd = m.light_curve(params)        #Calculates Fine-Grid Tranist model
    model=transit_mjd*p['f0'] * (p['p1']*phase + p['p2']*phase_hst + p['p3']*phase_hst**2. + p['p4']*phase_hst**3. + p['p5']*phase_hst**4. + 1.0)                  #Simple transit model is baseline flux X transit model
    return model

def model_fine(p):
    params.rp = p['rprs'].value                # Set Batman rprs to new fit rprs
    params.t0 = p['t0'].value                # Set Batman rprs to new fit rprs
    m = batman.TransitModel(params, t_fine)       #initializes model
    transit_mjd = m.light_curve(params)        #Calculates Fine-Grid Tranist model
    model_fine=transit_mjd*p['f0']
    return model_fine

# create Minimizer
mini = lmfit.Minimizer(residual, p, nan_policy='omit')
# first solve with Nelder-Mead
#out1 = mini.minimize(method='Nelder')
# then solve with Levenberg-Marquardt using the
# Nelder-Mead solution as a starting point
# https://lmfit.github.io/lmfit-py/fitting.html
result = mini.minimize(method='leastsq')
#result = mini.minimize(method='leastsq', params=out1.params)
print(dir(result))  # To View All Atributes of the
print("redchi",result.redchi)
print("chi2",result.chisqr)
print("nfree",result.nfree)
print("bic",result.bic)
print("aic",result.aic)
print("L-M FIT Variable")
print(lmfit.fit_report(result.params))
#print(p['f0'].value)
# file-output.py

#Update with best-fit parameters
p['rprs'].value=result.params['rprs'].value
p['t0'].value=result.params['t0'].value
p['f0'].value=result.params['f0'].value
p['p1'].value=result.params['p1'].value
p['p2'].value=result.params['p2'].value
p['p3'].value=result.params['p3'].value
p['p4'].value=result.params['p4'].value
p['p5'].value=result.params['p5'].value
# Re-calculate Bestfit Model
final_model=model(p)
final_model_fine=model_fine(p)
print("residual standard deviation",np.std((y-final_model)/p['f0'].value))
print("residual standard deviation (ppm)",1E6*np.std((y-final_model)/p['f0'].value))

#Plot data models
matplotlib.pyplot.scatter(mjd,total_counts, linewidth=1.5)
plt.errorbar(mjd,total_counts,yerr=total_error, fmt='o')
xlabel('MJD')
ylabel('Counts')

#overplot models
plot(t_fine,final_model_fine, linewidth=1.5)            #overplot fine-grid Transit model
matplotlib.pyplot.scatter(mjd,final_model, linewidth=2)  #overplot Transit model at data

show()
#Plot Residuals
matplotlib.pyplot.scatter(mjd,(y-final_model)/p['f0'].value, linewidth=2)  #overplot Transit model at data

show()


