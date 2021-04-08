

###############################################################
#*************************************************************#
####################### user inputs ###########################

p=2.2 # the power law index of the electron distribution
Ee=(3.0/7.0)*2.0e40 #1e42 # total energy in electrons
Eb=(4/3)*Ee # total energy in magnetic fields, ratio to Ee is (4/3) in equipartition
gamma10 = 1e1 # initial Lorentz factor of lowest energy electrons 
gamma20= 1e7 # initial Lorentz factor of highest energy electrons
r0= 10 * 0.08 * 2.99e10 #estimated from given (delta)T and v~0.1c # initial radius of blob when electron accelerated (e.g. 3e8 cm = 100 Schwarzschild radii for 10 M_sun BH)
vc= 0.1 #0.15 # expansion speed of blob in units of c (expansion here means expansion, not speed away from BH)
dkpc= 2.4 #8.0 # distance to source in kpc
#These can all be altered
n1ghz,n2ghz,n3ghz= 1.4, 5, 16 #  # frequencies (in GHz) to evaluate light curve at
l1="MeerKAT 1.4 GHz" # #Plot labels for those frequencies
l2="eMerlin 5 GHz"
l3="AMI-LA 16 GHz"

###############################################################
#*************************************************************#
###############################################################

print (' ______   ___   _  ____ ____ ___ __  __ ')
print ('/ ___\ \ / / \ | |/ ___/ ___|_ _|  \/  |')
print ('\__ _ \\ V /|  \| | |   \___ \| || |\/| |')
print (' ___) || | | |\  | |___ ___) | || |  | |')
print ('|____/ |_| |_| \_|\____|____/___|_|  |_| ___ raw')
print ('Version 1.0')
print ('Rob Fender 2020')

print ('Git desktop test')

import matplotlib.pyplot as plt

# the magic command which makes all the plots appear inline in the notebook
# %matplotlib inline 
import numpy as np
import math

################### physical constants #############################
pi=3.14159
e=4.803e-10 # electron charge cgs
c=2.998e10 # speed of light cgs
m=9.108e-28 # electron mass cgs
me=8.2e-7 # rest mass energy of electron in erg (511 keV)
mp=1.7e-24 # proton rest mass cgs

########## Exactly calculate Pacholczyk c1-c7 ##################

# constant c1 is just a combination of other constants
# 'check's are the actual formulae, you can print these to check if you like
c1=6.27e18
c1check=(3*e)/(4*pi*(m**3)*(c**5))

# same for c2
c2=2.37e-3
c2check=(2*e**4)/(3*(m**4)*(c**7))

# and c3
c3=1.87e-23
c3check=((3**0.5)*(e**3))/(4*pi*m*(c**2))

# c4 and c7 are combinations of these, and hence also constant
c4=2.63e26
c7=1.12e24

# Exact calculation of c5 and c6 rather than using tabulated versions
c5=((3**0.5)*(e**3)*(p+(7.0/3.0))*math.gamma((3*p-1)/12)*math.gamma((3*p+7)/12))/(16*pi*m*(c**2)*(p+1))
c6=((3**0.5)*pi*e*(m**5)*(c**10)*(p+(10.0/3.0))*math.gamma((3*p+2)/12)*math.gamma((3*p+10)/12))/(72.0)

#################################################################

# solid angle for a source of size s and distance d (small angle approximation)
def solidangle(r,d):
    omega=np.pi * (r**2.0) / (d**2.0)
    return omega # in steradians

# calculate nu1, frequency at which optical depth tau = 1
def nutau1 (r,p,N0,B):
    nu1= 2 * c1 * (r * c6)**(2/(p+4)) * N0**(2/(p+4)) * B**((p+2)/(p+4))
    return nu1 # in Hz

# Source function at nu1
def SFtau1 (B,nu1):
    SF=(c5 / c6 ) * B**(-0.5) * (nu1 / (2 * c1))**(2.5)
    return SF

# The following function will return the flux density as a function
# of frequency, nu1, source function and source angular size, in mJy
def F(freq,nu1,SF,omega):
    z=freq/nu1
    j=z**2.5 * (1 - np.exp(-z**(-(p+4)/2)))
    I=SF*j
    flux=I*omega*1.0e26
    return flux

# Determine initial conditions of new blob based upon input total electron energy,
# electron spectrum bounds, size, and deviation from equipartition

E10=gamma10*me # lower bound on electron energy converted to erg
E20=gamma20*me # upper bound on electron energy converted to erg

V=((4/3)*pi*(r0**3)) # volume

I2=(E20**(2-p)-E10**(2-p))/(2-p) # energy calculation integral
N00=Ee/(V * I2) # initial value of N0, the energy density density

#############################################################

# Integrated number of electrons, and corresponding mass
# for one proton for each electron

I1=(E20**(1-p)-E10**(1-p))/(1-p) # particle number calculation integral
Ntot=V * I1 * N00

Mtot=Ntot*mp # mass of the cold protons

#############################################################

# Determine B from ratio to equipartition, having established Eb

B0=(8 *np.pi * Eb / V)**(0.5)

# a constant for determining cyclotron frequencies
csynch=2.8e6*B0

#############################################################

print ('*****************************************************************************')
print ('Initial radius:','% .2e' % r0,'cm')
print ('Total energy in (electrons, magnetic field):','% .2e' % Ee,',','% .2e' % Eb,'erg')
print ('Upper and lower Lorentz factors','%.2e' % gamma10,',','%.2e' % gamma20)
print ('Corresponding initial number density density:','%.2e' % N00)
print ('Total number of electrons:','% .2e' % Ntot,'and mass in cold protons:','% .2e' % Mtot,'g')
print ('The magnetic field','% .2e' % B0, 'G')
print ('In this field the power law spectrum extends from','% .2e' % (gamma10*csynch),'to', '% .2e' % (gamma20*csynch), 'Hz')
print ('*****************************************************************************')

############################################
##                                        ##
##                 Part 2                 ##
##      Calculating the light curves      ##
##                                        ##
############################################

# Here we will calculate how the flux at a given frequency evolves with time
# as some synchrotron-emitting blob expands

# switch expansion speed, distance and observing freqs to real units
v=vc*c
d=dkpc*3.1e21 
n1=n1ghz*1.0e9
n2=n2ghz*1.0e9
n3=n3ghz*1.0e9

# here is a very simple function for radius as a function of time, 
# but you can make it as complicated as you like
def radius(t):   
# this is a power-law expansion. Linear expansion (e.g. vdL) gives a=1.0   
    a=1.0
    r=r0+(v*(t**a))
    return r
    
# similarly here is a very simple prescription for how B varies with radius
# we set a simple parameter b for power-law evolution
# in van der Laan it varies as R^{-2}, i.e. b=-2.0 which preserves equipartition ratio

def Br(r):
    b=-2.0 
    B=B0*((r/r0)**b)
    return B

# ... and how N0 varies (van der Laan is R^{-(2+p)})
def N0r(r):
    N0=N00*((r/r0)**(-(2+p)))
    return N0

# total mass in protons (should be conserved!!)
def MNr(N0,E10,E20,R):
    V=1.33 * np.pi * (R**3)
    I1=(E20**(1-p)-E10**(1-p))/(1-p)
    Ntot=V * I1 * N0
    Mtot=Ntot*mp
    return Mtot,Ntot


tmax= 100000#100000 # time range
tsteps=100000 # number of steps
tstep=tmax/tsteps # size of each step
t=np.linspace(1.0,tmax,tsteps)

rr=radius(t) # radius as a function of time # NB these functions make COPIES so t is not affected if they are changed
             # I make this clear because of course slices of nparrays are not copies but are just VIEWS
br=Br(rr) # magnetic fields as a function of time

# cyclotron freq, just for information (and as sanity check that cycl. is always self-absorbed)
cyc=2.8e6 * br

nr=N0r(rr)
# energy of upper and lower parts of electron spectrum, drops as blob expands
gamma1=gamma10*(r0/rr)
gamma2=gamma20*(r0/rr)
# do not allow Lorentz factor to be less than 1
gamma1[gamma1<1]=1
gamma2[gamma2<1]=1
# convert to physical units
E1=gamma1*me
E2=gamma2*me

# corresponding upper and lower synchrotron emission frequencies
nu1=2.8e6 * ((E1/me)**2)*br
nu2=2.8e6 * ((E2/me)**2)*br

# mass in, and number of, protons (for one proton for each electron)
M,N=MNr(nr,E1,E2,rr)

# solid angle
omega=pi * (rr**2.0) / (d**2.0)

# tau=1 frequency, and source function at that frequency
nu1r=nutau1(rr,p,nr,br)
sfr=SFtau1(br,nu1r)

# add here an optical depth term (for p=2.0, so slighty approximate for p not 2) [Pacholcyzk Equation 3.52]
# absorption coefficients
kappa1=c6 * nr * (br**2) * ((n1/(2 * c1))**(-3))
kappa2=c6 * nr * (br**2) * ((n2/(2 * c1))**(-3))
kappa3=c6 * nr * (br**2) * ((n3/(2 * c1))**(-3))
kappac=c6 * nr * (br**2) * ((cyc/(2 * c1))**(-3))

# optical depths
tau1=kappa1*rr
tau2=kappa2*rr
tau3=kappa3*rr
tauc=kappac*rr

# internal energy in magnetic fields
eb=((br**2)/(8*np.pi))*(1.33*pi*(rr**3))
# internal energy in particles
ee=nr*(1.33*pi*(rr**3))*(E2**(2.0-p)-E1**(2.0-p))/(2.0-p)

# maximum values
ebm=np.amax(eb)
eem=np.amax(ee)

# y are flux density curves
y1=F(n1,nu1r,sfr,omega)
y2=F(n2,nu1r,sfr,omega)
y3=F(n3,nu1r,sfr,omega)

lc1=np.vstack((t,y1)).T # this stacks time on top of flux then transposes array
lc2=np.vstack((t,y2)).T # which gives you nice set of (flux, time) pairs
lc3=np.vstack((t,y3)).T
lccombined=np.vstack((t,y1,y2,y3)).T # produces a combined array of time, fluxes

np.savetxt('lc1.csv', lc1, delimiter=',')
np.savetxt('lc2.csv', lc2, delimiter=',')
np.savetxt('lc3.csv', lc3, delimiter=',')
#np.savetxt('lccombined.csv', outputlc, delimiter=',') #NOTE NO OUTPUTLC DEFINED
# following this^ approach you can output any arrays in any way you like

print('')
print('Lightcurves saved in files lc[1-3].csv and lccombined.csv')
print('')

# We have calculated a LOT of different arrays here which we could potentially plot and output!

############################################
##                                        ##
##                 Part 3                 ##
##                 Plots                  ##
##                                        ##
############################################

# easy to calculate the light curves with these functions
# takes more lines to plot it nicely! in three different formats

plots=True # you can turn the plots on and off (they can take time)
if plots:

    # a row of three figures
    plt.figure(figsize=(14,4))
    plt.subplot(1,3,1)
    plt.plot(t,y1,'r', label=l1)
    plt.plot(t,y2,'g', label=l2)
    plt.plot(t,y3,'b', label=l3)
    plt.xlabel("Time")
    plt.ylabel("Flux density (mJy)")
    plt.title('Linear-Linear')
    plt.legend()

    plt.subplot(1,3,2)
    plt.yscale('log')
    plt.plot(t,y1,'r', label=l1)
    plt.plot(t,y2,'g', label=l2)
    plt.plot(t,y3,'b', label=l3)
    plt.xlabel("Time")
    plt.ylabel("Flux density (mJy)")
    plt.title('Log-Lin')
    plt.legend()

    plt.subplot(1,3,3)
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(t,y1,'r', label=l1)
    plt.plot(t,y2,'g', label=l2)
    plt.plot(t,y3,'b', label=l3)
    plt.xlabel("Time")
    plt.ylabel("Flux density (mJy)")
    plt.title('Log-Log')
    plt.legend()

    plt.tight_layout()
    plt.savefig("lightcurves.png")
    plt.show()

    
    

    # a(nother) row of three figures
    plt.figure(figsize=(14,4))
    plt.subplot(1,3,1)
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(t,rr,'r',label=r'Radius')
    plt.xlabel("Time")
    plt.ylabel("")
    plt.title('Evolution of blob size')
    plt.legend()

    plt.subplot(1,3,2)
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(t,tau1,'r',label=r'$\tau_{%.2f}$' % n1ghz) #Alex edited these
    plt.plot(t,tau2,'g',label=r'$\tau_{%.2f}$' % n2ghz)
    plt.plot(t,tau3,'b',label=r'$\tau_{%.2f}$' % n3ghz)
    plt.plot(t,tauc,'--',label=r'$\tau_{cycl.}$')
    plt.xlabel("Time")
    plt.ylabel("")
    plt.title("Optical depths")
    plt.legend()
    
    plt.subplot(1,3,3)
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(t,nu1r,'r',label=r'$\nu_1$')
    plt.xlabel("Time")
    plt.ylabel("")
    plt.title(r"Freq ($\tau=1$)")
    plt.legend()   

    plt.tight_layout()
    plt.savefig("anotherfigure.png")
    plt.show()

    # a(nother) row of three figures

    plt.figure(figsize=(14,4))
    plt.subplot(1,3,1)
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(t,gamma1,'r',label=r'$\gamma_{lower}$')
    plt.plot(t,gamma2,'g', label=r'$\gamma_{upper}$')
    plt.xlabel("Time")
    plt.ylabel("")
    plt.title('Evolution of electron spectrum')
    plt.legend()

    plt.subplot(1,3,2)
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(t,N,'r',label='N')
    plt.plot(t,M,'g', label='M')
    plt.xlabel("Time")
    plt.ylabel("")
    plt.title('Relativistic particle number')
    plt.legend()
    
    plt.subplot(1,3,3)
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(t,ee,'r',label='E(e)')
    plt.plot(t,eb,'g', label='E(B)')
    plt.xlabel("Time")
    plt.ylabel("")
    plt.title('Energy in e- and field')
    plt.legend()

    plt.tight_layout()
    plt.savefig("anotherfigure.png")
    plt.show()

############################################
##                                        ##
##                 Part 4                 ##
##      Analysis of the light curves      ##
##    Peaks and assoc. phys. parameters   ##
##                                        ##
############################################


# Analysis of the light curves
# np.argmax gives the index of the peak value of the light curve

def peaks(f,y):
    t=np.argmax(y) # index of peak 
    time=t*tstep
    peak=y[t] # flux density at peak 
    r=rr[t] # radius at peak
    b=br[t] # magnetic field at peak 
    nghz=f/1.e9 # freq in GHz
    etot=ee[t]+eb[t]
    print ('Time of peak at','% .2f' % nghz,'GHz:','% .2f' % time, 'sec, at which time peak flux density was','% .2e' % peak,' mJy')
    print ('Corresponding radius:','% .2e' % r,' cm, B:','% .2e' % b,'G, E(total):','% .2e' % etot,' erg.')
    print()
    return time,peak,nghz,r,b,etot
    
results=[]

# note that the 2nd column output is the time multiplier required to switch to seconds
# (some times are more readily recorded in hours/days etc)

time,peak,nghz,r,b,etot=peaks(n3,y3)
results.append([time,1.0,peak,nghz,r,b,etot])

time,peak,nghz,r,b,etot=peaks(n2,y2)
results.append([time,1.0,peak,nghz,r,b,etot])

time,peak,nghz,r,b,etot=peaks(n1,y1)
results.append([time,1.0,peak,nghz,r,b,etot])

# output data file which analysis package can work on
npr=np.array(results)
np.savetxt('sim_results.dat',npr,fmt='%.4e', delimiter=',')
#np.savetxt('sim_results1.dat',npr,fmt='%.4e', delimiter=',') # specify some particular path

print('')
print('Peaks analysis saved in sim_results.dat')
print('')
print('')
print('Fin.')
print('')
