'''
before you do anything, open terminal, type:
conda activate bb
cd Documents\VSCodeStuff\SNeDustFitting > NVM I changed the VScode explorer to have this folder as the default > it also auto runs activate bb
'''
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, FloatSlider
from scipy.constants import h, c, k

def B_nu(nu, T):
    L = ((2*h*nu**3)/c**2)/(np.exp((h*nu)/(k*T))-1)
    return L

def B_lam(lam, T):
    L = (2*h*c**2/lam**5)/(np.exp((h*c)/(lam*k*T))-1)
    return L

lam = np.logspace(-8, -2, 500) #reverses order of just defined lam
nu = c/lam

#fxn for betelgeuse expected flux, with known T, r, and d
def F_lam_betelgeuse(lam, T):
    L = (2*h*c**2/lam**5)/(np.exp((h*c)/(lam*k*T))-1) * np.pi * ((6.957e5*640)/(408*9.461e15)) #betelgeuse is 640 solar radii + 408 ly away
    return L

#fxn for sirius expected flux measurement across continuum
def F_lam_siriusa(lam, T):
    L = (2*h*c**2/lam**5)/(np.exp((h*c)/(lam*k*T))-1) * np.pi * ((6.957e5*1.714)/(8.61*9.461e15)) #sirius is 1.744 solar radii + 8.61 ly away
    return L

#fxn for sirius b expected flux, with known T, r, and d
def F_lam_siriusb(lam, T):
    L = (2*h*c**2/lam**5)/(np.exp((h*c)/(lam*k*T))-1) * np.pi * ((6.957e5*0.0081)/(8.61*9.461e15)) #sirius is 0.0081 solar radii + 8.61 ly away
    return L

#fxn for Proxima centauri expected flux, with known T, r, and d
def F_lam_pcent(lam, T):
    L = (2*h*c**2/lam**5)/(np.exp((h*c)/(lam*k*T))-1) * np.pi * ((6.957e5*0.154)/(9.461e15*4.25)) #pcent is 0.154 solar radii + 4.25 ly away
    return L

#fxn for Rigel expected flux, with known T, r, and d
def F_lam_rigel(lam, T):
    L = (2*h*c**2/lam**5)/(np.exp((h*c)/(lam*k*T))-1) * np.pi * ((6.957e5*74)/(9.461e15*850)) #rigel is 74 solar radii + 850 ly away
    return L

#fxn for alpha centauri a expected flux, with known T, r, and d
def F_lam_acenta(lam, T):
    L = (2*h*c**2/lam**5)/(np.exp((h*c)/(lam*k*T))-1) * np.pi * ((6.957e5*1.223)/(9.461e15*4.34)) #acenta is 1.223 solar radii + 4.34 ly away
    return L

#fxn for Aldebaran expected flux, with known T, r, and d
def F_lam_ald(lam, T):
    L = (2*h*c**2/lam**5)/(np.exp((h*c)/(lam*k*T))-1) * np.pi * ((6.957e5*45)/(9.461e15*67)) #aldebaran is 45 solar radii + 67 ly away
    return L



#plotting figure 1: ideal blackbodies
plt.figure(1)
for T in [100, 500, 1000, 5000, 10000]:
    plt.loglog(lam*1e6, B_lam(lam, T), label=f'{T} K') #f before smth allows variables via curly brackets
plt.xlabel(r"$\lambda$ in $\mu$m") #encompass backslash+greek letter spelledout in dollar signs. dont forget r, meaning 'raw' string > this is lateX stuff
plt.ylabel(r"$B_\lambda$ (W m$^{-2}$ sr$^{-1}$ m$^{-1}$)")
plt.ylim(1e-2, 1e16)
plt.title(r'Ideal Blackbody Spectral Radiance vs. $\lambda$')
plt.legend(title='Temps')

plt.figure(2)
plt.loglog(lam*1e6, F_lam_betelgeuse(lam, 3500), label='Betelgeuse (red supergiant, type M)')
plt.loglog(lam*1e6, F_lam_siriusa(lam, 9800), label='Sirius A (type A)')
plt.loglog(lam*1e6, F_lam_siriusb(lam, 25000), label='Sirius B (white dwarf, type DA)')
plt.loglog(lam*1e6, F_lam_rigel(lam, 12000), label='Rigel (blue supergiant, type B)')
plt.loglog(lam*1e6, F_lam_pcent(lam, 3000), label='Proxima Centauri (type M)')
plt.loglog(lam*1e6, F_lam_acenta(lam, 5700), label='Alpha Centauri A (type G)')
plt.loglog(lam*1e6, F_lam_ald(lam, 3900), label='Aldebaran (red giant, type K)')
plt.ylabel(r'$F_\lambda$ (W m$^{-2}$ m$^{-1}$ )')
plt.xlabel(r"$\lambda$ in $\mu$m")
plt.title(r'Expected stellar flux measurement across all $\lambda$ for known stars')
plt.ylim(1e-2, 1e6)
plt.legend()

plt.figure(3)


plt.show()
