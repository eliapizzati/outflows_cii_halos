"""This module contains natural constants in CGS units 

Translated from RADMC's IDL function problem_natconst.pro

List of natural constants:

====    =============
Name    Description
====    =============
gg      Gravitational constant
mp      Mass of proton [g]
me      Mass of electron [g]
kk      Bolzmann's constant [erg/K]
hh      Planck's constant [erg.s]
ee      Unit charge 
cc      Light speed [cm/s]
st      Thmpson cross-section [cm^2]
ss      Stefan-Boltzmann const [erg/cm^2/K^4/s]
aa      4 * ss / cc
muh2    Mean molecular weight (H2 + He + metals)
ev      Electronvolt [erg]
kev     Kilo electronvolt [erg]
micr    Micron [cm]
km      Kilometer [cm]
angs    Angstrom [cm]
ls      Solar luminosity [erg/s]
rs      Solar radius [cm]
ms      Solar mass [g]
ts      Solar effective temperature [K]
au      Astronomical unit [cm]
pc      Parsec [cm]
mea     Mass of Earth [g]
rea     Equatorila radius of Earth [cm]
mmo     Mass of Moon [g]
rmo     Radius of Moon [cm]
dmo     Distance earth-moon (center-to-center) [cm]
mju     Mass of Jupiter
rju     Equatorial radius of Jupiter [cm]
dju     Distance Jupiter-Sun [cm]
year    Year [s]
hour    Hour [s]
day     Day  [s]
====    =============
(MY STUFF AT THE END)
gamma   Adiabatic index
mus      Mean molecular weight
knorm   Normalization constant for temperature
A_C     Carbon abundance
A_H     Hydrogen abundance
Gamma_H_1000  H ionization rate at 1 kpc
Gamma_He_1000 He ionization rate at 1 kpc
Gamma_LW_1000 LW ionization rate at 1 kpc
Gamma_CI_1000 CI ionization rate at 1 kpc
Gamma_CI_EUV_1000 CI ionization rate (ionizing) at 1 kpc
Gamma_CI_FUV_1000 CI ionization rate (non ionizing) at 1 kpc
Gamma_CII_1000 CII ionization rate at 1 kpc
intensity_UV_1000  Intensity of the UV at 1 kpc
line_CII_rest_frame Rest frame wavelength of the CII line
line_UV_mixing_rest_frame Rest frame wavelength of the UV mixing line
CMB_temp  CMB temperature
T_star   CII line temperature
A_coeff_ul CII line upper level Einstein coefficient
A_coeff_ku CII line upper level Einstein coefficient for mixing
A_coeff_kl CII line lower level Einstein coefficient for mixing
"""
gg = 6.672e-8   # Gravitational constant
mp = 1.6726e-24  # Mass of proton [g]
me = 9.1095e-28  # Mass of electron [g]
kk = 1.3807e-16  # Bolzmann's constant [erg/K]
hh = 6.6262e-27  # Planck's constant [erg.s]
ee = 4.8032e-10  # Unit charge
cc = 2.99792458e10  # Light speed [cm/s]
st = 6.6524e-25  # Thmpson cross-section [cm^2]
ss = 5.6703e-5   # Stefan-Boltzmann const [erg/cm^2/K^4/s]
aa = 7.5657e-15  # 4 * ss / cc
#
# Gas constants
#
muh2 = 2.3000e0   # Mean molecular weight (H2 + He + metals)
#
# Alternative units
#
ev = 1.6022e-12  # Electronvolt [erg]
kev = 1.6022e-9   # Kilo electronvolt [erg]
micr = 1e-4        # Micron [cm]
km = 1e5         # Kilometer [cm]
angs = 1e-8        # Angstrom [cm]
#
# Astronomical constants
#
ls = 3.8525e33    # Solar luminosity [erg/s]
rs = 6.96e10      # Solar radius [cm]
ms = 1.99e33      # Solar mass [g]
ts = 5.780e3      # Solar effective temperature [K]
au = 1.496e13     # Astronomical unit [cm]
pc = 3.08572e18   # Parsec [cm]
mea = 5.9736e27   # Mass of Earth [g]
rea = 6.375e8     # Equatorila radius of Earth [cm]
mmo = 7.347e25    # Mass of Moon [g]
rmo = 1.738e8     # Radius of Moon [cm]
dmo = 3.844e10    # Distance earth-moon (center-to-center) [cm]
mju = 1.899e30    # Mass of Jupiter
rju = 7.1492e9    # Equatorial radius of Jupiter [cm]
dju = 7.78412e13  # Distance Jupiter-Sun [cm]
#
# Time units
#
year = 3.1536e7    # Year [s]
hour = 3.6000e3    # Hour [s]
day = 8.6400e4    # Day  [s]
#
#
#
# MY STUFF
#
#
#
# gas constants
gamma = 5./3
mus = 0.61 # mean molecular weight of the sun
knorm = kk/(mus*mp) 
A_C = 2.69e-4 # carbon percentage for Z=1
A_H = 0.76 # hydrogen percentage for Z=1
#
# photoionization rates (and UV intensity) at 1 kpc (from gal flux)
Gamma_H_1000 = 5.48031935502901e-09 # s^-1
Gamma_He_1000 = 1.7687762344020628e-09 # s^-1
Gamma_CI_1000 = 5.423253133640651e-08 # s^-1
Gamma_CI_FUV_1000 = 3.84971348468091e-08 # s^-1
Gamma_CI_EUV_1000 = 1.5249362754606713e-08 # s^-1
Gamma_CII_1000 = 9.56484853366846e-10 # s^-1
Gamma_LW_1000 = 1.4229125141877616e-08 # s^-1
intensity_UV_1000 = 1.3492085033628848e-17 # erg/s/cm^2/Hz/sr

#
# CMB suppression coefficients
line_CII_rest_frame = 1900*1e9 #Hz
line_UV_mixing_rest_frame = 22.4452673584*1e14 #Hz
CMB_temperature = 2.725 # K
T_star = 91.7 # K
A_coeff_ul = 2.36e-6 # s^-1
A_coeff_ku = 4.76e7 # s^-1
A_coeff_kl = 2.41e8 # s^-1
#
