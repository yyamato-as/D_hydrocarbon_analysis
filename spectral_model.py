import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as ac
import astropy.units as u
from scipy.interpolate import interp1d
import sqlite3

# constants
mH2 = 2.016*ac.u.to(u.g).value
h = ac.h.cgs.value
k_B = ac.k_B.cgs.value
c = ac.c.cgs.value
ckms = ac.c.to(u.km/u.s).value

##################### functions for spectroscopic data ######################
def get_partition_function(databasefile, entry):
    # execute sqlite command and get the partition function values
    conn = sqlite3.connect(databasefile)
    query_string = f'select * from partitionfunctions where PF_Name = "{entry}"'
    cursor = conn.cursor()
    cursor.execute(query_string)

    # fetch temperature info
    columns = list(map(lambda x: x[0], cursor.description))
    temperature = np.array([float("{}.{}".format(col.split("_")[1], col.split("_")[2])) for col in columns if any(chr.isdigit() for chr in col)])
    
    # fetch partition function info
    data = cursor.fetchone()
    pf = np.array([val for val, col in zip(data, columns) if any(chr.isdigit() for chr in col)])

    # close the connection
    cursor.close()
    conn.close()

    # interpolation
    Q = interp1d(temperature, pf, kind="cubic", fill_value="extrapolate")
    
    return Q

def get_spectroscopic_data(databasefile, entry, numin=75e9, numax=1000e9):
    # execute sqlite command and get the partition function values
    conn = sqlite3.connect(databasefile)

    query_string = f'select T_Frequency, T_EinsteinA, T_EnergyLower, T_UpperStateDegeneracy'
    query_string += f' from transitions where (T_Frequency >= {numin*1e-6} and T_Frequency <= {numax*1e-6} and (T_Name = "{entry}"))'
    query_string += f' order by T_Frequency'

    cursor = conn.cursor()
    cursor.execute(query_string)
    data = cursor.fetchall()

    # close
    cursor.close()
    conn.close()

    # post processes
    if len(data) == 0:
        return None
    else:
        nu0, Aul, El, gu = np.array(data).T
        nu0 *= 1e6
        El *= h * c / k_B
        return nu0, Aul, El, gu

##################### functions for spectral modeling ######################
def eta(source_size, beam):
    """calculate the beam filling factor. The unit of source_size and beam should be the same.

    Parameters
    ----------
    source_size : float or ndarray
        source size 
    beam : float or ndarray
        beam size

    Returns
    -------
    float or ndarray
        beam filling factor, eta
    """
    return source_size ** 2 / (beam ** 2 + source_size ** 2)

def calc_dust_optical_depth(nu, N_H2, kappa, beta, g2d=100., nu_ref=230e9):
    """calculate dust optical depth spectrum

    Parameters
    ----------
    nu : ndarray
        frequencies in Hz at which you want to calculate the optical depth of the dust
    N_H2 : float
        column density of H2 molecule in cm-2
    kappa : float
        dust opacity in cm2 g-1 at the reference frequency given by `nu_ref`
    beta : float
        spectral index, i.e., power-law index of the frequency dependence of the dust optical depth
    g2d : float, optional
        gas-to-dust mass ratio, by default 100
    nu_ref : float, optional
        reference frequency for the dust opacity in Hz, by default 230e9

    Returns
    -------
    ndarray
        dust optical depth spectrum
    """
    tau_d = (N_H2 * kappa * mH2 / g2d) * (nu / nu_ref) ** beta
    return tau_d

def line_profile_function(nu, nu0, sigma):
    """calculate line profile function. The returned ndarray shape will be (nu.size, nu0.size).

    Parameters
    ----------
    nu : ndarray
        frequencies in Hz at which you want to calculate phi
    nu0 : ndarray
        array of the central frequencies in Hz
    sigma : ndarray
        array of the frequency line width in Hz

    Returns
    -------
    ndarray with a shape of (nu.size, nu0.size)
        line profile function of each component
    """
    return 1. / ((2 * np.pi) ** 0.5 * sigma[None, :]) * np.exp(- 0.5 * (nu[:, None] - nu0[None, :]) ** 2 / sigma[None, :] ** 2)

def calc_line_optical_depth(nu, nu0, sigma, Tex, logN, Aul, gu, Eu, Q):
    """calculate line optical depth spectrum. The returned ndarray shape will be (nu.size, nu0.size).

    Parameters
    ----------
    nu : ndarray 
        Frequencies in Hz at which you want to calculate tau 
    nu0 : ndarray
        array of the central frequencies in Hz
    sigma : ndarray
        array of the frequency line width in Hz
    Tex : float or ndarray
        Excitation temperature in K of the molecule. Different Tex for each transition with ndarray is acceptable.
    logN : float or ndarray
        log10 of Column density in cm-2 of the molecule. Different Ntot for each transition with ndarray is acceptable.
    Aul : ndarray
        Einstein A coefficients in s-1 of the transitions
    gu : ndarray
        Upper state degeneracies of the transitions
    Eu : ndarray
        upper state enegies in K of the transitions
    Q : function
        partition function

    Returns
    -------
    1d array
        line optical depth profile
    """
    phi = line_profile_function(nu, nu0, sigma)
    tau_l = c ** 2 / (8 * np.pi * nu[:, None] ** 2) * Aul[None, :] * 10 ** logN
    tau_l *= gu[None, :] * np.exp(- Eu[None, :] / Tex) / Q(Tex)
    # tau_l *= (1 - np.exp(- h * nu0[None, :] / (k_B * Tex))) * phi
    tau_l *= (np.exp(h * nu0[None, :] / (k_B * Tex)) - 1) * phi
    return np.sum(tau_l, axis=1)

def Inu(nu, nu0, sigma, Tex, logN, size, beam, Aul, gu, El, Q, Tdust, N_H2, kappa, beta, g2d=100., nu_ref=230e9):
    """Intensity

    Parameters
    ----------
    nu : ndarray 
        Frequencies in Hz at which you want to calculate tau 
    nu0 : ndarray
        array of the central frequencies in Hz
    sigma : ndarray
        array of the frequency line width in Hz
    Tex : float or ndarray
        Excitation temperature in K of the molecule. Different Tex for each transition with ndarray is acceptable.
    logN : float or ndarray
        Column density in cm-2 of the molecule. Different Ntot for each transition with ndarray is acceptable.
    Aul : ndarray
        Einstein A coefficients in s-1 of the transitions
    gu : ndarray
        Upper state degeneracies of the transitions
    El : ndarray
        Lower state enegies in K of the transitions
    Q : function
        partition function
    Tdust : float 
        dust temperature in K
    N_H2 : float
        column density of H2 molecule in cm-2
    kappa : float
        dust opacity in cm2 g-1 at the reference frequency given by `nu_ref`
    beta : float
        spectral index, i.e., power-law index of the frequency dependence of the dust optical depth
    g2d : float, optional
        gas-to-dust mass ratio, by default 100
    nu_ref : float, optional
        reference frequency for the dust opacity in Hz, by default 230e9

    Returns
    -------
    1d array
        source function
    """
    # optical depths
    tau_d = calc_dust_optical_depth(nu, N_H2, kappa, beta, g2d=g2d, nu_ref=nu_ref)
    tau_l = calc_line_optical_depth(nu, nu0, sigma, Tex, logN, Aul, gu, El, Q)
    tau = tau_l + tau_d

    # source function
    Snu = (tau_l * Bnu(nu, Tex) + tau_d * Bnu(nu, Tdust)) / (tau_l + tau_d) 

    # return peak intensity
    return eta(size, beam) * (Snu - Bnu_CMB(nu)) * (1 - np.exp(-tau))

def Bnu(nu, T):
    """Planck function for blackbody radiation

    Parameters
    ----------
    nu : float or ndarray
        frequencies in Hz at which you wnat to calculate the blackbody
    T : float
        temperature in K

    Returns
    -------
    float or ndarray
        blackbody spectrum
    """
    return 2 * h * nu ** 3 / c ** 2 / (np.exp(h * nu / (k_B * T)) - 1.) * 1e23 # in Jy / sr

def Jnu(nu, T):
    """R-J brightness temperature for blackbody radiation

    Parameters
    ----------
    nu : float or ndarray
        frequencies in Hz at which you want to calculate the blackbody
    T : float
        temperature in K

    Returns
    -------
    float or ndarray
        R-J brightness temperature spectrum
    """
    return h * nu / k_B / (np.exp(h * nu / (k_B * T)) - 1.) # in K

def Bnu_CMB(nu, T_CMB=2.73):
    """Planck function for CMB radiation"""
    return Bnu(nu, T_CMB)

def Jnu_CMB(nu, T_CMB=2.73):
    """R-J brightness temperature for CMB radiation"""
    return Jnu(nu, T_CMB) 