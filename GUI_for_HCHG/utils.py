import numpy as np
from scipy import constants

def photonEnergyGrid(pulse):
    """
    Convert the grid corresponding to the complex electric field in frequency
    domain from angular frequency to photon energy in eV.
    
    Photon energy:
    E[J] = hbar * w = h * v
    E[eV] = E[J] / e
    
    Input(s):
    pulse: a pulse instance on which the conversion is performed.
    Output(s):
    W_eV: photon energy grid [eV].
    """
    
    W_eV = constants.h * pulse.F_mks / constants.e
    
    return W_eV

def dB(num):
    return 10 * np.log10(np.abs(num)**2)

def super_gaussian(x, x0, w, P, A=1):
    """
    Generate a Super-Gaussian function.

    Input(s):
    x (numpy.ndarray): The array over which the Super-Gaussian function is evaluated.
    x0 (float): The center of the Super-Gaussian peak.
    w (float): The Full Width at Half Maximum (FWHM) of the peak.
    P (int): The order of the Super-Gaussian.
    A (float): The amplitude of the Super-Gaussian function. Default is 1.

    Output(s):
    super_gaussian (numpy.ndarray): The Super-Gaussian function evaluated at each point in x.
    """
    
    super_gaussian = A * np.exp(-np.log(2) * ((4 * (x - x0)**2) / w**2)**P)
    return super_gaussian

def center_window(window, width, height):
    screen_width = window.winfo_screenwidth()
    screen_height = window.winfo_screenheight()
    position_top = int(screen_height / 2 - height / 2)
    position_right = int(screen_width / 2 - width / 2)
    window.geometry(f"{width}x{height}+{position_right}+{position_top}")