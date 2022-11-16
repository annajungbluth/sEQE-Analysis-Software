import math

import numpy as np
from numpy import exp
from scipy.interpolate import interp1d

from source.compilation import compile_EQE
from source.utils import R_squared


# -----------------------------------------------------------------------------------------------------------

# Function to calculate gaussian absorption

def calculate_gaussian_absorption(x, f, l, E, T):
    """Function to calculate gaussian absorption

    Parameters
    ----------
    x : list, required
        List of energy values [eV]
    f : float, required
        Oscillator strength [eV^2]
    l : float, required
        Reorganization energy [eV]
    E : float, required
        Peak energy [eV]
    T : float or int, required
        Temperature [K]
        
    Returns
    -------
    EQE : float
        Calculated EQE value
    """

    # Define variables
    k = 8.617 * math.pow(10, -5)  # [ev/K]

    return (f / (x * math.sqrt(4 * math.pi * l * T * k))) * exp(-(E + l - x) ** 2 / (4 * l * k * T))

# -----------------------------------------------------------------------------------------------------------

# Function to calculate gaussian absorption including disorder

def calculate_gaussian_disorder_absorption(x, f, l, E, sig, T):
    """Function to calculate gaussian absorption including disorder

    Parameters
    ----------
    x : list, required
        List of energy values [eV]
    f : float, required
        Oscillator strength [eV^2]
    l : float, required
        Reorganization energy [eV]
    E : float, required
        Peak energy [eV]
    sig : float, required
        Peak disorder parameter [eV]
    T : float or int, required
        Temperature [K]
        
    Returns
    -------
    EQE : float
        Calculated EQE value
    """

    # Define variables
    k = 8.617 * math.pow(10, -5)  # [ev/K]

    return (f / (x * math.sqrt(2 * math.pi * (2 * l * T * k + sig ** 2))) * exp(
        -(E + l - x) ** 2 / (4 * l * k * T + 2 * sig ** 2)))

# -----------------------------------------------------------------------------------------------------------

    # Function to calculate absorption using MLJ theory

def calculate_MLJ_absorption(x, f, l, E, T, S, hbarw):
    """Function to calculate absorption using MLJ equation

    Parameters
    ----------
    x : list, required
        List of energy values [eV]
    f : float, required
        Oscillator strength [eV^2]
    l : float, required
        Reorganization energy [eV]
    E : float, required
        Peak energy [eV]
    T : float or int, required
        Temperature [K]
    S : float, required
        Huang-Rhys parameter
    hbarw : float, required
        Vibrational energy [eV]
        
    Returns
    -------
    EQE : float
        Calculated EQE value
    """

    # Define variables
    k = 8.617 * math.pow(10, -5)  # [ev/K]

    EQE = 0
    for n in range(0, 6):
        EQE_n = (f / (x * math.sqrt(4 * math.pi * l * T * k))) \
                * (math.exp(-S) * S ** n / math.factorial(n)) \
                * exp(-(E + l - x + n * hbarw) ** 2 \
                      / (4 * l * k * T))
        EQE += EQE_n
    return EQE

# -----------------------------------------------------------------------------------------------------------

    # Function to calculate absorption including disorder using MLJ theory

def calculate_MLJ_disorder_absorption(x, f, l, E, T, sig, S, hbarw): 
    """Function to calculate absorption using MLJ equation including disorder

    Parameters
    ----------
    x : list, required
        List of energy values [eV]
    f : float, required
        Oscillator strength [eV^2]
    l : float, required
        Reorganization energy [eV]
    E : float, required
        Peak energy [eV]
    T : float or int, required
        Temperature [K]
    sig : float, required
        Peak disorder parameter [eV]
    S : float, required
        Huang-Rhys parameter
    hbarw : float, required
        Vibrational energy [eV]
        
    Returns
    -------
    EQE : float
        Calculated EQE value
    """

    # Define variables
    k = 8.617 * math.pow(10, -5)  # [ev/K]

    EQE = 0
    for n in range(0, 6):
        EQE_n = (f / (x * math.sqrt(2 * math.pi * (2 * l * T * k + sig ** 2))) \
                 * (math.exp(-S) * S ** n / math.factorial(n)) \
                 * exp(-(E + l - x + n * hbarw) ** 2 \
                       / (4 * l * k * T + 2 * sig ** 2)))
        EQE += EQE_n
    return EQE

# -----------------------------------------------------------------------------------------------------------

# Function to calculate parameters for double peak fit

def calculate_combined_fit(eqe,
                           stopE,
                           best_vals_Opt,
                           best_vals_CT,
                           T,
                           R2_Opt=None,
                           R2_CT=None,
                           include_disorder=False,
                           bias=False,
                           tolerance=0,
                           range=1.05
                           ):
    """Function to calculate parameters for double peak fit

    Parameters
    ----------
    eqe : list, required
        List of input EQE values [eV]
    stopE : float, required
        Stop energy of fit [eV]
    best_vals_Opt : list, required
        Optical gap/S1 peak fit values
    best_vals_CT : list, required
        CT state fit values
    T : float, required
        Temperature [K]
    R2_Opt : float, optional
        R squared of optical peak fit
    R2_CT : float, optional
        R squared of CT state fit
    include_disorder : bool, optional
        Boolean value specifying whether to include CT state disorder
    bias : bool, optional
        Boolean value specifying whether to bias fits above the data
    range : float, optional
        Defines upper bound of R2 calculation
        
    Returns
    -------
    result_dict : dict
        Dictionary with fit results
        Dict keys:
                R2_Combined: R squared of sum of CT and Opt fit [float]
                R2_CT: R squared of CT fit [float]
                R2_Opt: R squared of Opt fit [float]
                R2_Average: Average R squared of CT / Opt / Combined fit [float]
                Combined_Fit: Sum of CT and Opt fit [list]
                Opt_fit: Opt fit values [list]
                CT_fit: CT fit values [list]
                Energy: Energy values [list]
                EQE: Original EQE data [list]
    """

    wave_data, energy_data, eqe_data, log_eqe_data = compile_EQE(eqe,
                                                                 min(eqe['Energy']),
                                                                 stopE * range,  # NOTE: Increase stop energy to expand fit!
                                                                 1)
    # # Optional code to add interpolation to the data
    # int_func = interp1d(eqe['Energy'], eqe['EQE'])
    # energy_data = np.arange(min(eqe['Energy']), stopE * range, 1)
    # eqe_data = int_func(energy_data)

    if sum(best_vals_Opt) != 0 and sum(best_vals_CT) != 0:
        Opt_fit = np.array([calculate_gaussian_absorption(e,
                                                          best_vals_Opt[0],
                                                          best_vals_Opt[1],
                                                          best_vals_Opt[2],
                                                          T)
                            for e in energy_data])
        if R2_Opt is None:
            R2_Opt = R_squared(y_data=eqe_data,
                               yfit_data=Opt_fit.tolist(),
                               bias=bias,
                               tolerance=tolerance)
        if include_disorder:
            CT_fit = np.array([calculate_gaussian_disorder_absorption(e,
                                                                      best_vals_CT[0],
                                                                      best_vals_CT[1],
                                                                      best_vals_CT[2],
                                                                      best_vals_CT[3],
                                                                      T)
                               for e in energy_data])

        else:
            CT_fit = np.array([calculate_gaussian_absorption(e,
                                                             best_vals_CT[0],
                                                             best_vals_CT[1],
                                                             best_vals_CT[2],
                                                             T)
                               for e in energy_data])
        if R2_CT is None:
            R2_CT = R_squared(y_data=eqe_data,
                              yfit_data=CT_fit.tolist(),
                              bias=bias,
                              tolerance=tolerance)

        if len(Opt_fit) == len(CT_fit):
            combined_Fit = Opt_fit + CT_fit
            combined_R_Squared = R_squared(y_data=eqe_data,
                                           yfit_data=combined_Fit.tolist(),
                                           bias=bias,
                                           tolerance=tolerance)

    else:  # if any of the fits were unsuccessful
        Opt_fit = 0
        CT_fit = 0
        combined_Fit = 0
        combined_R_Squared = 0

    average_R_Squared = (R2_CT + R2_Opt + combined_R_Squared) / 3

    result_dict = {'R2_Combined': combined_R_Squared,
                   'R2_CT': R2_CT,
                   'R2_Opt': R2_Opt,
                   'R2_Average': average_R_Squared,
                   'Combined_Fit': combined_Fit,
                   'Opt_Fit': Opt_fit,
                   'CT_Fit': CT_fit,
                   'Energy': energy_data,
                   'EQE': eqe_data
                   }

    return result_dict

# -----------------------------------------------------------------------------------------------------------

# Function to calculate parameters for double peak MLJ fit

def calculate_combined_fit_MLJ(eqe,
                               stopE,
                               best_vals_Opt,
                               best_vals_CT,
                               T,
                               S,
                               hbarw,
                               R2_Opt=None,
                               R2_CT=None,
                               include_disorder=False,
                               bias=False,
                               tolerance=0,
                               range=1.05
                               ):
    """Function to calculate parameters for double peak MLJ fit

    Parameters
    ----------
    eqe : list, required
        List of input EQE values [eV]
    stopE : float, required
        Stop energy of fit [eV]
    best_vals_Opt : list, required
        Optical gap/S1 peak fit values
    best_vals_CT : list, required
        CT state fit values
    T : float, required
        Temperature [K]
    S : float, required
        Huang-Rhys parameter
    hbarw : float, required
        Vibrational energy [eV]
    R2_Opt : float, optional
        R squared of optical peak fit
    R2_CT : float, optional
        R squared of CT state fit
    include_disorder : bool, optional
        Boolean value specifying whether to include CT state disorder
    bias : bool, optional
        Boolean value specifying whether to bias fits above the data
    tolerance : float, optional
        Tolerance for fit above the data
    range : float, optional
        Defines upper bound of R2 calculation
        
    Returns
    -------
    result_dict : dict
        Dictionary with fit results
        Dict keys:
                R2_Combined: R squared of sum of CT and Opt fit [float]
                R2_CT: R squared of CT fit [float]
                R2_Opt: R squared of Opt fit [float]
                R2_Average: Average R squared of CT / Opt / Combined fit [float]
                Combined_Fit: Sum of CT and Opt fit [list]
                Opt_fit: Opt fit values [list]
                CT_fit: CT fit values [list]
                Energy: Energy values [list]
                EQE: Original EQE data [list]
    """

    wave_data, energy_data, eqe_data, log_eqe_data = compile_EQE(eqe,
                                                                 min(eqe['Energy']),
                                                                 stopE * range,  # NOTE: Increase stop energy to expand fit!
                                                                 1)
    # # Optional code to add interpolation to the data
    # int_func = interp1d(eqe['Energy'], eqe['EQE'])
    # energy_data = np.arange(min(eqe['Energy']), stopE * range, 1)
    # eqe_data = int_func(energy_data)

    if sum(best_vals_Opt) != 0 and sum(best_vals_CT) != 0:
        Opt_fit = np.array([calculate_gaussian_absorption(e,
                                                          best_vals_Opt[0],
                                                          best_vals_Opt[1],
                                                          best_vals_Opt[2],
                                                          T)
                            for e in energy_data])
        if R2_Opt is None:
            R2_Opt = R_squared(y_data=eqe_data,
                               yfit_data=Opt_fit.tolist(),
                               bias=bias,
                               tolerance=tolerance)
        if include_disorder:
            CT_fit = np.array([calculate_MLJ_disorder_absorption(x=e,
                                                                 f=best_vals_CT[0],
                                                                 l=best_vals_CT[1],
                                                                 E=best_vals_CT[2],
                                                                 sig=best_vals_CT[3],
                                                                 T=T,
                                                                 S=S,
                                                                 hbarw=hbarw)
                               for e in energy_data])

        else:
            CT_fit = np.array([calculate_MLJ_absorption(x=e,
                                                        f=best_vals_CT[0],
                                                        l=best_vals_CT[1],
                                                        E=best_vals_CT[2],
                                                        T=T,
                                                        S=S,
                                                        hbarw=hbarw)
                               for e in energy_data])
        if R2_CT is None:
            R2_CT = R_squared(y_data=eqe_data,
                              yfit_data=CT_fit.tolist(),
                              bias=bias,
                              tolerance=tolerance)

        if len(Opt_fit) == len(CT_fit):
            combined_Fit = Opt_fit + CT_fit
            combined_R_Squared = R_squared(y_data=eqe_data,
                                           yfit_data=combined_Fit.tolist(),
                                           bias=bias,
                                           tolerance=tolerance)

    else:  # if any of the fits were unsuccessful
        Opt_fit = 0
        CT_fit = 0
        combined_Fit = 0
        combined_R_Squared = 0

    average_R_Squared = (R2_CT + R2_Opt + combined_R_Squared) / 3

    result_dict = {'R2_Combined': combined_R_Squared,
                   'R2_CT': R2_CT,
                   'R2_Opt': R2_Opt,
                   'R2_Average': average_R_Squared,
                   'Combined_Fit': combined_Fit,
                   'Opt_Fit': Opt_fit,
                   'CT_Fit': CT_fit,
                   'Energy': energy_data,
                   'EQE': eqe_data
                   }

    return result_dict

# -----------------------------------------------------------------------------------------------------------
