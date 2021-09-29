import math

import numpy as np
from numpy import exp
from scipy.interpolate import interp1d

from source.compilation import compile_EQE
from source.utils import R_squared


# -----------------------------------------------------------------------------------------------------------

# Function to calculate gaussian absorption

def calculate_gaussian_absorption(x, f, l, E, T):
    """
    :param x: List of energy values [list]
    :param f: Oscillator strength [float]
    :param l: Reorganization Energy [float]
    :param E: Peak Energy [float]
    :param T: Temperature [float or int]
    :return: EQE value [float]
    """

    # Define variables
    k = 8.617 * math.pow(10, -5)  # [ev/K]

    return (f / (x * math.sqrt(4 * math.pi * l * T * k))) * exp(-(E + l - x) ** 2 / (4 * l * k * T))


# -----------------------------------------------------------------------------------------------------------

# Function to calculate gaussian absorption including disorder

def calculate_gaussian_disorder_absorption(x, f, l, E, sig, T):
    """
    :param x: List of energy values [list]
    :param f: Oscillator strength [float]
    :param l: Reorganization Energy [float]
    :param E: Peak Energy [float]
    :param sig: Peak disorder [float]
    :param T: Temperature [float or int]
    :return: EQE value [float]
    """

    # Define variables
    k = 8.617 * math.pow(10, -5)  # [ev/K]

    return (f / (x * math.sqrt(2 * math.pi * (2 * l * T * k + sig ** 2))) * exp(
        -(E + l - x) ** 2 / (4 * l * k * T + 2 * sig ** 2)))


# -----------------------------------------------------------------------------------------------------------

# Function to calculate parameters for double peak fit

def calculate_combined_fit(stopE,
                           best_vals_Opt,
                           best_vals_CT,
                           R2_Opt,
                           R2_CT,
                           eqe, T,
                           bias=False,
                           tolerance=0,
                           range=1.05,
                           include_disorder=False
                           ):
    """
    Function to compile combined fit for S1 and CT peak absorption after single peak fits
    :param stopE: stop energy of fit [float]
    :param best_vals_Opt: Opt fit values [list]
    :param best_vals_CT: CT fit values [list]
    :param R2_Opt: Opt fit R2 [float]
    :param R2_CT: CT fit R2 [float]
    :param eqe: EQE values [list]
    :param T: Temperature [float]
    :param bias: bias fit below data [boolean]
    :param tolerance: tolerance accepted of fit above data [float]
    :param range: defines upper bound of R2 calculation [float]
    :param include_disorder: boolean value to see whether to include disorder [bool]
    :return: result_dict : dictionary with fit results [dict]
                Dict keys:
                R2_Combined: R2 of sum of CT and Opt fit [float]
                R2_CT: R2 of CT fit [float]
                R2_Opt: R2 of Opt fit [float]
                R2_Average: average R2 of CT / Opt / Combined fit [float]
                Combined_Fit: sum of CT and Opt fit [list]
                Opt_fit: Opt fit values [list]
                CT_fit: CT fit values [list]
                Energy: Energy values [list]
                EQE: original EQE data [list]
    """

    wave_data, energy_data, eqe_data, log_eqe_data = compile_EQE(eqe,
                                                                 min(eqe['Energy']),
                                                                 stopE * range,  # Increase stop energy to expand fit!
                                                                 1)
    # # Optional code to add interpolation to the data
    # int_func = interp1d(eqe['Energy'], eqe['EQE'])
    # energy_data = np.arange(min(eqe['Energy']), stopE * range, 1)
    # eqe_data = int_func(energy_data)

    if R2_Opt != 0 and R2_CT != 0:
        Opt_fit = np.array([calculate_gaussian_absorption(e,
                                                          best_vals_Opt[0],
                                                          best_vals_Opt[1],
                                                          best_vals_Opt[2],
                                                          T)
                            for e in energy_data])
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

        if len(Opt_fit) == len(CT_fit):
            combined_Fit = Opt_fit + CT_fit
            combined_R_Squared = R_squared(eqe_data,
                                           combined_Fit.tolist(),
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

    # # Old code to return a list instead of dictionary
    # return [combined_R_Squared, combined_Fit, Opt_fit, CT_fit, energy_data, eqe_data]

# -----------------------------------------------------------------------------------------------------------
