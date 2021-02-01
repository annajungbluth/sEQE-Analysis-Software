import numpy as np
from scipy.optimize import curve_fit

from source.compilation import compile_EQE
from source.gaussian import calculate_gaussian_absorption
from source.utils import R_squared

# -----------------------------------------------------------------------------------------------------------

# Function to subtract optical fit from eqe

def subtract_Opt(eqe, best_vals, T):
    """
    :param eqe: EQE data [list]
    :param best_vals: Gaussian fit values [list]
    :param T: Temperature [float]
    :return: subtracted EQE [list]
    """

    eqe = eqe.copy()

    Opt_fit = np.array(
        [calculate_gaussian_absorption(e, best_vals[0], best_vals[1], best_vals[2], T) for e in
         eqe['Energy']])
    EQE_data = np.array(eqe['EQE'])

    subtracted_EQE = EQE_data - Opt_fit

    assert len(Opt_fit) == len(EQE_data)
    assert len(Opt_fit) == len(subtracted_EQE)

    eqe['EQE'] = subtracted_EQE

    return eqe

# -----------------------------------------------------------------------------------------------------------

# Function to perform fit with guess range

def guess_fit(eqe, startE, stopE, guessRange, function):
    """
    :param eqe: EQE data [list]
    :param startE: fit start energy value [float]
    :param stopE: fit stop energy value [float]
    :param guessRange: CT state energy initial values [list]
    :param function: function to fit [function]
    :return: best_vals: fit result [list]
             r_squared: R2 of fit [float]
    """

    if len(eqe) != 0:

        wave_fit, energy_fit, eqe_fit, log_eqe_fit = compile_EQE(eqe, startE, stopE, 1)

        # Attempt peak fit:
        p0 = None

        for E_guess in guessRange:
            try:
                # if include_Disorder:
                #     # Fit gaussian with disorder
                # else:
                #     # Fit gaussian without disorder

                best_vals, covar, y_fit, r_squared = fit_function(function, energy_fit, eqe_fit, p0=p0)
                if r_squared > 0:
                    return best_vals, r_squared
                else:
                    raise ArithmeticError
            except:
                p0 = [0.001, 0.1, E_guess]
                if E_guess == guessRange[-1]:
                    best_vals = [0, 0, 0]
                    r_squared = 0

                    return best_vals, r_squared

# -----------------------------------------------------------------------------------------------------------

# Mappable function to calculate guess fit

def calculate_guess_fit(x, df, eqe, guessRange, function):
    """
    :param x: row of the dataFrame [int]
    :param df: results dataFrame
    :param eqe: EQE values [list]
    :param guessRange: CT state energy initial values [list]
    :param function: function to fit [function]
    :return: best_vals: fit result [list]
             r_squared: R2 of fit [float]
             df['Start'][x]: Start value of the fit [float]
             df['Stop'][x]: Stop value of the fit [float]
    """

    best_vals, r_squared = guess_fit(eqe=eqe,
                                     startE=df['Start'][x],
                                     stopE=df['Stop'][x],
                                     guessRange=guessRange,
                                     function=function)

    return [best_vals, r_squared, df['Start'][x], df['Stop'][x]]

# -----------------------------------------------------------------------------------------------------------

# Wrapper function to perform curve fit

def fit_function(function, energy_fit, eqe_fit, p0=None):
    """
    :param function: function to fit against (i.e. gaussian, gaussian_disorder etc.)
    :param energy_fit: energy values to fit against [list or array]
    :param eqe_fit: EQE values to fit against [list or array]
    :param p0: list of initial guesses for curve_fit function [list]
    :return: best_vals: list of best fit parameters [list]
             covar: covariance matrix of fit
             y_fit: calculated EQE values of the fit [list]
             r_squared: R^2 of the fit [float]
    """
    best_vals, covar = curve_fit(function, energy_fit, eqe_fit, p0=p0)
    y_fit = [function(x, best_vals[0], best_vals[1], best_vals[2]) for x in energy_fit]
    r_squared = R_squared(eqe_fit, y_fit)

    return best_vals, covar, y_fit, r_squared

# -----------------------------------------------------------------------------------------------------------

# Mappable function to calculate combined fit

def calculate_combined_fit_df(x, df, eqe, T, bias = False, tolerance = 0):
    """
    :param x: row of the dataFrame [int]
    :param df: results dataFrame
    :param eqe: EQE values [list]
    :param T: Temperature [float]
    :param bias: bias fit below data [boolean]
    :param tolerance: tolerance accepted of fit above data [float]
    :return: list of :
             combined_R_Squared: R2 of sum of CT and Opt fit [float]
             combined_Fit: sum of CT and Opt fit [list]
             Opt_fit: Opt fit values [list]
             CT_fit: CT fit values [list]
             energy_data: Energy values [list]
             eqe_data: original EQE data [list]
    """

    wave_data, energy_data, eqe_data, log_eqe_data = compile_EQE(eqe, min(eqe['Energy']), df['Stop_Opt'][x] * 1.05, 1) # Increase the stop energy if you want to expand the fit!

    if df['R2_Opt'][x] != 0 and df['R2_CT'][x] != 0:

        Opt_fit = np.array([calculate_gaussian_absorption(e,
                                                          df['Fit_Opt'][x][0],
                                                          df['Fit_Opt'][x][1],
                                                          df['Fit_Opt'][x][2],
                                                          T)
                            for e in energy_data])
        CT_fit = np.array([calculate_gaussian_absorption(e,
                                                         df['Fit_CT'][x][0],
                                                         df['Fit_CT'][x][1],
                                                         df['Fit_CT'][x][2],
                                                         T)
                           for e in energy_data])

        if len(Opt_fit) == len(CT_fit):
                combined_Fit = Opt_fit + CT_fit
                combined_R_Squared = R_squared(eqe_data, combined_Fit.tolist(), bias=bias, tolerance=tolerance)

    else: # if any of the fits were unsuccessful
        Opt_fit = 0
        CT_fit = 0
        combined_Fit = 0
        combined_R_Squared = 0

    return [combined_R_Squared, combined_Fit, Opt_fit, CT_fit, energy_data, eqe_data]


def calculate_combined_fit(self, stopE, best_vals_Opt, best_vals_CT, R2_Opt, R2_CT, eqe):

    wave_data, energy_data, eqe_data, log_eqe_data = compile_EQE(eqe, min(eqe['Energy']), stopE * 1.05, 1) # Increase the stop energy if you want to expand the fit!

    if R2_Opt != 0 and R2_CT != 0:

        Opt_fit = np.array([calculate_gaussian_absorption(e,
                                                          best_vals_Opt[0],
                                                          best_vals_Opt[1],
                                                          best_vals_Opt[2],
                                                          self.T_double)
                            for e in energy_data])
        CT_fit = np.array([calculate_gaussian_absorption(e,
                                                         best_vals_CT[0],
                                                         best_vals_CT[1],
                                                         best_vals_CT[2],
                                                         self.T_double)
                           for e in energy_data])

        if len(Opt_fit) == len(CT_fit):
                combined_Fit = Opt_fit + CT_fit
                combined_R_Squared = R_squared(eqe_data, combined_Fit.tolist(), bias=self.bias, tolerance=self.tolerance)

    else: # if any of the fits were unsuccessful
        Opt_fit = 0
        CT_fit = 0
        combined_Fit = 0
        combined_R_Squared = 0

    return [combined_R_Squared, combined_Fit, Opt_fit, CT_fit, energy_data, eqe_data]