import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit import Model
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

from source.compilation import compile_EQE
from source.gaussian import calculate_gaussian_absorption, calculate_gaussian_disorder_absorption
from source.utils import R_squared
from source.utils import sep_list

# -----------------------------------------------------------------------------------------------------------

# Wrapper function to perform curve fit

def fit_function(function, energy_fit, eqe_fit, p0=None, bounds=None, include_disorder = False, double=False):
    """
    :param function: function to fit against (i.e. gaussian, gaussian_disorder etc.)
    :param energy_fit: energy values to fit against [list or array]
    :param eqe_fit: EQE values to fit against [list or array]
    :param p0: list of initial guesses for curve_fit function [list]
    :param bounds: tuple of bound values [tuple]
    :param include_disorder: boolean
    :param double: boolean
    :return: best_vals: list of best fit parameters [list]
             covar: covariance matrix of fit
             y_fit: calculated EQE values of the fit [list]
             r_squared: R^2 of the fit [float]
    """
    if bounds is not None:
        best_vals, covar = curve_fit(function, energy_fit, eqe_fit, p0=p0, bounds=bounds)
    else:
        best_vals, covar = curve_fit(function, energy_fit, eqe_fit, p0=p0)
    if double:
        if include_disorder:
            y_fit = [function(
                x,
                best_vals[0],
                best_vals[1],
                best_vals[2],
                best_vals[3],
                best_vals[4],
                best_vals[5],
                best_vals[6]
            ) for x in energy_fit]
        else:
            y_fit = [function(
                x,
                best_vals[0],
                best_vals[1],
                best_vals[2],
                best_vals[3],
                best_vals[4],
                best_vals[5]
            ) for x in energy_fit]
    else:
        if include_disorder:
            y_fit = [function(
                x,
                best_vals[0],
                best_vals[1],
                best_vals[2],
                best_vals[3]
            ) for x in energy_fit]
        else:
            y_fit = [function(
                x,
                best_vals[0],
                best_vals[1],
                best_vals[2]
            ) for x in energy_fit]
    r_squared = R_squared(eqe_fit, y_fit)

    return best_vals, covar, y_fit, r_squared

# -----------------------------------------------------------------------------------------------------------

# Function to perform fit with guess range

def guess_fit(eqe, startE, stopE, guessRange, function, guessRange_sig=None, include_disorder=False):
    """
    :param eqe: EQE data [list]
    :param startE: fit start energy value [float]
    :param stopE: fit stop energy value [float]
    :param guessRange: CT state energy initial values [list]
    :param function: function to fit [function]
    :param guessRange_sig: sigma initial values [list]
    :param include_disorder: boolean value to see whether to include disorder [bool]
    :return: best_vals: fit result [list]
             r_squared: R2 of fit [float]
    """

    if len(eqe) != 0:

        wave_fit, energy_fit, eqe_fit, log_eqe_fit = compile_EQE(eqe, startE, stopE, 1)

        # Attempt peak fit:
        p0 = None

        if include_disorder:
            best_guess_df = pd.DataFrame()
            p0_list = []
            R2_list = []
            for E_guess in guessRange:
                for sig_guess in guessRange_sig:
                    try:
                        best_vals, covar, y_fit, r_squared = fit_model(function, energy_fit, eqe_fit, p0=p0, include_disorder=True)
                        if r_squared > 0:
                            p0_list.append(p0)
                            R2_list.append(r_squared)
                        else:
                            raise Exception('Wrong fit determined.')
                        p0 = [0.001, 0.1, round(E_guess, 3), round(sig_guess, 3)]
                    except:
                        p0 = [0.001, 0.1, round(E_guess, 3), round(sig_guess, 3)]
                    # except Exception as e:
                    #     p0 = [0.001, 0.1, round(E_guess, 3), round(sig_guess, 3)]
                    #     print(e)

            best_guess_df['p0'] = p0_list
            best_guess_df['R2'] = R2_list

            best_R2 = max(best_guess_df['R2'])
            best_p0 = best_guess_df['p0'][best_guess_df['R2'] == best_R2].values[0]  # Find best initial guess

            # Determine fit values of fit with best intial guess
            best_vals, covar, y_fit, r_squared = fit_model(function, energy_fit, eqe_fit, p0=best_p0, include_disorder=True)

        else:
            for E_guess in guessRange:
                try:
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

# Function to perform simultaneous double peak fit with guess range

def guess_fit_sim(eqe, startE, stopE, function, guessRange_CT, guessRange_Opt, guessRange_sig=None, include_disorder=False, **kwargs):
    """
    :param eqe: EQE data [list]
    :param startE: fit start energy value [float]
    :param stopE: fit stop energy value [float]
    :param function: function to fit [function]
    :param guessRange_CT: CT state energy initial values [list]
    :param guessRange_Opt: Opt state energy initial values [list]
    :param guessRange_sig: sigma initial values [list]
    :param include_disorder: boolean value to see whether to include disorder [bool]
    :return: best_vals: fit result [list]
             r_squared: R2 of fit [float]
    """

    if len(eqe) != 0:
        int_func = interp1d(eqe['Energy'], eqe['EQE'])

        wave_fit, energy_fit_2, eqe_fit_2, log_eqe_fit = compile_EQE(eqe, startE, stopE, 1)

        energy_fit = np.arange(startE, stopE, 0.001)
        eqe_fit = int_func(energy_fit)

        # Attempt peak fit:
        p0 = None
        if include_disorder:
        #     bounds = ([0.00001, 0.01, 1.2, 0.00001, 0.01, 1.4, 0.01],  # fCT, lCT, ECT, fopt, lopt, Eopt, sig
        #               [0.04, 0.4, 1.5, 0.1, 0.05, 1.75, 0.5])
        # else:
        #     bounds = ([0.00001, 0.01, 1.2, 0.00001, 0.01, 1.4],  # fCT, lCT, ECT, fopt, lopt, Eopt
        #               [0.04, 0.4, 1.5, 0.1, 0.05, 1.75])
            bounds = ([kwargs['start_fCT'], kwargs['start_lCT'], kwargs['start_ECT'], kwargs['start_fopt'], kwargs['start_lopt'], kwargs['start_Eopt'], kwargs['start_sig']],
                      [kwargs['stop_fCT'], kwargs['stop_lCT'], kwargs['stop_ECT'], kwargs['stop_fopt'], kwargs['stop_lopt'], kwargs['stop_Eopt'], kwargs['stop_sig']])
        else:
            bounds = ([kwargs['start_fCT'], kwargs['start_lCT'], kwargs['start_ECT'], kwargs['start_fopt'], kwargs['start_lopt'], kwargs['start_Eopt']],
                      [kwargs['stop_fCT'], kwargs['stop_lCT'], kwargs['stop_ECT'], kwargs['stop_fopt'], kwargs['stop_lopt'], kwargs['stop_Eopt']])

        p0_list = []
        R2_list = []
        best_vals_list = []
        best_guess_df = pd.DataFrame()

        if include_disorder:
            for CT_guess in guessRange_CT:
                for Opt_guess in guessRange_Opt:
                    for sig_guess in guessRange_sig:
                        try:
                            best_vals, covar, y_fit, r_squared = fit_function(
                                function,
                                energy_fit,
                                eqe_fit,
                                p0=p0,
                                bounds=bounds,
                                include_disorder=include_disorder,
                                double=True
                            )
                            if r_squared > 0:
                                # return best_vals, r_squared # break loops if fit is ok
                                p0_list.append(p0)
                                R2_list.append(r_squared)
                                best_vals_list.append(best_vals)
                            else:
                                raise Exception('Wrong fit determined.')
                            p0 = [0.001, 0.1, round(CT_guess, 3), 0.01, 0.15, round(Opt_guess, 3), round(sig_guess, 3)]
                        except:
                            p0 = [0.001, 0.1, round(CT_guess, 3), 0.01, 0.15, round(Opt_guess, 3), round(sig_guess, 3)]
                            if CT_guess == guessRange_CT[-1]:
                                best_vals = [0, 0, 0, 0, 0, 0, 0]
                                r_squared = 0
                                p0_list.append(p0)
                                R2_list.append(r_squared)
                                best_vals_list.append(best_vals)

            best_guess_df['p0'] = p0_list
            best_guess_df['R2'] = R2_list
            best_guess_df['best_vals'] = best_vals_list

            best_R2 = max(best_guess_df['R2'])
            best_p0 = best_guess_df['p0'][best_guess_df['R2'] == best_R2].values[0]  # Find best initial guess

            # Determine fit values of fit with best intial guess
            best_vals, covar, y_fit, r_squared = fit_function(function, energy_fit, eqe_fit, p0=best_p0, bounds=bounds, include_disorder=True, double=True)

        else:
            for CT_guess in guessRange_CT:
                for Opt_guess in guessRange_Opt:
                    try:
                        best_vals, covar, y_fit, r_squared = fit_function(
                            function,
                            energy_fit,
                            eqe_fit,
                            p0=p0,
                            bounds=bounds,
                            include_disorder=include_disorder,
                            double=True
                        )
                        if r_squared > 0:
                            return best_vals, r_squared


                            p0_list.append(p0)
                            R2_list.append(r_squared)
                            best_vals_list.append(best_vals)
                        else:
                            raise Exception('Wrong fit determined.')
                        p0 = [0.001, 0.1, round(CT_guess, 3), 0.01, 0.15, round(Opt_guess, 3)]
                    except:
                        p0 = [0.001, 0.1, round(CT_guess, 3), 0.01, 0.15, round(Opt_guess, 3)]
                        if CT_guess == guessRange_CT[-1]:
                            best_vals = [0, 0, 0, 0, 0, 0, 0]
                            r_squared = 0
                            return best_vals, r_squared


                            p0_list.append(p0)
                            R2_list.append(r_squared)
                            best_vals_list.append(best_vals)

            best_guess_df['p0'] = p0_list
            best_guess_df['R2'] = R2_list
            best_guess_df['best_vals'] = best_vals_list

            best_guess_df.to_csv('~/Desktop/Test.csv')
            print('saved')

            best_R2 = max(best_guess_df['R2'])
            best_p0 = best_guess_df['p0'][best_guess_df['R2'] == best_R2].values[0]  # Find best initial guess

            # Determine fit values of fit with best intial guess
            best_vals, covar, y_fit, r_squared = fit_function(function, energy_fit, eqe_fit, p0=best_p0, bounds=bounds, include_disorder=False, double=True)

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

# Function to calculate parameters for double peak fit

def calculate_combined_fit(stopE, best_vals_Opt, best_vals_CT, R2_Opt, R2_CT, eqe, T, bias = False, tolerance = 0, range = 1.05, include_disorder=False):
    """
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
    :return: list of :
             combined_R_Squared: R2 of sum of CT and Opt fit [float]
             combined_Fit: sum of CT and Opt fit [list]
             Opt_fit: Opt fit values [list]
             CT_fit: CT fit values [list]
             energy_data: Energy values [list]
             eqe_data: original EQE data [list]
    """

    wave_data, energy_data, eqe_data, log_eqe_data = compile_EQE(eqe, min(eqe['Energy']), stopE * range, 1) # (+ 0.15) Increase the stop energy if you want to expand the fit!

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
                combined_R_Squared = R_squared(eqe_data, combined_Fit.tolist(), bias=bias, tolerance=tolerance)

    else: # if any of the fits were unsuccessful
        Opt_fit = 0
        CT_fit = 0
        combined_Fit = 0
        combined_R_Squared = 0

    return [combined_R_Squared, combined_Fit, Opt_fit, CT_fit, energy_data, eqe_data]

# -----------------------------------------------------------------------------------------------------------

# Function to calculate parameters for simultaneous double peak fit

def calculate_combined_fit_sim(stopE, best_vals_Opt, best_vals_CT, eqe, T, bias = False, tolerance = 0, range = 1.05, include_disorder=False):
    """
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
    :return: list of :
             combined_R_Squared: R2 of sum of CT and Opt fit [float]
             combined_Fit: sum of CT and Opt fit [list]
             Opt_fit: Opt fit values [list]
             CT_fit: CT fit values [list]
             energy_data: Energy values [list]
             eqe_data: original EQE data [list]
    """

    # wave_data, energy_data, eqe_data, log_eqe_data = compile_EQE(eqe, min(eqe['Energy']), stopE * range, 1) # (+ 0.15) Increase the stop energy if you want to expand the fit!

    int_func = interp1d(eqe['Energy'], eqe['EQE'])

    energy_data = np.arange(1.31, 2.00, 0.001)
    eqe_data = int_func(energy_data)



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

    R2_Opt = R_squared(eqe_data, Opt_fit, bias=bias, tolerance=tolerance)
    R2_CT = R_squared(eqe_data, CT_fit, bias=bias, tolerance=tolerance)

    if len(Opt_fit) == len(CT_fit):
            combined_Fit = Opt_fit + CT_fit
            R2_sum = R_squared(eqe_data, combined_Fit.tolist(), bias=bias, tolerance=tolerance)

    R2_comp = (R2_CT + R2_Opt + R2_sum)/3

    return [R2_sum, combined_Fit, Opt_fit, CT_fit, energy_data, eqe_data, R2_Opt, R2_CT, R2_comp]

# -----------------------------------------------------------------------------------------------------------

# Mappable function to determine individual / combined fits

def map_fit(x, df_Opt, df_CT, eqe, guessRange_CT, function, T, bias=False, tolerance=0, sub_fit = 1):
    """
    :param x: row of the dataFrame [int]
    :param df_Opt: dataFrame of Opt fit parameters [dataFrame]
    :param df_CT: dataFrame of CT fit parameters [dataFrame]
    :param eqe: EQE values [list]
    :param guessRange_CT: CT state energy initial values [list]
    :param function: function to fit [function]
    :param T: Temperature [float]
    :param bias: bias fit below data [boolean]
    :param tolerance: tolerance accepted of fit above data [float]
    :param sub_fit: indicates whether Opt fit is subtracted [0 - no subtract, 1 - subtract]
    :return: best_vals_Opt: fit values of Opt fit [list]
             R2_Opt: R2 of Opt fit [float]
             best_vals_CT: fit values of CT fit [list]
             R2_CT: R2 of CT fit [float]
             start_Opt_list: start energies of Opt fit [list]
             stop_Opt_list: stop energies of Opt fit [list]
             start_CT_list: start energies of CT fit [list]
             stop_CT_list: stop energies of CT fit [list]
             combined_R2_list: R2 of sum of Opt and CT fit [float]
             combined_Fit_list: sum of Opt and CT fit [list]
             Opt_Fit_list: Opt fit values [list]
             CT_Fit_list: CT fit values [list]
             Energy_list: Energy values [list]
             EQE_list: Original EQE data [list]
    """

    if sub_fit == 1:

        if df_Opt['R2'][x] > 0:  # Check that the fit was successful
            new_eqe = subtract_Opt(eqe=eqe, best_vals=df_Opt['Fit'][x], T=T)

            cal_vals_CT = list(map(lambda x: calculate_guess_fit(x, df_CT, new_eqe, guessRange_CT, function), range(len(df_CT))))

            best_vals_CT = list(map(lambda list_: sep_list(list_, 0), cal_vals_CT))
            R2_CT = list(map(lambda list_: sep_list(list_, 1), cal_vals_CT))

        else:
            best_vals_CT = [0, 0, 0]
            best_vals_CT = [best_vals_CT] * len(df_CT)

            R2_CT = [0]
            R2_CT = R2_CT * len(df_CT)

    elif sub_fit == 0 :

        cal_vals_CT = list(map(lambda x: calculate_guess_fit(x, df_CT, eqe, guessRange_CT, function), range(len(df_CT))))

        best_vals_CT = list(map(lambda list_: sep_list(list_, 0), cal_vals_CT))
        R2_CT = list(map(lambda list_: sep_list(list_, 1), cal_vals_CT))

    start_Opt_list = [df_Opt['Start'][x]]
    start_Opt_list = start_Opt_list * len(df_CT)

    stop_Opt_list = [df_Opt['Stop'][x]]
    stop_Opt_list = stop_Opt_list * len(df_CT)

    best_vals_Opt = df_Opt['Fit'][x]
    best_vals_Opt = [best_vals_Opt] * len(df_CT)

    R2_Opt = [df_Opt['R2'][x]]
    R2_Opt = R2_Opt * len(df_CT)

    start_CT_list = list(df_CT['Start'])
    stop_CT_list = list(df_CT['Stop'])

    # Calculate combined fit here
    parameter_list = list(map(lambda y: calculate_combined_fit(stopE = df_Opt['Stop'][x],
                                                               best_vals_Opt = df_Opt['Fit'][x],
                                                               best_vals_CT = best_vals_CT[y],
                                                               R2_Opt = df_Opt['R2'][x],
                                                               R2_CT = R2_CT[y],
                                                               eqe = eqe,
                                                               T = T,
                                                               bias = bias,
                                                               tolerance = tolerance),
                              range(len(df_CT))))

    combined_R2_list = list(map(lambda list_: sep_list(list_, 0), parameter_list))
    combined_Fit_list = list(map(lambda list_: sep_list(list_, 1), parameter_list))
    Opt_Fit_list = list(map(lambda list_: sep_list(list_, 2), parameter_list))
    CT_Fit_list = list(map(lambda list_: sep_list(list_, 3), parameter_list))
    Energy_list = list(map(lambda list_: sep_list(list_, 4), parameter_list))
    EQE_list = list(map(lambda list_: sep_list(list_, 5), parameter_list))

    return best_vals_Opt, R2_Opt, best_vals_CT, R2_CT, start_Opt_list, stop_Opt_list, start_CT_list, stop_CT_list, combined_R2_list, combined_Fit_list, Opt_Fit_list, CT_Fit_list, Energy_list, EQE_list

# -----------------------------------------------------------------------------------------------------------

# Wrapper function to perform curve fit using lmfit.Model ### This could be changed to use scipy.curve_fit instead

def fit_model(function, energy_fit, eqe_fit, p0=None, include_disorder=False):
    """
    :param function: function to fit against (i.e. gaussian, gaussian_disorder etc.)
    :param energy_fit: energy values to fit against [list or array]
    :param eqe_fit: EQE values to fit against [list or array]
    :param p0: list of initial guesses for curve_fit function [list]
    :param include_disorder: boolean value
    :return: best_vals: list of best fit parameters [list]
             covar: covariance matrix of fit
             y_fit: calculated EQE values of the fit [list]
             r_squared: R^2 of the fit [float]
    """
    gmodel = Model(function)

    gmodel.set_param_hint('Ect', min=0.8, max=1.6)
    gmodel.set_param_hint('l', min=0.01, max=0.3)
    gmodel.set_param_hint('f', min=0.001, max=0.4)

    if include_disorder:
        gmodel.set_param_hint('sig', min=0.002, max=0.2)

        result = gmodel.fit(eqe_fit, f=p0[0], l=p0[1], Ect=p0[2], sig=p0[3], E=energy_fit)

        f = float(result.params['f'].value)
        l = float(result.params['l'].value)
        Ect = float(result.params['Ect'].value)
        sig = float(result.params['sig'].value)

        best_vals = [f, l, Ect, sig]

        covar = result.covar
        if covar is None:
            covar = np.zeros((4,4))

        y_fit = gmodel.eval(E=energy_fit, f=f, l=l, Ect=Ect, sig=sig)
    else:
        result = gmodel.fit(eqe_fit, f=p0[0], l=p0[1], Ect=p0[2], E=energy_fit)

        f = float(result.params['f'].value)
        l = float(result.params['l'].value)
        Ect = float(result.params['Ect'].value)

        best_vals = [f, l, Ect]

        covar = result.covar
        if covar is None:
            covar = np.zeros((4,4))

        y_fit = gmodel.eval(E=energy_fit, f=f, l=l, Ect=Ect)

    r_squared = R_squared(eqe_fit, y_fit)

    return best_vals, covar, y_fit, r_squared

# -----------------------------------------------------------------------------------------------------------
