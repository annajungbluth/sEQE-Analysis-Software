import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit import Model
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

from source.compilation import compile_EQE
from source.gaussian import calculate_combined_fit, calculate_gaussian_absorption, calculate_gaussian_disorder_absorption
from source.utils import R_squared
from source.utils import sep_list
from source.add_subtract import subtract_Opt
from source.plot import set_up_plot


# -----------------------------------------------------------------------------------------------------------

# Define parameters for initial guesses
# NOTE: Modify initial guesses if fit is unsuccessful
# NOTE: Some of these parameters are repeated/hard coded in sEQE_Analysis.py

# These values are used for single and simultaneous double fitting
f_guess = 0.001
l_guess = 0.150

# These values are used for simultaneous double fitting
fopt_guess = 0.01
lopt_guess = 0.150

# Set floating point precision
precision = 8  # decimal places

# -----------------------------------------------------------------------------------------------------------

# Function to perform curve fit

def fit_function(function,
                 energy_fit,
                 eqe_fit,
                 p0=None,
                 bounds=None,
                 include_disorder=False,
                 double=False
                 ):
    """
    Function to perform fitting using curve_fit
    This function is capable of performing single and double peak fits.
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
        best_vals, covar = curve_fit(function,
                                     energy_fit,
                                     eqe_fit,
                                     p0=p0,
                                     bounds=bounds
                                     )
    else:
        best_vals, covar = curve_fit(function,
                                     energy_fit,
                                     eqe_fit,
                                     p0=p0
                                     )
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

# Function to perform curve fit using lmfit.Model

def fit_model(function,
              energy_fit,
              eqe_fit,
              p0=None,
              include_disorder=False
              ):
    """
    Function to perform curve fit using lmfit
    This function is used for standard / disorder single peak fits.
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

    if p0 is None: # if no guess is given, initialize with all ones. This is consistent with curve_fit.
        if include_disorder:
            p0 = [1, 1, 1, 1]
        else:
            p0 = [1, 1, 1]

    gmodel = Model(function)

    gmodel.set_param_hint('Ect', min=0, max=1.6)
    gmodel.set_param_hint('l', min=0, max=0.6) # changed from 0.4
    gmodel.set_param_hint('f', min=0, max=0.1) # changed from 0.4

    if include_disorder:
        gmodel.set_param_hint('sig', min=0, max=0.2)

        result = gmodel.fit(eqe_fit,
                            f=p0[0],
                            l=p0[1],
                            Ect=p0[2],
                            sig=p0[3],
                            E=energy_fit
                            )

        f = float(result.params['f'].value)
        l = float(result.params['l'].value)
        Ect = float(result.params['Ect'].value)
        sig = float(result.params['sig'].value)

        best_vals = [f, l, Ect, sig]

        covar = result.covar
        if covar is None:
            covar = np.ones((4, 4))

        y_fit = gmodel.eval(E=np.array(energy_fit),
                            f=f,
                            l=l,
                            Ect=Ect,
                            sig=sig
                            )
    else:
        result = gmodel.fit(eqe_fit,
                            f=p0[0],
                            l=p0[1],
                            Ect=p0[2],
                            E=energy_fit
                            )

        f = float(result.params['f'].value)
        l = float(result.params['l'].value)
        Ect = float(result.params['Ect'].value)

        best_vals = [f, l, Ect]

        covar = result.covar
        if covar is None:
            covar = np.ones((4, 4))

        y_fit = gmodel.eval(E=np.array(energy_fit),
                            f=f,
                            l=l,
                            Ect=Ect
                            )

    r_squared = R_squared(eqe_fit, y_fit)

    return best_vals, covar, y_fit, r_squared


# -----------------------------------------------------------------------------------------------------------

# Function to perform curve fit using lmfit.Model

def fit_model_double(function,
                     energy_fit,
                     eqe_fit,
                     bound_dict,
                     p0=None,
                     include_disorder=False,
                     print_report=False
                     ):
    """
    Function to perform curve fit using lmfit
    This function is used for simultaneous double peak fitting
    :param function: function to fit against (i.e. gaussian, gaussian_disorder etc.)
    :param energy_fit: energy values to fit against [list or array]
    :param eqe_fit: EQE values to fit against [list or array]
    :param bound_dict: dictionary of boundary values [dict]
                       dict keys:
                       start_ECT, stop_ECT
                       start_lCT, stop_lCT
                       start_fCT, stop_fCT
                       start_Eopt, stop_Eopt
                       start_lopt, stop_lopt
                       start_fopt, stop_fopt
                       start_sig, stop_sig
    :param p0: list of initial guesses for curve_fit function [list]
    :param include_disorder: boolean value to specify whether to include disorder [bool]
    :param print_report: boolean value to specify whether to print fit report [bool]
    :return: best_vals: list of best fit parameters [list]
             covar: covariance matrix of fit
             y_fit: calculated EQE values of the fit [list]
             r_squared: R^2 of the fit [float]
    """
    if p0 is None: # if no guess is given, initialize with all ones. This is consistent with curve_fit.
        if include_disorder:
            p0 = [1, 1, 1, 1, 1, 1, 1]
        else:
            p0 = [1, 1, 1, 1, 1, 1]

    gmodel = Model(function)

    gmodel.set_param_hint('ECT', min=bound_dict['start_ECT'], max=bound_dict['stop_ECT'])
    gmodel.set_param_hint('lCT', min=bound_dict['start_lCT'], max=bound_dict['stop_lCT'])
    gmodel.set_param_hint('fCT', min=bound_dict['start_fCT'], max=bound_dict['stop_fCT'])
    gmodel.set_param_hint('Eopt', min=bound_dict['start_Eopt'], max=bound_dict['stop_Eopt'])
    gmodel.set_param_hint('lopt', min=bound_dict['start_lopt'], max=bound_dict['stop_lopt'])
    gmodel.set_param_hint('fopt', min=bound_dict['start_fopt'], max=bound_dict['stop_fopt'])

    if include_disorder:
        gmodel.set_param_hint('sig', min=bound_dict['start_sig'], max=bound_dict['stop_sig'])

        result = gmodel.fit(eqe_fit,
                            E=energy_fit,
                            fCT=p0[0],
                            lCT=p0[1],
                            ECT=p0[2],
                            fopt=p0[3],
                            lopt=p0[4],
                            Eopt=p0[5],
                            sig=p0[6]
                            )

        if print_report:
            print(result.fit_report())

        fCT = float(result.params['fCT'].value)
        lCT = float(result.params['lCT'].value)
        ECT = float(result.params['ECT'].value)
        fopt = float(result.params['fopt'].value)
        lopt = float(result.params['lopt'].value)
        Eopt = float(result.params['Eopt'].value)
        sig = float(result.params['sig'].value)

        best_vals = [fCT, lCT, ECT, fopt, lopt, Eopt, sig]

        covar = result.covar
        if covar is None:
            covar = np.ones((7, 7))

        y_fit = gmodel.eval(E=np.array(energy_fit),
                            fCT=fCT,
                            lCT=lCT,
                            ECT=ECT,
                            fopt=fopt,
                            lopt=lopt,
                            Eopt=Eopt,
                            sig=sig
                            )
    else:
        result = gmodel.fit(eqe_fit,
                            E=energy_fit,
                            fCT=p0[0],
                            lCT=p0[1],
                            ECT=p0[2],
                            fopt=p0[3],
                            lopt=p0[4],
                            Eopt=p0[5]
                            )

        if print_report:
            print(result.fit_report())

        fCT = float(result.params['fCT'].value)
        lCT = float(result.params['lCT'].value)
        ECT = float(result.params['ECT'].value)
        fopt = float(result.params['fopt'].value)
        lopt = float(result.params['lopt'].value)
        Eopt = float(result.params['Eopt'].value)

        best_vals = [fCT, lCT, ECT, fopt, lopt, Eopt]

        covar = result.covar
        if covar is None:
            covar = np.ones((6, 6))

        y_fit = gmodel.eval(E=np.array(energy_fit),
                            fCT=fCT,
                            lCT=lCT,
                            ECT=ECT,
                            fopt=fopt,
                            lopt=lopt,
                            Eopt=Eopt
                            )

    r_squared = R_squared(eqe_fit, y_fit)

    return best_vals, covar, y_fit, r_squared


# -----------------------------------------------------------------------------------------------------------

# Function to perform fit with guess range

def guess_fit(eqe,
              startE,
              stopE,
              function,
              guessRange,
              guessRange_opt=None,
              guessRange_sig=None,
              include_disorder=False,
              simultaneous_double=False,
              bounds=None
              ):
    """
    Function to loop through guesses and determine best fit using lmfit-based fit_model function
    This function is used for both standard / disorder single and simultaneous double peak fitting.
    :param eqe: EQE data [list]
    :param startE: fit start energy value [float]
    :param stopE: fit stop energy value [float]
    :param function: function to fit [function]
    :param guessRange: primary peak initial values [list].
                       For single peak fitting, this can be for the CT or Opt peak.
                       For simultaneous double peak fitting, this will be for the CT peak only.
    :param guessRange_opt: Opt peak initial values for simultaneous double fits [list]
    :param guessRange_sig: sigma initial values [list]
    :param include_disorder: boolean value specifying whether to include disorder [bool]
    :param simultaneous_double: boolean value specifying whether to perform a simultaneous double fit [bool]
    :param bounds: dictionary of boundary values for simultaneous double fitting [dict]
                   Escalates the use of lmfit.Model for single peak fitting.
    :return: best_vals: fit result [list]
             r_squared: R2 of fit [float]
    """

    if len(eqe) != 0:

        # Compile EQE to fit
        wave_fit, energy_fit, eqe_fit, log_eqe_fit = compile_EQE(eqe, startE, stopE, 1) # precision 8

        # Attempt peak fit:
        p0 = None

        p0_list = []

        # Compile guesses
        if simultaneous_double:  # Simultaneous double fitting
            if include_disorder:  # Simultaneous double fitting including disorder
                for CT_guess in guessRange:
                    for Opt_guess in guessRange_opt:
                        for sig_guess in guessRange_sig:
                            p0_list.append([f_guess,
                                            l_guess,
                                            round(CT_guess, 3),
                                            fopt_guess,
                                            lopt_guess,
                                            round(Opt_guess),
                                            round(sig_guess, 3)
                                            ])
            else: # Standard simultaneous double fitting
                for CT_guess in guessRange:
                    for Opt_guess in guessRange_opt:
                        p0_list.append([f_guess,
                                        l_guess,
                                        round(CT_guess, 3),
                                        fopt_guess,
                                        lopt_guess,
                                        round(Opt_guess, 3)
                                        ])
        else:  # Single peak fitting
            if include_disorder:  # Single peak fitting including disorder
                for E_guess in guessRange:
                    for sig_guess in guessRange_sig:
                        p0_list.append([f_guess,
                                        l_guess,
                                        round(E_guess, 3),
                                        round(sig_guess, 3)
                                        ])
            else:  # Standard single peak fitting
                for E_guess in guessRange:
                    p0_list.append([f_guess,
                                    l_guess,
                                    round(E_guess, 3)
                                    ])

        # Simultaneous double peak fitting
        if simultaneous_double:
            for p0 in p0_list: # Start with initial guesses, rather than p0 = None
                try:
                    best_vals, covar, y_fit, r_squared = fit_model_double(function=function,
                                                                          energy_fit=energy_fit,
                                                                          eqe_fit=eqe_fit,
                                                                          bound_dict=bounds,
                                                                          p0=p0,
                                                                          include_disorder=include_disorder
                                                                          )
                    if r_squared > 0:
                        return best_vals, covar, p0, r_squared
                    else:
                        raise ArithmeticError
                except:
                    if p0 == p0_list[-1]:
                        if include_disorder:
                            best_vals = [0, 0, 0, 0, 0, 0, 0]
                            covar = np.ones((7, 7))
                        else:
                            best_vals = [0, 0, 0, 0, 0, 0]
                            covar = np.ones((6, 6))
                        r_squared = 0
        # Single peak fitting
        else:
            for p0_guess in p0_list:
                try:
                    # TODO: Replace permanently with fit_model?
                    # NOTE: fit_function works better with single peak Marcus fitting
                    if bounds is None:
                        best_vals, covar, y_fit, r_squared = fit_function(function=function,
                                                                          energy_fit=energy_fit,
                                                                          eqe_fit=eqe_fit,
                                                                          p0=p0,
                                                                          include_disorder=include_disorder
                                                                          )
                    else:
                        best_vals, covar, y_fit, r_squared = fit_model(function=function,
                                                                       energy_fit=energy_fit,
                                                                       eqe_fit=eqe_fit,
                                                                       p0=p0,
                                                                       include_disorder=include_disorder
                                                                       )
                    if r_squared > 0:
                        return best_vals, covar, p0, r_squared
                    else:
                        raise ArithmeticError
                except Exception as e:
                    print(e)
                    p0 = p0_guess
                    if p0_guess == p0_list[-1]:
                        if include_disorder:
                            best_vals = [0, 0, 0, 0]
                            covar = np.ones((4, 4))
                        else:
                            best_vals = [0, 0, 0]
                            covar = np.ones((3, 3))
                        r_squared = 0

        return best_vals, covar, p0, r_squared

        # NOTE: Old code to loop through all initial guesses and determine best fit
        # if include_disorder:
        #     best_guess_df = pd.DataFrame()
        #     p0_list = []
        #     R2_list = []
        #     for CT_guess in guessRange_CT:
        #         for sig_guess in guessRange_sig:
        #             try:
        #                 best_vals, covar, y_fit, r_squared = fit_model(function,
        #                                                                energy_fit,
        #                                                                eqe_fit,
        #                                                                p0=p0,
        #                                                                include_disorder=True
        #                                                                )
        #                 if r_squared > 0:
        #                     p0_list.append(p0)
        #                     R2_list.append(r_squared)
        #                 else:
        #                     raise Exception('Wrong fit determined.')
        #                 p0 = [0.001, 0.1, round(CT_guess, 3), round(sig_guess, 3)]
        #             except:
        #                 p0 = [0.001, 0.1, round(CT_guess, 3), round(sig_guess, 3)]
        #             # except Exception as e:
        #             #     p0 = [0.001, 0.1, round(E_guess, 3), round(sig_guess, 3)]
        #             #     print(e)
        #
        #     best_guess_df['p0'] = p0_list
        #     best_guess_df['R2'] = R2_list
        #
        #     best_R2 = max(best_guess_df['R2'])
        #     best_p0 = best_guess_df['p0'][best_guess_df['R2'] == best_R2].values[0]  # Find best initial guess
        #
        #     # Determine fit values of fit with best intial guess
        #     best_vals, covar, y_fit, r_squared = fit_model(function,
        #                                                    energy_fit,
        #                                                    eqe_fit,
        #                                                    p0=best_p0,
        #                                                    include_disorder=True
        #                                                    )
        #
        # else:
        #     for CT_guess in guessRange_CT:
        #         try:
        #             best_vals, covar, y_fit, r_squared = fit_function(function,
        #                                                               energy_fit,
        #                                                               eqe_fit,
        #                                                               p0=p0
        #                                                               )
        #             if r_squared > 0:
        #                 return best_vals, r_squared
        #             else:
        #                 raise ArithmeticError
        #         except:
        #             p0 = [0.001, 0.1, CT_guess]
        #             if CT_guess == guessRange_CT[-1]:
        #                 best_vals = [0, 0, 0]
        #                 r_squared = 0
        #
        # return best_vals, r_squared


# -----------------------------------------------------------------------------------------------------------

# Mappable function to calculate guess fit

def calculate_guess_fit(x,
                        df,
                        eqe,
                        function,
                        guessRange,
                        guessRange_opt=None,
                        guessRange_sig=None,
                        include_disorder=False,
                        simultaneous_double=False,
                        bounds=None
                        ):
    """
    Mappable wrapper function to loop through initial guesses
    This function is used for standard / disorder single and simultaneous double peak fits.
    :param x: row of the dataFrame [int]
    :param df: results dataFrame
    :param eqe: EQE values [list]
    :param function: function to fit [function]
    :param guessRange: primary peak (either CT or Opt) initial values [list]
    :param guessRange_opt: Opt peak initial values [list]
    :param guessRange_sig: sigma initial values [list]
    :param bounds: dictionary of boundary values [dict]
    :param include_disorder: boolean value specifying whether to include disorder [bool]
    :param simultaneous_double: boolean value specifying whether to perform a simultaneous double fit [bool]
    :return: best_vals: fit result [list]
             r_squared: R2 of fit [float]
             df['Start'][x]: Start value of the fit [float]
             df['Stop'][x]: Stop value of the fit [float]
    """

    best_vals, covar, p0, r_squared = guess_fit(eqe=eqe,
                                                startE=df['Start'][x],
                                                stopE=df['Stop'][x],
                                                function=function,
                                                guessRange=guessRange,
                                                guessRange_opt=guessRange_opt,
                                                guessRange_sig=guessRange_sig,
                                                include_disorder=include_disorder,
                                                simultaneous_double=simultaneous_double,
                                                bounds=bounds
                                                )

    return [best_vals, covar, r_squared, df['Start'][x], df['Stop'][x]]


# -----------------------------------------------------------------------------------------------------------

# Mappable function to determine individual / combined fits

def map_fit(x,
            df_Opt,
            df_CT,
            eqe,
            guessRange_CT,
            function,
            T,
            bias=False,
            tolerance=0,
            sub_fit=1
            ):
    """
    Mappable function to determine individual / combined fits
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

    if sub_fit == 1:  # to subtract optical fit

        if df_Opt['R2'][x] > 0:  # Check that the fit was successful
            new_eqe = subtract_Opt(eqe=eqe, best_vals=df_Opt['Fit'][x], T=T)

            cal_vals_CT = list(
                map(lambda x: calculate_guess_fit(x,
                                                  df_CT,
                                                  new_eqe,
                                                  guessRange_CT,
                                                  function
                                                  ), range(len(df_CT))))

            best_vals_CT = list(map(lambda list_: sep_list(list_, 0), cal_vals_CT))
            covar_CT = list(map(lambda list_: sep_list(list_, 1), cal_vals_CT))
            R2_CT = list(map(lambda list_: sep_list(list_, 2), cal_vals_CT))

        else:
            best_vals_CT = [0, 0, 0]
            best_vals_CT = [best_vals_CT] * len(df_CT)

            R2_CT = [0]
            R2_CT = R2_CT * len(df_CT)

    elif sub_fit == 0:  # to not subtract optical fit

        cal_vals_CT = list(
            map(lambda x: calculate_guess_fit(x,
                                              df_CT,
                                              eqe,
                                              guessRange_CT,
                                              function
                                              ), range(len(df_CT))))

        best_vals_CT = list(map(lambda list_: sep_list(list_, 0), cal_vals_CT))
        covar_CT = list(map(lambda list_: sep_list(list_, 1), cal_vals_CT))
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
    parameter_list = list(map(lambda y: calculate_combined_fit(stopE=df_Opt['Stop'][x],
                                                               best_vals_Opt=df_Opt['Fit'][x],
                                                               best_vals_CT=best_vals_CT[y],
                                                               R2_Opt=df_Opt['R2'][x],
                                                               R2_CT=R2_CT[y],
                                                               eqe=eqe,
                                                               T=T,
                                                               bias=bias,
                                                               tolerance=tolerance
                                                               ),
                              range(len(df_CT))))

    combined_R2_list = list(map(lambda list_: sep_list(list_, 0), parameter_list))
    combined_Fit_list = list(map(lambda list_: sep_list(list_, 1), parameter_list))
    Opt_Fit_list = list(map(lambda list_: sep_list(list_, 2), parameter_list))
    CT_Fit_list = list(map(lambda list_: sep_list(list_, 3), parameter_list))
    Energy_list = list(map(lambda list_: sep_list(list_, 4), parameter_list))
    EQE_list = list(map(lambda list_: sep_list(list_, 5), parameter_list))

    return best_vals_Opt, R2_Opt, best_vals_CT, R2_CT, start_Opt_list, stop_Opt_list, start_CT_list, stop_CT_list, \
           combined_R2_list, combined_Fit_list, Opt_Fit_list, CT_Fit_list, Energy_list, EQE_list


# -----------------------------------------------------------------------------------------------------------

    # Function to determine the best fit for separate double peak fitting

def find_best_fit(df_both,
                  eqe,
                  T,
                  label,
                  n_fit=0,
                  include_disorder=False,
                  simultaneous_double=False,
                  ext_factor=1.2
                  ):
    """
    Function to find best EQE fits
    :param df_both: dataFrame with final fit results [dataFrame]
    :param eqe: original EQE data [dataFrame]
    :param T: Temperature [float]
    :param label: label to use in plot [string]
    :param n_fit: fit number [int]
    :param include_disorder: boolean of whether to include disorder [bool]
    :param simultaneous_double: boolean value to specify whether double peaks were fit simultaneously [bool]
    :param ext_factor: multiplication factor to extend data for R2 calculation [float]
    :return: df_copy: copy of df_both [dataFrame]
    """

    if len(df_both) != 0:

        # self.logger.info('Determining Best Fit ...')

        # Determine best fit
        max_index = df_both[df_both['Total_R2'] == max(df_both['Total_R2'])].index.values[0]

        if simultaneous_double: # Adjusts some of the print statements
            wave_plot, energy_plot, eqe_plot, log_eqe_plot = compile_EQE(eqe,
                                                                         min(eqe['Energy']),
                                                                         df_both['Stop'][max_index] * ext_factor,
                                                                         1)
            Opt_fit_plot = np.array([calculate_gaussian_absorption(e,
                                                                   df_both['Fit_Opt'][max_index][0],
                                                                   df_both['Fit_Opt'][max_index][1],
                                                                   df_both['Fit_Opt'][max_index][2],
                                                                   T)
                                     for e in energy_plot])

            print('-' * 35)
            print('R_Squared : ', format(df_both['Total_R2'][max_index], '.6f'))
            print('Fit Range (eV): ', df_both['Start'][max_index], ' - ', df_both['Stop'][max_index])
            print('-' * 35)
            print('f_Opt (eV**2) : ', format(df_both['Fit_Opt'][max_index][0], '.6f'),
                  '+/-', format(math.sqrt(df_both['Covar'][max_index][3, 3]), '.6f'))
            print('l_Opt (eV) : ', format(df_both['Fit_Opt'][max_index][1], '.6f'),
                  '+/-', format(math.sqrt(df_both['Covar'][max_index][4, 4]), '.6f'))
            print('E_Opt (eV) : ', format(df_both['Fit_Opt'][max_index][2], '.6f'),
                  '+/-', format(math.sqrt(df_both['Covar'][max_index][5, 5]), '.6f'))
            print('-' * 35)
            print('f_CT (eV**2) : ', format(df_both['Fit_CT'][max_index][0], '.6f'),
                  '+/-', format(math.sqrt(df_both['Covar'][max_index][0, 0]), '.6f'))
            print('l_CT (eV) : ', format(df_both['Fit_CT'][max_index][1], '.6f'),
                  '+/-', format(math.sqrt(df_both['Covar'][max_index][1, 1]), '.6f'))
            print('E_CT (eV) : ', format(df_both['Fit_CT'][max_index][2], '.6f'),
                  '+/-', format(math.sqrt(df_both['Covar'][max_index][2, 2]), '.6f'))

            if include_disorder:
                print('Sigma (eV) : ', format(df_both['Fit_CT'][max_index][3], '.6f'),
                      '+/-', format(math.sqrt(df_both['Covar'][max_index][6, 6]), '.6f'))

        else:
            wave_plot, energy_plot, eqe_plot, log_eqe_plot = compile_EQE(eqe, min(eqe['Energy']),
                                                                         df_both['Stop_Opt'][max_index] * ext_factor, 1)
            Opt_fit_plot = np.array([calculate_gaussian_absorption(e,
                                                                   df_both['Fit_Opt'][max_index][0],
                                                                   df_both['Fit_Opt'][max_index][1],
                                                                   df_both['Fit_Opt'][max_index][2],
                                                                   T)
                                     for e in energy_plot])

            # print('-' * 80)
            # print(('Combined Best Fit:').format(n_fit))
            # print('-' * 25)

            print('-' * 35)
            print('R_Squared : ', format(df_both['Total_R2'][max_index], '.6f'))
            print('-' * 35)
            print('Opt Fit Range (eV): ', df_both['Start_Opt'][max_index], ' - ', df_both['Stop_Opt'][max_index])
            print('f_Opt (eV**2) : ', format(df_both['Fit_Opt'][max_index][0], '.6f'),
                  '+/-', format(math.sqrt(df_both['Covar_Opt'][max_index][0, 0]), '.6f'))
            print('l_Opt (eV) : ', format(df_both['Fit_Opt'][max_index][1], '.6f'),
                  '+/-', format(math.sqrt(df_both['Covar_Opt'][max_index][1, 1]), '.6f'))
            print('E_Opt (eV) : ', format(df_both['Fit_Opt'][max_index][2], '.6f'),
                  '+/-', format(math.sqrt(df_both['Covar_Opt'][max_index][2, 2]), '.6f'))
            print('-' * 35)
            print('CT Fit Range (eV): ', df_both['Start_CT'][max_index], ' - ', df_both['Stop_CT'][max_index])
            print('f_CT (eV**2) : ', format(df_both['Fit_CT'][max_index][0], '.6f'),
                  '+/-', format(math.sqrt(df_both['Covar_CT'][max_index][0, 0]), '.6f'))
            print('l_CT (eV) : ', format(df_both['Fit_CT'][max_index][1], '.6f'),
                  '+/-', format(math.sqrt(df_both['Covar_CT'][max_index][1, 1]), '.6f'))
            print('E_CT (eV) : ', format(df_both['Fit_CT'][max_index][2], '.6f'),
                  '+/-', format(math.sqrt(df_both['Covar_CT'][max_index][2, 2]), '.6f'))

            if include_disorder:
                print('Sigma (eV) : ', format(df_both['Fit_CT'][max_index][3], '.6f'),
                      '+/-', format(math.sqrt(df_both['Covar_CT'][max_index][3, 3]), '.6f'))

            # print('Temperature [T] (K) : ', T)
            # print('-' * 80)

        axDouble_1, axDouble_2 = set_up_plot(flag='Energy')

        axDouble_1.plot(eqe['Energy'],
                        eqe['EQE'],
                        linewidth=2,
                        linestyle='-',
                        label=label,
                        color='black'
                        )
        axDouble_1.plot(energy_plot,
                        Opt_fit_plot,
                        linewidth=2,
                        linestyle='dotted',
                        label='Optical Peak Fit'
                        )
        axDouble_1.plot(df_both['Energy'][max_index],
                        df_both['CT_Fit'][max_index],
                        linewidth=2,
                        linestyle='--',
                        label='CT State Fit'
                        )
        axDouble_1.plot(df_both['Energy'][max_index],
                        df_both['Total_Fit'][max_index],
                        linewidth=2,
                        linestyle='dashdot',
                        label='Total Fit'
                        )
        axDouble_1.set_xlim(min(eqe['Energy']), 2.5)
        axDouble_1.set_title(('Range No. {}').format(n_fit))
        axDouble_1.legend()

        axDouble_2.plot(eqe['Energy'],
                        eqe['EQE'],
                        linewidth=2,
                        linestyle='-',
                        label=label,
                        color='black'
                        )
        axDouble_2.plot(energy_plot,
                        Opt_fit_plot,
                        linewidth=2,
                        linestyle='--',
                        label='Optical Peak Fit'
                        )
        axDouble_2.plot(df_both['Energy'][max_index],
                        df_both['CT_Fit'][max_index],
                        linewidth=2,
                        linestyle='--',
                        label='CT State Fit'
                        )
        axDouble_2.plot(df_both['Energy'][max_index],
                        df_both['Total_Fit'][max_index],
                        linewidth=2,
                        linestyle='dashdot',
                        label='Total Fit'
                        )
        axDouble_2.set_xlim(min(eqe['Energy']), 2.5)
        axDouble_2.set_ylim([10 ** (-7), max(eqe['EQE']) * 1.4])
        # axDouble_2.set_title(('Range No. {}').format(n_fit))
        axDouble_2.legend()

        df_copy = df_both.copy()
        df_copy['Total_R2'][max_index] = 0

    return df_copy


