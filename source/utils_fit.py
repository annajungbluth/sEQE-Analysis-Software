import math
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit import Model
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from tkinter import filedialog

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

k = 8.617 * math.pow(10, -5)  # [ev/K]

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
    """Function to perform curve fit

    Parameters
    ----------
    function : function, required
        Function to perform fit with (i.e. gaussian, gaussian_disorder etc.)
    energy_fit : list or array, required
        Energy values to fit against
    eqe_fit : list or array, required
        EQE values to fit again
    p0 : list, optional
        List of initial guesses for curve_fit function
    bounds : tuple, optional
        Tuple of fit bound values
    include_disorder : bool, optional
        Boolean value specifying whether to include CT state disorder
    double : bool, optional
        Boolean value specifying whether to perform double peak fit

    Returns
    -------
    best_vals : list
        List of best fit parameters
    covar : array
        Covariance matrix of fit
    y_fit : list
        Calculated EQE values of fit
    r_squared : float
        R squared of fit
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
    """Function to perform curve fit using lmfit

    Parameters
    ----------
    function : function, required
        Function to perform fit with (i.e. gaussian, gaussian_disorder etc.)
    energy_fit : list or array, required
        Energy values to fit against
    eqe_fit : list or array, required
        EQE values to fit again
    p0 : list, optional
        List of initial guesses for curve_fit function
    include_disorder : bool, optional
        Boolean value specifying whether to include CT state disorder

    Returns
    -------
    best_vals : list
        List of best fit parameters
    covar : array
        Covariance matrix of fit
    y_fit : list
        Calculated EQE values of fit
    r_squared : float
        R squared of fit
    """

    if p0 is None: # if no guess is given, initialize with ones. This is consistent with curve_fit.
        if include_disorder:
            p0 = [1, 1, 1, 1]
        else:
            p0 = [1, 1, 1]

    gmodel = Model(function)

    gmodel.set_param_hint('Ect', min=0, max=1.6)
    # gmodel.set_param_hint('Ect', min=0, max=1.7)
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
            covar = np.zeros((4, 4))

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
            covar = np.zeros((4, 4))

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
    """Function to perform double peak curve fit using lmfit.Model

    Parameters
    ----------
    function : function, required
        Function to perform fit with (i.e. gaussian, gaussian_disorder etc.)
    energy_fit : list or array, required
        Energy values to fit against
    eqe_fit : list or array, required
        EQE values to fit again
    bound_dict : dict
        Dictionary of boundary values
        Dict keys:
            start_ECT, stop_ECT
            start_lCT, stop_lCT
            start_fCT, stop_fCT
            start_Eopt, stop_Eopt
            start_lopt, stop_lopt
            start_fopt, stop_fopt
            start_sig, stop_sig
    p0 : list, optional
        List of initial guesses for curve_fit function
    include_disorder : bool, optional
        Boolean value specifying whether to include CT state disorder
    print_report : bool, optional
        Boolean value specifying whether to print fit report
    
    Returns
    -------
    best_vals : list
        List of best fit parameters
    covar : array
        Covariance matrix of fit
    y_fit : list
        Calculated EQE values of fit
    r_squared : float
        R squared of fit
    """

    if p0 is None: # if no guess is given, initialize with ones. This is consistent with curve_fit.
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
            covar = np.zeros((7, 7))

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
            covar = np.zeros((6, 6))

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
    """Function to loop through guesses and determine best fit using lmfit-based fit_model function
    This function is used for both standard / disorder single and simultaneous double peak fitting.

    Parameters
    ----------
    eqe : list, required
        Input EQE data
    startE : float, required
        Fit start energy value [eV]
    stopE : float, required
        Fit stop energy value [eV]
    function : function, required
        Function to fit
    guessRange : list, required
        Primary peak initial values
        For single peak fitting, this can be for the CT or Opt peak
        For simultaneous double peak fitting, this will be for the CT peak only
    guessRange_opt : list, optional
        Optical (S1) peak initial values for simultaneous double peak fits
    guessRange_sig : list, optional 
        Disorder parameter initial values 
    include_disorder : bool, optional
        Boolean value specifying whether to include peak disorder
    simultaneous_double : bool, optional
        Boolean value specifying whether to perform simultaneous double peak fitting
    bounds : dict, optional
        Dictionary of boundary values for simultaneous double fitting
        Escalates the use of lmfit.Model for single peak fitting
    
    Returns
    -------
    best_vals : list
        List of best fit parameters
    r_squared : float
        R squared of fit
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
                            covar = np.zeros((7, 7))
                        else:
                            best_vals = [0, 0, 0, 0, 0, 0]
                            covar = np.zeros((6, 6))
                        r_squared = 0
        # Single peak fitting
        else:
            for p0_guess in p0_list:
                try:
                    # NOTE: Replace permanently with fit_model?
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
                            covar = np.zeros((4, 4))
                        else:
                            best_vals = [0, 0, 0]
                            covar = np.zeros((3, 3))
                        r_squared = 0

        return best_vals, covar, p0, r_squared

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
    """Mappable wrapper function to loop through initial guesses
    This function is used for standard / disorder single and simultaneous double peak fits.

    Parameters
    ----------
    x : int, required
        Number specifying row of dataFrame
    df : dataFrame, required
        DataFrame of fit results
    eqe : list, required
        List of EQE values
    function : function, required
        Function to fit
    guessRange : list, required
        Primary peak (either CT or Optical (S1) peak) initial values
    guessRange_opt : list, optional
        Optical (S1) peak initial values
    guessRange_sig : list, optional 
        Disorder parameter initial values 
    include_disorder : bool, optional
        Boolean value specifying whether to include peak disorder
    simultaneous_double : bool, optional
        Boolean value specifying whether to perform simultaneous double peak fitting
    
    Returns
    -------
    best_vals : list
        List of best fit parameters
    r_squared : float
        R squared of fit
    df['Start'][x] : float, required
        Start value of the fit
    df['Stop'][x] : float, 
        Stop value of the fit
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
    """Mappable function to determine individual / combined fits

    Parameters
    ----------
    x : int, required
        Number specifying row of dataFrame
    df_Opt : dataFrame, required
        DataFrame of Optical (S1) fit results
    df_CT : dataFrame, required
        DataFrame of CT state fit results
    eqe : list, required
        List of EQE values
    guessRange_CT : list, optional
        CT state initial values
    function : function, required
        Function to fit
    T : float, required
        Temperature [K]
    bias : bool, optional
        Boolean value specifying whether to bias fits above the data
    tolerance : float, optional
        Tolerance for fit above the data
    sub_fit : int, optional
        Number indicating whether Optical (S1) peak fit is substracted
        0 - Don't subtract
        1 - Subtract

    Returns
    -------
    best_vals_Opt : list
        List of Optical (S1) peak fit values
    R2_Opt : list
        R squared of Optical (S1) peak fit
    best_vals_CT : list
        List of CT state fit values
    R2_CT : list
        R squared of CT state fit
    start_Opt_list : list
        Start energies of Optical (S1) peak fit
    stop_Opt_list : list
        Stop energies of Optical (S1) peak fit
    start_CT_list : list
        Start energies of CT state fit
    stop_CT_list : list
        Stop energies of CT state fit
    combined_R2_list: float
        R squared of sum of Optical (S1) and CT state fit
    combined_Fit_list: float
        Sum of Optical (S1) and CT state fit
    Opt_Fit_list: list
        Optical (S1) fit values
    CT_Fit_list: list
        CT state fit values
    Energy_list: list
        Energy values
    EQE_list: list
        Original EQE data

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
                  ext_factor=1.2,
                  save_fit=False,
                  save_fit_file=None
                  ):
    """Function to determine the best fit for separate double peak fitting

    Parameters
    ----------
    df_both : dataFrame, required
        DataFrame with results of separate double peak fitting
    eqe : dataFrame, required
        Input EQE data
    T : float, required
        Temperature [K]
    label : str, optional
        Plot label
    n_fit : int, optional
        Fit number 
    include_disorder : bool, optional
        Boolean value specifying whether to include peak disorder
    simultaneous_double : bool, optional
        Boolean value specifying whether simultaneous double peak fitting was performed
    ext_factor : float, optional
        Multiplication factor to extend data for R squared calculation
    save_fit : bool, optional
        Boolean value specifying whether to save fits
    save_fit_file : str, optional
        Directory or folder name to save data to

    Returns
    -------
    df_copy : dataFrame
        Copy of dataFrame summarizing Optical (S1) and CT state fits

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

            print('-' * 80)
            print('R2 : ', format(df_both['Total_R2'][max_index], '.6f'))
            print('Fit Range (eV): ', df_both['Start'][max_index], ' - ', df_both['Stop'][max_index])
            print('-' * 35)
            print('f_Opt (eV^2) : ', format(df_both['Fit_Opt'][max_index][0], '.6f'),
                  '+/-', format(math.sqrt(abs(df_both['Covar'][max_index][3, 3])), '.6f'))
            print('l_Opt (eV) : ', format(df_both['Fit_Opt'][max_index][1], '.6f'),
                  '+/-', format(math.sqrt(abs(df_both['Covar'][max_index][4, 4])), '.6f'))
            print('E_Opt (eV) : ', format(df_both['Fit_Opt'][max_index][2], '.6f'),
                  '+/-', format(math.sqrt(abs(df_both['Covar'][max_index][5, 5])), '.6f'))
            print('-' * 35)
            print('f_CT (eV^2) : ', format(df_both['Fit_CT'][max_index][0], '.6f'),
                  '+/-', format(math.sqrt(abs(df_both['Covar'][max_index][0, 0])), '.6f'))
            print('l_CT (eV) : ', format(df_both['Fit_CT'][max_index][1], '.6f'),
                  '+/-', format(math.sqrt(abs(df_both['Covar'][max_index][1, 1])), '.6f'))
            print('E_CT (eV) : ', format(df_both['Fit_CT'][max_index][2], '.6f'),
                  '+/-', format(math.sqrt(abs(df_both['Covar'][max_index][2, 2])), '.6f'))

            if include_disorder:
                print('Sigma (eV) : ', format(df_both['Fit_CT'][max_index][3], '.6f'),
                      '+/-', format(math.sqrt(abs(df_both['Covar'][max_index][6, 6])), '.6f'))
                W = df_both['Fit_CT'][max_index][1] * T + (df_both['Fit_CT'][max_index][3] ** 2) / (2 * k)
                print('Gaussian Variance [W] (eV K) : ', format(W, '.2f'))

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

            print('-' * 80)
            print('R2 : ', format(df_both['Total_R2'][max_index], '.6f'))
            print('-' * 35)
            print('Opt Fit Range (eV): ', df_both['Start_Opt'][max_index], ' - ', df_both['Stop_Opt'][max_index])
            print('f_Opt (eV^2) : ', format(df_both['Fit_Opt'][max_index][0], '.6f'),
                  '+/-', format(math.sqrt(abs(df_both['Covar_Opt'][max_index][0, 0])), '.6f'))
            print('l_Opt (eV) : ', format(df_both['Fit_Opt'][max_index][1], '.6f'),
                  '+/-', format(math.sqrt(abs(df_both['Covar_Opt'][max_index][1, 1])), '.6f'))
            print('E_Opt (eV) : ', format(df_both['Fit_Opt'][max_index][2], '.6f'),
                  '+/-', format(math.sqrt(abs(df_both['Covar_Opt'][max_index][2, 2])), '.6f'))
            print('-' * 35)
            print('CT Fit Range (eV): ', df_both['Start_CT'][max_index], ' - ', df_both['Stop_CT'][max_index])
            print('f_CT (eV^2) : ', format(df_both['Fit_CT'][max_index][0], '.6f'),
                  '+/-', format(math.sqrt(abs(df_both['Covar_CT'][max_index][0, 0])), '.6f'))
            print('l_CT (eV) : ', format(df_both['Fit_CT'][max_index][1], '.6f'),
                  '+/-', format(math.sqrt(abs(df_both['Covar_CT'][max_index][1, 1])), '.6f'))
            print('E_CT (eV) : ', format(df_both['Fit_CT'][max_index][2], '.6f'),
                  '+/-', format(math.sqrt(abs(df_both['Covar_CT'][max_index][2, 2])), '.6f'))

            if include_disorder:
                print('Sigma (eV) : ', format(df_both['Fit_CT'][max_index][3], '.6f'),
                      '+/-', format(math.sqrt(abs(df_both['Covar_CT'][max_index][3, 3])), '.6f'))
                W = df_both['Fit_CT'][max_index][1] * T + (df_both['Fit_CT'][max_index][3] ** 2) / (2 * k)
                print('Gaussian Variance [W] (eV K) : ', format(W, '.2f'))

            # print('Temperature [T] (K) : ', T)
            # print('-' * 80)

        # Save fit data
        if save_fit:
            opt_file = pd.DataFrame()
            opt_file['Energy'] = df_both['Energy'][max_index]
            opt_file['Signal'] = df_both['Opt_Fit'][max_index]
            opt_file['Temperature'] = T
            opt_file['Oscillator Strength (eV**2)'] = df_both['Fit_Opt'][max_index][0]
            opt_file['Reorganization Energy (eV)'] = df_both['Fit_Opt'][max_index][1]
            opt_file['Optical Peak Energy (eV)'] = df_both['Fit_Opt'][max_index][2]

            CT_file = pd.DataFrame()
            CT_file['Energy'] = df_both['Energy'][max_index]
            CT_file['Signal'] = df_both['CT_Fit'][max_index]
            CT_file['Temperature'] = T
            CT_file['Oscillator Strength (eV**2)'] = df_both['Fit_CT'][max_index][0]
            CT_file['Reorganization Energy (eV)'] = df_both['Fit_CT'][max_index][1]
            CT_file['CT State Energy (eV)'] = df_both['Fit_CT'][max_index][2]

            if include_disorder:
                CT_file['Sigma (eV)'] = df_both['Fit_CT'][max_index][3]

            save_fit_path, save_fit_filename = os.path.split(save_fit_file)
            if len(save_fit_path) != 0:  # Check if the user actually selected a path
                os.chdir(save_fit_path)  # Change the working directory
                if simultaneous_double:
                    opt_file.to_csv(f'{save_fit_filename}_Fit_simDouble_Opt_Range{n_fit}')  # Save data to csv
                    CT_file.to_csv(f'{save_fit_filename}_Fit_simDouble_CT_Range{n_fit}')  # Save data to csv
                else:
                    opt_file.to_csv(f'{save_fit_filename}_Fit_Opt_Range{n_fit}')  # Save data to csv
                    CT_file.to_csv(f'{save_fit_filename}_Fit_subOpt_CT_Range{n_fit}')  # Save data to csv

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

        if save_fit:
            if len(save_fit_path) != 0:  # Check if the user actually selected a path
                os.chdir(save_fit_path)  # Change the working directory
                if simultaneous_double:
                    plt.savefig((f'{save_fit_filename}_Fit_simDouble_Range{n_fit}.png'))
                else:
                    plt.savefig((f'{save_fit_filename}_Fit_Range{n_fit}.png'))


        df_copy = df_both.copy()
        df_copy['Total_R2'][max_index] = -10000

    return df_copy

# -----------------------------------------------------------------------------------------------------------
