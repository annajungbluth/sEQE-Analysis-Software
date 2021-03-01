import logging
import sys
import numpy as np
from scipy.interpolate import interp1d

# -----------------------------------------------------------------------------------------------------------

# Function to interpolate values

def interpolate(num, x, y):
    """
    :param num: value to interpolate at [float]
    :param x: x data [list or array]
    :param y: y data [list or array]
    :return: interpolated y-value at x-value of num [float]
    """
    f = interp1d(x, y)
    return f(num)


# -----------------------------------------------------------------------------------------------------------

# Function to calculate R Squared of fit

def R_squared(y_data, yfit_data, bias = False, tolerance = None):
    """
    :param y_data: original y data [list or array]
    :param yfit_data: y data of fit [list or array]
    :param bias: add bias to R-Squared if True [bool]
    "param tolerance: allowed mean percent deviated above the data [float]
    :return: r_squared: R Squared of fit to data [float]
             r_squared_log: R Squared of the log of the fit to the data [float]
    """

    if len(y_data) == len(yfit_data):

        y_data = np.array(y_data)
        yfit_data = np.array(yfit_data)

        residuals = y_data - yfit_data
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)

        # Calculate the R Squared of the log of the data

        y_data_log = np.array(np.log(y_data))
        yfit_data_log = np.array(np.log(yfit_data))

        residuals_log = y_data_log - yfit_data_log

        ss_res_log = np.sum(residuals_log ** 2)
        ss_tot_log = np.sum((y_data_log - np.mean(y_data_log)) ** 2)
        r_squared_log = 1 - (ss_res_log / ss_tot_log)

        if bias == True and tolerance is not None:
            mean_percent = abs(np.mean(residuals[residuals<0]/y_data[residuals<0]))
            if mean_percent > tolerance:
                r_squared = 0.0001
                r_squared_log = 0.0001 # These values are artificially set very low to discard fit

        return (r_squared + r_squared_log)/2
    else:
        logger.error('Error Code 1: Length mismatch.')

# -----------------------------------------------------------------------------------------------------------

# Function to separate out elements from a list

def sep_list(list_, n):
    """
    :param list_: list of values [list]
    :param n: number of element to extract [int]
    :return: nth element of list
    """
    return list_[n]

# -----------------------------------------------------------------------------------------------------------

# Function to separate a list of lists

def sep_list_list(list_list):
    """
    :param list_list: list of lists [list]
    :return: one_list: list of all individual elements in list_list [list]
    """
    one_list = []

    for list in list_list:
        for x in list:
            one_list.append(x)

    return one_list

# -----------------------------------------------------------------------------------------------------------

# Function to set up logger

def get_logger():
    """
    Return a logger for current module
    Returns
    -------
    logger : logger instance
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(fmt="%(asctime)s %(levelname)s %(name)s: %(message)s",
                                  datefmt="%Y-%m-%d - %H:%M:%S")
    if logger.hasHandlers():
        logger.handlers.clear()

    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    console.setFormatter(formatter)

    logger.addHandler(console)

    return logger

logger = get_logger()