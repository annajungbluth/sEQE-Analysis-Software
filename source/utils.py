import logging
import sys
import numpy as np
from scipy.interpolate import interp1d


# -----------------------------------------------------------------------------------------------------------

# Function to interpolate values

def interpolate(num, x, y):
    """Function to interpolate values

    Parameters
    ----------
    num : float, required
        Value to interpolate
    x : list or array, required
        x data for interpolation
    y : list or array, required
        y data for interpolation
    
    Returns
    -------
    f(num) : float
        Interpolated y value
    """

    f = interp1d(x, y)
    return f(num)


# -----------------------------------------------------------------------------------------------------------

# Function to calculate R Squared of fit

def R_squared(y_data,
              yfit_data,
              bias=False,
              tolerance=None
              ):
    """Function to calculate R squared of fit

    Parameters
    ----------
    y_data : list or array, required
        Input y data
    y_fit_data : list or array, required
        y data of the fit
    bias : bool, optional
        Boolean value specifying whether to bias fits above the data
    tolerance : float, optional
        Tolerance (mean percent) allowed for fit above the data
    
    Returns
    -------
    r_squared: float
        R squared of fit to data
    r_squared_log: float
        R squared of the logarithm of the fit to the data
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

        if bias and tolerance is not None:
            if len(residuals[residuals < 0]) != 0:  # Check that there are fit values that fall above the data
                mean_percent = abs(np.nanmean(residuals[residuals < 0] / y_data[residuals < 0]))
                # mean_percent = abs(np.nanmean(residuals/y_data)) # To account for positive and negative deviations
                if mean_percent > tolerance:
                    r_squared = 0.0001
                    r_squared_log = 0.0001  # These values are artificially set very low to discard fit

        return (r_squared + r_squared_log) / 2
    else:
        logger.error('Error Code 1: Length mismatch.')


# -----------------------------------------------------------------------------------------------------------

# Function to separate out elements from a list

def sep_list(list_, n):
    """Function to separate out elements from a list

    Parameters
    ----------
    list_ : list, required
        List of values to separate
    n : int, required
        Number of element to extract

    Returns
    -------
        nth element of list
    """
    return list_[n]

# -----------------------------------------------------------------------------------------------------------

# Function to separate a list of lists

def sep_list_list(list_list):
    """Function to separate a list of lists

    Parameters
    ----------
    list_list : list, required
        List of lists

    Returns
    -------
    one_list : list
        List of all individual elements in list_list
    """

    one_list = []

    for list in list_list:
        for x in list:
            one_list.append(x)

    return one_list


# -----------------------------------------------------------------------------------------------------------

# Function to set up logger

def get_logger():
    """Function to set up logger

    Parameters
    ----------
    None

    Returns
    -------
    logger : logger instance
        Logger to log terminal outputs to
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

# -----------------------------------------------------------------------------------------------------------
