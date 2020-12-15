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

def R_squared(y_data, yfit_data):
    """
    :param y_data: original y data [list or array]
    :param yfit_data: y data of fit [list or array]
    :return: r_squared: R Squared of fit to data [float]
    """

    if len(y_data) == len(yfit_data):

        y_data = np.array(y_data)
        yfit_data = np.array(yfit_data)

        residuals = y_data - yfit_data
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)

        return r_squared
    else:
        print('Error Code 1: Length mismatch.')

# -----------------------------------------------------------------------------------------------------------
