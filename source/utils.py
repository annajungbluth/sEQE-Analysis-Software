from scipy.interpolate import interp1d

# Function to interpolate values

def interpolate(num, x, y):
    """
    :param num: value to interpolate at
    :param x: x data
    :param y: y data
    :return: interpolated y-value at x-value of num
    """
    f = interp1d(x, y)
    return f(num)