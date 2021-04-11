import math

from numpy import exp

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

    return (f / (x * math.sqrt(2 * math.pi * (2 * l * T * k + sig**2))) * exp(-(E + l - x) ** 2 / (4 * l * k * T + 2 * sig**2)))

# -----------------------------------------------------------------------------------------------------------
