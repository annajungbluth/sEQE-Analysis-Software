import math

# -----------------------------------------------------------------------------------------------------------

### Function to caluclate gaussian absorption

def calculate_gaussian_absorption(x, f, l, E, T):
    """
    :param x: List of energy values
    :param f: Oscillator strength
    :param l: Reorganization Energy
    :param E: Peak Energy
    :param T: Temperature
    :return: EQE value
    """

    # Define variables
    k = 8.617 * math.pow(10, -5)  # [ev/K]

    return (f / (x * math.sqrt(4 * math.pi * l * T * k))) * exp(-(E + l - x) ** 2 / (4 * l * k * T))

# -----------------------------------------------------------------------------------------------------------