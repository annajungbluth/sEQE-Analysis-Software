import math

# Function to calculate the black body spectrum

def bb_spectrum(E_list, T_EL):
    """
    :param E_list: list of energy values in eV
    :param T_EL: Temperature of EL measurement
    :return: phi_bb_dict: dictionary of calculated black body spectrum
    """

    h_2 = 4.136 * math.pow(10, -15)  # [eV s]
    c = 2.998 * math.pow(10, 8)  # [m/s]
    k = 8.617 * math.pow(10, -5)  # [ev/K]

    phi_bb_dict = {}

    for energy in E_list:
        phi_bb = 2 * math.pi * (energy)**2 * math.exp(-1 * energy / (k * T_EL)) / (h_2**3 * c**2) # -1) - without -1 as an approximation
        phi_bb_dict[energy] = phi_bb

    return phi_bb_dict #[s/kg m^4]