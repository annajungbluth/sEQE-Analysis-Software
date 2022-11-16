import math


# -----------------------------------------------------------------------------------------------------------

# Function to calculate the black body spectrum

def bb_spectrum(E_list,
                T_EL
                ):
    """Function to calculate black body spectrum

    Parameters
    ----------
    E_list : list, required
        List of input energy values [eV]
    T_EL : float, required
        Temperature of EL measurement [K]
        
    Returns
    -------
    phi_bb_dict : dict
        Dictionary of calculated black body spectrum
    """

    h_2 = 4.136 * math.pow(10, -15)  # [eV s]
    c = 2.998 * math.pow(10, 8)  # [m/s]
    k = 8.617 * math.pow(10, -5)  # [ev/K]

    phi_bb_dict = {}

    for energy in E_list:
        # phi_bb = 2 * math.pi * (energy) ** 2 * math.exp(-1 * energy / (k * T_EL)) / (
        #             h_2 ** 3 * c ** 2)  # -1) - without -1 as an approximation

        # this equation is confirmed in Thomas Kirchartz book chapter on EL
        phi_bb = (2 * math.pi * energy ** 2) / (h_2 ** 3 * c ** 2) * (1 / (math.exp(energy / (k * T_EL)) - 1))

        phi_bb_dict[energy] = phi_bb

    return phi_bb_dict  # [s/kg m^4]

# -----------------------------------------------------------------------------------------------------------
