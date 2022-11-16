import numpy as np

from source.gaussian import calculate_gaussian_absorption

# -----------------------------------------------------------------------------------------------------------

# Function to subtract optical fit from EQE

def subtract_Opt(eqe,
                 best_vals,
                 T
                 ):
    """Function to subtract optical peak fit from EQE data

    Parameters
    ----------
    eqe : list, required
        EQE input data
    T : float, required
        EQE Measurement Temperature [K]
        
    Returns
    -------
    eqe : list
        EQE with subtracted optical fit
    """

    eqe = eqe.copy()

    Opt_fit = np.array(
        [calculate_gaussian_absorption(e,
                                       best_vals[0],
                                       best_vals[1],
                                       best_vals[2],
                                       T
                                       ) for e in eqe['Energy']])
    EQE_data = np.array(eqe['EQE'])

    subtracted_EQE = EQE_data - Opt_fit

    assert len(Opt_fit) == len(EQE_data)
    assert len(Opt_fit) == len(subtracted_EQE)

    eqe['EQE'] = subtracted_EQE

    return eqe

# -----------------------------------------------------------------------------------------------------------
