from source.utils import get_logger, interpolate

logger = get_logger()

# -----------------------------------------------------------------------------------------------------------

# Function to normalize EQE data

def normalize_EQE(eqe_df, startNM, stopNM, normNM):
    """
    :param eqe_df: dataFrame of EQE values with columns ['Wavelength', ' Energy', 'EQE', 'Log_EQE']
    :param startNM: start wavelength [float or int]
    :param stopNM: stop wavelength [float or int]
    :param normNM: normalization wavelength [float or int]
    :return: Wavelength: list of wavelength values [list]
             Energy: list of energy values [list]
             EQE: list of normalized EQE values [list]
             log_EQE: list of normalized log EQE values [list]
    """

    # Define empty lists
    Wavelength = []
    Energy = []
    EQE = []
    log_EQE = []

    norm_EQE = interpolate(normNM, eqe_df['Wavelength'], eqe_df['EQE'])
    norm_log_EQE = interpolate(normNM, eqe_df['Wavelength'], eqe_df['Log_EQE'])

    for y in range(len(eqe_df['Wavelength'])):  # Iterate through columns of EQE file
        if startNM <= eqe_df['Wavelength'][y] <= stopNM:  # Compile EQE only if start <= wavelength <= stop, otherwise ignore
            Wavelength.append(eqe_df['Wavelength'][y])
            Energy.append(eqe_df['Energy'][y])
            EQE.append(eqe_df['EQE'][y] / norm_EQE)
            log_EQE.append(eqe_df['Log_EQE'][y] / norm_log_EQE)

    if len(Wavelength) == len(EQE) and len(Energy) == len(log_EQE):  # Check that the lengths are the same
        return Wavelength, Energy, EQE, log_EQE

    else:
        logger.error('Error Code 1: Length mismatch.')

# -----------------------------------------------------------------------------------------------------------
