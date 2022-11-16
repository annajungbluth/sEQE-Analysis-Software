from source.utils import get_logger, interpolate

logger = get_logger()

# -----------------------------------------------------------------------------------------------------------

# Function to normalize EQE data

def normalize_EQE(eqe_df,
                  startNM,
                  stopNM,
                  normNM
                  ):
    """Function to normalize EQE data

    Parameters
    ----------
    eqe_df : dataFrame, required
        DataFrame of EQE values with columns ['Wavelength', ' Energy', 'EQE', 'Log_EQE']
    startNM : float or int, required
        Start wavelength [nm]
    stopNM : float or int, required
        Stop wavelength [nm]
    normNM : float or int, required
        Normalization Wavelength [nm]
        
    Returns
    -------
    Wavelength : list
        List of wavelength values corresponding to normalized EQE
    Energy : list
        List of energy values corresponding to normalized EQE
    EQE : list
        List of normalized EQE values
    log_EQE : list
        List of normalized EQE values on a logarithmic scale
    """

    Wavelength = []
    Energy = []
    EQE = []
    log_EQE = []

    norm_EQE = interpolate(normNM, eqe_df['Wavelength'], eqe_df['EQE'])
    norm_log_EQE = interpolate(normNM, eqe_df['Wavelength'], eqe_df['Log_EQE'])

    for y in range(len(eqe_df['Wavelength'])):  # Iterate through columns of EQE file
        if startNM <= eqe_df['Wavelength'][y] <= stopNM:  # Compile EQE if start <= wavelength <= stop, otherwise ignore
            Wavelength.append(eqe_df['Wavelength'][y])
            Energy.append(eqe_df['Energy'][y])
            EQE.append(eqe_df['EQE'][y] / norm_EQE)
            log_EQE.append(eqe_df['Log_EQE'][y] / norm_log_EQE)

    if len(Wavelength) == len(EQE) and len(Energy) == len(log_EQE):  # Check that the lengths are the same
        return Wavelength, Energy, EQE, log_EQE

    else:
        logger.error('Error Code 1: Length mismatch.')

# -----------------------------------------------------------------------------------------------------------
