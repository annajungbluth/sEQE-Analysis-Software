from source.utils import get_logger

logger = get_logger()


# -----------------------------------------------------------------------------------------------------------

# Function to check if reference & data files are non-empty and within wavelength range

def Ref_Data_is_valid(ref_df,
                      data_df,
                      startNM,
                      stopNM,
                      range_no
                      ):
    """Function to check if data files are non-empty and within wavelength range

    Parameters
    ----------
    ref_df : dataFrame, required
        DataFrame of reference diode measurements
    data_df : dataFrame, required
        DataFrame of sample measurements
    startNM : float or int, required
        Start wavelength [nm]
    stopNM : float or int, required
        Stop wavelength [nm]
    range_no : int
        Number of measurement range to reference correct and plot

    Returns
    -------
    X : bool
        Boolean value specifying whether range is valid
    """

    if len(ref_df) != 0 and len(data_df) != 0:

        if startNM >= data_df['Wavelength'][0] and \
                startNM >= ref_df['Wavelength'][0] and \
                stopNM <= data_df['Wavelength'][int(len(data_df['Wavelength'])) - 1] and \
                stopNM <= ref_df['Wavelength'][int(len(ref_df['Wavelength'])) - 1]:
            return True

        elif startNM < data_df['Wavelength'][0] or startNM < ref_df['Wavelength'][0]:
            logger.error('Please select a valid start wavelength for Range %s.' % str(range_no))
            return False

        elif stopNM > data_df['Wavelength'][int(len(data_df['Wavelength'])) - 1] or stopNM > ref_df['Wavelength'][
            int(len(ref_df['Wavelength'])) - 1]:
            logger.error('Please select a valid stop wavelength for Range %s.' % str(range_no))
            return False

        else:
            logger.error('Please select a valid wavelength range for Range %s.' % str(range_no))
            return False

    elif len(ref_df) == 0 and len(data_df) != 0:  # If reference file is empty / hasn't been selected
        logger.error('Please import a valid reference file for Range %s.' % str(range_no))
        return False

    elif len(ref_df) != 0 and len(data_df) == 0:  # If data file is empty / hasn't been selected
        logger.error('Please import a valid data file for Range %s.' % str(range_no))
        return False

    else:  # If reference and data files are empty / haven't been selected
        logger.error('Please import valid reference and data files for Range %s.' % str(range_no))
        return False


# -----------------------------------------------------------------------------------------------------------

# Function to check if EQE files are non-empty and within wavelength range

def EQE_is_valid(eqe_df,
                 startNM,
                 stopNM,
                 EQE_no
                 ):
    """Function to check if EQE files are non-empty and within wavelength range

    Parameters
    ----------
    eqe_df : dataFrame, required
        DataFrame of EQE values with columns ['Wavelength', ' Energy', 'EQE', 'Log_EQE']
    startNM : float or int, required
        Start wavelength [nm]
    stopNM : float or int, required
        Stop wavelength [nm]
    EQE_no : int
        Number of EQE file to plot

    Returns
    -------
    X : bool
        Boolean value specifying whether EQE range is valid
    """

    if len(eqe_df) != 0:

        if startNM >= eqe_df['Wavelength'][0] and \
                stopNM <= eqe_df['Wavelength'][int(len(eqe_df['Wavelength'])) - 1]:

            return True

        elif startNM < eqe_df['Wavelength'][0] and stopNM <= eqe_df['Wavelength'][int(len(eqe_df['Wavelength'])) - 1]:
            logger.error('Please select a valid start wavelength for EQE File %s.' % str(EQE_no))
            return False

        elif startNM >= eqe_df['Wavelength'][0] and stopNM > eqe_df['Wavelength'][int(len(eqe_df['Wavelength'])) - 1]:
            logger.error('Please select a valid stop wavelength for EQE File %s.' % str(EQE_no))
            return False

        else:
            logger.error('Please select a valid wavelength range for EQE File %s.' % str(EQE_no))
            return False

    else:  # If EQE file is empty / hasn't been selected
        logger.error('Please import a valid file for EQE File %s.' % str(EQE_no))
        return False


# -----------------------------------------------------------------------------------------------------------

# Function to check if data files are non-empty and within energy range

def Data_is_valid(df,
                  startE,
                  stopE
                  ):
    """Function to check if data files are non-empty and within wavelength range

    Parameters
    ----------
    df : dataFrame, required
        DataFrame of EQE values with columns ['Wavelength', ' Energy', 'EQE', 'Log_EQE']
    startE : float or int, required
        Start energy [eV]
    stopE : float or int, required
        Stop energy [eV]
    
    Returns
    -------
    X : bool
        Boolean value specifying whether data range is valid
    """

    if len(df) != 0:

        if startE <= max(df['Energy']) and stopE >= min(df['Energy']):
            return True

        elif startE > max(df['Energy']) and stopE >= min(df['Energy']):
            logger.error('Please select a valid start energy.')
            return False

        elif startE <= max(df['Energy']) and stopE < min(df['Energy']):
            logger.error('Please select a valid stop energy.')
            return False

        else:
            logger.error('Please select a valid energy range.')
            return False

    else:  # If file is empty / hasn't been selected
        logger.error('Please import a valid file.')
        return False


# -----------------------------------------------------------------------------------------------------------

# Function to check if normalization wavelength is within wavelength range

def Normalization_is_valid(eqe_df,
                           normNM,
                           EQE_no
                           ):
    """Function to check if normalization wavelength is within wavelength range

    Parameters
    ----------
    eqe_df : dataFrame, required
        DataFrame of EQE values with columns ['Wavelength', ' Energy', 'EQE', 'Log_EQE']
    normNM : float or int, required
        Normalization wavelength [nm]
    EQE_no : int
        Number of EQE file to plot

    Returns
    -------
    X : bool
        Boolean value specifying whether normalization is valid
    """

    if len(eqe_df) != 0:

        min_wave = int(eqe_df['Wavelength'].min())
        max_wave = int(eqe_df['Wavelength'].max())

        if min_wave <= int(normNM) <= max_wave:
            return True

        else:
            logger.error('Please select a valid normalization wavelength for EQE file %s.' % str(EQE_no))
            return False


# -----------------------------------------------------------------------------------------------------------

# Function to check if EQE files are non-empty and within fit energy range

def Fit_is_valid(eqe_df,
                 startE,
                 stopE,
                 startFitE,
                 stopFitE,
                 EQE_no
                 ):
    """Function to check if EQE files are non-empty and within fit energy range

    Parameters
    ----------
    eqe_df : dataFrame, required
        DataFrame of EQE measurements
    startE : float or int, required
        Start energy [eV]
    stopE : float or int, required
        Stop energy [eV]
    startFitE : float
        Start energy of fit
    stopFitE : float
        Stop energy of fit
    EQE_no : int
        Number of EQE file to fit

    Returns
    -------
    X : bool
        Boolean value specifying whether fit is valid
    """

    if len(eqe_df) != 0:

        if startE <= eqe_df['Energy'][0] and \
                startFitE <= eqe_df['Energy'][0] and \
                stopE >= eqe_df['Energy'][int(len(eqe_df['Energy'])) - 1] and \
                stopFitE >= eqe_df['Energy'][int(len(eqe_df['Energy'])) - 1]:
            return True

        elif startE > eqe_df['Energy'][0] and stopE >= eqe_df['Energy'][int(len(eqe_df['Energy'])) - 1]:
            logger.error('Please select a valid start energy for EQE File %s.' % str(EQE_no))
            return False

        elif startE <= eqe_df['Energy'][0] and stopE < eqe_df['Energy'][int(len(eqe_df['Energy'])) - 1]:
            logger.error('Please select a valid stop energy for EQE File %s.' % str(EQE_no))
            return False

        else:
            logger.error('Please select a valid energy range for EQE File %s.' % str(EQE_no))
            return False

    else:  # If EQE file is empty / hasn't been selected
        logger.error('Please import a valid file for EQE File %s.' % str(EQE_no))
        return False


# -----------------------------------------------------------------------------------------------------------

# Function to check if start fit energy is smaller than stop fit energy

def StartStop_is_valid(start, stop):
    """Function to check if fit energy is smaller than stop fit energy

    Parameters
    ----------
    start : float, required
        Start value (i.e. start energy or fit start)
    stop : float, required
        Stop value (i.e. stop energy or fit stop)

    Returns
    -------
    X : bool
        Boolean value specifying whether start-stop combination is valid
    """

    if start < stop:
        return True
    else:
        logger.info('Please select valid start and stop energies.')
        return False

# -----------------------------------------------------------------------------------------------------------
