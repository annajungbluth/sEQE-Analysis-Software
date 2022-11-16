import math

from source.utils import get_logger

logger = get_logger()


# -----------------------------------------------------------------------------------------------------------

# Function to compile EQE data

def compile_EQE(eqe_df,
                start,
                stop,
                number,
                precision=8
                ):
    """Function to compile EQE data based on start/stop values

    Parameters
    ----------
    eqe_df : dataFrame, required
        Dataframe of EQE values with columns ['Wavelength', ' Energy', 'EQE', 'Log_EQE']
    start : float, required
        Start wavelength or energy [eV/nm]
    stop : float, required
        Stop wavelength or energy [eV/nm]
    number : int, required
        Number indicating wavelength or energy compilation
        number = 0 => compile wavelength
        number = 1 => compile energy
    precision : int, optional
        Decimal point precision to compile data
        
    Returns
    -------
    Wavelength : list
        Wavelength values corresponding to compiled EQE [nm]
    Energy : list
        Energy values corresponding to compiled EQE [eV]
    EQE : list
        EQE values of compiled EQE
    log_EQE : list
        Logarithmic EQE values of compiled EQE spectra
    """

    # Define variables
    h = 6.626 * math.pow(10, -34)  # [m^2 kg/s]
    c = 2.998 * math.pow(10, 8)  # [m/s]
    q = 1.602 * math.pow(10, -19)  # [C]

    # Define empty lists
    Wavelength = []
    Energy = []
    EQE = []
    log_EQE = []

    if number == 0:  # If a wavelength range is given
        startNM = start
        stopNM = stop

    elif number == 1:  # If an energy range is given
        startNM = (h * c * math.pow(10, 9)) / (
                stop * q)  # The start wavelength corresponds to the high energy stop value
        stopNM = (h * c * math.pow(10, 9)) / (
                start * q)  # The stop wavelength corresponds to the low energy start value

    for y in range(len(eqe_df['Wavelength'])):  # Iterate through columns of EQE file
        if startNM <= eqe_df['Wavelength'][y] <= stopNM:  # Compile EQE if start <= wavelength <= stop, otherwise ignore
            Wavelength.append(round(eqe_df['Wavelength'][y], precision))
            Energy.append(round(eqe_df['Energy'][y], precision))
            EQE.append(round(eqe_df['EQE'][y], precision))  # * eqe_df['Energy'][y] - Eventually this might need to be added for fitting
            log_EQE.append(round(eqe_df['Log_EQE'][y], precision))

    if len(Wavelength) == len(EQE) and len(Energy) == len(log_EQE):  # Check that the lengths are the same
        return Wavelength, Energy, EQE, log_EQE

    else:
        logger.error('Error Code 1: Length mismatch.')


# -----------------------------------------------------------------------------------------------------------

# Function to compile EL data

def compile_EL(el_df,
               start,
               stop,
               number,
               precision=8
               ):
    """Function to compile EL data based on start/stop values

    Parameters
    ----------
    el_df : dataFrame, required
        Dataframe of EL values 
    start : float, required
        Start wavelength or energy [eV/nm]
    stop : float, required
        Stop wavelength or energy [eV/nm]
    number : int, required
        Number indicating wavelength or energy compilation
        number = 0 => compile wavelength
        number = 1 => compile energy
    precision : int, optional
        Decimal point precision to compile data
        
    Returns
    -------
    Wavelength : list
        Wavelength values corresponding to compiled EL [nm]
    Energy : list
        Energy values corresponding to compiled EL [eV]
    EL : list
        EL values of compiled EL spectra
    """

    # Define variables
    h = 6.626 * math.pow(10, -34)  # [m^2 kg/s]
    c = 2.998 * math.pow(10, 8)  # [m/s]
    q = 1.602 * math.pow(10, -19)  # [C]

    # Define empty lists
    Wavelength = []
    Energy = []
    EL = []

    if number == 0:  # If a wavelength range is given
        startNM = start
        stopNM = stop

    elif number == 1:  # If an energy range is given
        startNM = (h * c * math.pow(10, 9)) / (
                stop * q)  # The start wavelength corresponds to the high energy stop value
        stopNM = (h * c * math.pow(10, 9)) / (
                start * q)  # The stop wavelength corresponds to the low energy start value

    for y in range(len(el_df['Wavelength'])):  # Iterate through columns of EL file
        if startNM <= el_df['Wavelength'][y] <= stopNM:  # Compile EL if start <= wavelength <= stop, otherwise ignore
            Wavelength.append(round(el_df['Wavelength'][y], precision))
            Energy.append(round(el_df['Energy'][y], precision))
            EL.append(round(el_df['Signal'][y], precision))

    if len(Wavelength) == len(EL) and len(Energy) == len(EL):  # Check that the lengths are the same
        return Wavelength, Energy, EL

    else:
        logger.error('Error Code 1: Length mismatch.')


# -----------------------------------------------------------------------------------------------------------

# Function to compile any data

def compile_Data(energy,
                 y,
                 startE,
                 stopE,
                 precision=8
                 ):
    """Function to compile any data based on start/stop values

    Parameters
    ----------
    energy : list, required
        List of input energy values [eV]
    y : list, required
        List of input y values (e.g. EQE or EL)
    startE : float, required
        Start energy [eV]
    stopE : float, required
        Stop energy [eV]
    precision : int, optional
        Decimal point precision to compile data
        
    Returns
    -------
    Energy_comp : list
        List of compiled energy values [eV]
    y_comp : list
        List of compiled y values
    """

    # Define empty lists
    Energy_comp = []
    y_comp = []

    for x in range(len(energy)):
        if startE <= energy[x] <= stopE:  # Compile data only if start <= energy <= stop, otherwise ignore
            Energy_comp.append(round(energy[x], precision))
            y_comp.append(round(y[x], precision))

    if len(Energy_comp) == len(y_comp):  # Check that the lengths are the same
        return Energy_comp, y_comp

    else:
        logger.error('Error Code 1: Length mismatch.')

# -----------------------------------------------------------------------------------------------------------
