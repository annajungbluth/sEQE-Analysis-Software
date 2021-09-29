import math

from source.utils import get_logger

logger = get_logger()


# -----------------------------------------------------------------------------------------------------------

# Function to compile EQE data

def compile_EQE(eqe_df,
                start,
                stop,
                number
                ):
    """
    :param eqe_df: dataFrame of EQE values with columns ['Wavelength', ' Energy', 'EQE', 'Log_EQE']
    :param start: start wavelength of energy [float or int]
    :param stop: stop wavelength or energy [float or int]
    :param number: number indicating wavelength or energy compilation [int]
                   number = 0 => compile wavelength
                   number = 1 => compile energy
    :return: Wavelength: list of wavelength values [list]
             Energy: list of energy values [list]
             EQE: list of compiled EQE values [list]
             log_EQE: list of compiled log EQE values [list]
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
            Wavelength.append(eqe_df['Wavelength'][y])
            Energy.append(eqe_df['Energy'][y])
            EQE.append(eqe_df['EQE'][y])  # * eqe_df['Energy'][y] - Eventually this might need to be added for fitting
            log_EQE.append(eqe_df['Log_EQE'][y])
            # (math.log10(eqe_df['EQE'][y] * eqe_df['Energy'][y]) - Eventually this might need to be added for fitting

    if len(Wavelength) == len(EQE) and len(Energy) == len(log_EQE):  # Check that the lengths are the same
        return Wavelength, Energy, EQE, log_EQE

    else:
        logger.error('Error Code 1: Length mismatch.')


# -----------------------------------------------------------------------------------------------------------

# Function to compile EL data

def compile_EL(el_df,
               start,
               stop,
               number
               ):
    """
    :param el_df: dataFrame of EL values
    :param start: start wavelength or energy [float or int]
    :param stop: stop wavelength or energy [float or int]
    :param number: number indicating wavelength or energy compilation [int]
                   number = 0 => compile wavelength
                   number = 1 => compile energy
    :return: Wavelength: list of compiled wavelength values [list]
             Energy: list of compiled energy values [list]
             EL: list of compiled EL values [list]
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
            Wavelength.append(el_df['Wavelength'][y])
            Energy.append(el_df['Energy'][y])
            EL.append(el_df['Signal'][y])

    if len(Wavelength) == len(EL) and len(Energy) == len(EL):  # Check that the lengths are the same
        return Wavelength, Energy, EL

    else:
        logger.error('Error Code 1: Length mismatch.')


# -----------------------------------------------------------------------------------------------------------

# Function to compile any data

def compile_Data(energy,
                 y,
                 startE,
                 stopE
                 ):
    """
    :param energy: list of energy values [list]
    :param y: list of y values [list]
    :param startE: start energy [float]
    :param stopE: stop energy [float]
    :return: Energy_comp: list of compiled energy values [list]
             y_comp: list of compiled y values [float]
    """

    # Define empty lists
    Energy_comp = []
    y_comp = []

    for x in range(len(energy)):
        if startE <= energy[x] <= stopE:  # Compile data only if start <= energy <= stop, otherwise ignore
            Energy_comp.append(energy[x])
            y_comp.append(y[x])

    if len(Energy_comp) == len(y_comp):  # Check that the lengths are the same
        return Energy_comp, y_comp

    else:
        logger.error('Error Code 1: Length mismatch.')

# -----------------------------------------------------------------------------------------------------------
