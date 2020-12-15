import math

# -----------------------------------------------------------------------------------------------------------

### Function to compile EL data

def compile_EL(el_df, start, stop, number):
    """
    :param el_df: dataFrame of EL values
    :param start: start wavelength or energy [float or int]
    :param stop: stop wavelength or energy [float or int]
    :param number: file number [int]
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
        startNM = (h * c * math.pow(10, 9)) / (stop * q)  # The start wavelength corresponds to the high energy stop value
        stopNM = (h * c * math.pow(10, 9)) / (start * q)  # The stop wavelength corresponds to the low energy start value

    for y in range(len(el_df['Wavelength'])):  # Iterate through columns of EL file
        if startNM <= el_df['Wavelength'][y] <= stopNM:  # Compile EL only if start <= wavelength <= stop, otherwise ignore
            Wavelength.append(el_df['Wavelength'][y])
            Energy.append(el_df['Energy'][y])
            EL.append(el_df['Signal'][y])

    if len(Wavelength) == len(EL) and len(Energy) == len(EL):  # Check that the lengths are the same
        return Wavelength, Energy, EL

    else:
        print('Error Code 1: Length mismatch.')

# -----------------------------------------------------------------------------------------------------------

### Function to compile any data

def compile_Data(energy, y, startE, stopE):
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
        if startE <= energy[x] <= stopE :  # Compile data only if start <= energy <= stop, otherwise ignore
            Energy_comp.append(energy[x])
            y_comp.append(y[x])

    if len(Energy_comp) == len(y_comp):  # Check that the lengths are the same
        return Energy_comp, y_comp

    else:
        print('Error Code 1: Length mismatch.')

# -----------------------------------------------------------------------------------------------------------