# -----------------------------------------------------------------------------------------------------------

# Function to check if reference & data files are non-empty and within wavelength range

def Ref_Data_is_valid(ref_df, data_df, startNM, stopNM, range_no):
    """
    :param ref_df: dataFrame of reference diode
    :param data_df: dataFrame of measurement data
    :param startNM: start wavelength [float or int]
    :param stopNM: stop wavelength [float or int]
    :param range_no: number of measurment range to plot / reference correct [int]
    :return: True or False
    """

    if len(ref_df) != 0 and len(data_df) != 0:

        if startNM >= data_df['Wavelength'][0] and \
                startNM >= ref_df['Wavelength'][0] and \
                stopNM <= data_df['Wavelength'][int(len(data_df['Wavelength'])) - 1] and \
                stopNM <= ref_df['Wavelength'][int(len(ref_df['Wavelength'])) - 1]:
            return True

        elif startNM < data_df['Wavelength'][0] or startNM < ref_df['Wavelength'][0]:
            print('Please select a valid start wavelength for Range %s.' % str(range_no))
            return False

        elif stopNM > data_df['Wavelength'][int(len(data_df['Wavelength'])) - 1] or stopNM > ref_df['Wavelength'][
            int(len(ref_df['Wavelength'])) - 1]:
            print('Please select a valid stop wavelength for Range %s.' % str(range_no))
            return False

        else:
            print('Please select a valid wavelength range for Range %s.' % str(range_no))
            return False

    elif len(ref_df) == 0 and len(data_df) != 0:  # If reference file is empty / hasn't been selected
        print('Please import a valid reference file for Range %s.' % str(range_no))
        return False

    elif len(ref_df) != 0 and len(data_df) == 0:  # If data file is empty / hasn't been selected
        print('Please import a valid data file for Range %s.' % str(range_no))
        return False

    else:  # If reference and data files are empty / haven't been selected
        print('Please import valid reference and data files for Range %s.' % str(range_no))
        return False


# -----------------------------------------------------------------------------------------------------------

# Function to check if EQE files are non-empty and within wavelength range

def EQE_is_valid(eqe_df, startNM, stopNM, EQE_no):
    """
    :param eqe_df: dataFrame of EQE values with columns ['Wavelength', ' Energy', 'EQE', 'Log_EQE']
    :param startNM: start wavelength [float or int]
    :param stopNM: stop wavelength [float or int]
    :param EQE_no: number of EQE file to plot [int]
    :return: True or False
    """

    if len(eqe_df) != 0:

        if startNM >= eqe_df['Wavelength'][0] and \
                stopNM <= eqe_df['Wavelength'][int(len(eqe_df['Wavelength'])) - 1]:

            return True

        elif startNM < eqe_df['Wavelength'][0] and stopNM <= eqe_df['Wavelength'][int(len(eqe_df['Wavelength'])) - 1]:
            print('Please select a valid start wavelength for EQE File %s.' % str(EQE_no))
            return False

        elif startNM >= eqe_df['Wavelength'][0] and stopNM > eqe_df['Wavelength'][int(len(eqe_df['Wavelength'])) - 1]:
            print('Please select a valid stop wavelength for EQE File %s.' % str(EQE_no))
            return False

        else:
            print('Please select a valid wavelength range for EQE File %s.' % str(EQE_no))
            return False

    else:  # If EQE file is empty / hasn't been selected
        print('Please import a valid file for EQE File %s.' % str(EQE_no))
        return False


# -----------------------------------------------------------------------------------------------------------

# Function to check if data files are non-empty and within energy range

def Data_is_valid(df, startE, stopE):
    """
    :param df: dataFrame of EQE values with columns ['Wavelength', ' Energy', 'EQE', 'Log_EQE']
    :param startE: start energy [float]
    :param stopE: stop energy [float]
    :return: True or False
    """

    if len(df) != 0:

        if startE <= max(df['Energy']) and stopE >= min(df['Energy']):
            return True

        elif startE > max(df['Energy']) and stopE >= min(df['Energy']):
            print('Please select a valid start energy.')
            return False

        elif startE <= max(df['Energy']) and stopE < min(df['Energy']):
            print('Please select a valid stop energy.')
            return False

        else:
            print('Please select a valid energy range.')
            return False

    else:  # If file is empty / hasn't been selected
        print('Please import a valid file.')
        return False


# -----------------------------------------------------------------------------------------------------------

# Function to check if normalization wavelength is within wavelength range

def Normalization_is_valid(eqe_df, normNM, EQE_no):
    """
    :param eqe_df: dataFrame of EQE values with columns ['Wavelength', ' Energy', 'EQE', 'Log_EQE']
    :param normNM: normalization wavelength [float or int]
    :param EQE_no: number of EQE file to normalize [int]
    :return: True or False
    """

    if len(eqe_df) != 0:

        min_wave = int(eqe_df['Wavelength'].min())
        max_wave = int(eqe_df['Wavelength'].max())

        if min_wave <= int(normNM) <= max_wave:
            return True

        else:
            print('Please select a valid normalization wavelength for EQE file %s.' % str(EQE_no))
            return False


# -----------------------------------------------------------------------------------------------------------

# Function to check if EQE files are non-empty and within fit energy range

def Fit_is_valid(eqe_df, startE, stopE, startFitE, stopFitE, EQE_no):
    """
    :param eqe_df: dataFrame of EQE values with columns ['Wavelength', ' Energy', 'EQE', 'Log_EQE']
    :param startE: EQE start energy [float]
    :param stopE: EQE stop energy [float]
    :param startFitE: start energy of fit [float]
    :param stopFitE: stop energy of fit [float]
    :param EQE_no: number of EQE file to fit [int]
    :return: True or False
    """

    if len(eqe_df) != 0:

        if startE <= eqe_df['Energy'][0] and \
                startFitE <= eqe_df['Energy'][0] and \
                stopE >= eqe_df['Energy'][int(len(eqe_df['Energy'])) - 1] and \
                stopFitE >= eqe_df['Energy'][int(len(eqe_df['Energy'])) - 1]:
            return True

        elif startE > eqe_df['Energy'][0] and stopE >= eqe_df['Energy'][int(len(eqe_df['Energy'])) - 1]:
            print('Please select a valid start energy for EQE File %s.' % str(EQE_no))
            return False

        elif startE <= eqe_df['Energy'][0] and stopE < eqe_df['Energy'][int(len(eqe_df['Energy'])) - 1]:
            print('Please select a valid stop energy for EQE File %s.' % str(EQE_no))
            return False

        else:
            print('Please select a valid energy range for EQE File %s.' % str(EQE_no))
            return False

    else:  # If EQE file is empty / hasn't been selected
        print('Please import a valid file for EQE File %s.' % str(EQE_no))
        return False


# -----------------------------------------------------------------------------------------------------------

# Function to check if start fit energy is smaller than stop fit energy

def StartStop_is_valid(start, stop):
    """
    :param start: start value (i.e. start energy or fit start) [float]
    :param stop: stop value (i.e. stop energy or fit stop) [float]
    :return: True or False
    """
    if start < stop:
        return True
    else:
        print('Please select valid start and stop energies.')
        return False

# -----------------------------------------------------------------------------------------------------------
