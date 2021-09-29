from source.utils import interpolate


# -----------------------------------------------------------------------------------------------------------

# Function to calculate the reference power

def calculate_Power(ref_df,
                    cal_df
                    ):
    """
    :param ref_df: dataFrame of measured reference diode values including columns ['Wavelength', 'Mean Current']
    :param cal_df: dataFrame of diode calibration values including columns ['Wavelength [nm]', 'Responsivity [A/W]']
    :return: ref_df['Power']: column of power values calculated in reference dataFrame
    """

    # Define empty dictionaries and lists
    cal_wave_dict = {}
    power = []

    for x in range(len(cal_df['Wavelength [nm]'])):  # Iterate through columns of calibration file
        cal_wave_dict[cal_df['Wavelength [nm]'][x]] = cal_df['Responsivity [A/W]'][
            x]  # Add wavelength and corresponding responsivity to dictionary

    for y in range(len(ref_df['Wavelength'])):  # Iterate through columns of reference file
        if ref_df['Wavelength'][y] in cal_wave_dict.keys():  # Check if reference wavelength is in calibration file
            power.append(float(ref_df['Mean Current'][y]) / float(
                cal_wave_dict[ref_df['Wavelength'][y]]))  # Add power to the list
        else:  # If reference wavelength is not in calibration file
            resp_int = interpolate(ref_df['Wavelength'][y], cal_df['Wavelength [nm]'],
                                   cal_df['Responsivity [A/W]'])  # Interpolate responsivity
            power.append(float(ref_df['Mean Current'][y]) / float(resp_int))  # Add power to the list

    ref_df['Power'] = power  # Create new column in reference file

    return ref_df['Power']

# -----------------------------------------------------------------------------------------------------------
