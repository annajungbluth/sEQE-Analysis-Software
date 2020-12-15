import random

from colour import Color
from numpy import random


# -----------------------------------------------------------------------------------------------------------

# Function to check if input is a color

def is_Colour(colour):
    """
    :param colour: colour input [string]
    :return: True if colour is valide, False if it is not
    """

    try:
        Color(colour)
        return True
    except:
        return False


# -----------------------------------------------------------------------------------------------------------

# Function to pick plot color

def pick_EQE_Color(colour_Box, file_no):
    """
    :param colour_Box: textbox with color information [from ui]
    :param file_no: indicator of associated file [string or int
    :return: HEX code of random color [string]
    """

    colour = colour_Box.toPlainText()
    colour = colour.replace(" ", "")  # If the user inputs "Sky Blue" instead of "SkyBlue" etc.

    color_comp = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']
    color_choice = [random.choice(color_comp) for j in range(6)]
    random_colour = '#' + ''.join(color_choice)

    if len(colour) != 0:
        if is_Colour(colour):
            return colour
        else:
            if file_no == 100:
                print('Please name a valid colour.')
            else:
                print('Please name a valid colour for EQE File %s.' % str(file_no))
            return random_colour
    else:
        return random_colour


# -----------------------------------------------------------------------------------------------------------

# Function to pick EQE plot label

def pick_EQE_Label(label_Box, filename_Box):
    """
    :param label_Box: textbox with label information [from ui]
    :param filename_Box: textbox with file name [from ui]
    :return: filename: string of filename / label to plot [string]
    """

    label = label_Box.toPlainText()
    filename = filename_Box.toPlainText()

    if len(label) != 0:
        return label
    else:  # We don't need to check that there is a filename, as the "pick_label" function is only called after checking "EQE_is_valid"
        return filename


# -----------------------------------------------------------------------------------------------------------

# Function to pick "Calculate EQE" plot label

def pick_Label(range_no, startNM, stopNM):
    """
    :param range_no: number of range to plot [int]
    :param startNM: start wavelength [float]
    :param stopNM: stop wavelength [float
    :return: label: string of the lable to plot [string]
    """

    label = str("Range") + str(range_no) + "_" + str(int(startNM)) + "nm" + "-" + str(
        int(stopNM)) + "nm"  # Example: Range1_360nm_800nm
    return label

# -----------------------------------------------------------------------------------------------------------
