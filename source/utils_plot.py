import random
import matplotlib.pyplot as plt

from colour import Color
from numpy import random

from source.utils import get_logger

logger = get_logger()


# -----------------------------------------------------------------------------------------------------------

# Function to check if input is a colour

def is_Colour(colour):
    """Function to check if input is a colour

    Parameters
    ----------
    colour : str, required
        Colour input

    Returns
    -------
    Colour : bool
        Boolean value specifying whether colour is valid

    """

    try:
        Color(colour)
        return True
    except:
        return False


# -----------------------------------------------------------------------------------------------------------

# Function to pick plot colour

def pick_EQE_Color(colour_Box, file_no):
    """Function to pick plot colour

    Parameters
    ----------
    colour_Box : gui object, required
        GUI textbox with colour information
    file_no : str or int, required
        Indicator of associated file 

    Returns
    -------
    Colour : HEX
        HEX code of randomly generated colour

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
                logger.error('Please name a valid colour.')
            else:
                logger.error('Please name a valid colour for EQE File %s.' % str(file_no))
            return random_colour
    else:
        return random_colour


# -----------------------------------------------------------------------------------------------------------

# Function to pick EQE plot label

def pick_EQE_Label(label_Box, filename_Box):
    """Function to pick EQE plot label

    Parameters
    ----------
    label_Box : gui object, required
        GUI textbox with label information
    filename_Box : gui object, required
        GUI textbox with filename information

    Returns
    -------
    label/filename : str
        Label or filename to plot
    """

    label = label_Box.toPlainText()
    filename = filename_Box.toPlainText()

    if len(label) != 0:
        return label
    else:  # Don't need to check that filename exists, "pick_label" function is called after checking "EQE_is_valid"
        return filename

# -----------------------------------------------------------------------------------------------------------

# Function to pick "Calculate EQE" plot label

def pick_Label(range_no, startNM, stopNM):
    """Function to pick "Calculate EQE" plot label

    Parameters
    ----------
    range_no : int, required
        Number of data range to plot
    startNM : float, required
        Start wavelength to plot [nm]
    stopNM : float, required
        Stop wavelength to plot [nm]

    Returns
    -------
    label : str
        Label to plot
    """

    label = str("Range") + str(range_no) + "_" + str(int(startNM)) + "nm" + "-" + str(
        int(stopNM)) + "nm"  # Example: Range1_360nm_800nm
    return label

# -----------------------------------------------------------------------------------------------------------
