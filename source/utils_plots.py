import random
from colour import Color
from numpy import random

# -----------------------------------------------------------------------------------------------------------

### Function to check if input is a color

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

### Function to pick plot color

def pick_EQE_Color(colour_Box, file_no):
    """
    :param colour_Box: textbox with colour information [from ui]
    :param file_no: indicator of associated file [string or int
    :return: HEX code of random colour [string]
    """

    colour = colour_Box.toPlainText()
    colour = colour.replace(" ", "") # If the user inputs "Sky Blue" instead of "SkyBlue" etc.

    color_comp = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']
    color_choice = [random.choice(color_comp) for j in range(6)]
    random_colour = '#'+''.join(color_choice)

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