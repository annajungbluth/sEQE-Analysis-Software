from colour import Color

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
