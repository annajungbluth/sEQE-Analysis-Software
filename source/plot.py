import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


# -----------------------------------------------------------------------------------------------------------

# Function to plot any data

def plot(ax1,
         ax2,
         x,
         y,
         label_,
         color_
         ):
    """Function to plot any data

    Parameters
    ----------
    ax_1 : axis object, required
        matplotlib axis object
    ax_2 : axis object, required
        matplotlib axis object
    x : list or array, required
        x input data
    y : list or array, required
        y input data
    label_ : str, required
        Plot label
    color_ : str, required
        Plot color
        
    Returns
    -------
    None
    """

    ax1.plot(x, y, linewidth=3, label=label_, color=color_)
    ax1.legend()

    ax2.plot(x, y, linewidth=3, label=label_, color=color_)
    ax2.legend()

    plt.draw()


# -----------------------------------------------------------------------------------------------------------

# Function to set up "Calculate EQE" plot

def set_up_plot(flag='Wavelength'):
    """Function to set up "Calculate EQE" plot

    Parameters
    ----------
    flag : str, optional
        Flag to specify whether to plot wavelength or energy on x-axis.
        Options: "Wavelength", "Energy"
        
    Returns
    -------
    ax_1 : axis object
        matplotlib axis object to plot EQE on linear scale
    ax_2 : axis object
        matplotlib axis object to plot EQE on logarithmic scale
    """

    if flag == 'Wavelength':

        # style.use('ggplot')
        fig1 = plt.figure()

        ax1 = fig1.add_subplot(2, 1, 1)
        plt.ylabel('EQE', fontsize=17, fontweight='medium')
        plt.grid(False)
        # plt.box()
        plt.rcParams['figure.facecolor'] = 'xkcd:white'
        plt.rcParams['figure.edgecolor'] = 'xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
        plt.minorticks_on()

        ax2 = fig1.add_subplot(2, 1, 2)
        ax2.set_yscale('log')
        plt.xlabel('Wavelength (nm)', fontsize=17, fontweight='medium')
        plt.ylabel('EQE', fontsize=17, fontweight='medium')
        plt.grid(False)
        # plt.box()
        plt.rcParams['figure.facecolor'] = 'xkcd:white'
        plt.rcParams['figure.edgecolor'] = 'xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
        plt.minorticks_on()
        plt.show()

    elif flag == 'Energy':

        fig1 = plt.figure(figsize=(7, 10))

        ax1 = fig1.add_subplot(2, 1, 1)
        plt.ylabel('EQE', fontsize=17, fontweight='medium')
        plt.grid(False)
        # plt.box()
        plt.rcParams['figure.facecolor'] = 'xkcd:white'
        plt.rcParams['figure.edgecolor'] = 'xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
        plt.minorticks_on()

        ax2 = fig1.add_subplot(2, 1, 2)
        ax2.set_yscale('log')
        plt.xlabel('Energy (eV)', fontsize=17, fontweight='medium')
        plt.ylabel('EQE', fontsize=17, fontweight='medium')
        plt.grid(False)
        # plt.box()
        plt.rcParams['figure.facecolor'] = 'xkcd:white'
        plt.rcParams['figure.edgecolor'] = 'xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
        plt.minorticks_on()
        plt.show()

    return ax1, ax2


# -----------------------------------------------------------------------------------------------------------

# Function to set up EQE plot

def set_up_EQE_plot(number=None,
                    norm_num=None
                    ):
    """Function to set up EQE plot

    Parameters
    ----------
    number : int, optional
        Number indicating whether wavelength or energy are plotted
        number = 0 => plot wavelength
        number = 1 => plot energy
        number = None => plot energy
    norm_num : int, optional
        Number indicating whether raw, normalized, or reduced data are plotted [int]
        norm_num = 0 => plot raw EQE
        norm_num = 1 => plot normalized EQE
        norm_num = None => plot raw EQE
        
    Returns
    -------
    ax_1 : axis object
        matplotlib axis object to plot EQE on linear scale
    ax_2 : axis object
        matplotlib axis object to plot EQE on logarithmic scale
    """

    fontsize = 15
    # fontsize = 17
    plt.ion()

    fig_1, ax_1 = plt.subplots()

    if number == 0:  # number determines whether the x-axis is in wavelength or energy
        plt.xlabel('Wavelength (nm)', fontsize=fontsize, fontweight='medium')
    elif number == 1:
        plt.xlabel('Energy (eV)', fontsize=fontsize, fontweight='medium')
    elif number is None:
        plt.xlabel('Energy (eV)', fontsize=fontsize, fontweight='medium')

    if norm_num == 0:  # norm_num determines whether the y-axis is "EQE" or "Normalized EQE"
        plt.ylabel('EQE', fontsize=fontsize, fontweight='medium')
    elif norm_num == 1:
        plt.ylabel('Normalized EQE', fontsize=fontsize, fontweight='medium')
    elif norm_num is None:
        plt.ylabel('EQE', fontsize=fontsize, fontweight='medium')

    plt.rcParams['figure.facecolor'] = 'xkcd:white'
    plt.rcParams['figure.edgecolor'] = 'xkcd:white'
    plt.tick_params(labelsize=fontsize, direction='in', axis='both', which='major', length=8, width=2)
    plt.tick_params(labelsize=fontsize, direction='in', axis='both', which='minor', length=4, width=2)
    # plt.tick_params(labelsize=fontsize-2, direction='in', axis='both', which='major', length=8, width=2)
    # plt.tick_params(labelsize=fontsize-2, direction='in', axis='both', which='minor', length=4, width=2)
    plt.minorticks_on()
    plt.show()

    fig_2, ax_2 = plt.subplots()
    ax_2.set_yscale('log')  # To generate log scale axis

    if number == 0:
        plt.xlabel('Wavelength (nm)', fontsize=fontsize, fontweight='medium')
    elif number == 1:
        plt.xlabel('Energy (eV)', fontsize=fontsize, fontweight='medium')
    elif number is None:
        plt.xlabel('Energy (eV)', fontsize=fontsize, fontweight='medium')

    if norm_num == 0:
        plt.ylabel('EQE', fontsize=fontsize, fontweight='medium')
    elif norm_num == 1:
        plt.ylabel('Normalized EQE', fontsize=fontsize, fontweight='medium')
    elif norm_num is None:
        plt.ylabel('EQE', fontsize=fontsize, fontweight='medium')

    plt.rcParams['figure.facecolor'] = 'xkcd:white'
    plt.rcParams['figure.edgecolor'] = 'xkcd:white'
    plt.tick_params(labelsize=fontsize, direction='in', axis='both', which='major', length=8, width=2)
    plt.tick_params(labelsize=fontsize, direction='in', axis='both', which='minor', length=4, width=2)
    # plt.tick_params(labelsize=fontsize-2, direction='in', axis='both', which='major', length=8, width=2)
    # plt.tick_params(labelsize=fontsize-2, direction='in', axis='both', which='minor', length=4, width=2)
    plt.minorticks_on()
    plt.show()

    return ax_1, ax_2


# -----------------------------------------------------------------------------------------------------------

# Function to set up EL and EQE plot

def set_up_EL_plot():
    """Function to set up EL and eEQE plot

    Parameters
    ----------
    None
        
    Returns
    -------
    ax_1 : axis object
        matplotlib axis object to plot EQE on linear scale
    ax_2 : axis object
        matplotlib axis object to plot EQE on logarithmic scale
    """

    plt.ion()

    fig_1, ax_1 = plt.subplots()
    plt.ylabel('Red. EL (1/eV), Red. EQE (eV)', fontsize=17, fontweight='medium')
    plt.xlabel('Energy (eV)', fontsize=17, fontweight='medium')
    plt.grid(False)
    plt.rcParams['figure.facecolor'] = 'xkcd:white'
    plt.rcParams['figure.edgecolor'] = 'xkcd:white'
    plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
    plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
    plt.minorticks_on()
    plt.show()

    fig_2, ax_2 = plt.subplots()
    ax_2.set_yscale('log')
    plt.ylabel('Red. EL (1/eV), Red. EQE (eV)', fontsize=17, fontweight='medium')
    plt.xlabel('Energy (eV)', fontsize=17, fontweight='medium')
    plt.grid(False)
    plt.rcParams['figure.facecolor'] = 'xkcd:white'
    plt.rcParams['figure.edgecolor'] = 'xkcd:white'
    plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
    plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
    plt.minorticks_on()
    plt.show()

    return ax_1, ax_2

# -----------------------------------------------------------------------------------------------------------
