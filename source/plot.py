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
    """
    :param ax1: plt axis object
    :param ax2: plt axis object
    :param x: x data [list or array]
    :param y: y data [list or array]
    :param label_: label to plot [string]
    :param color_: color to plot [string]
    :return:
    """
    ax1.plot(x, y, linewidth=3, label=label_, color=color_)
    ax1.legend()

    ax2.plot(x, y, linewidth=3, label=label_, color=color_)
    ax2.legend()

    plt.draw()


# -----------------------------------------------------------------------------------------------------------

# Function to set up "Calculate EQE" plot

def set_up_plot(flag='Wavelength'):
    """
    :param flag: Flag to specify if wavelength or energy are plotted [string]
    :return: ax_1: plt axis object to plot EQE on linear scale
             ax_2: plt axis objects to plot EQE on log scale
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
    """
    :param: number: number indicating whether wavelength or energy are plotted [int]
            number = 0 => plot wavelength
            number = 1 => plot energy
            number = None => plot energy
    :param: norm_num: number indicating whether raw, normalized, or reduced data are plotted [int]
            norm_num = 0 => plot raw EQE
            norm_num = 1 => plot normalized EQE
            norm_num = None => plot raw EQE
    :return: ax_1: plt axis object to plot EQE on linear scale
             ax_2: plt axis objects to plot EQE on log scale
    """

    # with mpl.rc_context({'axes.linewidth': 2}):
    #
    #     fontsize = 15
    #     # fontsize = 17
    #     plt.ion()
    #
    #     fig_1, ax_1 = plt.subplots(figsize=(6, 6))
    #
    #     if number == 0:  # number determines whether the x-axis is in wavelength or energy
    #         plt.xlabel('Wavelength / nm', fontsize=fontsize, fontweight='medium')
    #     elif number == 1:
    #         plt.xlabel('Energy / eV', fontsize=fontsize, fontweight='medium')
    #     elif number is None:
    #         plt.xlabel('Energy / eV', fontsize=fontsize, fontweight='medium')
    #
    #     if norm_num == 0:  # norm_num determines whether the y-axis is "EQE" or "Normalized EQE"
    #         plt.ylabel('EQE / a.u.', fontsize=fontsize, fontweight='medium')
    #     elif norm_num == 1:
    #         plt.ylabel('Normalized EQE / a.u.', fontsize=fontsize, fontweight='medium')
    #     elif norm_num is None:
    #         plt.ylabel('EQE / a.u.', fontsize=fontsize, fontweight='medium')
    #
    #     plt.xticks(np.arange(0.5, 2.5, 0.1))
    #     plt.xlim(1.1, 1.7)
    #
    #     plt.rcParams['figure.facecolor'] = 'xkcd:white'
    #     plt.rcParams['figure.edgecolor'] = 'xkcd:white'
    #     plt.tick_params(labelsize=fontsize, direction='in', axis='both', which='major', length=6, width=1, top=True, right=False, left=True)
    #     plt.tick_params(labelsize=fontsize, direction='in', axis='both', which='minor', length=3, width=1, left=True, bottom=True, top=True)
    #     plt.minorticks_on()
    #     plt.show()
    #
    #     fig_2, ax_2 = plt.subplots(figsize=(6, 6))
    #     ax_2.set_yscale('log')  # To generate log scale axis
    #
    #     if number == 0:
    #         plt.xlabel('Wavelength / nm', fontsize=fontsize, fontweight='medium')
    #     elif number == 1:
    #         plt.xlabel('Energy / eV', fontsize=fontsize, fontweight='medium')
    #     elif number is None:
    #         plt.xlabel('Energy / eV', fontsize=fontsize, fontweight='medium')
    #
    #     if norm_num == 0:
    #         plt.ylabel('EQE / a.u.', fontsize=fontsize, fontweight='medium')
    #     elif norm_num == 1:
    #         plt.ylabel('Normalized EQE / a.u.', fontsize=fontsize, fontweight='medium')
    #     elif norm_num is None:
    #         plt.ylabel('EQE / a.u.', fontsize=fontsize, fontweight='medium')
    #
    #     plt.ylim(5*10E-6, 10E0)
    #     plt.xticks(np.arange(0.5, 2.5, 0.1))
    #     plt.xlim(1.1, 1.7)
    #
    #     plt.rcParams['figure.facecolor'] = 'xkcd:white'
    #     plt.rcParams['figure.edgecolor'] = 'xkcd:white'
    #     plt.tick_params(labelsize=fontsize, direction='in', axis='both', which='major', length=6, width=1, top=True, right=False, left=True)
    #     plt.tick_params(labelsize=fontsize, direction='in', axis='both', which='minor', length=3, width=1, left=True, bottom=True, top=True)
    #     plt.minorticks_on()
    #     plt.show()

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
    """
    :return: ax_1: plt axis object to plot EQE on linear scale
             ax_2: plt axis objects to plot EQE on log scale
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
