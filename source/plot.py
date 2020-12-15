import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------------------------------------

### Function to set up generic EQE plot

def set_up_EQE_plot():

    plt.ion()

    fig_1, ax_1 = plt.subplots()
    plt.xlabel('Energy (eV)', fontsize=17, fontweight='medium')
    plt.ylabel('EQE', fontsize=17, fontweight='medium')
    plt.rcParams['figure.facecolor'] = 'xkcd:white'
    plt.rcParams['figure.edgecolor'] = 'xkcd:white'
    plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
    plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
    plt.minorticks_on()
    plt.show()

    fig_2, ax_2 = plt.subplots()
    ax_2.set_yscale('log')  # To generate log scale axis
    plt.xlabel('Energy (eV)', fontsize=17, fontweight='medium')
    plt.ylabel('EQE', fontsize=17, fontweight='medium')
    plt.rcParams['figure.facecolor'] = 'xkcd:white'
    plt.rcParams['figure.edgecolor'] = 'xkcd:white'
    plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
    plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
    plt.minorticks_on()
    plt.show()

    print('ok')

    return ax_1, ax_2

# -----------------------------------------------------------------------------------------------------------