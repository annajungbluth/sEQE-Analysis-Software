#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 11:59:40 2018

@author: jungbluth
"""

import itertools
import math
import os
import random
import sys
import tkinter as tk
from tkinter import filedialog

import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sEQE_Analysis_template
import seaborn
import serial
import xlrd

# for the gui
from PyQt5 import QtCore, QtGui, QtWidgets
from colour import Color
from matplotlib import style
from numpy import exp, linspace, random
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):

        QtWidgets.QMainWindow.__init__(self)


        # Set up the user interface from Designer

        self.ui = sEQE_Analysis_template.Ui_MainWindow()
        self.ui.setupUi(self)


        # Tkinter

        root = tk.Tk()
        root.withdraw()


        ##### Page 1 - Calculate EQE

        self.ref_1 = []
        self.ref_2 = []
        self.ref_3 = []
        self.ref_4 = []
        self.ref_5 = []

        self.data_1 = []
        self.data_2 = []
        self.data_3 = []
        self.data_4 = []
        self.data_5 = []

        # Handle Browse Buttons

        self.ui.browseButton_1.clicked.connect(lambda: self.writeText(self.ui.textBox_1, 1))
        self.ui.browseButton_2.clicked.connect(lambda: self.writeText(self.ui.textBox_2, 2))
        self.ui.browseButton_3.clicked.connect(lambda: self.writeText(self.ui.textBox_3, 3))
        self.ui.browseButton_4.clicked.connect(lambda: self.writeText(self.ui.textBox_4, 4))
        self.ui.browseButton_5.clicked.connect(lambda: self.writeText(self.ui.textBox_5, 5))
        self.ui.browseButton_6.clicked.connect(lambda: self.writeText(self.ui.textBox_6, 6))
        self.ui.browseButton_7.clicked.connect(lambda: self.writeText(self.ui.textBox_7, 7))
        self.ui.browseButton_8.clicked.connect(lambda: self.writeText(self.ui.textBox_8, 8))
        self.ui.browseButton_9.clicked.connect(lambda: self.writeText(self.ui.textBox_9, 9))
        self.ui.browseButton_10.clicked.connect(lambda: self.writeText(self.ui.textBox_10, 10))

        # Handle Calculate Buttons

        self.ui.calculateButton_1.clicked.connect(lambda: self.pre_EQE(self.ref_1, self.data_1, self.ui.startNM_1, self.ui.stopNM_1, 1))
        self.ui.calculateButton_2.clicked.connect(lambda: self.pre_EQE(self.ref_2, self.data_2, self.ui.startNM_2, self.ui.stopNM_2, 2))
        self.ui.calculateButton_3.clicked.connect(lambda: self.pre_EQE(self.ref_3, self.data_3, self.ui.startNM_3, self.ui.stopNM_3, 3))
        self.ui.calculateButton_4.clicked.connect(lambda: self.pre_EQE(self.ref_4, self.data_4, self.ui.startNM_4, self.ui.stopNM_4, 4))
        self.ui.calculateButton_5.clicked.connect(lambda: self.pre_EQE(self.ref_5, self.data_5, self.ui.startNM_5, self.ui.stopNM_5, 5))

        # Handle Export All Data Button

        self.ui.exportButton.clicked.connect(self.export_EQE)

        # Handle Clear Plot Button

        self.ui.clearButton.clicked.connect(self.clear_plot)


        ##### Page 2 - Plot EQE  

        self.EQE_1 = []
        self.EQE_2 = []
        self.EQE_3 = []
        self.EQE_4 = []
        self.EQE_5 = []
        self.EQE_6 = []
        self.EQE_7 = []
        self.EQE_8 = []
        self.EQE_9 = []
        self.EQE_10 = []

        # Handle Import EQE Buttons

        self.ui.browseEQEButton_1.clicked.connect(lambda: self.writeText(self.ui.textBox_11, 11))
        self.ui.browseEQEButton_2.clicked.connect(lambda: self.writeText(self.ui.textBox_14, 14))
        self.ui.browseEQEButton_3.clicked.connect(lambda: self.writeText(self.ui.textBox_17, 17))
        self.ui.browseEQEButton_4.clicked.connect(lambda: self.writeText(self.ui.textBox_20, 20))
        self.ui.browseEQEButton_5.clicked.connect(lambda: self.writeText(self.ui.textBox_23, 23))
        self.ui.browseEQEButton_6.clicked.connect(lambda: self.writeText(self.ui.textBox_26, 26))
        self.ui.browseEQEButton_7.clicked.connect(lambda: self.writeText(self.ui.textBox_29, 29))
        self.ui.browseEQEButton_8.clicked.connect(lambda: self.writeText(self.ui.textBox_32, 32))
        self.ui.browseEQEButton_9.clicked.connect(lambda: self.writeText(self.ui.textBox_35, 35))
        self.ui.browseEQEButton_10.clicked.connect(lambda: self.writeText(self.ui.textBox_38, 38))

        # Handle Plot EQE Buttons

        self.ui.plotEQEButton_Wavelength.clicked.connect(lambda: self.pre_plot_EQE(0))
        self.ui.plotEQEButton_Energy.clicked.connect(lambda: self.pre_plot_EQE(1))


        ##### Page 3 - Fit CT States

        self.data_fit_1 = []
        self.data_fit_2 = []
        self.data_fit_3 = []
        self.data_fit_4 = []
        self.data_fit_5 = []

        # Handle Import EQE Buttons

        self.ui.browseFitButton_1.clicked.connect(lambda: self.writeText(self.ui.textBox_f1, 'f1'))
        self.ui.browseFitButton_2.clicked.connect(lambda: self.writeText(self.ui.textBox_f4, 'f4'))
        self.ui.browseFitButton_3.clicked.connect(lambda: self.writeText(self.ui.textBox_f7, 'f7'))
        self.ui.browseFitButton_4.clicked.connect(lambda: self.writeText(self.ui.textBox_f10, 'f10'))
        self.ui.browseFitButton_5.clicked.connect(lambda: self.writeText(self.ui.textBox_f13, 'f13'))

        # Handle Gaussian Fit Buttons

        self.ui.gaussianFit_1.clicked.connect(lambda: self.pre_fit_EQE(self.data_fit_1, self.ui.startPlot_1, self.ui.stopPlot_1, self.ui.startFit_1, self.ui.stopFit_1, self.ui.startFitPlot_1, self.ui.stopFitPlot_1, self.ui.textBox_f1, self.ui.textBox_f2, self.ui.textBox_f3, 1))
        self.ui.gaussianFit_2.clicked.connect(lambda: self.pre_fit_EQE(self.data_fit_2, self.ui.startPlot_2, self.ui.stopPlot_2, self.ui.startFit_2, self.ui.stopFit_2, self.ui.startFitPlot_2, self.ui.stopFitPlot_2, self.ui.textBox_f4, self.ui.textBox_f5, self.ui.textBox_f6, 2))
        self.ui.gaussianFit_3.clicked.connect(lambda: self.pre_fit_EQE(self.data_fit_3, self.ui.startPlot_3, self.ui.stopPlot_3, self.ui.startFit_3, self.ui.stopFit_3, self.ui.startFitPlot_3, self.ui.stopFitPlot_3, self.ui.textBox_f7, self.ui.textBox_f8, self.ui.textBox_f9, 3))
        self.ui.gaussianFit_4.clicked.connect(lambda: self.pre_fit_EQE(self.data_fit_4, self.ui.startPlot_4, self.ui.stopPlot_4, self.ui.startFit_4, self.ui.stopFit_4, self.ui.startFitPlot_4, self.ui.stopFitPlot_4, self.ui.textBox_f10, self.ui.textBox_f11, self.ui.textBox_f12, 4))
        self.ui.gaussianFit_5.clicked.connect(lambda: self.pre_fit_EQE(self.data_fit_5, self.ui.startPlot_5, self.ui.stopPlot_5, self.ui.startFit_5, self.ui.stopFit_5, self.ui.startFitPlot_5, self.ui.stopFitPlot_5, self.ui.textBox_f13, self.ui.textBox_f14, self.ui.textBox_f15, 5))


        # Handle Heat Map Buttons

        self.ui.heatButton_1.clicked.connect(lambda: self.heatMap(self.data_fit_1, self.ui.startPlot_1, self.ui.stopPlot_1, self.ui.startStart_1, self.ui.startStop_1, self.ui.stopStart_1, self.ui.stopStop_1, self.ui.textBox_f1, self.ui.textBox_f2, self.ui.textBox_f3, 1))
        self.ui.heatButton_2.clicked.connect(lambda: self.heatMap(self.data_fit_2, self.ui.startPlot_2, self.ui.stopPlot_2, self.ui.startStart_2, self.ui.startStop_2, self.ui.stopStart_2, self.ui.stopStop_2, self.ui.textBox_f4, self.ui.textBox_f5, self.ui.textBox_f6, 2))
        self.ui.heatButton_3.clicked.connect(lambda: self.heatMap(self.data_fit_3, self.ui.startPlot_3, self.ui.stopPlot_3, self.ui.startStart_3, self.ui.startStop_3, self.ui.stopStart_3, self.ui.stopStop_3, self.ui.textBox_f7, self.ui.textBox_f8, self.ui.textBox_f9, 3))
        self.ui.heatButton_4.clicked.connect(lambda: self.heatMap(self.data_fit_4, self.ui.startPlot_4, self.ui.stopPlot_4, self.ui.startStart_4, self.ui.startStop_4, self.ui.stopStart_4, self.ui.stopStop_4, self.ui.textBox_f10, self.ui.textBox_f11, self.ui.textBox_f12, 4))
        self.ui.heatButton_5.clicked.connect(lambda: self.heatMap(self.data_fit_5, self.ui.startPlot_5, self.ui.stopPlot_5, self.ui.startStart_5, self.ui.startStop_5, self.ui.stopStart_5, self.ui.stopStop_5, self.ui.textBox_f13, self.ui.textBox_f14, self.ui.textBox_f15, 5))

        self.ui.clearButton_2.clicked.connect(self.clear_EQE_fit_plot)


        ##### Page 4 - Extended Fits - MLJ Theory

        self.data_xFit_1 = []

        self.S_i = self.ui.Huang_Rhys.value()
        self.hbarw_i = self.ui.vib_Energy.value()

        # Handle Import EQE Buttons

        self.ui.browseButton_extraFit.clicked.connect(lambda: self.writeText(self.ui.textBox_xF1, 'xF1'))

        # Handle Fit Button

        self.ui.extraFitButton.clicked.connect(lambda: self.pre_fit_EQE(self.data_xFit_1, self.ui.startExtraPlot, self.ui.stopExtraPlot, self.ui.startExtraFit, self.ui.stopExtraFit, self.ui.startExtraFitPlot, self.ui.stopExtraFitPlot, self.ui.textBox_xF1, self.ui.textBox_xF2, self.ui.textBox_xF3, 'x1'))

        # Handle Heat Map Button

        self.ui.extraHeatButton.clicked.connect(lambda: self.heatMap(self.data_xFit_1, self.ui.startExtraPlot, self.ui.stopExtraPlot, self.ui.extraStartStart, self.ui.extraStartStop, self.ui.extraStopStart, self.ui.extraStopStop, self.ui.textBox_xF1, self.ui.textBox_xF2, self.ui.textBox_xF3, 'x1'))

        # Handle Clear Extra Fit Button

        self.ui.clearButton_extraFit.clicked.connect(self.clear_EQE_fit_plot)


        ##### Page 5 - Fit EL and EQE

        self.EL = []
        self.EL_EQE = []
        self.Red_EL_cal = []
        self.Red_EL_meas = []
        self.Red_EQE_meas = []

        # Handle Import EL and EQE Data Buttons

        self.ui.browseELButton_1.clicked.connect(lambda: self.writeText(self.ui.textBox_EL1, 'el1'))
        self.ui.browseELButton_2.clicked.connect(lambda: self.writeText(self.ui.textBox_EL4, 'el2'))

        # Handle EL Plot Buttons

        self.ui.plotELButton_1.clicked.connect(lambda: self.pre_plot_EL_EQE(self.EL, self.ui.startPlot_EL1, self.ui.stopPlot_EL1, 0)) # Plot EL
        self.ui.plotELButton_2.clicked.connect(lambda: self.pre_plot_EL_EQE(self.EL, self.ui.startPlot_EL2, self.ui.stopPlot_EL2, 1)) # Plot Abs
        self.ui.plotELButton_3.clicked.connect(lambda: self.pre_plot_EL_EQE(self.EL_EQE, self.ui.startPlot_EQE, self.ui.stopPlot_EQE, 2)) # Plot EQE

        # Handle Fit Buttons

        self.ui.fitButton_EL1.clicked.connect(lambda: self.pre_plot_EL_EQE(self.EL, self.ui.startPlot_EL1, self.ui.stopPlot_EL1, 0, fit=True)) # Fit EL
        self.ui.fitButton_EL2.clicked.connect(lambda: self.pre_plot_EL_EQE(self.EL, self.ui.startPlot_EL2, self.ui.stopPlot_EL2, 1, fit=True)) # Fit Abs
        self.ui.fitButton_EL3.clicked.connect(lambda: self.pre_plot_EL_EQE(self.EL_EQE, self.ui.startPlot_EQE, self.ui.stopPlot_EQE, 2, fit=True)) # Fit EQE

        # Handle Intersection Button

        #self.ui.intersectionButton.clicked.connect(lambda: self.intersection()) ##### Ability to determine intersection between EQE & EL Data. Currently removed.

        # Handle Clear EL Plot Button

        self.ui.clearButton_EL.clicked.connect(lambda: self.clear_EL_plot())


        # Import Photodiode Calibration Files

        Si_file = pd.ExcelFile("FDS100-CAL.xlsx") # The files are in the sEQE Analysis folder, just as this program
#        print(Si_file.sheet_names)
        self.Si_cal = Si_file.parse('Sheet1')

        InGaAs_file = pd.ExcelFile("FGA21-CAL.xlsx")
        self.InGaAs_cal = InGaAs_file.parse('Sheet1')


        # Data Directory

        self.data_dir = '/home/jungbluth/Desktop'   ### CHANGE THIS IF NEEDED LATER


        # Define Variables

        self.h = 6.626 * math.pow(10,-34) # [m^2 kg/s]
        self.h_2 = 4.136 * math.pow(10, -15) # [eV s]
        self.c = 2.998 * math.pow(10,8) # [m/s]
        self.q = 1.602 * math.pow(10,-19) # [C]
        self.k = 8.617 * math.pow(10, -5) # [ev/K]
        #self.T = 300 # [K]


        # To Export Calculated EQE Files

        self.export = False
        self.do_plot = True
        self.fit_plot = True
        self.do_plot_EL = True

# -----------------------------------------------------------------------------------------------------------

    #### Functions to read in file and update text box

# -----------------------------------------------------------------------------------------------------------

    def writeText(self, text_Box, textBox_no):
        self.change_dir(self.data_dir)
        file_ = filedialog.askopenfilename()
        if len(file_) != 0:
            path_, filename_ = os.path.split(file_)
#            print("Path : ", path_)
#            print("Filename : ", filename_)

            text_Box.clear() # Clear the text box in case sth has been uploaded already
            text_Box.insertPlainText(filename_) # Insert filename into text box

            # For reference files:

            if textBox_no == 1:
                self.ref_1 = pd.read_csv(file_) # Turn file into dataFrame
#                print(self.ref_1)
#                print(len(self.ref_1))

            elif textBox_no == 3:
                self.ref_2 = pd.read_csv(file_)

            elif textBox_no == 5:
                self.ref_3 = pd.read_csv(file_)

            elif textBox_no == 7:
                self.ref_4 = pd.read_csv(file_)

            elif textBox_no == 9:
                self.ref_5 = pd.read_csv(file_)

            # For data files:

            elif textBox_no == 2:
                self.data_1 = pd.read_csv(file_)

            elif textBox_no == 4:
                self.data_2 = pd.read_csv(file_)

            elif textBox_no == 6:
                self.data_3 = pd.read_csv(file_)

            elif textBox_no == 8:
                self.data_4 = pd.read_csv(file_)

            elif textBox_no == 10:
                self.data_5 = pd.read_csv(file_)

            # For EQE files:

            elif textBox_no == 11:
                self.EQE_1 = pd.read_csv(file_)

            elif textBox_no == 14:
                self.EQE_2 = pd.read_csv(file_)

            elif textBox_no == 17:
                self.EQE_3 = pd.read_csv(file_)

            elif textBox_no == 20:
                self.EQE_4 = pd.read_csv(file_)

            elif textBox_no == 23:
                self.EQE_5 = pd.read_csv(file_)

            elif textBox_no == 26:
                self.EQE_6 = pd.read_csv(file_)

            elif textBox_no == 29:
                self.EQE_7 = pd.read_csv(file_)

            elif textBox_no == 32:
                self.EQE_8 = pd.read_csv(file_)

            elif textBox_no == 35:
                self.EQE_9 = pd.read_csv(file_)

            elif textBox_no == 38:
                self.EQE_10 = pd.read_csv(file_)


            # For fit files

            elif textBox_no == 'f1':
                self.data_fit_1 = pd.read_csv(file_)

            elif textBox_no == 'f4':
                self.data_fit_2 = pd.read_csv(file_)

            elif textBox_no == 'f7':
                self.data_fit_3 = pd.read_csv(file_)

            elif textBox_no == 'f10':
                self.data_fit_4 = pd.read_csv(file_)

            elif textBox_no == 'f13':
                self.data_fit_5 = pd.read_csv(file_)

            # For El files

            elif textBox_no == 'el1':
                self.EL = pd.read_table(file_, sep= ',', index_col=0)

            elif textBox_no == 'el2':
                self.EL_EQE = pd.read_csv((file_))


            # For extra fit files

            elif textBox_no == 'xF1':
                self.data_xFit_1 = pd.read_csv(file_)


# -----------------------------------------------------------------------------------------------------------

    #### Functions to calculate EQE

# -----------------------------------------------------------------------------------------------------------

    # Function to calculate the reference power

    def calculatePower(self, ref_df, cal_df):
        cal_wave_dict = {} # Create an empty dictionary
        power = [] # Create an empty list

        for x in range(len(cal_df['Wavelength [nm]'])): # Iterate through columns of calibration file
            cal_wave_dict[cal_df['Wavelength [nm]'][x]] = cal_df['Responsivity [A/W]'][x] # Add wavelength and corresponding responsivity to dictionary


        for y in range(len(ref_df['Wavelength'])): # Iterate through columns of reference file
#            print(ref_df['Wavelength'][y])
            if ref_df['Wavelength'][y] in cal_wave_dict.keys(): # Check if reference wavelength is in calibraton file
                power.append(float(ref_df['Mean Current'][y]) / float(cal_wave_dict[ref_df['Wavelength'][y]])) # Add power to the list
            else: # If reference wavelength is not in calibration file
                resp_int = self.interpolate(ref_df['Wavelength'][y], cal_df['Wavelength [nm]'], cal_df['Responsivity [A/W]']) # Interpolate responsivity
                power.append(float(ref_df['Mean Current'][y]) / float(resp_int)) # Add power to the list

        ref_df['Power'] = power # Create new column in reference file
#        print(ref_df['Power'])

        return ref_df['Power']

# -----------------------------------------------------------------------------------------------------------

    ### Function to select data and reference file   

    def pre_EQE(self, ref_df, data_df, start, stop, range_no):
        startNM = start.value()
        stopNM = stop.value()
        if self.Ref_Data_is_valid(ref_df, data_df, startNM, stopNM, range_no):
            self.calculate_EQE(ref_df, data_df, startNM, stopNM, range_no)

# -----------------------------------------------------------------------------------------------------------

    ### Function to calculate EQE

    def calculate_EQE (self, ref_df, data_df, startNM, stopNM, range_no):

        power_dict = {}
        Wavelength = []
        Energy = []
        EQE = []
        log_EQE = []

        if 'Power' not in ref_df.columns:

            print('Calculating power values.')

            if range_no == 1:
                if self.ui.Range1_Si_button.isChecked() and not self.ui.Range1_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = self.calculatePower(ref_df, self.Si_cal)
                    except:
                        print('Please select a valid reference diode.')
                elif self.ui.Range1_InGaAs_button.isChecked() and not self.ui.Range1_Si_button.isChecked():
                    try:
                        ref_df['Power'] = self.calculatePower(ref_df, self.InGaAs_cal)
                    except:
                        print('Please select a valid reference diode.')
                else:
                    print('Please select a valid reference diode.')

            elif range_no == 2:
                if self.ui.Range2_Si_button.isChecked() and not self.ui.Range2_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = self.calculatePower(ref_df, self.Si_cal)
                    except:
                        print('Please select a valid reference diode.')
                elif self.ui.Range2_InGaAs_button.isChecked() and not self.ui.Range2_Si_button.isChecked():
                    try:
                        ref_df['Power'] = self.calculatePower(ref_df, self.InGaAs_cal)
                    except:
                        print('Please select a valid reference diode.')
                else:
                    print('Please select a valid reference diode.')

            elif range_no == 3:
                if self.ui.Range3_Si_button.isChecked() and not self.ui.Range3_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = self.calculatePower(ref_df, self.Si_cal)
                    except:
                        print('Please select a valid reference diode.')
                elif self.ui.Range3_InGaAs_button.isChecked() and not self.ui.Range3_Si_button.isChecked():
                    try:
                        ref_df['Power'] = self.calculatePower(ref_df, self.InGaAs_cal)
                    except:
                        print('Please select a valid reference diode.')
                else:
                    print('Please select a valid reference diode.')

            elif range_no == 4:
                if self.ui.Range4_Si_button.isChecked() and not self.ui.Range4_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = self.calculatePower(ref_df, self.Si_cal)
                    except:
                        print('Please select a valid reference diode.')
                elif self.ui.Range4_InGaAs_button.isChecked() and not self.ui.Range4_Si_button.isChecked():
                    try:
                        ref_df['Power'] = self.calculatePower(ref_df, self.InGaAs_cal)
                    except:
                        print('Please select a valid reference diode.')
                else:
                    print('Please select a valid reference diode.')

            elif range_no == 5:
                if self.ui.Range5_Si_button.isChecked() and not self.ui.Range5_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = self.calculatePower(ref_df, self.Si_cal)
                    except:
                        print('Please select a valid reference diode.')
                elif self.ui.Range5_InGaAs_button.isChecked() and not self.ui.Range5_Si_button.isChecked():
                    try:
                        ref_df['Power'] = self.calculatePower(ref_df, self.InGaAs_cal)
                    except:
                        print('Please select a valid reference diode.')
                else:
                    print('Please select a valid reference diode.')

        if 'Power' in ref_df.columns:  # Check if the power has been calculated already

            for x in range(len(ref_df['Wavelength'])):  # Iterate through columns of reference file
                power_dict[ref_df['Wavelength'][x]] = ref_df['Power'][x]  # Add wavelength and corresponding power to dictionary

            for y in range(len(data_df['Wavelength'])): # Iterate through columns of data file
                if startNM <= data_df['Wavelength'][y] <= stopNM: # Calculate EQE only if start <= wavelength <= stop, otherwise ignore
    #                print(data_df['Wavelength'][y])

                    if data_df['Wavelength'][y] in power_dict.keys():  # Check if data wavelength is in reference file
                        Wavelength.append(data_df['Wavelength'][y])
                        Energy_val = (self.h * self.c) / (data_df['Wavelength'][y] * math.pow(10,-9) * self.q) # Calculate energy in eV
                        Energy.append(Energy_val)
                        EQE_val = (data_df['Mean Current'][y] * Energy_val) / (power_dict[data_df['Wavelength'][y]]) # Easier way to calculate EQE
    #                    EQE_val = ((data_df['Mean Current'][y] * self.h * self.c) / (data_df['Wavelength'][y] * math.pow(10,-9) * power_dict[data_df['Wavelength'][y]] * self.q))
    #                    EQE_val = (100 * data_df['Mean Current'][y] * Energy_val) / (power_dict[data_df['Wavelength'][y]]) # *100 to turn into percent
                        EQE.append(EQE_val)
                        log_EQE.append(math.log10(EQE_val))

                    else: # If data wavelength is not in reference file
                        Wavelength.append(data_df['Wavelength'][y])
                        Energy_val = (self.h * self.c) / (data_df['Wavelength'][y] * math.pow(10,-9) * self.q)
                        Energy.append(Energy_val)
                        Power_int = self.interpolate(data_df['Wavelength'][y], ref_df['Wavelength'], ref_df['Power']) # Interpolate power
                        EQE_int = (data_df['Mean Current'][y] * Energy_val) / (Power_int)
    #                    EQE_int = ((data_df['Mean Current'][y] * self.h * self.c) / (data_df['Wavelength'][y] * math.pow(10,-9) * Power_int * self.q))
                        EQE.append(EQE_int)
                        log_EQE.append(math.log10(EQE_int))


        if len(Wavelength) == len(EQE) and len(Energy) == len(log_EQE): # Check if the lists have the same length

            if self.export: # If the "Export Data" button has been clicked
                return (Wavelength, Energy, EQE, log_EQE)

            else: # If the "Calculate EQE" button has been clicked
                if self.do_plot: # This is set to true during setup of the program
                    self.set_up_plot()
                    self.do_plot = False # Set self.do_plot to False to plot on the same graph

                label_ = self.pick_Label(range_no, startNM, stopNM)
                self.ax1.plot(Wavelength, EQE, linewidth = 3, label = label_)
#                self.ax2.plot(Wavelength, log_EQE, linewidth = 3)
                self.ax2.semilogy(Wavelength, EQE, linewidth = 3) # Equivalent to the line above but with proper log scale axes

#                self.ax1.plot(Energy, EQE, linewidth = 3)
#                self.ax2.semilogy(Energy, EQE, linewidth = 3)

                self.ax1.legend()
                plt.draw()

        else:
            print('Error Code 1: Length mismatch.')

# -----------------------------------------------------------------------------------------------------------

    #### Functions to export EQE data

# -----------------------------------------------------------------------------------------------------------

    def export_EQE(self):

        self.export = True

        ok_1 = True # Create boolean variable to use for "is_valid" function
        ok_2 = True
        ok_3 = True
        ok_4 = True
        ok_5 = True

        Wave_1 = [] # Create empty lists for all the data
        Wave_2 = []
        Wave_3 = []
        Wave_4 = []
        Wave_5 = []

        Energy_1 = []
        Energy_2 = []
        Energy_3 = []
        Energy_4 = []
        Energy_5 = []

        EQE_1 = []
        EQE_2 = []
        EQE_3 = []
        EQE_4 = []
        EQE_5 = []

        log_EQE_1 = []
        log_EQE_2 = []
        log_EQE_3 = []
        log_EQE_4 = []
        log_EQE_5 = []

        columns = ['Wavelength', 'Energy', 'EQE', 'Log_EQE'] # Columns for dataFrame

        export_file = pd.DataFrame(columns = columns) # Create empty dataFrame

        wave_inc = {} # create empty dictionary

        if self.ui.exportBox_1.isChecked(): # If the checkBox is checked
            startNM1 = self.ui.startNM_1.value() # Pick start wavelength
            stopNM1 = self.ui.stopNM_1.value() # Pick stop wavelength
            if self.Ref_Data_is_valid(self.ref_1, self.data_1, startNM1, stopNM1, 1): # Check that files are non-empty and within wavelength range
                Wave_1, Energy_1, EQE_1, log_EQE_1 = self.calculate_EQE(self.ref_1, self.data_1, startNM1, stopNM1, 1) # Extract data
                export_1 = pd.DataFrame({'Wavelength': Wave_1, 'Energy': Energy_1, 'EQE': EQE_1, 'Log_EQE': log_EQE_1}) # Create dataFrame with EQE data
                wave_inc['1'] = Wave_1[0] # Add the first wavelength value to the wave_inc list
            else:
                ok_1 = False # Set variable to False if calculation is invalid

        if self.ui.exportBox_2.isChecked():
            startNM2 = self.ui.startNM_2.value()
            stopNM2 = self.ui.stopNM_2.value()
            if self.Ref_Data_is_valid(self.ref_2, self.data_2, startNM2, stopNM2, 2):
                Wave_2, Energy_2, EQE_2, log_EQE_2 = self.calculate_EQE(self.ref_2, self.data_2, startNM2, stopNM2, 2)
                export_2 = pd.DataFrame({'Wavelength': Wave_2, 'Energy': Energy_2, 'EQE': EQE_2, 'Log_EQE': log_EQE_2})
                wave_inc['2'] = Wave_2[0]
            else:
                ok_2 = False

        if self.ui.exportBox_3.isChecked():
            startNM3 = self.ui.startNM_3.value()
            stopNM3 = self.ui.stopNM_3.value()
            if self.Ref_Data_is_valid(self.ref_3, self.data_3, startNM3, stopNM3, 3):
                Wave_3, Energy_3, EQE_3, log_EQE_3 = self.calculate_EQE(self.ref_3, self.data_3, startNM3, stopNM3, 3)
                export_3 = pd.DataFrame({'Wavelength': Wave_3, 'Energy': Energy_3, 'EQE': EQE_3, 'Log_EQE': log_EQE_3})
                wave_inc['3'] = Wave_3[0]
            else:
                ok_3 = False

        if self.ui.exportBox_4.isChecked():
            startNM4 = self.ui.startNM_4.value()
            stopNM4 = self.ui.stopNM_4.value()
            if self.Ref_Data_is_valid(self.ref_4, self.data_4, startNM4, stopNM4, 4):
                Wave_4, Energy_4, EQE_4, log_EQE_4 = self.calculate_EQE(self.ref_4, self.data_4, startNM4, stopNM4, 4)
                export_4 = pd.DataFrame({'Wavelength': Wave_4, 'Energy': Energy_4, 'EQE': EQE_4, 'Log_EQE': log_EQE_4})
                wave_inc['4'] = Wave_4[0]
            else:
                ok_4 = False

        if self.ui.exportBox_5.isChecked():
            startNM5 = self.ui.startNM_5.value()
            stopNM5 = self.ui.stopNM_5.value()
            if self.Ref_Data_is_valid(self.ref_5, self.data_5, startNM5, stopNM5, 5):
                Wave_5, Energy_5, EQE_5, log_EQE_5 = self.calculate_EQE(self.ref_5, self.data_5, startNM5, stopNM5, 5)
                export_5 = pd.DataFrame({'Wavelength': Wave_5, 'Energy': Energy_5, 'EQE': EQE_5, 'Log_EQE': log_EQE_5})
                wave_inc['5'] = Wave_5[0]
            else:
                ok_5 = False


        if ok_1 and ok_2 and ok_3 and ok_4 and ok_5: # Check if all operations are ok or if fields are left empty

            for x in range(len(wave_inc.keys())): # Iterate through wave_inc list
#                print(x)
#                print(wave_inc)
                minimum = min(wave_inc, key = wave_inc.get) # Find key corresponding to minimum value
#                print(minimum)
                if minimum == '1': # Append correct dataFrame in order of decending wavelength
                    export_file = export_file.append(export_1, ignore_index = True)
#                    print(export_file)
                elif minimum == '2':
                    export_file = export_file.append(export_2, ignore_index = True)
                elif minimum == '3':
                    export_file = export_file.append(export_3, ignore_index = True)
                elif minimum == '4':
                    export_file = export_file.append(export_4, ignore_index = True)
                elif minimum == '5':
                    export_file = export_file.append(export_5, ignore_index = True)

                del wave_inc[minimum] # Delete recently appended value
#            print(export_file)

            EQE_file = filedialog.asksaveasfilename() # Prompt the user to pick a folder & name to save data to
            export_path, export_filename = os.path.split(EQE_file)
            if len(export_path) != 0: # Check if the user actually selected a path
                self.change_dir(export_path) # Change the working directory
                export_file.to_csv(export_filename) # Save data to csv
                print('Saving data to: %s' % str(EQE_file))
                self.change_dir(self.data_dir) # Change the directory back

        self.export = False

# -----------------------------------------------------------------------------------------------------------

    #### Functions to plot EQE data

# -----------------------------------------------------------------------------------------------------------

    ### Function to select EQE data  

    def pre_plot_EQE(self, number):

        ok_EQE_1 = True
        ok_EQE_2 = True
        ok_EQE_3 = True
        ok_EQE_4 = True
        ok_EQE_5 = True
        ok_EQE_6 = True
        ok_EQE_7 = True
        ok_EQE_8 = True
        ok_EQE_9 = True
        ok_EQE_10 = True

        if self.ui.normalizeBox.isChecked():
            norm_num = 1
        else:
            norm_num = 0

        self.set_up_EQE_plot(number, norm_num)

        if self.ui.plotBox_1.isChecked():
            ok_EQE_1 = self.plot_EQE(self.EQE_1, self.ui.startEQE_1, self.ui.stopEQE_1, self.ui.textBox_11, self.ui.textBox_12, self.ui.textBox_13, 1, number)

        if self.ui.plotBox_2.isChecked():
            ok_EQE_2 = self.plot_EQE(self.EQE_2, self.ui.startEQE_2, self.ui.stopEQE_2, self.ui.textBox_14, self.ui.textBox_15, self.ui.textBox_16, 2, number)

        if self.ui.plotBox_3.isChecked():
            ok_EQE_3 = self.plot_EQE(self.EQE_3, self.ui.startEQE_3, self.ui.stopEQE_3, self.ui.textBox_17, self.ui.textBox_18, self.ui.textBox_19, 3, number)

        if self.ui.plotBox_4.isChecked():
            ok_EQE_4 = self.plot_EQE(self.EQE_4, self.ui.startEQE_4, self.ui.stopEQE_4, self.ui.textBox_20, self.ui.textBox_21, self.ui.textBox_22, 4, number)

        if self.ui.plotBox_5.isChecked():
            ok_EQE_5 = self.plot_EQE(self.EQE_5, self.ui.startEQE_5, self.ui.stopEQE_5, self.ui.textBox_23, self.ui.textBox_24, self.ui.textBox_25, 5, number)

        if self.ui.plotBox_6.isChecked():
            ok_EQE_6 = self.plot_EQE(self.EQE_6, self.ui.startEQE_6, self.ui.stopEQE_6, self.ui.textBox_26, self.ui.textBox_27, self.ui.textBox_28, 6, number)

        if self.ui.plotBox_7.isChecked():
            ok_EQE_7 = self.plot_EQE(self.EQE_7, self.ui.startEQE_7, self.ui.stopEQE_7, self.ui.textBox_29, self.ui.textBox_30, self.ui.textBox_31, 7, number)

        if self.ui.plotBox_8.isChecked():
            ok_EQE_8 = self.plot_EQE(self.EQE_8, self.ui.startEQE_8, self.ui.stopEQE_8, self.ui.textBox_32, self.ui.textBox_33, self.ui.textBox_33, 8, number)

        if self.ui.plotBox_9.isChecked():
            ok_EQE_9 = self.plot_EQE(self.EQE_9, self.ui.startEQE_9, self.ui.stopEQE_9, self.ui.textBox_34, self.ui.textBox_35, self.ui.textBox_36, 9, number)

        if self.ui.plotBox_10.isChecked():
            ok_EQE_10 = self.plot_EQE(self.EQE_10, self.ui.startEQE_10, self.ui.stopEQE_10, self.ui.textBox_37, self.ui.textBox_38, self.ui.textBox_40, 10, number)


        if ok_EQE_1 and ok_EQE_2 and ok_EQE_3 and ok_EQE_4 and ok_EQE_5 and ok_EQE_6 and ok_EQE_7 and ok_EQE_8 and ok_EQE_9 and ok_EQE_10:
            self.ax3.legend()
            self.ax4.legend()
            plt.show()
        else:
            plt.close()
            plt.close()

# -----------------------------------------------------------------------------------------------------------

    def plot_EQE(self, eqe_df, startNM, stopNM, filename_Box, label_Box, color_Box, file_no, number):

       startNM = startNM.value() # Pick start wavelength
       stopNM = stopNM.value() # Pick stop wavelength

       if self.EQE_is_valid(eqe_df, startNM, stopNM, file_no): # Check that files are non-empty and within wavelength range

           if self.ui.normalizeBox.isChecked():
               normNM = self.ui.normalizeNM.value()
               if self.Normalization_is_valid(eqe_df, normNM, file_no):
                   wave, energy, eqe, log_eqe = self.normalize_EQE(eqe_df, startNM, stopNM, normNM)
               else:
                   return False
           else:
               wave, energy, eqe, log_eqe = self.compile_EQE(eqe_df, startNM, stopNM, 0)

           label_ = self.pick_EQE_Label(label_Box, filename_Box)
           color_ = self.pick_EQE_Color(color_Box, file_no)

           if number == 0:
               self.ax3.plot(wave, eqe, linewidth = 3, label = label_, color = color_)
               self.ax4.semilogy(wave, eqe, linewidth = 3, label = label_, color = color_)
           elif number == 1:
               self.ax3.plot(energy, eqe, linewidth = 3, label = label_, color = color_)
               self.ax4.semilogy(energy, eqe, linewidth = 3, label = label_, color = color_)

           return True

       else:
           return False


# -----------------------------------------------------------------------------------------------------------

    #### Functions to fit EQE data

# -----------------------------------------------------------------------------------------------------------

    ### Function to select EQE data  

    def pre_fit_EQE(self, eqe_df, startE, stopE, startFit, stopFit, startPlotFit, stopPlotFit, filename_Box, label_Box, color_Box, file_no):

        # Change this if you want to plot on the same graph 
        if self.fit_plot:
            self.set_up_EQE_fit_plot() # This sets up a plot with x-axis "Energy" and y-label of "EQE" / "Log EQE"
            self.fit_plot = False

        ok_plot_Fit = self.plot_fit_EQE(eqe_df, startE, stopE, startFit, stopFit, startPlotFit, stopPlotFit, filename_Box, label_Box, color_Box, file_no)

        if ok_plot_Fit:
            self.ax6.legend()
            plt.show()
#        else:
#            plt.close()

# -----------------------------------------------------------------------------------------------------------

    ### Function to plot EQE data and gaussian fit 

    def plot_fit_EQE(self, eqe_df, startE, stopE, startFit, stopFit, startPlotFit, stopPlotFit, filename_Box, label_Box, color_Box, file_no):

       include_Disorder = False

       startE = startE.value() # Pick start energy
       stopE = stopE.value() # Pick stop energy

       startFit = startFit.value() # Pick start fit energy
       stopFit = stopFit.value() # Pick stop fit energy

       startPlotFit = startPlotFit.value()
       stopPlotFit = stopPlotFit.value()

       if self.Fit_is_valid(eqe_df, startE, stopE, startFit, stopFit, file_no): # Check that files are non-empty and within energy range
           wave, energy, eqe, log_eqe = self.compile_EQE(eqe_df, startE, stopE, 1) # Compile EQE file
           wave_fit, energy_fit, eqe_fit, log_eqe_fit = self.compile_EQE(eqe_df, startFit, stopFit, 1) # Compile fit range of EQE file

           label_ = self.pick_EQE_Label(label_Box, filename_Box)
           color_ = self.pick_EQE_Color(color_Box, file_no)


           if (str(file_no)).isnumeric(): # For Marcus theory fitting

               if file_no == 1:
                   self.T_CT = self.ui.Temperature_1.value()
                   if self.ui.static_Disorder_1.isChecked():
                       include_Disorder = True
                       self.sig = self.ui.Disorder_1.value()

               elif file_no ==2:
                   self.T_CT = self.ui.Temperature_2.value()
                   if self.ui.static_Disorder_2.isChecked():
                       include_Disorder = True
                       self.sig = self.ui.Disorder_2.value()

               elif file_no == 3:
                   self.T_CT = self.ui.Temperature_3.value()
                   if self.ui.static_Disorder_3.isChecked():
                       include_Disorder = True
                       self.sig = self.ui.Disorder_3.value()

               elif file_no ==4:
                   self.T_CT = self.ui.Temperature_4.value()
                   if self.ui.static_Disorder_4.isChecked():
                       include_Disorder = True
                       self.sig = self.ui.Disorder_4.value()

               elif file_no ==5:
                   self.T_CT = self.ui.Temperature_5.value()
                   if self.ui.static_Disorder_5.isChecked():
                       include_Disorder = True
                       self.sig = self.ui.Disorder_5.value()

               try:
                   x_gaussian = linspace(startPlotFit, stopPlotFit,50)  # Create more x values to perform the fit on. This is useful to plot more of the gaussian.
                   y_gaussian = []

                   if include_Disorder:
                       best_vals, covar = curve_fit(self.gaussian_disorder, energy_fit, eqe_fit) # , p0 = init_values
                       for value in x_gaussian:
                           y_gaussian.append(self.gaussian_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                   else:
                       best_vals, covar = curve_fit(self.gaussian, energy_fit, eqe_fit)
                       for value in x_gaussian:
                           y_gaussian.append(self.gaussian(value, best_vals[0], best_vals[1], best_vals[2]))

                   print('-'*80)
                   print('Temperature [T] (K) : ', self.T_CT)
                   print('Oscillator Strength [f] (eV**2) : ', format(best_vals[0], '.6f'), '+/-', format(math.sqrt(covar[0,0]), '.6f'))
                   print('Reorganization Energy [l] (eV) : ',  format(best_vals[1], '.6f'), '+/-', format(math.sqrt(covar[1,1]), '.6f'))
                   print('CT State Energy [ECT] (eV) : ',  format(best_vals[2], '.6f'), '+/-', format(math.sqrt(covar[2,2]), '.6f'))
                   print('-'*80)

                   self.ax6.semilogy(energy, eqe, linewidth = 3, label = label_, color = color_)
                   plt.draw()
                   if include_Disorder:
                       self.ax6.plot(x_gaussian, y_gaussian, linewidth=2, label='Gaussian Fit + Disorder', color='#000000', linestyle='--')
                   else:
                       self.ax6.plot(x_gaussian, y_gaussian, linewidth=2, label='Gaussian Fit', color='#000000', linestyle='--')
                   plt.draw()

               except:
                   print('Optimal parameters not found.')
                   return False

               return True

           elif file_no == 'x1': # For extended fitting with MLJ Theory

               self.S_i = self.ui.Huang_Rhys.value()
               self.hbarw_i = self.ui.vib_Energy.value()
               self.T_x = self.ui.extra_Temperature.value()

               if self.ui.extra_static_Disorder.isChecked():
                   include_Disorder = True
                   self.sig_x = self.ui.extra_Disorder.value()

               try:
                   x_MLJ_theory = linspace(startPlotFit, stopPlotFit, 50)
                   y_MLJ_theory = []

                   if include_Disorder:
                       best_vals, covar = curve_fit(self.MLJ_gaussian_disorder, energy_fit, eqe_fit)
                       for value in x_MLJ_theory:
                           y_MLJ_theory.append(self.MLJ_gaussian_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                   else:
                       best_vals, covar = curve_fit(self.MLJ_gaussian, energy_fit, eqe_fit)
                       for value in x_MLJ_theory:
                           y_MLJ_theory.append(self.MLJ_gaussian(value, best_vals[0], best_vals[1], best_vals[2]))

                   print('-'*80)
                   print('Temperature [T] (K) : ', self.T_x)
                   print('Oscillator Strength [f] (eV**2) : ', format(best_vals[0], '.6f'), '+/-', format(math.sqrt(covar[0,0]), '.6f'))
                   print('Reorganization Energy [l] (eV) : ',  format(best_vals[1], '.6f'), '+/-', format(math.sqrt(covar[1,1]), '.6f'))
                   print('CT State Energy [ECT] (eV) : ',  format(best_vals[2], '.6f'), '+/-', format(math.sqrt(covar[2,2]), '.6f'))
                   print('-'*80)

                   self.ax6.semilogy(energy, eqe, linewidth=3, label=label_, color=color_)
                   plt.draw()
                   if include_Disorder:
                       self.ax6.plot(x_MLJ_theory, y_MLJ_theory, linewidth=2, label='MLJ Fit + Disorder', color='#000000', linestyle='--')
                   else:
                       self.ax6.plot(x_MLJ_theory, y_MLJ_theory, linewidth=2, label='MLJ Fit', color='#000000', linestyle='--')
                   plt.draw()

               except:
                   print('Optimal parameters not found.')
                   return False

               return True

       else:
           return False


# -----------------------------------------------------------------------------------------------------------

    ### Function to generate heat map of fitting values 

    def heatMap(self, eqe_df, startE, stopE, startStartE, startStopE, stopStartE, stopStopE, filename_Box, label_Box, color_Box, file_no):

        include_Disorder = False

        startE = startE.value() # Pick start energy value
        stopE = stopE.value() # Pick stop energy value
        startStartE = startStartE.value()
        startStopE = startStopE.value()
        stopStartE = stopStartE.value()
        stopStopE = stopStopE.value()

        startStart_ok = self.StartStop_is_valid(startStartE, startStopE) # Check that start energy is lower than stop energy
        startStop_ok = self.StartStop_is_valid(startStopE, stopStartE)
        stopStop_ok = self.StartStop_is_valid(stopStartE, stopStopE)

        if startStart_ok and startStop_ok and stopStop_ok: # If all operations are valid, proceed with heat map calculations

            startEnergies = [] # Create empty lists for start and stop energies
            stopEnergies = []

            start_df = [] # Create empty lists for fit data collection
            stop_df = []
            f_df = []
            l_df = []
            Ect_df = []

            x = float(startStartE)
            y = float(stopStartE)

            while x <= float(startStopE): # As long as start range start value is smaller than start range stop value
                startEnergies.append(round(x, 3)) # Round to three digits
                x += 0.005 # Adjust this number to change resolution

            while y <= float(stopStopE):
                stopEnergies.append(round(y, 3))
                y += 0.005

            for start in startEnergies: # Iterate through start energies
                for stop in stopEnergies: # Iterature through stop energies

                    if self.Fit_is_valid(eqe_df, startE, stopE, start, stop, file_no): # If the fit is valid, perform the heat map calculations

                        wave, energy, eqe, log_eqe = self.compile_EQE(eqe_df, startE, stopE, 1)
                        wave_fit, energy_fit, eqe_fit, log_eqe_fit = self.compile_EQE(eqe_df, start, stop, 1)

                        if (str(file_no)).isnumeric():

                            if file_no == 1:
                                self.T_CT = self.ui.Temperature_1.value()
                                if self.ui.static_Disorder_1.isChecked():
                                    include_Disorder = True
                                    self.sig = self.ui.Disorder_1.value()

                            elif file_no == 2:
                                self.T_CT = self.ui.Temperature_2.value()
                                if self.ui.static_Disorder_2.isChecked():
                                    include_Disorder = True
                                    self.sig = self.ui.Disorder_2.value()

                            elif file_no == 3:
                                self.T_CT = self.ui.Temperature_3.value()
                                if self.ui.static_Disorder_3.isChecked():
                                    include_Disorder = True
                                    self.sig = self.ui.Disorder_3.value()

                            elif file_no == 4:
                                self.T_CT = self.ui.Temperature_4.value()
                                if self.ui.static_Disorder_4.isChecked():
                                    include_Disorder = True
                                    self.sig = self.ui.Disorder_4.value()

                            elif file_no == 5:
                                self.T_CT = self.ui.Temperature_5.value()
                                if self.ui.static_Disorder_5.isChecked():
                                    include_Disorder = True
                                    self.sig = self.ui.Disorder_5.value()

                            try:
                                if include_Disorder:
                                    best_vals, covar = curve_fit(self.gaussian_disorder, energy_fit,eqe_fit)  # , p0 = init_values
                                else:
                                    best_vals, covar = curve_fit(self.gaussian, energy_fit, eqe_fit)

                                start_df.append(start)
                                stop_df.append(stop)
                                f_df.append(best_vals[0])
                                l_df.append(best_vals[1])
                                Ect_df.append(best_vals[2])

                            except:
                               print('Optimal parameters not found.')

                        elif file_no == 'x1': # For extended fitting with MLJ Theory

                            self.S_i = self.ui.Huang_Rhys.value()
                            self.hbarw_i = self.ui.vib_Energy.value()
                            self.T_x = self.ui.extra_Temperature.value()

                            if self.ui.extra_static_Disorder.isChecked():
                                include_Disorder = True
                                self.sig_x = self.ui.extra_Disorder.value()

                            try:
                                if include_Disorder:
                                    best_vals, covar = curve_fit(self.MLJ_gaussian_disorder, energy_fit, eqe_fit)
                                else:
                                    best_vals, covar = curve_fit(self.MLJ_gaussian, energy_fit, eqe_fit)

                                start_df.append(start)
                                stop_df.append(stop)
                                f_df.append(best_vals[0])
                                l_df.append(best_vals[1])
                                Ect_df.append(best_vals[2])

                            except:
                                print('Optimal parameters not found.')


            parameter_df = pd.DataFrame({'Start': start_df, 'Stop': stop_df, 'f': f_df, 'l': l_df, 'Ect': Ect_df}) # Create a dataFrame with all results

            print('-'*80)

            if file_no == 'x1':
                print('Temperature [T] (K) : ', self.T_x)
            else:
                print('Temperature [T] (K) : ', self.T_CT)

            print('Average Oscillator Strength [f] (eV**2) : ', format(parameter_df['f'].mean(), '.6f'), ', Stdev : ', format(parameter_df['f'].std(), '.6f')) # Determine the average value and standard deviation
            print('Average Reorganization Energy [l] (eV) : ', format(parameter_df['l'].mean(), '.6f'), ', Stdev : ', format(parameter_df['l'].std(), '.6f'))
            print('Average CT State Energy [ECT] (eV) : ', format(parameter_df['Ect'].mean(), '.6f'), ', Stdev : ', format(parameter_df['Ect'].std(), '.6f'))
            print('-'*80)


            f_df = parameter_df.pivot('Stop', 'Start', 'f') # Pivot the dataFrame: x-value = Stop, y-value = Start, value = f
            l_df = parameter_df.pivot('Stop', 'Start', 'l')
            Ect_df = parameter_df.pivot('Stop', 'Start', 'Ect')

            #seaborn.set(font_scale=1.2)

            if len(f_df) != 0:
                plt.ion()
                plt.figure() # Create a new figure
                self.ax5 = seaborn.heatmap(f_df, xticklabels=3, yticklabels=3) # Create the heat map
                plt.xlabel('Inital Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.ylabel('Final Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.title('Oscillator Strength ($eV^2$)', fontsize=17, fontweight='medium')
                cbar = self.ax5.collections[0].colorbar
                cbar.ax.tick_params(labelsize=15)
                plt.yticks(rotation = 360)
                plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
                #plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
                plt.show()
            else:
                print('No fitting oscillator strengths determined.')

            if len(l_df) != 0:
                plt.ion()
                plt.figure()
                self.ax6 = seaborn.heatmap(l_df, xticklabels=3, yticklabels=3)
                plt.xlabel('Inital Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.ylabel('Final Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.title('Reorganization Energy (eV)', fontsize=17, fontweight='medium')
                cbar = self.ax6.collections[0].colorbar
                cbar.ax.tick_params(labelsize=15)
                plt.yticks(rotation = 360)
                plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
                #plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
                plt.show()
            else:
                print('No fitting reorganization energies determined.')

            if len(l_df) != 0:
                plt.ion()
                plt.figure()
                self.ax7 = seaborn.heatmap(Ect_df, xticklabels=3, yticklabels=3)
                plt.xlabel('Inital Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.ylabel('Final Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.title('CT State Energy (eV)', fontsize=17, fontweight='medium')
                cbar = self.ax7.collections[0].colorbar
                cbar.ax.tick_params(labelsize=15)
                plt.yticks(rotation = 360)
                plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
                #plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
                plt.show()
            else:
                print('No fitting CT state energies determined.')


# -----------------------------------------------------------------------------------------------------------

    ### Gaussian fitting function

    def gaussian(self, E, f, l, Ect):
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """
        return (f/(E * math.sqrt(4 * math.pi * l * self.T_CT * self.k))) * exp(-(Ect + l - E)**2 / (4 * l * self.k * self.T_CT))

    ### Gaussian fitting function including disorder

    def gaussian_disorder(self, E, f, l, Ect):
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """
        return (f/(E * math.sqrt(4 * math.pi * l * self.T_CT * self.k + 2 * self.sig**2))) * exp(-(Ect + l - E)**2 / (4 * l * self.k * self.T_CT + 2 * self.sig**2))

    # -----------------------------------------------------------------------------------------------------------

    ### MLJ Theory fitting function

    def MLJ_gaussian(self, E, f, l_o, Ect): # Double check if this equation is correct
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """
        EQE = 0
        for n in range(0, 6):
            EQE_n = (f/(E * math.sqrt(4 * math.pi * l_o * self.T_x * self.k))) \
                    * (math.exp(-self.S_i) * self.S_i**n / math.factorial(n)) \
                    * exp(-(Ect + l_o - E + n*self.hbarw_i)**2 \
                    / (4 * l_o * self.k * self.T_x))
            EQE += EQE_n
        return EQE

    ### MLJ Theory fitting function including disorder

    def MLJ_gaussian_disorder(self, E, f, l_o, Ect): # Double check if this equation is correct
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """
        EQE = 0
        for n in range(0, 6):
            EQE_n = (f/(E * math.sqrt(4 * math.pi * l_o * self.T_x * self.k + 2 * self.sig_x**2))) \
                    * (math.exp(-self.S_i) * self.S_i**n / math.factorial(n)) \
                    * exp(-(Ect + l_o - E + n*self.hbarw_i)**2 \
                    / (4 * l_o * self.k * self.T_x + 2 * self.sig_x**2))
            EQE += EQE_n
        return EQE

# -----------------------------------------------------------------------------------------------------------

    #### Functions to fit EL and EQE data

# -----------------------------------------------------------------------------------------------------------

    ### Function to scale and reduce EL and EQE data

    def pre_plot_EL_EQE(self, data_df, startE, stopE, data_no, fit = False):

        self.T_EL = self.ui.EL_Temperature.value()

        startE = startE.value()
        stopE = stopE.value()

        if self.do_plot_EL:
            self.set_up_EL_plot()
            self.do_plot_EL = False

        if data_no < 2: # EL data

            Energy = []
            scaleFactor = self.ui.scalePlot.value()

            if len(data_df) != 0: # Check that file is non-empty

                for y in range(len(data_df['Wavelength'])): # Calculate energy values

                    Energy_val = (self.h * self.c) / (data_df['Wavelength'][y] * math.pow(10, -9) * self.q)
                    Energy.append(Energy_val)

                data_df['Energy'] = Energy

                if self.Data_is_valid(data_df, startE, stopE) and self.StartStop_is_valid(startE, stopE):

                    EL_wave, EL_energy, EL_signal = self.compile_EL(data_df, startE, stopE, 1)
                    red_EL_scaled = [EL_signal[x] / (scaleFactor * EL_energy[x])  for x in range(len(EL_signal))] # Divide by energy to reduce

                    if data_no == 0: # EL Data

                        if not fit:
                            label_ = self.pick_EQE_Label(self.ui.textBox_EL2, self.ui.textBox_EL1)
                            #color_ = self.pick_EQE_Color(self.ui.textBox_EL3, 100) # not currently used
                            color_ = '#1f77b4' # Blue
                            self.plot_EL_EQE(EL_energy, red_EL_scaled, label_, color_)

                        elif fit:
                            self.fit_EL_EQE(EL_energy, red_EL_scaled, self.ui.startFit_EL1, self.ui.stopFit_EL1, 0)

                    elif data_no == 1: # Abs Data

                        bb_dict = self.bb_spectrum(EL_energy)
                        EL_abs = [EL_signal[x] / bb_dict[EL_energy[x]] for x in range(len(EL_energy))]
                        red_EL_abs_scaled = [EL_abs[x] / (scaleFactor * EL_energy[x]) for x in range(len(EL_abs))]

                        if not fit:
                            #label_ = self.pick_EQE_Label(self.ui.textBox_EL2, self.ui.textBox_EL1)
                            #color_ = self.pick_EQE_Color(self.ui.textBox_EL3, 100) # not currently used
                            label_ = '$\mathregular{Red. EQE_{cal}}$'
                            color_ = '#ff7716' # Orange
                            self.plot_EL_EQE(EL_energy, red_EL_abs_scaled, label_, color_)

                        elif fit:
                            self.fit_EL_EQE(EL_energy, red_EL_abs_scaled, self.ui.startFit_EL2, self.ui.stopFit_EL2, 1)

            else:
                print('Please select a valid EL file.')

        elif data_no == 2: # EQE Data

            if self.Data_is_valid(data_df, startE, stopE) and self.StartStop_is_valid(startE, stopE):

                self.Red_EQE_meas = pd.DataFrame() # For determining the intersect between abs and emission
                EQE_wave, EQE_energy, EQE, EQE_log = self.compile_EQE(data_df, startE, stopE, 1)
                red_EQE = [EQE[x] * EQE_energy[x] for x in range(len(EQE))]

                if not fit:
                    label_ = self.pick_EQE_Label(self.ui.textBox_EL5, self.ui.textBox_EL4)
                    #color_ = self.pick_EQE_Color(self.ui.textBox_EL6, 100)
                    color_ = '#000000' # Black
                    self.plot_EL_EQE(EQE_energy, red_EQE, label_, color_)

                elif fit:
                    self.fit_EL_EQE(EQE_energy, red_EQE, self.ui.startFit_EQE, self.ui.stopFit_EQE, 1)

    ### Function to plot any data

    def plot_EL_EQE(self, x, y, label_, color_):

        self.ax5.plot(x, y, linewidth = 3, label = label_, color = color_)
        plt.legend()
        plt.draw()

# -----------------------------------------------------------------------------------------------------------

    ### Function to fit reduced EL and EQE data

    def fit_EL_EQE(self, energy, y, startE, stopE, data_no):

        include_Disorder = False

        startFit = startE.value()
        stopFit = stopE.value()

        df = pd.DataFrame()
        df['Energy'] = energy

        if self.Data_is_valid(df, startFit, stopFit):

            energy_fit, y_fit = self.compile_Data(energy, y, startFit, stopFit)

            diff = stopFit - startFit  # Find difference between start and stop fit energy
            x_gaussian = linspace(startFit, stopFit + 0.5*diff, 50)  # Create more x values to perform the fit on. This is useful to plot more of the gaussian.
            y_gaussian = []

            if self.ui.static_Disorder_EL.isChecked():
                include_Disorder = True
                self.sig_EL = self.ui.EL_Disorder.value()

            try:

                if self.ui.Gaussian_EL_EQE.isChecked(): # Marcus Theory Fitting

                    if data_no == 0: # EL Data
                        if include_Disorder: # Include Disorder in Fit
                            #y_fit_smooth = savgol_filter(y_fit, 51, 3) # In case you need to smooth the data
                            #y_fit_smooth = [x for x in y_fit_smooth]
                            #log_y_fit = [math.log(x) for x in y_fit]
                            #self.plot_EL_EQE(energy_fit, y_fit_smooth, 'Smoothed Data', '#330000')
                            best_vals, covar = curve_fit(self.gaussian_EL_disorder, energy_fit, y_fit, p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(self.gaussian_EL_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                        else: # Without Disorder
                            best_vals, covar = curve_fit(self.gaussian_EL, energy_fit, y_fit, p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(self.gaussian_EL(value, best_vals[0], best_vals[1], best_vals[2]))

                    elif data_no == 1: # EQE / Abs Data
                        x_gaussian = linspace(startFit, stopFit + 2 * diff, 50)
                        if include_Disorder:
                            best_vals, covar = curve_fit(self.gaussian_EQE_disorder, energy_fit, y_fit, p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(self.gaussian_EQE_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                        else:
                            best_vals, covar = curve_fit(self.gaussian_EQE, energy_fit, y_fit, p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(self.gaussian_EQE(value, best_vals[0], best_vals[1], best_vals[2]))

                    print('-'*80)
                    print('Temperature [T] (K): ', self.T_EL)
                    print('Oscillator Strength [f] (eV**2) : ', format(best_vals[0], '.6f'), '+/-', format(math.sqrt(covar[0, 0]), '.6f'))
                    print('Reorganization Energy [l] (eV) : ', format(best_vals[1], '.6f'), '+/-', format(math.sqrt(covar[1, 1]), '.6f'))
                    print('CT State Energy [ECT] (eV) : ', format(best_vals[2], '.6f'), '+/-', format(math.sqrt(covar[2, 2]), '.6f'))
                    print('-'*80)

                    if include_Disorder:
                        self.ax5.plot(x_gaussian, y_gaussian, linewidth=2, label='Gaussian Fit + Disorder', color='#000000', linestyle='--')
                    else:
                        self.ax5.plot(x_gaussian, y_gaussian, linewidth=2, label='Gaussian Fit', color='#000000', linestyle='--')
                    plt.legend()
                    plt.draw()

                elif self.ui.MLJ_Gaussian_EL_EQE.isChecked():  # MLJ Theory Fitting

                    self.S_i_EL = self.ui.EL_Huang_Rhys.value()
                    self.hbarw_i_EL = self.ui.EL_vib_Energy.value()

                    x_gaussian = linspace(1.18, stopFit + 0.5 * diff, 50)

                    if data_no == 0:  # EL Data
                        if include_Disorder:
                            best_vals, covar = curve_fit(self.MLJ_gaussian_EL_disorder, energy_fit, y_fit, p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(
                                    self.MLJ_gaussian_EL_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                        else:
                            best_vals, covar = curve_fit(self.MLJ_gaussian_EL, energy_fit, y_fit, p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(self.MLJ_gaussian_EL(value, best_vals[0], best_vals[1], best_vals[2]))

                    elif data_no == 1:  # EQE / Abs Data
                        x_gaussian = linspace(startFit, stopFit + 2 * diff, 50)
                        if include_Disorder:
                            best_vals, covar = curve_fit(self.MLJ_gaussian_EQE_disorder, energy_fit, y_fit, p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(
                                    self.MLJ_gaussian_EQE_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                        else:
                            best_vals, covar = curve_fit(self.MLJ_gaussian_EQE, energy_fit, y_fit, p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(self.MLJ_gaussian_EQE(value, best_vals[0], best_vals[1], best_vals[2]))

                    print('-' * 80)
                    print('Temperature [T] (K): ', self.T_EL)
                    print('Oscillator Strength [f] (eV**2) : ', format(best_vals[0], '.6f'), '+/-', format(math.sqrt(covar[0, 0]), '.6f'))
                    print('Reorganization Energy [l] (eV) : ', format(best_vals[1], '.6f'), '+/-', format(math.sqrt(covar[1, 1]), '.6f'))
                    print('CT State Energy [ECT] (eV) : ', format(best_vals[2], '.6f'), '+/-',format(math.sqrt(covar[2, 2]), '.6f'))
                    print('-' * 80)

                    if include_Disorder:
                        self.ax5.plot(x_gaussian, y_gaussian, linewidth=2, label='MLJ Fit + Disorder', color='#000000', linestyle='--')
                    else:
                        self.ax5.plot(x_gaussian, y_gaussian, linewidth=2, label='MLJ Fit', color='#000000',linestyle='--')
                    plt.legend()
                    plt.draw()

            except:
                print('Optimal parameters not found.')

# -----------------------------------------------------------------------------------------------------------

    ### Gaussian fitting function for reduced EL data

    def gaussian_EL(self, E, f, l, Ect):
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: reduced EL value
        """

        return (f / (math.sqrt(4 * math.pi * l * self.T_EL * self.k))) * exp(-(Ect - l - E) ** 2 / (4 * l * self.k * self.T_EL))

    ### Gaussian fitting function inlcuding disorder for reduced EL data

    def gaussian_EL_disorder(self, E, f, l, Ect):
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: reduced EL value
        """

        return (f / (math.sqrt(4 * math.pi * l * self.T_EL * self.k + 2 * self.sig_EL**2))) * exp(-(Ect - l - E) ** 2 / (4 * l * self.k * self.T_EL + 2 * self.sig_EL**2))

    # -----------------------------------------------------------------------------------------------------------

    ### MLJ Theory fitting function for reduced EL data

    def MLJ_gaussian_EL(self, E, f, l_o, Ect): # Double check if this equation is correct
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EL value
        """
        EL = 0
        for n in range(0, 6):
            EL_n = (f/(math.sqrt(4 * math.pi * l_o * self.T_EL * self.k))) \
                    * (math.exp(-self.S_i_EL) * self.S_i_EL**n / math.factorial(n)) \
                    * exp(-(Ect - E - l_o - n*self.hbarw_i_EL)**2 \
                    / (4 * l_o * self.k * self.T_EL))
            EL += EL_n
        return EL

    ### MLJ Theory fitting function including disorder for reduced EL data

    def MLJ_gaussian_EL_disorder(self, E, f, l_o, Ect): # Double check if this equation is correct
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """

        EL = 0
        for n in range(0, 6):
            EL_n = (f/(math.sqrt(4 * math.pi * l_o * self.T_EL * self.k + 2 * self.sig_EL**2))) \
                    * (math.exp(-self.S_i_EL) * self.S_i_EL**n / math.factorial(n)) \
                    * exp(-(Ect - E - l_o - n*self.hbarw_i_EL)**2 \
                    / (4 * l_o * self.k * self.T_EL + 2 * self.sig_EL**2))
            EL += EL_n
        return EL


# -----------------------------------------------------------------------------------------------------------

    ### Gaussian fitting function for reduced EQE data

    def gaussian_EQE(self, E, f, l, Ect):
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: reduced EL value
        """

        return (f / (math.sqrt(4 * math.pi * l * self.T_EL * self.k))) * exp(-(Ect + l - E) ** 2 / (4 * l * self.k * self.T_EL))

    ### Gaussian fitting function for reduced EQE data

    def gaussian_EQE_disorder(self, E, f, l, Ect):
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: reduced EL value
        """

        return (f / (math.sqrt(4 * math.pi * l * self.T_EL * self.k + 2 * self.sig_EL**2))) * exp(-(Ect + l - E) ** 2 / (4 * l * self.k * self.T_EL + 2 * self.sig_EL**2))

    # -----------------------------------------------------------------------------------------------------------

    ### MLJ Theory fitting function

    def MLJ_gaussian_EQE(self, E, f, l_o, Ect): # Double check if this equation is correct
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """
        EQE = 0
        for n in range(0, 6):
            EQE_n = (f/(math.sqrt(4 * math.pi * l_o * self.T_EL * self.k))) \
                    * (math.exp(-self.S_i_EL) * self.S_i_EL**n / math.factorial(n)) \
                    * exp(-(Ect - E + l_o + n*self.hbarw_i_EL)**2 \
                    / (4 * l_o * self.k * self.T_EL))
            EQE += EQE_n
        return EQE

    ### MLJ Theory fitting function including disorder

    def MLJ_gaussian_EQE_disorder(self, E, f, l_o, Ect): # Double check if this equation is correct
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """

        EQE = 0
        for n in range(0, 6):
            EQE_n = (f/(math.sqrt(4 * math.pi * l_o * self.T_EL * self.k + 2 * self.sig_EL**2))) \
                    * (math.exp(-self.S_i_EL) * self.S_i_EL**n / math.factorial(n)) \
                    * exp(-(Ect - E + l_o + n*self.hbarw_i_EL)**2 \
                    / (4 * l_o * self.k * self.T_EL + 2 * self.sig_EL**2))
            EQE += EQE_n
        return EQE


# -----------------------------------------------------------------------------------------------------------

    #### Functions to compile & normalize EQE data

# -----------------------------------------------------------------------------------------------------------

    ### Function to compile EQE data  

    def compile_EQE(self, eqe_df, start, stop, number):

        Wavelength = []
        Energy = []
        EQE = []
        log_EQE = []

        if number == 0: # If a wavelength range is given
            startNM = start
            stopNM = stop

        elif number == 1: # If an energy range is given
            startNM = (self.h * self.c * math.pow(10,9)) / (stop * self.q) # The start wavelength corresponds to the high energy stop value
            stopNM = (self.h * self.c * math.pow(10,9)) / (start * self.q) # The stop wavelength corresponds to the low energy start value


        for y in range(len(eqe_df['Wavelength'])): # Iterate through columns of EQE file
            if startNM <= eqe_df['Wavelength'][y] <= stopNM: # Compile EQE only if start <= wavelength <= stop, otherwise ignore
                Wavelength.append(eqe_df['Wavelength'][y])
                Energy.append(eqe_df['Energy'][y])
                EQE.append(eqe_df['EQE'][y]) #  * eqe_df['Energy'][y] - Eventually this might need to be added for fitting
                log_EQE.append(eqe_df['Log_EQE'][y]) #  (math.log10(eqe_df['EQE'][y] * eqe_df['Energy'][y]) - Eventually this might need to be added for fitting

        if len(Wavelength) == len(EQE) and len(Energy) == len(log_EQE): # Check that the lengths are the same
            return Wavelength, Energy, EQE, log_EQE

        else:
            print('Error Code 1: Length mismatch.')

# -----------------------------------------------------------------------------------------------------------

    ### Function to normalize EQE data  

    def normalize_EQE(self, eqe_df, startNM, stopNM, normNM):

        Wavelength = []
        Energy = []
        EQE = []
        log_EQE = []


        norm_EQE = self.interpolate(normNM, eqe_df['Wavelength'], eqe_df['EQE'])
        norm_log_EQE = self.interpolate(normNM, eqe_df['Wavelength'], eqe_df['Log_EQE'])

        print(norm_EQE)
        print(norm_log_EQE)

        for y in range(len(eqe_df['Wavelength'])): # Iterate through columns of EQE file
            if startNM <= eqe_df['Wavelength'][y] <= stopNM: # Compile EQE only if start <= wavelength <= stop, otherwise ignore
                Wavelength.append(eqe_df['Wavelength'][y])
                Energy.append(eqe_df['Energy'][y])
                EQE.append(eqe_df['EQE'][y] / norm_EQE)
                log_EQE.append(eqe_df['Log_EQE'][y] / norm_log_EQE)

        if len(Wavelength) == len(EQE) and len(Energy) == len(log_EQE): # Check that the lengths are the same
            return Wavelength, Energy, EQE, log_EQE

        else:
            print('Error Code 1: Length mismatch.')


# -----------------------------------------------------------------------------------------------------------

    #### Functions to compile & scale EL data

# -----------------------------------------------------------------------------------------------------------

    ### Function to compile EL data

    def compile_EL(self, el_df, start, stop, number):

        Wavelength = []
        Energy = []
        EL = []

        if number == 0:  # If a wavelength range is given
            startNM = start
            stopNM = stop

        elif number == 1:  # If an energy range is given
            startNM = (self.h * self.c * math.pow(10, 9)) / (stop * self.q)  # The start wavelength corresponds to the high energy stop value
            stopNM = (self.h * self.c * math.pow(10, 9)) / (start * self.q)  # The stop wavelength corresponds to the low energy start value

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

    def compile_Data(self, energy, y, startE, stopE):

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

    #### Validity functions

# -----------------------------------------------------------------------------------------------------------

    ### Function to check if reference & data files are non-empty and within wavelength range              

    def Ref_Data_is_valid(self, ref_df, data_df, startNM, stopNM, range_no):

        if len(ref_df) !=0 and len(data_df) != 0:

            if startNM >= data_df['Wavelength'][0] and \
               startNM >= ref_df['Wavelength'][0] and \
               stopNM <= data_df['Wavelength'][int(len(data_df['Wavelength']))-1] and \
               stopNM <= ref_df['Wavelength'][int(len(ref_df['Wavelength']))-1]:
                return True

            elif startNM < data_df['Wavelength'][0] or startNM < ref_df['Wavelength'][0]:
                print('Please select a valid start wavelength for Range %s.' % str(range_no))
                return False

            elif stopNM > data_df['Wavelength'][int(len(data_df['Wavelength']))-1] or stopNM > ref_df['Wavelength'][int(len(ref_df['Wavelength']))-1]:
                print('Please select a valid stop wavelength for Range %s.' % str(range_no))
                return False

            else:
                print('Please select a valid wavelength range for Range %s.' % str(range_no))
                return False


        elif len(ref_df) == 0 and len(data_df) != 0: # If reference file is empty / hasn't been selected
            print('Please import a valid reference file for Range %s.' % str(range_no))
            return False

        elif len(ref_df) != 0 and len(data_df) == 0: # If data file is empty / hasn't been selected
            print('Please import a valid data file for Range %s.' % str(range_no))
            return False

        else: # If reference and data files are empty / haven't been selected
            print('Please import valid reference and data files for Range %s.' % str(range_no))
            return False

# -----------------------------------------------------------------------------------------------------------

    ### Function to check if EQE files are non-empty and within wavelength range              

    def EQE_is_valid(self, eqe_df, startNM, stopNM, EQE_no):

        if len(eqe_df) != 0:

            if startNM >= eqe_df['Wavelength'][0] and \
               stopNM <= eqe_df['Wavelength'][int(len(eqe_df['Wavelength']))-1]:

                   return True

            elif startNM < eqe_df['Wavelength'][0] and stopNM <= eqe_df['Wavelength'][int(len(eqe_df['Wavelength']))-1]:
                print('Please select a valid start wavelength for EQE File %s.' % str(EQE_no))
                return False

            elif startNM >= eqe_df['Wavelength'][0] and stopNM > eqe_df['Wavelength'][int(len(eqe_df['Wavelength']))-1]:
                print('Please select a valid stop wavelength for EQE File %s.' % str(EQE_no))
                return False

            else:
                print('Please select a valid wavelength range for EQE File %s.' % str(EQE_no))
                return False

        else: # If EQE file is empty / hasn't been selected
            print('Please import a valid file for EQE File %s.' % str(EQE_no))
            return False

# -----------------------------------------------------------------------------------------------------------

    ### Function to check if any data file is non-empty and within energy range

    # def Data_is_valid(self, df, startE, stopE):
    #
    #     if len(df) != 0:
    #
    #         if startE >= min(df['Energy']) and stopE <= max(df['Energy']):
    #             return True
    #
    #         elif startE < min(df['Energy']) and stopE <= max(df['Energy']):
    #             print('Please select a valid start energy.')
    #             return False
    #
    #         elif startE >= min(df['Energy']) and stopE > max(df['Energy']):
    #             print('Please select a valid stop energy.')
    #             return False
    #
    #         else:
    #             print('Please select a valid energy range.')
    #             return False
    #
    #     else: # If file is empty / hasn't been selected
    #         print('Please import a valid file.')
    #         return False

    def Data_is_valid(self, df, startE, stopE):

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

        else: # If file is empty / hasn't been selected
            print('Please import a valid file.')
            return False

# -----------------------------------------------------------------------------------------------------------

    ### Function to check if normalization wavelength is within wavelength range   

    def Normalization_is_valid(self, eqe_df, normNM, EQE_no):

        if len(eqe_df) != 0:

            min_wave = int(eqe_df['Wavelength'].min())
            max_wave = int(eqe_df['Wavelength'].max())

            if min_wave <= int(normNM) <= max_wave:
                return True

            else:
                print('Please select a valid normalization wavelength for EQE file %s.' % str(EQE_no))
                return False

# -----------------------------------------------------------------------------------------------------------

    ### Function to check if EQE files are non-empty and within energy range              

    def Fit_is_valid(self, eqe_df, startE, stopE, startFitE, stopFitE, EQE_no):

        if len(eqe_df) != 0:

            if startE <= eqe_df['Energy'][0] and \
               startFitE <= eqe_df['Energy'][0] and \
               stopE >= eqe_df['Energy'][int(len(eqe_df['Energy']))-1] and \
               stopFitE >= eqe_df['Energy'][int(len(eqe_df['Energy']))-1]:
                return True

            elif startE > eqe_df['Energy'][0] and stopE >= eqe_df['Energy'][int(len(eqe_df['Energy']))-1]:
                print('Please select a valid start energy for EQE File %s.' % str(EQE_no))
                return False

            elif startE <= eqe_df['Energy'][0] and stopE <eqe_df['Energy'][int(len(eqe_df['Energy']))-1]:
                print('Please select a valid stop energy for EQE File %s.' % str(EQE_no))
                return False

            else:
                print('Please select a valid energy range for EQE File %s.' % str(EQE_no))
                return False

        else: # If EQE file is empty / hasn't been selected
            print('Please import a valid file for EQE File %s.' % str(EQE_no))
            return False

# -----------------------------------------------------------------------------------------------------------

    ### Function to check if start fit energy is smaller than stop fit energy

    def StartStop_is_valid(self, startE, stopE):
        if startE < stopE:
            return True
        else:
            print('Please select valid start and stop energies.')
            return False


# -----------------------------------------------------------------------------------------------------------

    #### Funtions to set up plot

# -----------------------------------------------------------------------------------------------------------

    ### Function to set up "Calculate EQE" plot

    def set_up_plot(self):

#        style.use('ggplot')
        fig1 = plt.figure()

        self.ax1 = fig1.add_subplot(2,1,1)
        plt.ylabel('EQE', fontsize=17, fontweight='medium')
        plt.grid(False)
#        plt.box()
        plt.rcParams['figure.facecolor']='xkcd:white'
        plt.rcParams['figure.edgecolor']='xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
        plt.minorticks_on()

        self.ax2 = fig1.add_subplot(2,1,2)
        self.ax2.set_yscale('log') # To generate log scale axis
        plt.xlabel('Wavelength (nm)', fontsize=17, fontweight='medium')
        plt.ylabel('EQE', fontsize=17, fontweight='medium')
        plt.grid(False)
#        plt.box()
        plt.rcParams['figure.facecolor']='xkcd:white'
        plt.rcParams['figure.edgecolor']='xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
        plt.minorticks_on()
        plt.show()

# -----------------------------------------------------------------------------------------------------------

    ### Function to set up EQE plot

    def set_up_EQE_plot(self, number, norm_num):

#        style.use('ggplot')

        plt.ion()

        self.fig3, self.ax3 = plt.subplots()

        if number == 0: # number determines whether the x-axis is in wavelength or energy
            plt.xlabel('Wavelength (nm)', fontsize=17, fontweight='medium')
        elif number == 1:
            plt.xlabel('Energy (eV)', fontsize=17, fontweight='medium')

        if norm_num == 0: # norm_num determines whether the y-axis is "EQE" or "Normalized EQE"
            plt.ylabel('EQE', fontsize=17, fontweight='medium')
        elif norm_num == 1:
            plt.ylabel('Normalized EQE', fontsize=17, fontweight='medium')

        plt.rcParams['figure.facecolor']='xkcd:white'
        plt.rcParams['figure.edgecolor']='xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
        plt.minorticks_on()
        plt.show()


        self.fig4, self.ax4 = plt.subplots()
        self.ax4.set_yscale('log')  # To generate log scale axis

        if number == 0:
            plt.xlabel('Wavelength (nm)', fontsize=17, fontweight='medium')
        elif number == 1:
            plt.xlabel('Energy (eV)', fontsize=17, fontweight='medium')

        if norm_num == 0:
            plt.ylabel('EQE', fontsize=17, fontweight='medium')
        elif norm_num == 1:
            plt.ylabel('Normalized EQE', fontsize=17, fontweight='medium')
        elif norm_num == 2:
            plt.ylabel('Reduced EQE', fontsize=17, fontweight='medium')

        plt.rcParams['figure.facecolor'] = 'xkcd:white'
        plt.rcParams['figure.edgecolor'] = 'xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
        plt.minorticks_on()
        plt.show()


# -----------------------------------------------------------------------------------------------------------

    ### Function to set up EQE fit plot

    def set_up_EQE_fit_plot(self):

        plt.ion()

        self.fig6, self.ax6 = plt.subplots()
        self.ax6.set_yscale('log')  # To generate log scale axis

        plt.xlabel('Energy (eV)', fontsize=17, fontweight='medium')
        plt.ylabel('EQE', fontsize=17, fontweight='medium')

        plt.rcParams['figure.facecolor'] = 'xkcd:white'
        plt.rcParams['figure.edgecolor'] = 'xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
        plt.minorticks_on()
        plt.show()

# -----------------------------------------------------------------------------------------------------------

### Function to set up EL and EQE plot

    def set_up_EL_plot(self):

        plt.ion()

        self.fig5, self.ax5 = plt.subplots()
        self.ax5.set_yscale('log')

        plt.ylabel('Red. EL (1/eV), Red. EQE (eV)', fontsize=17, fontweight='medium')
        plt.xlabel('Energy (eV)', fontsize=17, fontweight='medium')
        plt.grid(False)
        plt.rcParams['figure.facecolor'] = 'xkcd:white'
        plt.rcParams['figure.edgecolor'] = 'xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
        plt.minorticks_on()
        plt.show()

# -----------------------------------------------------------------------------------------------------------

    ### Function to clear plot        

    def clear_plot(self):

        plt.close() # Close the current plot
        self.set_up_plot() # Set up a new plot, this is preferred over plt.clf() in case the plot window was closed
#        self.do_plot = True

# -----------------------------------------------------------------------------------------------------------

    ### Function to clear fit plot        

    def clear_EQE_fit_plot(self):

        plt.close() # Close the current plot
        self.set_up_EQE_fit_plot() # Set up a new plot, this is preferred over plt.clf() in case the plot window was closed

# -----------------------------------------------------------------------------------------------------------

    ### Function to clear EL plot

    def clear_EL_plot(self):

        plt.close()
        self.set_up_EL_plot()

# -----------------------------------------------------------------------------------------------------------

    #### Plotting helper functions

# -----------------------------------------------------------------------------------------------------------

    ### Function to pick "Calculate EQE" plot label

    def pick_Label(self, range_no, startNM, stopNM):

        label = str("Range") + str(range_no) + "_" + str(int(startNM)) + "nm" + "-" + str(int(stopNM)) + "nm" # Range1_360nm_800nm
        return label

# -----------------------------------------------------------------------------------------------------------

    ### Function to pick EQE plot label

    def pick_EQE_Label(self, label_Box, filename_Box):

        label = label_Box.toPlainText()
        filename = filename_Box.toPlainText()

        if len(label) != 0:
            return label
        else: # We don't need to check that there is a filename, as the "pick_label" function is only called after checking "EQE_is_valid"
            return filename

# -----------------------------------------------------------------------------------------------------------

    ### Function to pick plot color

    def pick_EQE_Color(self, colour_Box, file_no):

        colour = colour_Box.toPlainText()
        colour = colour.replace(" ", "") # If the user inputs "Sky Blue" instead of "SkyBlue" etc.

        color_list = ['#808080', '#000000', '#FF0000', '#800000', '#FFFF00', '#808000', '#00FF00', '#008000', '#00FFFF', '#008080', '#0000FF', '#000080', '#FF00FF', '#800080']
#        random_colour = random.choice(color_list)

        color_comp = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']
        color_choice = [random.choice(color_comp) for j in range(6)]
        random_colour = '#'+''.join(color_choice)

#        random_colour = '#'+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
#        print(random_colour)

        if len(colour) != 0:
            if self.is_Colour(colour):
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

    ### Function to check if input is a color

    def is_Colour(self, colour):

        try:
            Color(colour)
            return True
        except:
            return False

# -----------------------------------------------------------------------------------------------------------

    #### Other helper functions

# -----------------------------------------------------------------------------------------------------------

    ### Function to interpolate values    

    def interpolate(self, num, x, y):

        f = interp1d(x, y)
        return f(num)

# -----------------------------------------------------------------------------------------------------------

    ### Function to change working directory

    def change_dir(self, directory):

        os.chdir(directory)

# -----------------------------------------------------------------------------------------------------------

    #### Additional function

# -----------------------------------------------------------------------------------------------------------

    def bb_spectrum(self, E_list):

        phi_bb_dict = {}

        for energy in E_list:
            phi_bb = 2 * math.pi * (energy)**2 * math.exp(-1 * energy / (self.k * self.T_EL)) / (self.h_2**3 * self.c**2) # -1) - without -1 as an approximation
            phi_bb_dict[energy] = phi_bb

        return phi_bb_dict #[s/kg m^4]

# -----------------------------------------------------------------------------------------------------------


def main():

  app = QtWidgets.QApplication(sys.argv)
  monoUI = MainWindow()
  monoUI.show()
  sys.exit(app.exec_())

if __name__ == '__main__':
    try:
        main()
    except:
        app = QtWidgets.QApplication(sys.argv)
