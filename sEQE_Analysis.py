#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 11:59:40 2018

@author: jungbluth
"""

import math
import os
import sys
import tkinter as tk
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn
# for the gui
from PyQt5 import QtWidgets
from numpy import exp, linspace
from scipy.optimize import curve_fit
from tqdm import tqdm

import sEQE_Analysis_template
from source.compilation import compile_EQE, compile_EL, compile_Data
from source.electroluminescence import bb_spectrum
from source.gaussian import calculate_gaussian_absorption
from source.normalization import normalize_EQE
from source.plot import plot, set_up_plot, set_up_EQE_plot, set_up_EL_plot
from source.reference_correction import calculate_Power
from source.utils import interpolate, R_squared
from source.utils_plot import is_Colour, pick_EQE_Color, pick_EQE_Label, pick_Label
from source.validity import Ref_Data_is_valid, EQE_is_valid, Data_is_valid, Normalization_is_valid, Fit_is_valid, \
    StartStop_is_valid


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):

        QtWidgets.QMainWindow.__init__(self)

        # Set up the user interface from Designer
        self.ui = sEQE_Analysis_template.Ui_MainWindow()
        self.ui.setupUi(self)

        # Tkinter
        root = tk.Tk()
        root.withdraw()

        # Page 1 - Calculate EQE

        self.ref_1 = []
        self.ref_2 = []
        self.ref_3 = []
        self.ref_4 = []
        self.ref_5 = []
        self.ref_6 = []

        self.data_1 = []
        self.data_2 = []
        self.data_3 = []
        self.data_4 = []
        self.data_5 = []
        self.data_6 = []

        # Handle Browse Buttons
        self.ui.browseButton_1.clicked.connect(lambda: self.writeText(self.ui.textBox_p1_1, 1))
        self.ui.browseButton_2.clicked.connect(lambda: self.writeText(self.ui.textBox_p1_2, 2))
        self.ui.browseButton_3.clicked.connect(lambda: self.writeText(self.ui.textBox_p1_3, 3))
        self.ui.browseButton_4.clicked.connect(lambda: self.writeText(self.ui.textBox_p1_4, 4))
        self.ui.browseButton_5.clicked.connect(lambda: self.writeText(self.ui.textBox_p1_5, 5))
        self.ui.browseButton_6.clicked.connect(lambda: self.writeText(self.ui.textBox_p1_6, 6))
        self.ui.browseButton_7.clicked.connect(lambda: self.writeText(self.ui.textBox_p1_7, 7))
        self.ui.browseButton_8.clicked.connect(lambda: self.writeText(self.ui.textBox_p1_8, 8))
        self.ui.browseButton_9.clicked.connect(lambda: self.writeText(self.ui.textBox_p1_9, 9))
        self.ui.browseButton_10.clicked.connect(lambda: self.writeText(self.ui.textBox_p1_10, 10))
        self.ui.browseButton_11.clicked.connect(lambda: self.writeText(self.ui.textBox_p1_11, 11))
        self.ui.browseButton_12.clicked.connect(lambda: self.writeText(self.ui.textBox_p1_12, 12))

        # Handle Calculate Buttons
        self.ui.calculateButton_1.clicked.connect(
            lambda: self.pre_EQE(self.ref_1, self.data_1, self.ui.startNM_1, self.ui.stopNM_1, 1))
        self.ui.calculateButton_2.clicked.connect(
            lambda: self.pre_EQE(self.ref_2, self.data_2, self.ui.startNM_2, self.ui.stopNM_2, 2))
        self.ui.calculateButton_3.clicked.connect(
            lambda: self.pre_EQE(self.ref_3, self.data_3, self.ui.startNM_3, self.ui.stopNM_3, 3))
        self.ui.calculateButton_4.clicked.connect(
            lambda: self.pre_EQE(self.ref_4, self.data_4, self.ui.startNM_4, self.ui.stopNM_4, 4))
        self.ui.calculateButton_5.clicked.connect(
            lambda: self.pre_EQE(self.ref_5, self.data_5, self.ui.startNM_5, self.ui.stopNM_5, 5))
        self.ui.calculateButton_6.clicked.connect(
            lambda: self.pre_EQE(self.ref_6, self.data_6, self.ui.startNM_6, self.ui.stopNM_6, 6))

        # Handle Export All Data Button
        self.ui.exportButton.clicked.connect(self.export_EQE)

        # Handle Clear Plot Button
        self.ui.clearButton.clicked.connect(self.clear_plot)

        # Page 2 - Plot EQE

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
        self.ui.browseEQEButton_1.clicked.connect(lambda: self.writeText(self.ui.textBox_p2_1, 'p1'))
        self.ui.browseEQEButton_2.clicked.connect(lambda: self.writeText(self.ui.textBox_p2_4, 'p4'))
        self.ui.browseEQEButton_3.clicked.connect(lambda: self.writeText(self.ui.textBox_p2_7, 'p7'))
        self.ui.browseEQEButton_4.clicked.connect(lambda: self.writeText(self.ui.textBox_p2_10, 'p10'))
        self.ui.browseEQEButton_5.clicked.connect(lambda: self.writeText(self.ui.textBox_p2_13, 'p13'))
        self.ui.browseEQEButton_6.clicked.connect(lambda: self.writeText(self.ui.textBox_p2_16, 'p16'))
        self.ui.browseEQEButton_7.clicked.connect(lambda: self.writeText(self.ui.textBox_p2_19, 'p19'))
        self.ui.browseEQEButton_8.clicked.connect(lambda: self.writeText(self.ui.textBox_p2_22, 'p22'))
        self.ui.browseEQEButton_9.clicked.connect(lambda: self.writeText(self.ui.textBox_p2_25, 'p25'))
        self.ui.browseEQEButton_10.clicked.connect(lambda: self.writeText(self.ui.textBox_p2_28, 'p28'))

        # Handle Plot EQE Buttons
        self.ui.plotEQEButton_Wavelength.clicked.connect(lambda: self.pre_plot_EQE(0))
        self.ui.plotEQEButton_Energy.clicked.connect(lambda: self.pre_plot_EQE(1))

        # Page 3 - Fit EQE (Marcus Theory)

        self.data_fit_1 = []
        self.data_fit_2 = []

        # Handle Import EQE Buttons
        self.ui.browseFitButton_1.clicked.connect(lambda: self.writeText(self.ui.textBox_f1, 'f1'))
        self.ui.browseFitButton_2.clicked.connect(lambda: self.writeText(self.ui.textBox_f4, 'f4'))

        # Handle Gaussian Fit Buttons
        self.ui.gaussianFit_1.clicked.connect(
            lambda: self.pre_fit_EQE(self.data_fit_1, self.ui.startPlot_1, self.ui.stopPlot_1, self.ui.startFit_1,
                                     self.ui.stopFit_1, self.ui.startFitPlot_1, self.ui.stopFitPlot_1,
                                     self.ui.textBox_f1, self.ui.textBox_f2, self.ui.textBox_f3, 1))
        self.ui.gaussianFit_2.clicked.connect(
            lambda: self.pre_fit_EQE(self.data_fit_2, self.ui.startPlot_2, self.ui.stopPlot_2, self.ui.startFit_2,
                                     self.ui.stopFit_2, self.ui.startFitPlot_2, self.ui.stopFitPlot_2,
                                     self.ui.textBox_f4, self.ui.textBox_f5, self.ui.textBox_f6, 2))

        # Handle Heat Map Buttons
        self.ui.heatButton_1.clicked.connect(
            lambda: self.heatMap(self.data_fit_1, self.ui.startStart_1, self.ui.startStop_1,
                                 self.ui.stopStart_1, self.ui.stopStop_1, self.ui.textBox_f1,
                                 self.ui.textBox_f2, self.ui.textBox_f3, 1))
        self.ui.heatButton_2.clicked.connect(
            lambda: self.heatMap(self.data_fit_2, self.ui.startStart_2, self.ui.startStop_2,
                                 self.ui.stopStart_2, self.ui.stopStop_2, self.ui.textBox_f4,
                                 self.ui.textBox_f5, self.ui.textBox_f6, 2))

        self.ui.clearButton_2.clicked.connect(self.clear_EQE_plot)

        # Double Fits

        self.data_double = []

        # Handle Import Data Button
        self.ui.browseDoubleFitButton.clicked.connect(lambda: self.writeText(self.ui.textBox_dF1, 'double1'))

        # Handle Double Fit Button
        self.ui.doubleFitButton.clicked.connect(lambda: self.pre_double_fit())

        # Page 4 - Fit EQE (MLJ Theory)

        self.data_xFit_1 = []

        self.S_i = self.ui.Huang_Rhys.value()
        self.hbarw_i = self.ui.vib_Energy.value()

        # Handle Import EQE Buttons
        self.ui.browseButton_extraFit.clicked.connect(lambda: self.writeText(self.ui.textBox_xF1, 'xF1'))

        # Handle Fit Button
        self.ui.extraFitButton.clicked.connect(
            lambda: self.pre_fit_EQE(self.data_xFit_1, self.ui.startExtraPlot, self.ui.stopExtraPlot,
                                     self.ui.startExtraFit, self.ui.stopExtraFit, self.ui.startExtraFitPlot,
                                     self.ui.stopExtraFitPlot, self.ui.textBox_xF1, self.ui.textBox_xF2,
                                     self.ui.textBox_xF3, 'x1'))

        # Handle Heat Map Button
        self.ui.extraHeatButton.clicked.connect(
            lambda: self.heatMap(self.data_xFit_1, self.ui.extraStartStart, self.ui.extraStartStop,
                                 self.ui.extraStopStart, self.ui.extraStopStop, self.ui.textBox_xF1,
                                 self.ui.textBox_xF2, self.ui.textBox_xF3, 'x1'))

        # Handle Clear Extra Fit Button
        self.ui.clearButton_extraFit.clicked.connect(self.clear_EQE_plot)

        # Page 5 - Fit EL and EQE

        self.EL = []
        self.EL_EQE = []
        self.Red_EL_cal = []
        self.Red_EL_meas = []
        self.Red_EQE_meas = []

        # Handle Import EL and EQE Data Buttons
        self.ui.browseELButton_1.clicked.connect(lambda: self.writeText(self.ui.textBox_EL1, 'el1'))
        self.ui.browseELButton_2.clicked.connect(lambda: self.writeText(self.ui.textBox_EL4, 'el2'))

        # Handle EL Plot Buttons
        self.ui.plotELButton_1.clicked.connect(
            lambda: self.pre_plot_EL_EQE(self.EL, self.ui.startPlot_EL1, self.ui.stopPlot_EL1, 0))  # Plot EL
        self.ui.plotELButton_2.clicked.connect(
            lambda: self.pre_plot_EL_EQE(self.EL, self.ui.startPlot_EL2, self.ui.stopPlot_EL2, 1))  # Plot Abs
        self.ui.plotELButton_3.clicked.connect(
            lambda: self.pre_plot_EL_EQE(self.EL_EQE, self.ui.startPlot_EQE, self.ui.stopPlot_EQE, 2))  # Plot EQE

        # Handle Fit Buttons
        self.ui.fitButton_EL1.clicked.connect(
            lambda: self.pre_plot_EL_EQE(self.EL, self.ui.startPlot_EL1, self.ui.stopPlot_EL1, 0, fit=True))  # Fit EL
        self.ui.fitButton_EL2.clicked.connect(
            lambda: self.pre_plot_EL_EQE(self.EL, self.ui.startPlot_EL2, self.ui.stopPlot_EL2, 1, fit=True))  # Fit Abs
        self.ui.fitButton_EL3.clicked.connect(
            lambda: self.pre_plot_EL_EQE(self.EL_EQE, self.ui.startPlot_EQE, self.ui.stopPlot_EQE, 2,
                                         fit=True))  # Fit EQE

        # Handle Clear EL Plot Button
        self.ui.clearButton_EL.clicked.connect(lambda: self.clear_EL_plot())

        # Page 6 - Subtract and Add Peak Fits

        # Subtract Peak Fits

        self.data_subFit = []
        self.data_subEQE = []

        # Handle Import Fit and EQE Data Buttons
        self.ui.browseSubButton_Fit.clicked.connect(lambda: self.writeText(self.ui.textBox_p6_1, 'sub1'))
        self.ui.browseSubButton_EQE.clicked.connect(lambda: self.writeText(self.ui.textBox_p6_4, 'sub2'))

        # Handle Subtract Fit Data Button
        self.ui.subButton.clicked.connect(
            lambda: self.subtract_Fit(self.data_subFit, self.data_subEQE, self.ui.textBox_p6_2, self.ui.textBox_p6_5,
                                      self.ui.textBox_p6_3, self.ui.textBox_p6_6))

        # Add Peak Fits

        self.data_addOptFit = []
        self.data_addCTFit = []
        self.data_addEQE = []

        # Handle Import Fit and EQE Data Buttons
        self.ui.browseAddButton_optFit.clicked.connect(lambda: self.writeText(self.ui.textBox_p7_1, 'add1'))
        self.ui.browseAddButton_CTFit.clicked.connect(lambda: self.writeText(self.ui.textBox_p7_4, 'add2'))
        self.ui.browseAddButton_EQE.clicked.connect(lambda: self.writeText(self.ui.textBox_p7_7, 'add3'))

        # Handle Plot Fit Button
        self.ui.plotAddButton.clicked.connect(
            lambda: self.add_Fits(self.data_addOptFit, self.data_addCTFit, self.data_addEQE))

        # Import Photodiode Calibration Files

        Si_file = pd.ExcelFile(
            "calibration_files/FDS100-CAL.xlsx")  # The files are in the sEQE Analysis folder, just as this program
        self.Si_cal = Si_file.parse('Sheet1')

        InGaAs_file = pd.ExcelFile("calibration_files/FGA21-CAL.xlsx")
        self.InGaAs_cal = InGaAs_file.parse('Sheet1')

        # Define Variables

        self.data_dir = '/home/jungbluth/Desktop'  # FIX: change this if needed

        self.h = 6.626 * math.pow(10, -34)  # [m^2 kg/s]
        self.h_2 = 4.136 * math.pow(10, -15)  # [eV s]
        self.c = 2.998 * math.pow(10, 8)  # [m/s]
        self.q = 1.602 * math.pow(10, -19)  # [C]
        self.k = 8.617 * math.pow(10, -5)  # [ev/K]

        self.export = False
        self.do_plot = True
        self.fit_plot = True
        self.do_plot_EL = True

    # -----------------------------------------------------------------------------------------------------------

    # Functions to read file and update textbox

    # -----------------------------------------------------------------------------------------------------------

    def writeText(self, text_Box, textBox_no):

        os.chdir(self.data_dir)
        file_ = filedialog.askopenfilename()

        if len(file_) != 0:
            path_, filename_ = os.path.split(file_)

            text_Box.clear()  # Clear the text box in case sth has been uploaded already
            text_Box.insertPlainText(filename_)  # Insert filename into text box

            # Page 1 - Calculate EQE

            # Reference files:

            if textBox_no == 1:
                self.ref_1 = pd.read_csv(file_)  # Turn file into dataFrame

            elif textBox_no == 3:
                self.ref_2 = pd.read_csv(file_)

            elif textBox_no == 5:
                self.ref_3 = pd.read_csv(file_)

            elif textBox_no == 7:
                self.ref_4 = pd.read_csv(file_)

            elif textBox_no == 9:
                self.ref_5 = pd.read_csv(file_)

            elif textBox_no == 11:
                self.ref_6 = pd.read_csv(file_)

            # Data files:

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

            elif textBox_no == 12:
                self.data_6 = pd.read_csv(file_)

            # Page 2 - Plot EQE

            elif textBox_no == 'p1':
                self.EQE_1 = pd.read_csv(file_)

            elif textBox_no == 'p4':
                self.EQE_2 = pd.read_csv(file_)

            elif textBox_no == 'p7':
                self.EQE_3 = pd.read_csv(file_)

            elif textBox_no == 'p10':
                self.EQE_4 = pd.read_csv(file_)

            elif textBox_no == 'p13':
                self.EQE_5 = pd.read_csv(file_)

            elif textBox_no == 'p16':
                self.EQE_6 = pd.read_csv(file_)

            elif textBox_no == 'p19':
                self.EQE_7 = pd.read_csv(file_)

            elif textBox_no == 'p22':
                self.EQE_8 = pd.read_csv(file_)

            elif textBox_no == 'p25':
                self.EQE_9 = pd.read_csv(file_)

            elif textBox_no == 'p28':
                self.EQE_10 = pd.read_csv(file_)

            # Page 3 - Fit EQE (Marcus Theory)

            elif textBox_no == 'f1':
                self.data_fit_1 = pd.read_csv(file_)

            elif textBox_no == 'f4':
                self.data_fit_2 = pd.read_csv(file_)

            # For double fits

            elif textBox_no == 'double1':
                self.data_double = pd.read_csv(file_)

            # Page 4 - Fit EQE (MLJ Theory)

            elif textBox_no == 'xF1':
                self.data_xFit_1 = pd.read_csv(file_)

            # Page 5 - Fit EL and EQE

            elif textBox_no == 'el1':
                self.EL = pd.read_table(file_, sep=',', index_col=0)

            elif textBox_no == 'el2':
                self.EL_EQE = pd.read_csv(file_)

            # Page 6 - Subtract and Add Peak Fits

            # Subtract Peak Fits

            elif textBox_no == 'sub1':
                self.data_subFit = pd.read_csv(file_)

            elif textBox_no == 'sub2':
                self.data_subEQE = pd.read_csv(file_)

            # Add Peak Fits

            elif textBox_no == 'add1':
                self.data_addOptFit = pd.read_csv(file_)

            elif textBox_no == 'add2':
                self.data_addCTFit = pd.read_csv(file_)

            elif textBox_no == 'add3':
                self.data_addEQE = pd.read_csv(file_)

    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------

    # Page 1 - Calculate EQE

    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------

    # Function to select data and reference files

    def pre_EQE(self, ref_df, data_df, start, stop, range_no):

        startNM = start.value()
        stopNM = stop.value()

        if Ref_Data_is_valid(ref_df, data_df, startNM, stopNM, range_no):
            self.calculate_EQE(ref_df, data_df, startNM, stopNM, range_no)

    # -----------------------------------------------------------------------------------------------------------

    # Function to calculate EQE

    def calculate_EQE(self, ref_df, data_df, startNM, stopNM, range_no):

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
                        ref_df['Power'] = calculate_Power(ref_df, self.Si_cal)
                    except:
                        print('Please select a valid reference diode.')
                elif self.ui.Range1_InGaAs_button.isChecked() and not self.ui.Range1_Si_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.InGaAs_cal)
                    except:
                        print('Please select a valid reference diode.')
                else:
                    print('Please select a valid reference diode.')

            elif range_no == 2:
                if self.ui.Range2_Si_button.isChecked() and not self.ui.Range2_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.Si_cal)
                    except:
                        print('Please select a valid reference diode.')
                elif self.ui.Range2_InGaAs_button.isChecked() and not self.ui.Range2_Si_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.InGaAs_cal)
                    except:
                        print('Please select a valid reference diode.')
                else:
                    print('Please select a valid reference diode.')

            elif range_no == 3:
                if self.ui.Range3_Si_button.isChecked() and not self.ui.Range3_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.Si_cal)
                    except:
                        print('Please select a valid reference diode.')
                elif self.ui.Range3_InGaAs_button.isChecked() and not self.ui.Range3_Si_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.InGaAs_cal)
                    except:
                        print('Please select a valid reference diode.')
                else:
                    print('Please select a valid reference diode.')

            elif range_no == 4:
                if self.ui.Range4_Si_button.isChecked() and not self.ui.Range4_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.Si_cal)
                    except:
                        print('Please select a valid reference diode.')
                elif self.ui.Range4_InGaAs_button.isChecked() and not self.ui.Range4_Si_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.InGaAs_cal)
                    except:
                        print('Please select a valid reference diode.')
                else:
                    print('Please select a valid reference diode.')

            elif range_no == 5:
                if self.ui.Range5_Si_button.isChecked() and not self.ui.Range5_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.Si_cal)
                    except:
                        print('Please select a valid reference diode.')
                elif self.ui.Range5_InGaAs_button.isChecked() and not self.ui.Range5_Si_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.InGaAs_cal)
                    except:
                        print('Please select a valid reference diode.')
                else:
                    print('Please select a valid reference diode.')

            elif range_no == 6:
                if self.ui.Range6_Si_button.isChecked() and not self.ui.Range6_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.Si_cal)
                    except:
                        print('Please select a valid reference diode.')
                elif self.ui.Range6_InGaAs_button.isChecked() and not self.ui.Range6_Si_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.InGaAs_cal)
                    except:
                        print('Please select a valid reference diode.')
                else:
                    print('Please select a valid reference diode.')

        if 'Power' in ref_df.columns:  # Check if the power has been calculated already

            for x in range(len(ref_df['Wavelength'])):  # Iterate through columns of reference file
                power_dict[ref_df['Wavelength'][x]] = ref_df['Power'][
                    x]  # Add wavelength and corresponding power to dictionary

            for y in range(len(data_df['Wavelength'])):  # Iterate through columns of data file
                if startNM <= data_df['Wavelength'][
                    y] <= stopNM:  # Calculate EQE only if start <= wavelength <= stop, otherwise ignore

                    if data_df['Wavelength'][y] in power_dict.keys():  # Check if data wavelength is in reference file
                        Wavelength.append(data_df['Wavelength'][y])
                        Energy_val = (self.h * self.c) / (
                                data_df['Wavelength'][y] * math.pow(10, -9) * self.q)  # Calculate energy in eV
                        Energy.append(Energy_val)
                        EQE_val = (data_df['Mean Current'][y] * Energy_val) / (
                            power_dict[data_df['Wavelength'][y]])  # Easier way to calculate EQE
                        # EQE_val = ((data_df['Mean Current'][y] * self.h * self.c) / (data_df['Wavelength'][y] * math.pow(10,-9) * power_dict[data_df['Wavelength'][y]] * self.q))
                        # EQE_val = (100 * data_df['Mean Current'][y] * Energy_val) / (power_dict[data_df['Wavelength'][y]]) # *100 to turn into percent
                        EQE.append(EQE_val)
                        log_EQE.append(math.log10(EQE_val))

                    else:  # If data wavelength is not in reference file
                        Wavelength.append(data_df['Wavelength'][y])
                        Energy_val = (self.h * self.c) / (data_df['Wavelength'][y] * math.pow(10, -9) * self.q)
                        Energy.append(Energy_val)
                        Power_int = interpolate(data_df['Wavelength'][y], ref_df['Wavelength'],
                                                ref_df['Power'])  # Interpolate power
                        EQE_int = (data_df['Mean Current'][y] * Energy_val) / (Power_int)
                        # EQE_int = ((data_df['Mean Current'][y] * self.h * self.c) / (data_df['Wavelength'][y] * math.pow(10,-9) * Power_int * self.q))
                        EQE.append(EQE_int)
                        log_EQE.append(math.log10(EQE_int))

        if len(Wavelength) == len(EQE) and len(Energy) == len(log_EQE):  # Check if the lists have the same length

            if self.export:  # If the "Export Data" button has been clicked
                return (Wavelength, Energy, EQE, log_EQE)

            else:  # If the "Calculate EQE" button has been clicked
                if self.do_plot:  # This is set to true during setup of the program
                    self.ax1, self.ax2 = set_up_plot()
                    self.do_plot = False  # Set self.do_plot to False to plot on the same graph

                label_ = pick_Label(range_no, startNM, stopNM)

                self.ax1.plot(Wavelength, EQE, linewidth=3, label=label_)
                self.ax2.semilogy(Wavelength, EQE,
                                  linewidth=3)  # Equivalent to the line above but with proper log scale axes
                self.ax1.legend()
                plt.draw()

        else:
            print('Error Code 1: Length mismatch.')

    # -----------------------------------------------------------------------------------------------------------

    # Function to export EQE

    def export_EQE(self):

        self.export = True

        ok_1 = True  # Create boolean variable to use for "is_valid" function
        ok_2 = True
        ok_3 = True
        ok_4 = True
        ok_5 = True
        ok_6 = True

        columns = ['Wavelength', 'Energy', 'EQE', 'Log_EQE']  # Columns for dataFrame

        export_file = pd.DataFrame(columns=columns)  # Create empty dataFrame

        wave_inc = {}  # create empty dictionary

        if self.ui.exportBox_1.isChecked():  # If the checkBox is checked
            startNM1 = self.ui.startNM_1.value()  # Pick start wavelength
            stopNM1 = self.ui.stopNM_1.value()  # Pick stop wavelength
            if Ref_Data_is_valid(self.ref_1, self.data_1, startNM1, stopNM1,
                                 1):  # Check that files are non-empty and within wavelength range
                Wave_1, Energy_1, EQE_1, log_EQE_1 = self.calculate_EQE(self.ref_1, self.data_1, startNM1, stopNM1,
                                                                        1)  # Extract data
                export_1 = pd.DataFrame({'Wavelength': Wave_1, 'Energy': Energy_1, 'EQE': EQE_1,
                                         'Log_EQE': log_EQE_1})  # Create dataFrame with EQE data
                wave_inc['1'] = Wave_1[0]  # Add the first wavelength value to the wave_inc list
            else:
                ok_1 = False  # Set variable to False if calculation is invalid

        if self.ui.exportBox_2.isChecked():
            startNM2 = self.ui.startNM_2.value()
            stopNM2 = self.ui.stopNM_2.value()
            if Ref_Data_is_valid(self.ref_2, self.data_2, startNM2, stopNM2, 2):
                Wave_2, Energy_2, EQE_2, log_EQE_2 = self.calculate_EQE(self.ref_2, self.data_2, startNM2, stopNM2, 2)
                export_2 = pd.DataFrame({'Wavelength': Wave_2, 'Energy': Energy_2, 'EQE': EQE_2, 'Log_EQE': log_EQE_2})
                wave_inc['2'] = Wave_2[0]
            else:
                ok_2 = False

        if self.ui.exportBox_3.isChecked():
            startNM3 = self.ui.startNM_3.value()
            stopNM3 = self.ui.stopNM_3.value()
            if Ref_Data_is_valid(self.ref_3, self.data_3, startNM3, stopNM3, 3):
                Wave_3, Energy_3, EQE_3, log_EQE_3 = self.calculate_EQE(self.ref_3, self.data_3, startNM3, stopNM3, 3)
                export_3 = pd.DataFrame({'Wavelength': Wave_3, 'Energy': Energy_3, 'EQE': EQE_3, 'Log_EQE': log_EQE_3})
                wave_inc['3'] = Wave_3[0]
            else:
                ok_3 = False

        if self.ui.exportBox_4.isChecked():
            startNM4 = self.ui.startNM_4.value()
            stopNM4 = self.ui.stopNM_4.value()
            if Ref_Data_is_valid(self.ref_4, self.data_4, startNM4, stopNM4, 4):
                Wave_4, Energy_4, EQE_4, log_EQE_4 = self.calculate_EQE(self.ref_4, self.data_4, startNM4, stopNM4, 4)
                export_4 = pd.DataFrame({'Wavelength': Wave_4, 'Energy': Energy_4, 'EQE': EQE_4, 'Log_EQE': log_EQE_4})
                wave_inc['4'] = Wave_4[0]
            else:
                ok_4 = False

        if self.ui.exportBox_5.isChecked():
            startNM5 = self.ui.startNM_5.value()
            stopNM5 = self.ui.stopNM_5.value()
            if Ref_Data_is_valid(self.ref_5, self.data_5, startNM5, stopNM5, 5):
                Wave_5, Energy_5, EQE_5, log_EQE_5 = self.calculate_EQE(self.ref_5, self.data_5, startNM5, stopNM5, 5)
                export_5 = pd.DataFrame({'Wavelength': Wave_5, 'Energy': Energy_5, 'EQE': EQE_5, 'Log_EQE': log_EQE_5})
                wave_inc['5'] = Wave_5[0]
            else:
                ok_5 = False

        if self.ui.exportBox_6.isChecked():
            startNM6 = self.ui.startNM_6.value()
            stopNM6 = self.ui.stopNM_6.value()
            if Ref_Data_is_valid(self.ref_6, self.data_6, startNM6, stopNM6, 6):
                Wave_6, Energy_6, EQE_6, log_EQE_6 = self.calculate_EQE(self.ref_6, self.data_6, startNM6, stopNM6, 6)
                export_6 = pd.DataFrame({'Wavelength': Wave_6, 'Energy': Energy_6, 'EQE': EQE_6, 'Log_EQE': log_EQE_6})
                wave_inc['6'] = Wave_6[0]
            else:
                ok_6 = False

        if ok_1 and ok_2 and ok_3 and ok_4 and ok_5 and ok_6:  # Check if all operations are ok or if fields are left empty

            for x in range(len(wave_inc.keys())):  # Iterate through wave_inc list
                minimum = min(wave_inc, key=wave_inc.get)  # Find key corresponding to minimum value
                if minimum == '1':  # Append correct dataFrame in order of decending wavelength
                    export_file = export_file.append(export_1, ignore_index=True)
                elif minimum == '2':
                    export_file = export_file.append(export_2, ignore_index=True)
                elif minimum == '3':
                    export_file = export_file.append(export_3, ignore_index=True)
                elif minimum == '4':
                    export_file = export_file.append(export_4, ignore_index=True)
                elif minimum == '5':
                    export_file = export_file.append(export_5, ignore_index=True)
                elif minimum == '6':
                    export_file = export_file.append(export_6, ignore_index=True)

                del wave_inc[minimum]  # Delete recently appended value

            EQE_file = filedialog.asksaveasfilename()  # Prompt the user to pick a folder & name to save data to
            export_path, export_filename = os.path.split(EQE_file)

            if len(export_path) != 0:  # Check if the user actually selected a path
                os.chdir(export_path)  # Change the working directory
                export_file.to_csv(export_filename)  # Save data to csv
                print('Saving data to: %s' % str(EQE_file))
                os.chdir(self.data_dir)  # Change the directory back

        self.export = False

    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------

    # Page 2 - Calculate EQE

    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------

    # Function to select EQE for plotting

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

        self.axEQE_1, self.axEQE_2 = set_up_EQE_plot(number, norm_num)

        if self.ui.plotBox_1.isChecked():
            ok_EQE_1 = self.plot_EQE(self.EQE_1, self.ui.startEQE_1, self.ui.stopEQE_1, self.ui.textBox_p2_1,
                                     self.ui.textBox_p2_2, self.ui.textBox_p2_3, 1, number)

        if self.ui.plotBox_2.isChecked():
            ok_EQE_2 = self.plot_EQE(self.EQE_2, self.ui.startEQE_2, self.ui.stopEQE_2, self.ui.textBox_p2_4,
                                     self.ui.textBox_p2_5, self.ui.textBox_p2_6, 2, number)

        if self.ui.plotBox_3.isChecked():
            ok_EQE_3 = self.plot_EQE(self.EQE_3, self.ui.startEQE_3, self.ui.stopEQE_3, self.ui.textBox_p2_7,
                                     self.ui.textBox_p2_8, self.ui.textBox_p2_9, 3, number)

        if self.ui.plotBox_4.isChecked():
            ok_EQE_4 = self.plot_EQE(self.EQE_4, self.ui.startEQE_4, self.ui.stopEQE_4, self.ui.textBox_p2_10,
                                     self.ui.textBox_p2_11, self.ui.textBox_p2_12, 4, number)

        if self.ui.plotBox_5.isChecked():
            ok_EQE_5 = self.plot_EQE(self.EQE_5, self.ui.startEQE_5, self.ui.stopEQE_5, self.ui.textBox_p2_13,
                                     self.ui.textBox_p2_14, self.ui.textBox_p2_15, 5, number)

        if self.ui.plotBox_6.isChecked():
            ok_EQE_6 = self.plot_EQE(self.EQE_6, self.ui.startEQE_6, self.ui.stopEQE_6, self.ui.textBox_p2_16,
                                     self.ui.textBox_p2_17, self.ui.textBox_p2_18, 6, number)

        if self.ui.plotBox_7.isChecked():
            ok_EQE_7 = self.plot_EQE(self.EQE_7, self.ui.startEQE_7, self.ui.stopEQE_7, self.ui.textBox_p2_19,
                                     self.ui.textBox_p2_20, self.ui.textBox_p2_21, 7, number)

        if self.ui.plotBox_8.isChecked():
            ok_EQE_8 = self.plot_EQE(self.EQE_8, self.ui.startEQE_8, self.ui.stopEQE_8, self.ui.textBox_p2_22,
                                     self.ui.textBox_p2_23, self.ui.textBox_p2_24, 8, number)

        if self.ui.plotBox_9.isChecked():
            ok_EQE_9 = self.plot_EQE(self.EQE_9, self.ui.startEQE_9, self.ui.stopEQE_9, self.ui.textBox_p2_25,
                                     self.ui.textBox_p2_26, self.ui.textBox_p2_27, 9, number)

        if self.ui.plotBox_10.isChecked():
            ok_EQE_10 = self.plot_EQE(self.EQE_10, self.ui.startEQE_10, self.ui.stopEQE_10, self.ui.textBox_p2_28,
                                      self.ui.textBox_p2_29, self.ui.textBox_p2_30, 10, number)

        if ok_EQE_1 and ok_EQE_2 and ok_EQE_3 and ok_EQE_4 and ok_EQE_5 and ok_EQE_6 and ok_EQE_7 and ok_EQE_8 and ok_EQE_9 and ok_EQE_10:
            self.axEQE_1.legend()
            self.axEQE_2.legend()
            plt.show()
        else:
            plt.close()
            plt.close()

    # -----------------------------------------------------------------------------------------------------------

    # Function to plot EQE

    def plot_EQE(self, eqe_df, startNM, stopNM, filename_Box, label_Box, color_Box, file_no, number):

        startNM = startNM.value()  # Pick start wavelength
        stopNM = stopNM.value()  # Pick stop wavelength

        if EQE_is_valid(eqe_df, startNM, stopNM, file_no):  # Check that files are non-empty and within wavelength range

            if self.ui.normalizeBox.isChecked():
                normNM = self.ui.normalizeNM.value()
                if Normalization_is_valid(eqe_df, normNM, file_no):
                    wave, energy, eqe, log_eqe = normalize_EQE(eqe_df, startNM, stopNM, normNM)
                else:
                    return False
            else:
                wave, energy, eqe, log_eqe = compile_EQE(eqe_df, startNM, stopNM, 0)

            label_ = pick_EQE_Label(label_Box, filename_Box)
            color_ = pick_EQE_Color(color_Box, file_no)

            if number == 0:
                self.axEQE_1.plot(wave, eqe, linewidth=3, label=label_, color=color_)
                self.axEQE_2.semilogy(wave, eqe, linewidth=3, label=label_, color=color_)
            elif number == 1:
                self.axEQE_1.plot(energy, eqe, linewidth=3, label=label_, color=color_)
                self.axEQE_2.semilogy(energy, eqe, linewidth=3, label=label_, color=color_)

            return True

        else:
            return False

    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------

    # Page 3 - Fit EQE (Marcus Theory)
    # Page 4 - Fit EQE (MLJ Theory)

    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------

    # Function to select EQE to fit

    def pre_fit_EQE(self, eqe_df, startE, stopE, startFit, stopFit, startPlotFit, stopPlotFit, filename_Box, label_Box,
                    color_Box, file_no):

        # Change this if you want to plot on the same graph 
        if self.fit_plot:
            self.axFit_1, self.axFit_2 = set_up_EQE_plot()  # Sets up plot of energy vs EQE / log(EQE)
            self.fit_plot = False

        ok_plot_Fit = self.plot_fit_EQE(eqe_df, startE, stopE, startFit, stopFit, startPlotFit, stopPlotFit,
                                        filename_Box, label_Box, color_Box, file_no)

        if ok_plot_Fit:
            self.axFit_1.legend()
            self.axFit_2.legend()
            plt.show()

    # -----------------------------------------------------------------------------------------------------------

    # Function to fit and plot EQE

    def plot_fit_EQE(self, eqe_df, startE, stopE, startFit, stopFit, startPlotFit, stopPlotFit, filename_Box, label_Box,
                     color_Box, file_no):

        include_Disorder = False
        fit_opticalPeak = False

        startE = startE.value()  # Pick start energy
        stopE = stopE.value()  # Pick stop energy

        startFit = startFit.value()  # Pick start fit energy
        stopFit = stopFit.value()  # Pick stop fit energy

        startPlotFit = startPlotFit.value()
        stopPlotFit = stopPlotFit.value()

        if Fit_is_valid(eqe_df, startE, stopE, startFit, stopFit,
                        file_no):  # Check that files are non-empty and within energy range
            wave, energy, eqe, log_eqe = compile_EQE(eqe_df, startE, stopE, 1)  # Compile EQE file
            wave_fit, energy_fit, eqe_fit, log_eqe_fit = compile_EQE(eqe_df, startFit, stopFit,
                                                                     1)  # Compile fit range of EQE file

            label_ = pick_EQE_Label(label_Box, filename_Box)
            color_ = pick_EQE_Color(color_Box, file_no)

            # Fit EQE (Marcus Theory)
            if (str(file_no)).isnumeric():

                if file_no == 1:
                    self.T_CT = self.ui.Temperature_1.value()
                    guessStart = self.ui.guessStart_1.value()
                    guessStop = self.ui.guessStop_1.value()
                    if self.ui.static_Disorder_1.isChecked():
                        include_Disorder = True
                        self.sig = self.ui.Disorder_1.value()
                    if self.ui.OptButton_1.isChecked() and not self.ui.CTButton_1.isChecked():  # Check if CT state or optical gap are fitted.
                        fit_opticalPeak = True
                    elif self.ui.OptButton_1.isChecked() and self.ui.CTButton_1.isChecked():
                        print('Please select a valid peak to fit.')
                    elif not self.ui.OptButton_1.isChecked() and not self.ui.CTButton_1.isChecked():
                        print('Please select a valid peak to fit.')

                elif file_no == 2:
                    self.T_CT = self.ui.Temperature_2.value()
                    guessStart = self.ui.guessStart_2.value()
                    guessStop = self.ui.guessStop_2.value()
                    if self.ui.static_Disorder_2.isChecked():
                        include_Disorder = True
                        self.sig = self.ui.Disorder_2.value()
                    if self.ui.OptButton_2.isChecked() and not self.ui.CTButton_2.isChecked():
                        fit_opticalPeak = True
                    elif self.ui.OptButton_2.isChecked() and self.ui.CTButton_2.isChecked():
                        print('Please select a valid peak to fit.')
                    elif not self.ui.OptButton_2.isChecked() and not self.ui.CTButton_2.isChecked():
                        print('Please select a valid peak to fit.')

                # Attempt peak fit:
                x_gaussian = linspace(startPlotFit, stopPlotFit, 50)

                ECT_guess = np.arange(guessStart, guessStop + 0.1, 0.05)
                p0 = None

                for ECT in ECT_guess:
                    y_gaussian = []
                    try:
                        if include_Disorder:
                            best_vals, covar, y_fit, r_squared = self.fit_function(self.gaussian_disorder, energy_fit,
                                                                                   eqe_fit, p0=p0)
                            for value in x_gaussian:
                                y_gaussian.append(
                                    self.gaussian_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                        else:
                            best_vals, covar, y_fit, r_squared = self.fit_function(self.gaussian, energy_fit, eqe_fit,
                                                                                   p0=p0)
                            for value in x_gaussian:
                                y_gaussian.append(self.gaussian(value, best_vals[0], best_vals[1], best_vals[2]))

                        if r_squared > 0:
                            print('Initial Guess (eV) : ', p0)

                            print('-' * 80)
                            print('Temperature [T] (K) : ', self.T_CT)
                            print('Oscillator Strength [f] (eV**2) : ', format(best_vals[0], '.6f'), '+/-',
                                  format(math.sqrt(covar[0, 0]), '.6f'))
                            print('Reorganization Energy [l] (eV) : ', format(best_vals[1], '.6f'), '+/-',
                                  format(math.sqrt(covar[1, 1]), '.6f'))
                            if fit_opticalPeak:
                                print('Optical Peak Energy [E_Opt] (eV) : ', format(best_vals[2], '.6f'), '+/-',
                                      format(math.sqrt(covar[2, 2]), '.6f'))
                            else:
                                print('CT State Energy [ECT] (eV) : ', format(best_vals[2], '.6f'), '+/-',
                                      format(math.sqrt(covar[2, 2]), '.6f'))
                            print('R_Squared : ', format(r_squared, '.6f'))
                            print('-' * 80)

                            # Plot EQE data and CT fit
                            self.axFit_1.plot(energy, eqe, linewidth=3, label=label_, color=color_)
                            plt.draw()
                            if include_Disorder:
                                if fit_opticalPeak:
                                    self.axFit_1.plot(x_gaussian, y_gaussian, linewidth=2,
                                                      label='Gaussian Fit + Disorder', color='#000000',
                                                      linestyle='dotted')
                                else:
                                    self.axFit_1.plot(x_gaussian, y_gaussian, linewidth=2,
                                                      label='Gaussian Fit + Disorder', color='#000000', linestyle='--')
                            else:
                                if fit_opticalPeak:
                                    self.axFit_1.plot(x_gaussian, y_gaussian, linewidth=2, label='Gaussian Fit',
                                                      color='#000000', linestyle='dotted')
                                else:
                                    self.axFit_1.plot(x_gaussian, y_gaussian, linewidth=2, label='Gaussian Fit',
                                                      color='#000000', linestyle='--')
                            # plt.draw()

                            self.axFit_2.semilogy(energy, eqe, linewidth=3, label=label_, color=color_)
                            plt.draw()
                            if include_Disorder:
                                if fit_opticalPeak:
                                    self.axFit_2.plot(x_gaussian, y_gaussian, linewidth=2,
                                                      label='Gaussian Fit + Disorder', color='#000000',
                                                      linestyle='dotted')
                                else:
                                    self.axFit_2.plot(x_gaussian, y_gaussian, linewidth=2,
                                                      label='Gaussian Fit + Disorder', color='#000000', linestyle='--')
                            else:
                                if fit_opticalPeak:
                                    self.axFit_2.plot(x_gaussian, y_gaussian, linewidth=2, label='Gaussian Fit',
                                                      color='#000000', linestyle='dotted')
                                else:
                                    self.axFit_2.plot(x_gaussian, y_gaussian, linewidth=2, label='Gaussian Fit',
                                                      color='#000000', linestyle='--')
                            # plt.draw()

                            # Save fit data
                            if self.ui.save_gaussianFit.isChecked():
                                fit_file = pd.DataFrame()
                                fit_file['Energy'] = x_gaussian
                                fit_file['Signal'] = y_gaussian

                                # FIX to add as header rather than columns
                                fit_file['Temperature'] = self.T_CT
                                fit_file['Oscillator Strength (eV**2)'] = best_vals[0]
                                fit_file['Reorganization Energy (eV)'] = best_vals[1]
                                if fit_opticalPeak:
                                    fit_file['Optical Peak Energy (eV)'] = best_vals[2]
                                else:
                                    fit_file['CT State Energy (eV)'] = best_vals[2]
                                save_fit_file = filedialog.asksaveasfilename()  # User to pick folder & name to save to
                                save_fit_path, save_fit_filename = os.path.split(save_fit_file)
                                if len(save_fit_path) != 0:  # Check if the user actually selected a path
                                    os.chdir(save_fit_path)  # Change the working directory
                                    fit_file.to_csv(save_fit_filename)  # Save data to csv
                                    print('Saving fit data to: %s' % str(save_fit_file))
                                    os.chdir(self.data_dir)  # Change the directory back

                            return True  # This breaks the fore loop

                        else:
                            raise Exception('Wrong fit determined.')

                    except:
                        p0 = [0.001, 0.1, round(ECT, 3)]

                print('Optimal parameters not found.')
                return False

            # Fit EQE (MLJ Theory)
            elif file_no == 'x1':

                self.S_i = self.ui.Huang_Rhys.value()
                self.hbarw_i = self.ui.vib_Energy.value()
                self.T_x = self.ui.extra_Temperature.value()

                guessStart = self.ui.extraGuessStart.value()
                guessStop = self.ui.extraGuessStop.value()

                if self.ui.extra_static_Disorder.isChecked():
                    include_Disorder = True
                    self.sig_x = self.ui.extra_Disorder.value()

                # Attempt peak fit:
                x_MLJ_theory = linspace(startPlotFit, stopPlotFit, 50)

                ECT_guess = np.arange(guessStart, guessStop + 0.1, 0.05)
                p0 = None

                for ECT in ECT_guess:
                    y_MLJ_theory = []
                    try:
                        if include_Disorder:
                            best_vals, covar, y_fit, r_squared = self.fit_function(self.MLJ_gaussian_disorder,
                                                                                   energy_fit, eqe_fit, p0=p0)
                            for value in x_MLJ_theory:
                                y_MLJ_theory.append(
                                    self.MLJ_gaussian_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                        else:
                            best_vals, covar, y_fit, r_squared = self.fit_function(self.MLJ_gaussian, energy_fit,
                                                                                   eqe_fit, p0=p0)
                            for value in x_MLJ_theory:
                                y_MLJ_theory.append(self.MLJ_gaussian(value, best_vals[0], best_vals[1], best_vals[2]))

                        if r_squared > 0:
                            print('Initial Guess (eV) : ', p0)

                            print('-' * 80)
                            print('Temperature [T] (K) : ', self.T_x)
                            print('Oscillator Strength [f] (eV**2) : ', format(best_vals[0], '.6f'), '+/-',
                                  format(math.sqrt(covar[0, 0]), '.6f'))
                            print('Reorganization Energy [l] (eV) : ', format(best_vals[1], '.6f'), '+/-',
                                  format(math.sqrt(covar[1, 1]), '.6f'))
                            print('CT State Energy [ECT] (eV) : ', format(best_vals[2], '.6f'), '+/-',
                                  format(math.sqrt(covar[2, 2]), '.6f'))
                            print('R_Squared : ', format(r_squared, '.6f'))
                            print('-' * 80)

                            # Plot EQE data and CT fit
                            self.axFit_1.plot(energy, eqe, linewidth=3, label=label_, color=color_)
                            plt.draw()
                            if include_Disorder:
                                self.axFit_1.plot(x_MLJ_theory, y_MLJ_theory, linewidth=2, label='MLJ Fit + Disorder',
                                                  color='#000000', linestyle='--')
                            else:
                                self.axFit_1.plot(x_MLJ_theory, y_MLJ_theory, linewidth=2, label='MLJ Fit',
                                                  color='#000000', linestyle='--')
                            plt.draw()

                            self.axFit_2.semilogy(energy, eqe, linewidth=3, label=label_, color=color_)
                            plt.draw()
                            if include_Disorder:
                                self.axFit_2.plot(x_MLJ_theory, y_MLJ_theory, linewidth=2, label='MLJ Fit + Disorder',
                                                  color='#000000', linestyle='--')
                            else:
                                self.axFit_2.plot(x_MLJ_theory, y_MLJ_theory, linewidth=2, label='MLJ Fit',
                                                  color='#000000', linestyle='--')
                            plt.draw()

                            return True

                        else:
                            raise Exception('Wrong fit determined.')

                    except:
                        p0 = [0.001, 0.1, round(ECT, 3)]

                print('Optimal parameters not found.')
                return False

        else:
            return False

    # -----------------------------------------------------------------------------------------------------------

    # Function to generate heat map of fitting values

    def heatMap(self, eqe_df, startStartE, startStopE, stopStartE, stopStopE, filename_Box, label_Box,
                color_Box, file_no):

        include_Disorder = False
        fit_opticalPeak = False

        startStartE = startStartE.value()
        startStopE = startStopE.value()
        stopStartE = stopStartE.value()
        stopStopE = stopStopE.value()

        startStart_ok = StartStop_is_valid(startStartE, startStopE)  # Check that start energy is lower than stop energy
        startStop_ok = StartStop_is_valid(startStopE, stopStartE)
        stopStop_ok = StartStop_is_valid(stopStartE, stopStopE)

        if startStart_ok and startStop_ok and stopStop_ok:  # If all operations are valid, proceed with heat map calculations

            startEnergies = []  # Create empty lists for start and stop energies
            stopEnergies = []

            start_df = []  # Create empty lists for fit data collection
            stop_df = []
            f_df = []
            l_df = []
            Ect_df = []
            R_df = []

            x = float(startStartE)
            y = float(stopStartE)

            while x <= float(startStopE):  # As long as start range start value is smaller than start range stop value
                startEnergies.append(round(x, 3))  # Round to three digits
                x += 0.005  # Adjust this number to change resolution

            while y <= float(stopStopE):
                stopEnergies.append(round(y, 3))
                y += 0.005

            for start in startEnergies:  # Iterate through start energies
                for stop in stopEnergies:  # Iterature through stop energies

                    wave_fit, energy_fit, eqe_fit, log_eqe_fit = compile_EQE(eqe_df, start, stop, 1)

                    # Fit EQE (Marcus Theory)
                    if (str(file_no)).isnumeric():

                        if file_no == 1:
                            self.T_CT = self.ui.Temperature_1.value()
                            guessStart = self.ui.guessStart_1.value()
                            guessStop = self.ui.guessStop_1.value()
                            if self.ui.static_Disorder_1.isChecked():
                                include_Disorder = True
                                self.sig = self.ui.Disorder_1.value()

                            if self.ui.OptButton_1.isChecked() and not self.ui.CTButton_1.isChecked():  # Check if CT state or optical gap are fitted.
                                fit_opticalPeak = True
                            elif self.ui.OptButton_1.isChecked() and self.ui.CTButton_1.isChecked():
                                print('Please select a valid peak to fit.')
                            elif not self.ui.OptButton_1.isChecked() and not self.ui.CTButton_1.isChecked():
                                print('Please select a valid peak to fit.')

                        elif file_no == 2:
                            self.T_CT = self.ui.Temperature_2.value()
                            guessStart = self.ui.guessStart_2.value()
                            guessStop = self.ui.guessStop_2.value()
                            if self.ui.static_Disorder_2.isChecked():
                                include_Disorder = True
                                self.sig = self.ui.Disorder_2.value()

                            if self.ui.OptButton_2.isChecked() and not self.ui.CTButton_2.isChecked():
                                fit_opticalPeak = True
                            elif self.ui.OptButton_2.isChecked() and self.ui.CTButton_2.isChecked():
                                print('Please select a valid peak to fit.')
                            elif not self.ui.OptButton_2.isChecked() and not self.ui.CTButton_2.isChecked():
                                print('Please select a valid peak to fit.')

                        # Attempt peak fit:
                        ECT_guess = np.arange(guessStart, guessStop + 0.1, 0.05)
                        p0 = None

                        for ECT in ECT_guess:
                            try:
                                if include_Disorder:
                                    best_vals, covar, y_fit, r_squared = self.fit_function(self.gaussian_disorder,
                                                                                           energy_fit, eqe_fit,
                                                                                           p0=p0)
                                else:
                                    best_vals, covar, y_fit, r_squared = self.fit_function(self.gaussian,
                                                                                           energy_fit, eqe_fit,
                                                                                           p0=p0)
                                if r_squared > 0:
                                    start_df.append(start)
                                    stop_df.append(stop)
                                    f_df.append(best_vals[0])
                                    l_df.append(best_vals[1])
                                    Ect_df.append(best_vals[2])
                                    R_df.append(r_squared)
                                break
                            except:
                                p0 = [0.001, 0.1, ECT]
                                if ECT == ECT_guess[-1]:
                                    print('Optimal parameters not found.')

                    # Fit EQE (MLJ Theory)
                    elif file_no == 'x1':

                        self.S_i = self.ui.Huang_Rhys.value()
                        self.hbarw_i = self.ui.vib_Energy.value()
                        self.T_x = self.ui.extra_Temperature.value()

                        guessStart = self.ui.extraGuessStart.value()
                        guessStop = self.ui.extraGuessStop.value()

                        if self.ui.extra_static_Disorder.isChecked():
                            include_Disorder = True
                            self.sig_x = self.ui.extra_Disorder.value()

                        # Attempt peak fit:
                        ECT_guess = np.arange(guessStart, guessStop + 0.1, 0.05)
                        p0 = None

                        for ECT in ECT_guess:
                            try:
                                if include_Disorder:
                                    best_vals, covar, y_fit, r_squared = self.fit_function(
                                        self.MLJ_gaussian_disorder, energy_fit, eqe_fit, p0=p0)
                                else:
                                    best_vals, covar, y_fit, r_squared = self.fit_function(self.MLJ_gaussian,
                                                                                           energy_fit, eqe_fit,
                                                                                           p0=p0)
                                if r_squared > 0:
                                    start_df.append(start)
                                    stop_df.append(stop)
                                    f_df.append(best_vals[0])
                                    l_df.append(best_vals[1])
                                    Ect_df.append(best_vals[2])
                                    R_df.append(r_squared)
                                    break
                            except:
                                p0 = [0.001, 0.1, ECT]
                                if ECT == ECT_guess[-1]:
                                    print('Optimal parameters not found.')

            if len(R_df) != 0:
                parameter_df = pd.DataFrame({'Start': start_df, 'Stop': stop_df, 'f': f_df, 'l': l_df, 'Ect': Ect_df,
                                             'R_Squared': R_df})  # Create a dataFrame with all results
                max_index = parameter_df[parameter_df['R_Squared'] == max(parameter_df['R_Squared'])].index.values[0]

                print('-' * 80)
                if file_no == 'x1':
                    print('Temperature [T] (K) : ', self.T_x)
                else:
                    print('Temperature [T] (K) : ', self.T_CT)

                print('Average Oscillator Strength [f] (eV**2) : ', format(parameter_df['f'].mean(), '.6f'), '+/-',
                      format(parameter_df['f'].std(), '.6f'))  # Determine the average value and standard deviation
                print('Average Reorganization Energy [l] (eV) : ', format(parameter_df['l'].mean(), '.6f'), '+/-',
                      format(parameter_df['l'].std(), '.6f'))

                if fit_opticalPeak:
                    print('Average Optical Peak Energy [E_Opt] (eV) : ', format(parameter_df['Ect'].mean(), '.6f'),
                          '+/-', format(parameter_df['Ect'].std(), '.6f'))
                else:
                    print('Average CT State Energy [ECT] (eV) : ', format(parameter_df['Ect'].mean(), '.6f'), '+/-',
                          format(parameter_df['Ect'].std(), '.6f'))

                print('Average R_Squared : ', format(parameter_df['R_Squared'].mean(), '.6f'), '+/-',
                      format(parameter_df['R_Squared'].std(), '.6f'))

                print('-' * 80)

                if max(parameter_df['R_Squared']) == 1.0:
                    print('Max R_squared : ', format(max(parameter_df['R_Squared']), '.6f'))
                    print('Start Energies (eV) : ',
                          np.array(parameter_df['Start'][parameter_df['R_Squared'] == 1.0]).tolist())
                    print('Stop Energies (eV) : ',
                          np.array(parameter_df['Stop'][parameter_df['R_Squared'] == 1.0]).tolist())
                    print('Average Oscillator Strength [f] (eV**2) : ',
                          format(parameter_df['f'][parameter_df['R_Squared'] == 1.0].mean(), '.6f'), '+/-',
                          format(parameter_df['f'][parameter_df['R_Squared'] == 1.0].std(), '.6f'))
                    print('Average Reorganization Energy [l] (eV) : ',
                          format(parameter_df['l'][parameter_df['R_Squared'] == 1.0].mean(), '.6f'), '+/-',
                          format(parameter_df['l'][parameter_df['R_Squared'] == 1.0].std(), '.6f'))
                    if fit_opticalPeak:
                        print('Average Optical Peak Energy [E_Opt] (eV) : ',
                              format(parameter_df['Ect'][parameter_df['R_Squared'] == 1.0].mean(), '.6f'), '+/-',
                              format(parameter_df['Ect'][parameter_df['R_Squared'] == 1.0].std(), '.6f'))
                    else:
                        print('Average CT State Energy [ECT] (eV) : ',
                              format(parameter_df['Ect'][parameter_df['R_Squared'] == 1.0].mean(), '.6f'), '+/-',
                              format(parameter_df['Ect'][parameter_df['R_Squared'] == 1.0].std(), '.6f'))

                else:
                    print('Max R_squared : ', format(max(parameter_df['R_Squared']), '.6f'))
                    print('Start Energy (eV) : ', parameter_df['Start'][max_index])
                    print('Stop Energy (eV) : ', parameter_df['Stop'][max_index])
                print('-' * 80)

                f_df = parameter_df.pivot('Stop', 'Start',
                                          'f')  # Pivot the dataFrame: x-value = Stop, y-value = Start, value = f
                l_df = parameter_df.pivot('Stop', 'Start', 'l')
                Ect_df = parameter_df.pivot('Stop', 'Start', 'Ect')
                R_df = parameter_df.pivot('Stop', 'Start', 'R_Squared')

                plt.ion()
                plt.figure(figsize=(11, 9))  # Create a new figure
                self.heatmap_1 = seaborn.heatmap(f_df, xticklabels=3, yticklabels=3)  # Create the heat map
                plt.xlabel('Initial Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.ylabel('Final Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.title('Oscillator Strength ($eV^2$)', fontsize=17, fontweight='medium')
                cbar = self.heatmap_1.collections[0].colorbar
                cbar.ax.tick_params(labelsize=15)
                plt.yticks(rotation=360)
                plt.xticks(rotation=90)
                plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
                # plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
                plt.show()

                plt.ion()
                plt.figure(figsize=(11, 9))
                self.heatmap_2 = seaborn.heatmap(l_df, xticklabels=3, yticklabels=3)
                plt.xlabel('Initial Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.ylabel('Final Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.title('Reorganization Energy (eV)', fontsize=17, fontweight='medium')
                cbar = self.heatmap_2.collections[0].colorbar
                cbar.ax.tick_params(labelsize=15)
                plt.yticks(rotation=360)
                plt.xticks(rotation=90)
                plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
                plt.show()

                plt.ion()
                plt.figure(figsize=(11, 9))
                self.heatmap_3 = seaborn.heatmap(Ect_df, xticklabels=3, yticklabels=3)
                plt.xlabel('Initial Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.ylabel('Final Energy Value (eV)', fontsize=17, fontweight='medium')
                if fit_opticalPeak:
                    plt.title('Optical Peak Energy (eV)', fontsize=17, fontweight='medium')
                else:
                    plt.title('CT State Energy (eV)', fontsize=17, fontweight='medium')
                cbar = self.heatmap_3.collections[0].colorbar
                cbar.ax.tick_params(labelsize=15)
                plt.yticks(rotation=360)
                plt.xticks(rotation=90)
                plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
                plt.show()

                plt.ion()
                plt.figure(figsize=(11, 9))
                self.heatmap_4 = seaborn.heatmap(R_df, xticklabels=3, yticklabels=3)
                plt.xlabel('Initial Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.ylabel('Final Energy Value (eV)', fontsize=17, fontweight='medium')
                plt.title('$\mathregular{R^{2}}$', fontsize=17, fontweight='medium')
                cbar = self.heatmap_4.collections[0].colorbar
                cbar.ax.tick_params(labelsize=15)
                plt.yticks(rotation=360)
                plt.xticks(rotation=90)
                plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
                plt.show()

            else:
                print('No fits determined.')

    # -----------------------------------------------------------------------------------------------------------

    # Gaussian function

    def gaussian(self, E, f, l, Ect):
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """
        return (f / (E * math.sqrt(4 * math.pi * l * self.T_CT * self.k))) * exp(
            -(Ect + l - E) ** 2 / (4 * l * self.k * self.T_CT))

    # Gaussian function including disorder

    def gaussian_disorder(self, E, f, l, Ect):
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """
        return (f / (E * math.sqrt(4 * math.pi * l * self.T_CT * self.k + 2 * self.sig ** 2))) * exp(
            -(Ect + l - E) ** 2 / (4 * l * self.k * self.T_CT + 2 * self.sig ** 2))

    # -----------------------------------------------------------------------------------------------------------

    # MLJ function

    def MLJ_gaussian(self, E, f, l_o, Ect):  # Double check if this equation is correct
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """
        EQE = 0
        for n in range(0, 6):
            EQE_n = (f / (E * math.sqrt(4 * math.pi * l_o * self.T_x * self.k))) \
                    * (math.exp(-self.S_i) * self.S_i ** n / math.factorial(n)) \
                    * exp(-(Ect + l_o - E + n * self.hbarw_i) ** 2 \
                          / (4 * l_o * self.k * self.T_x))
            EQE += EQE_n
        return EQE

    # MLJ function including disorder

    def MLJ_gaussian_disorder(self, E, f, l_o, Ect):  # Double check if this equation is correct
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """
        EQE = 0
        for n in range(0, 6):
            EQE_n = (f / (E * math.sqrt(4 * math.pi * l_o * self.T_x * self.k + 2 * self.sig_x ** 2))) \
                    * (math.exp(-self.S_i) * self.S_i ** n / math.factorial(n)) \
                    * exp(-(Ect + l_o - E + n * self.hbarw_i) ** 2 \
                          / (4 * l_o * self.k * self.T_x + 2 * self.sig_x ** 2))
            EQE += EQE_n
        return EQE

    # -----------------------------------------------------------------------------------------------------------

    # Wrapper function to perform curve fit

    def fit_function(self, function, energy_fit, eqe_fit, p0=None):
        """
        :param function: function to fit against (i.e. gaussian, gaussian_disorder etc.)
        :param energy_fit: energy values to fit against (x values)
        :param eqe_fit: EQE values to fit against (y values
        :return: best_vals: list of best fit parameters,
                 covar: covariance matrix of fit
                 y_fit: calculated EQE values of the fit,
                 r_squared: R^2 of the fit
        """
        best_vals, covar = curve_fit(function, energy_fit, eqe_fit, p0=p0)
        y_fit = [function(x, best_vals[0], best_vals[1], best_vals[2]) for x in energy_fit]
        r_squared = R_squared(eqe_fit, y_fit)

        return best_vals, covar, y_fit, r_squared

    # -----------------------------------------------------------------------------------------------------------

    # Page 5 - Fit EL and EQE

    # -----------------------------------------------------------------------------------------------------------

    # Function to scale and reduce EL and EQE

    def pre_plot_EL_EQE(self, data_df, startE, stopE, data_no, fit=False):

        self.T_EL = self.ui.EL_Temperature.value()

        startE = startE.value()
        stopE = stopE.value()

        if self.do_plot_EL:
            self.axEL_1, self.axEL_2 = set_up_EL_plot()
            self.do_plot_EL = False

        if data_no < 2:  # EL data

            Energy = []
            scaleFactor = self.ui.scalePlot.value()

            if len(data_df) != 0:  # Check that file is non-empty

                for y in range(len(data_df['Wavelength'])):  # Calculate energy values

                    Energy_val = (self.h * self.c) / (data_df['Wavelength'][y] * math.pow(10, -9) * self.q)
                    Energy.append(Energy_val)

                data_df['Energy'] = Energy

                if Data_is_valid(data_df, startE, stopE) and StartStop_is_valid(startE, stopE):

                    EL_wave, EL_energy, EL_signal = compile_EL(data_df, startE, stopE, 1)
                    red_EL_scaled = [EL_signal[x] / (scaleFactor * EL_energy[x]) for x in
                                     range(len(EL_signal))]  # Divide by energy to reduce

                    if data_no == 0:  # EL Data

                        if not fit:
                            label_ = pick_EQE_Label(self.ui.textBox_EL2, self.ui.textBox_EL1)
                            # color_ = pick_EQE_Color(self.ui.textBox_EL3, 100) # not currently used
                            color_ = '#1f77b4'  # Blue
                            plot(self.axEL_1, self.axEL_2, EL_energy, red_EL_scaled, label_, color_)

                        elif fit:
                            self.fit_EL_EQE(EL_energy, red_EL_scaled, self.ui.startFit_EL1, self.ui.stopFit_EL1, 0)

                    elif data_no == 1:  # Abs Data

                        bb_dict = bb_spectrum(EL_energy, self.T_EL)
                        EL_abs = [EL_signal[x] / bb_dict[EL_energy[x]] for x in range(len(EL_energy))]
                        red_EL_abs_scaled = [EL_abs[x] / (scaleFactor * EL_energy[x]) for x in range(len(EL_abs))]

                        if not fit:
                            # label_ = pick_EQE_Label(self.ui.textBox_EL2, self.ui.textBox_EL1)
                            # color_ = pick_EQE_Color(self.ui.textBox_EL3, 100) # not currently used
                            label_ = '$\mathregular{Red. EQE_{cal}}$'
                            color_ = '#ff7716'  # Orange
                            plot(self.axEL_1, self.axEL_2, EL_energy, red_EL_abs_scaled, label_, color_)

                        elif fit:
                            self.fit_EL_EQE(EL_energy, red_EL_abs_scaled, self.ui.startFit_EL2, self.ui.stopFit_EL2, 1)

            else:
                print('Please select a valid EL file.')

        elif data_no == 2:  # EQE Data

            if Data_is_valid(data_df, startE, stopE) and StartStop_is_valid(startE, stopE):

                self.Red_EQE_meas = pd.DataFrame()  # For determining the intersect between abs and emission
                EQE_wave, EQE_energy, EQE, EQE_log = compile_EQE(data_df, startE, stopE, 1)
                red_EQE = [EQE[x] * EQE_energy[x] for x in range(len(EQE))]

                if not fit:
                    label_ = pick_EQE_Label(self.ui.textBox_EL5, self.ui.textBox_EL4)
                    # color_ = pick_EQE_Color(self.ui.textBox_EL6, 100)
                    color_ = '#000000'  # Black
                    plot(self.axEL_1, self.axEL_2, EQE_energy, red_EQE, label_, color_)

                elif fit:
                    self.fit_EL_EQE(EQE_energy, red_EQE, self.ui.startFit_EQE, self.ui.stopFit_EQE, 1)

    # -----------------------------------------------------------------------------------------------------------

    # Function to fit reduced EL and EQE

    def fit_EL_EQE(self, energy, y, startE, stopE, data_no):

        include_Disorder = False

        startFit = startE.value()
        stopFit = stopE.value()

        df = pd.DataFrame()
        df['Energy'] = energy

        if Data_is_valid(df, startFit, stopFit):

            energy_fit, y_fit = compile_Data(energy, y, startFit, stopFit)

            diff = stopFit - startFit  # Find difference between start and stop fit energy
            x_gaussian = linspace(startFit, stopFit + 0.5 * diff, 50)  # Create x values to plot fit
            y_gaussian = []

            if self.ui.static_Disorder_EL.isChecked():
                include_Disorder = True
                self.sig_EL = self.ui.EL_Disorder.value()

            try:

                if self.ui.Gaussian_EL_EQE.isChecked():  # Marcus Theory Fitting

                    if data_no == 0:  # EL Data
                        if include_Disorder:  # Include Disorder in Fit
                            # y_fit_smooth = savgol_filter(y_fit, 51, 3) # In case you need to smooth the data
                            # y_fit_smooth = [x for x in y_fit_smooth]
                            # log_y_fit = [math.log(x) for x in y_fit]
                            # plot(self.axEL_1, self.axEL_2, energy_fit, y_fit_smooth, 'Smoothed Data', '#330000')
                            best_vals, covar = curve_fit(self.gaussian_EL_disorder, energy_fit, y_fit,
                                                         p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(
                                    self.gaussian_EL_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                        else:  # Without Disorder
                            best_vals, covar = curve_fit(self.gaussian_EL, energy_fit, y_fit,
                                                         p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(self.gaussian_EL(value, best_vals[0], best_vals[1], best_vals[2]))

                    elif data_no == 1:  # EQE / Abs Data
                        x_gaussian = linspace(startFit, stopFit + 2 * diff, 50)
                        if include_Disorder:
                            best_vals, covar = curve_fit(self.gaussian_EQE_disorder, energy_fit, y_fit,
                                                         p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(
                                    self.gaussian_EQE_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                        else:
                            best_vals, covar = curve_fit(self.gaussian_EQE, energy_fit, y_fit,
                                                         p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(self.gaussian_EQE(value, best_vals[0], best_vals[1], best_vals[2]))

                    print('-' * 80)
                    print('Temperature [T] (K): ', self.T_EL)
                    print('Oscillator Strength [f] (eV**2) : ', format(best_vals[0], '.6f'), '+/-',
                          format(math.sqrt(covar[0, 0]), '.6f'))
                    print('Reorganization Energy [l] (eV) : ', format(best_vals[1], '.6f'), '+/-',
                          format(math.sqrt(covar[1, 1]), '.6f'))
                    print('CT State Energy [ECT] (eV) : ', format(best_vals[2], '.6f'), '+/-',
                          format(math.sqrt(covar[2, 2]), '.6f'))
                    print('-' * 80)

                    if include_Disorder:
                        self.axEL_1.plot(x_gaussian, y_gaussian, linewidth=2, label='Gaussian Fit + Disorder',
                                         color='#000000', linestyle='--')
                        self.axEL_2.plot(x_gaussian, y_gaussian, linewidth=2, label='Gaussian Fit + Disorder',
                                         color='#000000', linestyle='--')
                    else:
                        self.axEL_1.plot(x_gaussian, y_gaussian, linewidth=2, label='Gaussian Fit', color='#000000',
                                         linestyle='--')
                        self.axEL_2.plot(x_gaussian, y_gaussian, linewidth=2, label='Gaussian Fit', color='#000000',
                                         linestyle='--')
                    plt.legend()
                    plt.draw()

                elif self.ui.MLJ_Gaussian_EL_EQE.isChecked():  # MLJ Theory Fitting

                    self.S_i_EL = self.ui.EL_Huang_Rhys.value()
                    self.hbarw_i_EL = self.ui.EL_vib_Energy.value()

                    x_gaussian = linspace(1.18, stopFit + 0.5 * diff, 50)

                    if data_no == 0:  # EL Data
                        if include_Disorder:
                            best_vals, covar = curve_fit(self.MLJ_gaussian_EL_disorder, energy_fit, y_fit,
                                                         p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(
                                    self.MLJ_gaussian_EL_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                        else:
                            best_vals, covar = curve_fit(self.MLJ_gaussian_EL, energy_fit, y_fit,
                                                         p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(self.MLJ_gaussian_EL(value, best_vals[0], best_vals[1], best_vals[2]))

                    elif data_no == 1:  # EQE / Abs Data
                        x_gaussian = linspace(startFit, stopFit + 2 * diff, 50)
                        if include_Disorder:
                            best_vals, covar = curve_fit(self.MLJ_gaussian_EQE_disorder, energy_fit, y_fit,
                                                         p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(
                                    self.MLJ_gaussian_EQE_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                        else:
                            best_vals, covar = curve_fit(self.MLJ_gaussian_EQE, energy_fit, y_fit,
                                                         p0=[0.001, 0.1, self.ui.EL_CT_State.value()])
                            for value in x_gaussian:
                                y_gaussian.append(
                                    self.MLJ_gaussian_EQE(value, best_vals[0], best_vals[1], best_vals[2]))

                    print('-' * 80)
                    print('Temperature [T] (K): ', self.T_EL)
                    print('Oscillator Strength [f] (eV**2) : ', format(best_vals[0], '.6f'), '+/-',
                          format(math.sqrt(covar[0, 0]), '.6f'))
                    print('Reorganization Energy [l] (eV) : ', format(best_vals[1], '.6f'), '+/-',
                          format(math.sqrt(covar[1, 1]), '.6f'))
                    print('CT State Energy [ECT] (eV) : ', format(best_vals[2], '.6f'), '+/-',
                          format(math.sqrt(covar[2, 2]), '.6f'))
                    print('-' * 80)

                    if include_Disorder:
                        self.axEL_1.plot(x_gaussian, y_gaussian, linewidth=2, label='MLJ Fit + Disorder',
                                         color='#000000', linestyle='--')
                        self.axEL_2.plot(x_gaussian, y_gaussian, linewidth=2, label='MLJ Fit + Disorder',
                                         color='#000000', linestyle='--')
                    else:
                        self.axEL_1.plot(x_gaussian, y_gaussian, linewidth=2, label='MLJ Fit', color='#000000',
                                         linestyle='--')
                        self.axEL_2.plot(x_gaussian, y_gaussian, linewidth=2, label='MLJ Fit', color='#000000',
                                         linestyle='--')
                    plt.legend()
                    plt.draw()

            except:
                print('Optimal parameters not found.')

    # -----------------------------------------------------------------------------------------------------------

    # Gaussian function for reduced EL

    def gaussian_EL(self, E, f, l, Ect):
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: reduced EL value
        """

        return (f / (math.sqrt(4 * math.pi * l * self.T_EL * self.k))) * exp(
            -(Ect - l - E) ** 2 / (4 * l * self.k * self.T_EL))

    # Gaussian function for reduced EL including disorder

    def gaussian_EL_disorder(self, E, f, l, Ect):
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: reduced EL value
        """

        return (f / (math.sqrt(4 * math.pi * l * self.T_EL * self.k + 2 * self.sig_EL ** 2))) * exp(
            -(Ect - l - E) ** 2 / (4 * l * self.k * self.T_EL + 2 * self.sig_EL ** 2))

    # -----------------------------------------------------------------------------------------------------------

    # MLJ function for reduced EL

    def MLJ_gaussian_EL(self, E, f, l_o, Ect):  # Double check if this equation is correct
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EL value
        """
        EL = 0
        for n in range(0, 6):
            EL_n = (f / (math.sqrt(4 * math.pi * l_o * self.T_EL * self.k))) \
                   * (math.exp(-self.S_i_EL) * self.S_i_EL ** n / math.factorial(n)) \
                   * exp(-(Ect - E - l_o - n * self.hbarw_i_EL) ** 2 \
                         / (4 * l_o * self.k * self.T_EL))
            EL += EL_n
        return EL

    # MLJ function for reduced EL including disorder

    def MLJ_gaussian_EL_disorder(self, E, f, l_o, Ect):  # Double check if this equation is correct
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """

        EL = 0
        for n in range(0, 6):
            EL_n = (f / (math.sqrt(4 * math.pi * l_o * self.T_EL * self.k + 2 * self.sig_EL ** 2))) \
                   * (math.exp(-self.S_i_EL) * self.S_i_EL ** n / math.factorial(n)) \
                   * exp(-(Ect - E - l_o - n * self.hbarw_i_EL) ** 2 \
                         / (4 * l_o * self.k * self.T_EL + 2 * self.sig_EL ** 2))
            EL += EL_n
        return EL

    # -----------------------------------------------------------------------------------------------------------

    # Gaussian function for reduced EQE

    def gaussian_EQE(self, E, f, l, Ect):
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: reduced EL value
        """

        return (f / (math.sqrt(4 * math.pi * l * self.T_EL * self.k))) * exp(
            -(Ect + l - E) ** 2 / (4 * l * self.k * self.T_EL))

    # Gaussian function for reduced EQE including disorder

    def gaussian_EQE_disorder(self, E, f, l, Ect):
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: reduced EL value
        """

        return (f / (math.sqrt(4 * math.pi * l * self.T_EL * self.k + 2 * self.sig_EL ** 2))) * exp(
            -(Ect + l - E) ** 2 / (4 * l * self.k * self.T_EL + 2 * self.sig_EL ** 2))

    # -----------------------------------------------------------------------------------------------------------

    # MLJ function for reduced EQE

    def MLJ_gaussian_EQE(self, E, f, l_o, Ect):  # Double check if this equation is correct
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """
        EQE = 0
        for n in range(0, 6):
            EQE_n = (f / (math.sqrt(4 * math.pi * l_o * self.T_EL * self.k))) \
                    * (math.exp(-self.S_i_EL) * self.S_i_EL ** n / math.factorial(n)) \
                    * exp(-(Ect - E + l_o + n * self.hbarw_i_EL) ** 2 \
                          / (4 * l_o * self.k * self.T_EL))
            EQE += EQE_n
        return EQE

    # MLJ function for reduced EQE including disorder

    def MLJ_gaussian_EQE_disorder(self, E, f, l_o, Ect):  # Double check if this equation is correct
        """
        :param E: List of energy values
        :param f: Oscillator strength
        :param l_o: Reorganization Energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """

        EQE = 0
        for n in range(0, 6):
            EQE_n = (f / (math.sqrt(4 * math.pi * l_o * self.T_EL * self.k + 2 * self.sig_EL ** 2))) \
                    * (math.exp(-self.S_i_EL) * self.S_i_EL ** n / math.factorial(n)) \
                    * exp(-(Ect - E + l_o + n * self.hbarw_i_EL) ** 2 \
                          / (4 * l_o * self.k * self.T_EL + 2 * self.sig_EL ** 2))
            EQE += EQE_n
        return EQE

    # -----------------------------------------------------------------------------------------------------------

    # Page 6 - Subtract and Add Peak Fits

    # -----------------------------------------------------------------------------------------------------------

    # Function to subtract peak fit from EQE

    def subtract_Fit(self, data_Fit, data_EQE, label_Fit, label_EQE, color_Fit, color_EQE):

        E_fit = 0
        sub_EQE = []

        label_fit = pick_EQE_Label(label_Fit, self.ui.textBox_p6_1)
        label_eqe = pick_EQE_Label(label_EQE, self.ui.textBox_p6_4)

        color_fit = color_Fit.toPlainText()
        color_fit = color_fit.replace(" ", "")

        color_eqe = color_EQE.toPlainText()
        color_eqe = color_eqe.replace(" ", "")

        if len(color_fit) == 0 or not is_Colour(color_fit):
            color_fit = '#ff7716'

        if len(color_eqe) == 0 or not is_Colour(color_eqe):
            color_eqe = 'black'

        try:  # Check if fit for an optical peak was imported
            T_fit = data_Fit['Temperature'][0]
            f_fit = data_Fit['Oscillator Strength (eV**2)'][0]
            l_fit = data_Fit['Reorganization Energy (eV)'][0]
            E_fit = data_Fit['Optical Peak Energy (eV)'][0]
        except:
            try:  # Check if fit for a CT state was imported
                T_fit = data_Fit['Temperature'][0]
                f_fit = data_Fit['Oscillator Strength (eV**2)'][0]
                l_fit = data_Fit['Reorganization Energy (eV)'][0]
                E_fit = data_Fit['CT State Energy (eV)'][0]
            except:
                print('Please import a valid fit file.')

        if E_fit != 0:  # Only progress if a valid energy was imported

            for x in range(len(data_EQE['Energy'])):
                fit_value = calculate_gaussian_absorption(data_EQE['Energy'][x], f_fit, l_fit, E_fit, T_fit)
                sub_EQE.append(data_EQE['EQE'][x] - fit_value)

            # Save fit data
            if self.ui.save_subEQE.isChecked():

                sub_file = pd.DataFrame()
                sub_file['EQE'] = sub_EQE
                sub_file['Energy'] = data_EQE['Energy']
                sub_file['Log_EQE'] = np.log(sub_EQE)
                sub_file['Wavelength'] = data_EQE['Wavelength']

                save_sub_file = filedialog.asksaveasfilename()  # User to pick a folder & name to save data to
                save_sub_path, save_sub_filename = os.path.split(save_sub_file)
                if len(save_sub_path) != 0:  # Check if the user actually selected a path
                    os.chdir(save_sub_path)  # Change the working directory
                    sub_file.to_csv(save_sub_filename)  # Save data to csv
                    print('Saving fit data to: %s' % str(save_sub_file))
                    os.chdir(self.data_dir)  # Change the directory back

            self.axSub_1, self.axSub_2 = set_up_EQE_plot()

            self.axSub_1.plot(data_Fit['Energy'], data_Fit['Signal'], linewidth=2, linestyle='--', color=color_fit,
                              label=label_fit)
            self.axSub_1.plot(data_EQE['Energy'], data_EQE['EQE'], linewidth=2, linestyle='-', color=color_eqe,
                              label=label_eqe)
            self.axSub_1.plot(data_EQE['Energy'], sub_EQE, linewidth=2, linestyle='-', color='#1f77b4',
                              label='Subtracted EQE')
            self.axSub_1.legend()

            self.axSub_2.plot(data_Fit['Energy'], data_Fit['Signal'], linewidth=2, linestyle='--', color=color_fit,
                              label=label_fit)
            self.axSub_2.plot(data_EQE['Energy'], data_EQE['EQE'], linewidth=2, linestyle='-', color=color_eqe,
                              label=label_eqe)
            self.axSub_2.plot(data_EQE['Energy'], sub_EQE, linewidth=2, linestyle='-', color='#1f77b4',
                              label='Subtracted EQE')
            self.axSub_2.legend()

    # -----------------------------------------------------------------------------------------------------------

    # Function to add peak fits

    def add_Fits(self, data_OptFit, data_CTFit, data_EQE):

        add_Energy = []
        add_Fits = []

        label_OptFit = pick_EQE_Label(self.ui.textBox_p7_2, self.ui.textBox_p7_1)
        label_CTFit = pick_EQE_Label(self.ui.textBox_p7_5, self.ui.textBox_p7_4)
        label_EQE = pick_EQE_Label(self.ui.textBox_p7_8, self.ui.textBox_p7_7)

        color_OptFit = self.ui.textBox_p7_3.toPlainText()
        color_OptFit = color_OptFit.replace(" ", "")
        color_CTFit = self.ui.textBox_p7_6.toPlainText()
        color_CTFit = color_CTFit.replace(" ", "")
        color_EQE = self.ui.textBox_p7_9.toPlainText()
        color_EQE = color_EQE.replace(" ", "")

        if len(color_OptFit) == 0 or not is_Colour(color_OptFit):
            color_OptFit = '#ff7716'

        if len(color_CTFit) == 0 or not is_Colour(color_CTFit):
            color_CTFit = '#1f77b4'

        if len(color_EQE) == 0 or not is_Colour(color_EQE):
            color_EQE = 'black'

        try: # Check if fit for an optical peak was imported
            T_OptFit = data_OptFit['Temperature'][0]
            f_OptFit = data_OptFit['Oscillator Strength (eV**2)'][0]
            l_OptFit = data_OptFit['Reorganization Energy (eV)'][0]
            E_OptFit = data_OptFit['Optical Peak Energy (eV)'][0]
        except:
            print('No optical peak fit imported.')

        try: # Check if fit for a CT state fit was importated
            T_CTFit = data_CTFit['Temperature'][0]
            f_CTFit = data_CTFit['Oscillator Strength (eV**2)'][0]
            l_CTFit = data_CTFit['Reorganization Energy (eV)'][0]
            E_CTFit = data_CTFit['CT State Energy (eV)'][0]
        except:
            print('No CT state fit imported.')

        if E_OptFit != 0 and E_CTFit != 0:  # Only progress if a valid energy was imported

            if len(data_EQE) != 0:

                for x in range(len(data_EQE['Energy'])):
                    if data_EQE['Energy'][x] < max(data_OptFit['Energy']):
                        OptFit_value = calculate_gaussian_absorption(data_EQE['Energy'][x], f_OptFit, l_OptFit, E_OptFit,
                                                                     T_OptFit)
                        CTFit_value = calculate_gaussian_absorption(data_EQE['Energy'][x], f_CTFit, l_CTFit, E_CTFit,
                                                                    T_CTFit)
                        add_Energy.append(data_EQE['Energy'][x])
                        add_Fits.append(OptFit_value + CTFit_value)

                self.axAdd_1, self.axAdd_2 = set_up_EQE_plot()

                self.axAdd_1.plot(data_OptFit['Energy'], data_OptFit['Signal'], linewidth=2, linestyle='--',
                                  color=color_OptFit, label=label_OptFit)
                self.axAdd_1.plot(data_CTFit['Energy'], data_CTFit['Signal'], linewidth=2, linestyle='--',
                                  color=color_CTFit, label=label_CTFit)
                self.axAdd_1.plot(add_Energy, add_Fits, linewidth=2, linestyle='dotted', color='grey',
                                  label='Optical + CT Fit')
                self.axAdd_1.plot(data_EQE['Energy'], data_EQE['EQE'], linewidth=2, linestyle='-', color=color_EQE,
                                  label=label_EQE)
                self.axAdd_1.legend()

                self.axAdd_2.plot(data_OptFit['Energy'], data_OptFit['Signal'], linewidth=2, linestyle='--',
                                  color=color_OptFit, label=label_OptFit)
                self.axAdd_2.plot(data_CTFit['Energy'], data_CTFit['Signal'], linewidth=2, linestyle='--',
                                  color=color_CTFit, label=label_CTFit)
                self.axAdd_2.plot(add_Energy, add_Fits, linewidth=2, linestyle='dotted', color='grey',
                                  label='Optical + CT Fit')
                self.axAdd_2.plot(data_EQE['Energy'], data_EQE['EQE'], linewidth=2, linestyle='-', color=color_EQE,
                                  label=label_EQE)
                self.axAdd_2.legend()

            else:
                print('Please import a valid EQE file.')

    # -----------------------------------------------------------------------------------------------------------

    # Page 3 - Double Fit

    # -----------------------------------------------------------------------------------------------------------

    def pre_double_fit(self):

        # Import relevant parameters

        eqe = self.data_double
        self.T_double = self.ui.double_Temperature.value()

        startStart_Opt = float(self.ui.startStart_Opt.value())
        startStop_Opt = float(self.ui.startStop_Opt.value())
        stopStart_Opt = float(self.ui.stopStart_Opt.value())
        stopStop_Opt = float(self.ui.stopStop_Opt.value())

        startStart_CT = float(self.ui.startStart_CT.value())
        startStop_CT = float(self.ui.startStop_CT.value())
        stopStart_CT = float(self.ui.stopStart_CT.value())
        stopStop_CT = float(self.ui.stopStop_CT.value())

        startGuess_Opt = float(self.ui.guessStart_Opt.value())
        stopGuess_Opt = float(self.ui.guessStop_Opt.value())
        startGuess_CT = float(self.ui.guessStart_CT.value())
        stopGuess_CT = float(self.ui.guessStop_CT.value())

        # Check that all start and stop energies are valid

        startOpt_ok = StartStop_is_valid(startStart_Opt, startStop_Opt)
        stopOpt_ok = StartStop_is_valid(stopStart_Opt, stopStop_Opt)
        startStopOpt_ok = StartStop_is_valid(startStop_Opt, stopStart_Opt)

        startCT_ok = StartStop_is_valid(startStart_CT, startStop_CT)
        stopCT_ok = StartStop_is_valid(stopStart_CT, stopStop_CT)
        startStopCT_ok = StartStop_is_valid(startStop_CT, stopStart_CT)

        guessOpt_ok = StartStop_is_valid(startGuess_Opt, stopGuess_Opt)
        guessCT_ok = StartStop_is_valid(startGuess_CT, stopGuess_CT)

        # Compile all start / stop energies for Opt and CT fit

        if startOpt_ok and stopOpt_ok and startStopOpt_ok and guessOpt_ok and startCT_ok and stopCT_ok and startStopCT_ok and guessCT_ok:

            startRange_Opt = np.round(np.arange(startStart_Opt, startStop_Opt + 0.01, 0.01),
                                      3).tolist()  # Change increment to 0.05
            stopRange_Opt = np.round(np.arange(stopStart_Opt, stopStop_Opt + 0.01, 0.01), 3).tolist()

            startRange_CT = np.round(np.arange(startStart_CT, startStop_CT + 0.01, 0.01), 3).tolist()
            stopRange_CT = np.round(np.arange(stopStart_CT, stopStop_CT + 0.01, 0.01), 3).tolist()

            guessRange_Opt = np.round(np.arange(startGuess_Opt, stopGuess_Opt + 0.1, 0.05), 2).tolist()
            guessRange_CT = np.round(np.arange(startGuess_CT, stopGuess_CT + 0.1, 0.05), 2).tolist()

            # Compile a dataFrame with all combinations of start / stop values for Opt and CT fit

            startOpt_list = []
            stopOpt_list = []
            startCT_list = []
            stopCT_list = []

            Opt_df = pd.DataFrame()
            CT_df = pd.DataFrame()

            # print('Compiling fit ranges ...')

            for startOpt in startRange_Opt:
                for stopOpt in stopRange_Opt:
                    startOpt_list.append(startOpt)
                    stopOpt_list.append(stopOpt)

            for startCT in startRange_CT:
                for stopCT in stopRange_CT:
                    startCT_list.append(startCT)
                    stopCT_list.append(stopCT)

            Opt_df['Start'] = startOpt_list
            Opt_df['Stop'] = stopOpt_list

            CT_df['Start'] = startCT_list
            CT_df['Stop'] = stopCT_list

            f_Opt = []
            l_Opt = []
            E_Opt = []
            R_Opt = []

            f_CT = []
            l_CT = []
            E_CT = []
            R_CT = []

            if self.ui.subtract_DoubleFit.isChecked():

                sub_df = pd.DataFrame()

                sub_startOpt_list = []
                sub_stopOpt_list = []
                sub_startCT_list = []
                sub_stopCT_list = []

                print('Calculating fits ...')
                for x in tqdm(range(len(Opt_df))):
                    best_vals_opt, r_squared_opt = self.single_fit(eqe=eqe, startE=Opt_df['Start'][x],
                                                                   stopE=Opt_df['Stop'][x],
                                                                   guessRange=guessRange_Opt)

                    if r_squared_opt > 0:  # Check that the optical peak fit was successful
                        new_eqe = self.subtract_Opt(eqe, best_vals_opt)

                        for y in range(len(CT_df)):
                            best_vals_ct, r_squared_ct = self.single_fit(eqe=new_eqe, startE=CT_df['Start'][y],
                                                                         stopE=CT_df['Stop'][y],
                                                                         guessRange=guessRange_CT)

                            sub_startOpt_list.append(Opt_df['Start'][x])
                            sub_stopOpt_list.append(Opt_df['Stop'][x])

                            f_Opt.append(best_vals_opt[0])
                            l_Opt.append(best_vals_opt[1])
                            E_Opt.append(best_vals_opt[2])
                            R_Opt.append(r_squared_opt)

                            sub_startCT_list.append(CT_df['Start'][y])
                            sub_stopCT_list.append(CT_df['Stop'][y])

                            f_CT.append(best_vals_ct[0])
                            l_CT.append(best_vals_ct[1])
                            E_CT.append(best_vals_ct[2])
                            R_CT.append(r_squared_ct)

                    else:  # If optical peak fit was unsuccessful, add zeros to CT fit as well

                        for y in range(len(CT_df)):
                            sub_startOpt_list.append(Opt_df['Start'][x])
                            sub_stopOpt_list.append(Opt_df['Stop'][x])

                            f_Opt.append(best_vals_opt[0])
                            l_Opt.append(best_vals_opt[1])
                            E_Opt.append(best_vals_opt[2])
                            R_Opt.append(r_squared_opt)

                            sub_startCT_list.append(CT_df['Start'][y])
                            sub_stopCT_list.append(CT_df['Stop'][y])

                            f_CT.append(0)
                            l_CT.append(0)
                            E_CT.append(0)
                            R_CT.append(0)

                assert len(f_Opt) == len(f_CT)

                sub_df['Start_Opt'] = sub_startOpt_list
                sub_df['Stop_Opt'] = sub_stopOpt_list

                sub_df['f_Opt'] = f_Opt
                sub_df['l_Opt'] = l_Opt
                sub_df['E_Opt'] = E_Opt
                sub_df['R2_Opt'] = R_Opt

                sub_df['Start_CT'] = sub_startCT_list
                sub_df['Stop_CT'] = sub_stopCT_list

                sub_df['f_CT'] = f_CT
                sub_df['l_CT'] = l_CT
                sub_df['E_CT'] = E_CT
                sub_df['R2_CT'] = R_CT

                print(sub_df)

                total_df = self.add_sub_fits(eqe=eqe, sub_df=sub_df)

            else:

                print('Calculating optical peak fits ...')

                for x in tqdm(range(len(Opt_df))):
                    best_vals, r_squared = self.single_fit(eqe=eqe, startE=Opt_df['Start'][x],
                                                           stopE=Opt_df['Stop'][x],
                                                           guessRange=guessRange_Opt)
                    f_Opt.append(best_vals[0])
                    l_Opt.append(best_vals[1])
                    E_Opt.append(best_vals[2])
                    R_Opt.append(r_squared)

                print('Calculating CT state fits ...')

                for x in tqdm(range(len(CT_df))):
                    best_vals, r_squared = self.single_fit(eqe=eqe, startE=CT_df['Start'][x],
                                                           stopE=CT_df['Stop'][x],
                                                           guessRange=guessRange_CT)
                    f_CT.append(best_vals[0])
                    l_CT.append(best_vals[1])
                    E_CT.append(best_vals[2])
                    R_CT.append(r_squared)

                Opt_df['f'] = f_Opt
                Opt_df['l'] = l_Opt
                Opt_df['E'] = E_Opt
                Opt_df['R2'] = R_Opt

                CT_df['f'] = f_CT
                CT_df['l'] = l_CT
                CT_df['E'] = E_CT
                CT_df['R2'] = R_CT

                # Determine R Squared of double fit

                total_df = self.add_single_fits(eqe=eqe, Opt_df=Opt_df, CT_df=CT_df)

    def single_fit(self, eqe, startE, stopE, guessRange):

        if len(eqe) != 0:

            # include_Disorder = False

            y_fit = []

            wave_fit, energy_fit, eqe_fit, log_eqe_fit = compile_EQE(eqe, startE, stopE, 1)

            ### FIX and ADD Disorder ###
            # if self.ui.STATIC_DISORDER.isChecked():
            #     include_Disorder = True
            #     SIG = self.ui.DISORDER.value()

            # Attempt optical peak fit:
            p0 = None

            for E_guess in guessRange:
                try:
                    # if include_Disorder:
                    #     # Fit gaussian with disorder
                    # else:
                    #     # Fit gaussian without disorder

                    best_vals, covar, y_fit, r_squared = self.fit_function(self.gaussian_double, energy_fit,
                                                                           eqe_fit,
                                                                           p0=p0)
                    if r_squared > 0:
                        break
                except:
                    p0 = [0.001, 0.1, E_guess]
                    if E_guess == guessRange[-1]:
                        print('Optimal parameters not found.')

                        best_vals = [0, 0, 0]
                        r_squared = 0

            return best_vals, r_squared

    def add_single_fits(self, eqe, Opt_df, CT_df):

        total_df = pd.DataFrame()

        R_2 = []

        start_Opt = []
        stop_Opt = []

        start_CT = []
        stop_CT = []

        Energy_list = []
        EQE_list = []
        Opt_fit_list = []
        CT_fit_list = []
        Total_fit_list = []

        if len(Opt_df) != 0 and len(CT_df) != 0:

            print('Determining best fit combination ...')

            for x_opt in tqdm(range(len(Opt_df))):
                wave_data, energy_data, eqe_data, log_eqe_data = compile_EQE(eqe, min(eqe['Energy']),
                                                                             Opt_df['Stop'][x_opt], 1)

                Opt_fit = np.array([calculate_gaussian_absorption(e, Opt_df['f'][x_opt], Opt_df['l'][x_opt],
                                                                  Opt_df['E'][x_opt], self.T_double) for e in
                                    energy_data])

                for x_ct in range(len(CT_df)):
                    CT_fit = np.array([calculate_gaussian_absorption(e, CT_df['f'][x_ct], CT_df['l'][x_ct],
                                                                     CT_df['E'][x_ct], self.T_double) for e in
                                       energy_data])

                    total_Fit = Opt_fit + CT_fit

                    assert len(total_Fit) == len(CT_fit)
                    assert len(total_Fit) == len(Opt_fit)

                    total_R_Squared = R_squared(eqe_data, total_Fit.tolist())

                    start_Opt.append(Opt_df['Start'][x_opt])
                    stop_Opt.append(Opt_df['Stop'][x_opt])

                    start_CT.append(CT_df['Start'][x_ct])
                    stop_CT.append(CT_df['Stop'][x_ct])

                    R_2.append(total_R_Squared)

                    Opt_fit_list.append(Opt_fit)
                    CT_fit_list.append(CT_fit)
                    Energy_list.append(energy_data)
                    EQE_list.append(eqe_data)
                    Total_fit_list.append(total_Fit)

            total_df['Start_Opt'] = start_Opt
            total_df['Stop_Opt'] = stop_Opt
            total_df['Start_CT'] = start_CT
            total_df['Stop_CT'] = stop_CT
            total_df['R2'] = R_2

            total_df['Opt_Fit'] = Opt_fit_list
            total_df['CT_Fit'] = CT_fit_list
            total_df['Energy'] = Energy_list
            total_df['EQE'] = EQE_list
            total_df['Total_Fit'] = Total_fit_list

            max_index = total_df[total_df['R2'] == max(total_df['R2'])].index.values[0]
            print('Best fit:')
            print(total_df['R2'][max_index])
            print('Opt Fit : ', total_df['Start_Opt'][max_index], total_df['Stop_Opt'][max_index])
            print('CT Fit : ', total_df['Start_CT'][max_index], total_df['Stop_CT'][max_index])

            self.axDouble_1, self.axDouble_2 = set_up_EQE_plot()

            self.axDouble_1.plot(total_df['Energy'][max_index], total_df['Opt_Fit'][max_index], linewidth=2,
                                 linestyle='--', label='Optical Fit')
            self.axDouble_1.plot(total_df['Energy'][max_index], total_df['CT_Fit'][max_index], linewidth=2,
                                 linestyle='--', label='CT Fit')
            self.axDouble_1.plot(total_df['Energy'][max_index], total_df['EQE'][max_index], linewidth=2,
                                 linestyle='-',
                                 label='EQE')
            self.axDouble_1.plot(total_df['Energy'][max_index], total_df['Total_Fit'][max_index], linewidth=2,
                                 linestyle='dotted', color='grey', label='Optical + CT Fit')
            self.axDouble_1.legend()

            self.axDouble_2.plot(total_df['Energy'][max_index], total_df['Opt_Fit'][max_index], linewidth=2,
                                 linestyle='--', label='Optical Fit')
            self.axDouble_2.plot(total_df['Energy'][max_index], total_df['CT_Fit'][max_index], linewidth=2,
                                 linestyle='--', label='CT Fit')
            self.axDouble_2.plot(total_df['Energy'][max_index], total_df['EQE'][max_index], linewidth=2,
                                 linestyle='-',
                                 label='EQE')
            self.axDouble_2.plot(total_df['Energy'][max_index], total_df['Total_Fit'][max_index], linewidth=2,
                                 linestyle='dotted', color='grey', label='Optical + CT Fit')
            self.axDouble_2.legend()

            return total_df

    def add_sub_fits(self, eqe, sub_df):

        total_df = pd.DataFrame()

        R_2 = []

        start_Opt = []
        stop_Opt = []

        start_CT = []
        stop_CT = []

        Energy_list = []
        EQE_list = []
        Opt_fit_list = []
        CT_fit_list = []
        Total_fit_list = []

        if len(sub_df) != 0:

            print('Determining best fit combination ...')

            for x in tqdm(range(len(sub_df))):
                if sub_df['R2_Opt'][x] != 0:
                    wave_data, energy_data, eqe_data, log_eqe_data = compile_EQE(eqe, min(eqe['Energy']),
                                                                                 sub_df['Stop_Opt'][x], 1)

                    Opt_fit = np.array([calculate_gaussian_absorption(e, sub_df['f_Opt'][x], sub_df['l_Opt'][x],
                                                                      sub_df['E_Opt'][x], self.T_double) for e in
                                        energy_data])
                    CT_fit = np.array([calculate_gaussian_absorption(e, sub_df['f_CT'][x], sub_df['l_CT'][x],
                                                                     sub_df['E_CT'][x], self.T_double) for e in
                                       energy_data])

                    total_Fit = Opt_fit + CT_fit

                    assert len(total_Fit) == len(CT_fit)
                    assert len(total_Fit) == len(Opt_fit)

                    total_R_Squared = R_squared(eqe_data, total_Fit.tolist())

                    start_Opt.append(sub_df['Start_Opt'][x])
                    stop_Opt.append(sub_df['Stop_Opt'][x])

                    start_CT.append(sub_df['Start_CT'][x])
                    stop_CT.append(sub_df['Stop_CT'][x])

                    R_2.append(total_R_Squared)

                    Opt_fit_list.append(Opt_fit)
                    CT_fit_list.append(CT_fit)
                    Energy_list.append(energy_data)
                    EQE_list.append(eqe_data)
                    Total_fit_list.append(total_Fit)

            total_df['Start_Opt'] = start_Opt
            total_df['Stop_Opt'] = stop_Opt
            total_df['Start_CT'] = start_CT
            total_df['Stop_CT'] = stop_CT
            total_df['R2'] = R_2

            total_df['Opt_Fit'] = Opt_fit_list
            total_df['CT_Fit'] = CT_fit_list
            total_df['Energy'] = Energy_list
            total_df['EQE'] = EQE_list
            total_df['Total_Fit'] = Total_fit_list

            max_index = total_df[total_df['R2'] == max(total_df['R2'])].index.values[0]
            print('Best fit:')
            print(total_df['R2'][max_index])
            print('Opt Fit : ', total_df['Start_Opt'][max_index], total_df['Stop_Opt'][max_index])
            print('CT Fit : ', total_df['Start_CT'][max_index], total_df['Stop_CT'][max_index])

            self.axDouble_1, self.axDouble_2 = set_up_EQE_plot()

            self.axDouble_1.plot(total_df['Energy'][max_index], total_df['Opt_Fit'][max_index], linewidth=2,
                                 linestyle='--', label='Optical Fit')
            self.axDouble_1.plot(total_df['Energy'][max_index], total_df['CT_Fit'][max_index], linewidth=2,
                                 linestyle='--', label='CT Fit')
            self.axDouble_1.plot(eqe['Energy'], eqe['EQE'], linewidth=2, linestyle='-', label='EQE')
            self.axDouble_1.plot(total_df['Energy'][max_index], total_df['Total_Fit'][max_index], linewidth=2,
                                 linestyle='dotted', color='grey', label='Optical + CT Fit')
            self.axDouble_1.legend()

            self.axDouble_2.plot(total_df['Energy'][max_index], total_df['Opt_Fit'][max_index], linewidth=2,
                                 linestyle='--', label='Optical Fit')
            self.axDouble_2.plot(total_df['Energy'][max_index], total_df['CT_Fit'][max_index], linewidth=2,
                                 linestyle='--', label='CT Fit')
            self.axDouble_2.plot(eqe['Energy'], eqe['EQE'], linewidth=2, linestyle='-', label='EQE')
            self.axDouble_2.plot(total_df['Energy'][max_index], total_df['Total_Fit'][max_index], linewidth=2,
                                 linestyle='dotted', color='grey', label='Optical + CT Fit')
            self.axDouble_2.legend()

            return total_df

    def subtract_Opt(self, eqe, best_vals):

        eqe = eqe.copy()

        Opt_fit = np.array(
            [calculate_gaussian_absorption(e, best_vals[0], best_vals[1], best_vals[2], self.T_double) for e in
             eqe['Energy']])
        EQE_data = np.array(eqe['EQE'])

        subtracted_EQE = EQE_data - Opt_fit

        assert len(Opt_fit) == len(EQE_data)
        assert len(Opt_fit) == len(subtracted_EQE)

        eqe['EQE'] = subtracted_EQE

        return eqe

    # -----------------------------------------------------------------------------------------------------------

    # Gaussian fitting function for double fit

    def gaussian_double(self, x, f, l, E):
        """
        :param x: List of energy values
        :param f: Oscillator strength
        :param l: Reorganization Energy
        :param E: Peak Energy
        :return: EQE value
        """
        return (f / (x * math.sqrt(4 * math.pi * l * self.T_double * self.k))) * exp(
            -(E + l - x) ** 2 / (4 * l * self.k * self.T_double))

    # -----------------------------------------------------------------------------------------------------------

    # Funtions to clear plots

    # -----------------------------------------------------------------------------------------------------------

    # Function to clear "Calculate EQE" plot

    def clear_plot(self):

        plt.close()  # Close the current plot
        self.ax1, self.ax2 = set_up_plot()  # Set up a new plot, this is preferred over plt.clf() in case the plot window was closed

    # -----------------------------------------------------------------------------------------------------------

    # Function to clear EQE plot

    def clear_EQE_plot(self):

        plt.close()
        plt.close()
        self.axFit_1, self.axFit_2 = set_up_EQE_plot()

    # -----------------------------------------------------------------------------------------------------------

    # Function to clear EL plot

    def clear_EL_plot(self):

        plt.close()
        plt.close()
        self.axEl_1, self.axEL_2 = set_up_EL_plot()


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
