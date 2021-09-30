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
import warnings
from tkinter import filedialog
from lmfit import Model

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn
# for the gui
from PyQt5 import QtWidgets
from numpy import exp, linspace
from scipy.optimize import curve_fit
from tqdm import tqdm
from collections import defaultdict
from scipy.interpolate import interp1d

import sEQE_Analysis_template
from source.add_subtract import subtract_Opt
from source.compilation import compile_EQE, compile_EL, compile_Data
from source.electroluminescence import bb_spectrum
from source.gaussian import calculate_gaussian_absorption, calculate_gaussian_disorder_absorption, \
    calculate_combined_fit
from source.normalization import normalize_EQE
from source.plot import plot, set_up_plot, set_up_EQE_plot, set_up_EL_plot
from source.reference_correction import calculate_Power
from source.utils import interpolate, sep_list, get_logger
from source.utils_plot import is_Colour, pick_EQE_Color, pick_EQE_Label, pick_Label
from source.validity import Ref_Data_is_valid, EQE_is_valid, Data_is_valid, Normalization_is_valid, Fit_is_valid, \
    StartStop_is_valid
from source.utils_fit import guess_fit, fit_function, calculate_guess_fit, fit_model, fit_model_double, find_best_fit
from source.utils import R_squared

warnings.filterwarnings("ignore")
warnings.simplefilter('ignore', np.RankWarning)
warnings.simplefilter('ignore', np.ComplexWarning)
warnings.filterwarnings('ignore', "Intel MKL ERROR")


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):

        QtWidgets.QMainWindow.__init__(self)

        # Set up the user interface from Designer
        self.ui = sEQE_Analysis_template.Ui_MainWindow()
        self.ui.setupUi(self)

        # Tkinter
        root = tk.Tk()
        root.withdraw()

        # Logger
        self.logger = get_logger()

        ## Page 1 - Calculate EQE

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

        ## Page 2 - Plot EQE

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

        self.ref_label = '$C_{60}$'

        ## Page 3 - Fit EQE (Marcus Theory)

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

        ## Page 4 - Extended Fits (Marcus Theory)

        # Double Fits

        self.bias = False
        self.tolerance = None

        self.data_double = []

        # Handle Import Data Button
        self.ui.browseDoubleFitButton.clicked.connect(lambda: self.writeText(self.ui.textBox_dF1, 'double1'))

        # Handle Double Fit Button
        self.ui.doubleFitButton.clicked.connect(lambda: self.double_fit())

        # Simultaneous Peak Fitting

        self.bias_sim = False
        self.tolerance_sim = False

        self.data_sim = []

        # Handle Import Data Button
        self.ui.browseSimFitButton.clicked.connect(lambda: self.writeText(self.ui.textBox_simFit, 'sim'))

        # Handle Sim Fit Button
        self.ui.simDoubleFitButton_single.clicked.connect(lambda: self.sim_double_fit_single())
        self.ui.simDoubleFitButton.clicked.connect(lambda: self.sim_double_fit())

        ## Page 5 - Fit EQE (MLJ Theory)

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

        ## Page 6 - Fit EL and EQE

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

        ## Page 7 - Subtract and Add Peak Fits

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

        # NOTE: Modify path if switching from Linux to another operating system
        self.data_dir = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')

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

    def writeText(self,
                  text_Box,
                  textBox_no
                  ):
        """
        Function to load data and update text box in GUI
        :param text_Box: GUI text box to write filename into [ui object]
        :param textBox_no: Internal name of text box to specify which data variable to define [int or str]
        :return: None
        """

        os.chdir(self.data_dir)
        file_ = filedialog.askopenfilename()

        if len(file_) != 0:
            path_, filename_ = os.path.split(file_)

            text_Box.clear()  # Clear the text box in case sth has been uploaded already
            text_Box.insertPlainText(filename_)  # Insert filename into text box

            ## Page 1 - Calculate EQE

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

            ## Page 2 - Plot EQE

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

            ## Page 3 - Fit EQE (Marcus Theory)

            elif textBox_no == 'f1':
                self.data_fit_1 = pd.read_csv(file_)

            elif textBox_no == 'f4':
                self.data_fit_2 = pd.read_csv(file_)

            ## Page 4 - Extended Fits (Marcus Theory)

            # For Double Fits

            elif textBox_no == 'double1':
                self.data_double = pd.read_csv(file_)

            # For Simultaneous Fits

            elif textBox_no == 'sim':
                self.data_sim = pd.read_csv(file_)

            ## Page 5 - Fit EQE (MLJ Theory)

            elif textBox_no == 'xF1':
                self.data_xFit_1 = pd.read_csv(file_)

            ## Page 6 - Fit EL and EQE

            elif textBox_no == 'el1':
                self.EL = pd.read_table(file_, sep=',', index_col=0)

            elif textBox_no == 'el2':
                self.EL_EQE = pd.read_csv(file_)

            ## Page 7 - Subtract and Add Peak Fits

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

    def pre_EQE(self,
                ref_df,
                data_df,
                start,
                stop,
                range_no
                ):
        """
        Wrapper function to load variables and data for EQE calculation
        :param ref_df: Reference data [dataFrame]
        :param data_df: Measured signal [dataFrame]
        :param start: GUI field with start wavelength [ui object]
        :param stop: GUI field with stop wavelength [ui object]
        :param range_no: Number to specify which data range to compile [int]
        :return: None
        """

        startNM = start.value()
        stopNM = stop.value()

        if Ref_Data_is_valid(ref_df, data_df, startNM, stopNM, range_no):
            self.calculate_EQE(ref_df, data_df, startNM, stopNM, range_no)

    # -----------------------------------------------------------------------------------------------------------

    # Function to calculate EQE

    def calculate_EQE(self,
                      ref_df,
                      data_df,
                      startNM,
                      stopNM,
                      range_no
                      ):
        """
        Function to calculate EQE from signal and reference data
        :param ref_df: Reference data [dataFrame]
        :param data_df: Measured signal [dataFrame]
        :param startNM: Start wavelength [float]
        :param stopNM: Stop wavelength [float]
        :param range_no: Number to specify which data range to compile [int]
        :return: None
        """

        power_dict = {}
        Wavelength = []
        Energy = []
        EQE = []
        log_EQE = []

        if 'Power' not in ref_df.columns:

            self.logger.info('Calculating power values.')

            if range_no == 1:
                if self.ui.Range1_Si_button.isChecked() and not self.ui.Range1_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.Si_cal)
                    except:
                        self.logger.error('Please select a valid reference diode.')
                elif self.ui.Range1_InGaAs_button.isChecked() and not self.ui.Range1_Si_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.InGaAs_cal)
                    except:
                        self.logger.error('Please select a valid reference diode.')
                else:
                    self.logger.error('Please select a valid reference diode.')

            elif range_no == 2:
                if self.ui.Range2_Si_button.isChecked() and not self.ui.Range2_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.Si_cal)
                    except:
                        self.logger.error('Please select a valid reference diode.')
                elif self.ui.Range2_InGaAs_button.isChecked() and not self.ui.Range2_Si_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.InGaAs_cal)
                    except:
                        self.logger.error('Please select a valid reference diode.')
                else:
                    self.logger.error('Please select a valid reference diode.')

            elif range_no == 3:
                if self.ui.Range3_Si_button.isChecked() and not self.ui.Range3_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.Si_cal)
                    except:
                        self.logger.error('Please select a valid reference diode.')
                elif self.ui.Range3_InGaAs_button.isChecked() and not self.ui.Range3_Si_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.InGaAs_cal)
                    except:
                        self.logger.error('Please select a valid reference diode.')
                else:
                    self.logger.error('Please select a valid reference diode.')

            elif range_no == 4:
                if self.ui.Range4_Si_button.isChecked() and not self.ui.Range4_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.Si_cal)
                    except:
                        self.logger.error('Please select a valid reference diode.')
                elif self.ui.Range4_InGaAs_button.isChecked() and not self.ui.Range4_Si_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.InGaAs_cal)
                    except:
                        self.logger.error('Please select a valid reference diode.')
                else:
                    self.logger.error('Please select a valid reference diode.')

            elif range_no == 5:
                if self.ui.Range5_Si_button.isChecked() and not self.ui.Range5_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.Si_cal)
                    except:
                        self.logger.error('Please select a valid reference diode.')
                elif self.ui.Range5_InGaAs_button.isChecked() and not self.ui.Range5_Si_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.InGaAs_cal)
                    except:
                        self.logger.error('Please select a valid reference diode.')
                else:
                    self.logger.error('Please select a valid reference diode.')

            elif range_no == 6:
                if self.ui.Range6_Si_button.isChecked() and not self.ui.Range6_InGaAs_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.Si_cal)
                    except:
                        self.logger.error('Please select a valid reference diode.')
                elif self.ui.Range6_InGaAs_button.isChecked() and not self.ui.Range6_Si_button.isChecked():
                    try:
                        ref_df['Power'] = calculate_Power(ref_df, self.InGaAs_cal)
                    except:
                        self.logger.error('Please select a valid reference diode.')
                else:
                    self.logger.error('Please select a valid reference diode.')

        if 'Power' in ref_df.columns:  # Check if the power has been calculated already

            for x in range(len(ref_df['Wavelength'])):  # Iterate through columns of reference file
                power_dict[ref_df['Wavelength'][x]] = ref_df['Power'][
                    x]  # Add wavelength and corresponding power to dictionary

            for y in range(len(data_df['Wavelength'])):  # Iterate through columns of data file
                if startNM <= data_df['Wavelength'][y] <= stopNM:  # Calculate EQE if start <= wave <= stop, else ignore

                    if data_df['Wavelength'][y] in power_dict.keys():  # Check if data wavelength is in reference file
                        Wavelength.append(data_df['Wavelength'][y])
                        Energy_val = (self.h * self.c) / (
                                data_df['Wavelength'][y] * math.pow(10, -9) * self.q)  # Calculate energy in eV
                        Energy.append(Energy_val)
                        EQE_val = (data_df['Mean Current'][y] * Energy_val) / (
                            power_dict[data_df['Wavelength'][y]])  # Easier way to calculate EQE
                        # EQE_val = ((data_df['Mean Current'][y] * self.h * self.c) / (
                        #   data_df['Wavelength'][y] * math.pow(10,-9) * power_dict[data_df['Wavelength'][y]] * self.q))
                        # EQE_val = (100 * data_df['Mean Current'][y] * Energy_val) / (
                        #   power_dict[data_df['Wavelength'][y]]) # *100 to turn into percent
                        EQE.append(EQE_val)
                        log_EQE.append(math.log10(EQE_val))

                    else:  # If data wavelength is not in reference file
                        Wavelength.append(data_df['Wavelength'][y])
                        Energy_val = (self.h * self.c) / (data_df['Wavelength'][y] * math.pow(10, -9) * self.q)
                        Energy.append(Energy_val)
                        Power_int = interpolate(data_df['Wavelength'][y], ref_df['Wavelength'],
                                                ref_df['Power'])  # Interpolate power
                        EQE_int = (data_df['Mean Current'][y] * Energy_val) / (Power_int)
                        # EQE_int = ((data_df['Mean Current'][y] * self.h * self.c) / (
                        #         data_df['Wavelength'][y] * math.pow(10,-9) * Power_int * self.q))
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
            self.logger.error('Length mismatch.')

    # -----------------------------------------------------------------------------------------------------------

    # Function to export EQE

    def export_EQE(self):
        """
        Functionm to export EQE to csv file
        :return: None
        """

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

        if ok_1 and ok_2 and ok_3 and ok_4 and ok_5 and ok_6:  # Check if all operations are ok or if fields are empty

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
                self.logger.info('Saving data to: %s' % str(EQE_file))
                os.chdir(self.data_dir)  # Change the directory back

        self.export = False

    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------

    # Page 2 - Calculate EQE

    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------

    # Function to select EQE for plotting

    def pre_plot_EQE(self,
                     number
                     ):
        """
        Wrapper function to select EQE for plotting
        :param number: Number of the EQE file to plot [int]
        :return: None
        """

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

        if ok_EQE_1 and ok_EQE_2 and ok_EQE_3 and ok_EQE_4 and ok_EQE_5 and ok_EQE_6 and ok_EQE_7 and ok_EQE_8 and \
                ok_EQE_9 and ok_EQE_10:
            self.axEQE_1.legend()
            self.axEQE_2.legend()
            plt.show()
        else:
            plt.close()
            plt.close()

    # -----------------------------------------------------------------------------------------------------------

    # Function to plot EQE

    def plot_EQE(self,
                 eqe_df,
                 startNM,
                 stopNM,
                 filename_Box,
                 label_Box,
                 color_Box,
                 file_no,
                 number
                 ):
        """
        Function to plot EQE data
        :param eqe_df: EQE data [dataFrame]
        :param startNM: Start wavelength [float]
        :param stopNM: Stop wavelength [float]
        :param filename_Box: GUI text box with filename to use for plot labeling [ui object]
        :param label_Box: GUI text box with plot label [ui object]
        :param color_Box: GUI text box with plot color [ui object]
        :param file_no: Number of EQE file to plot [int]
        :param number: Number to specify whether to plot wavelength (0) or energy (1) [boolean int]
        :return: None
        """

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

            ls_ = '-'

            if label_ == self.ref_label:
                ls_ = '--'

            if number == 0:
                self.axEQE_1.plot(wave, eqe, linewidth=3, linestyle=ls_, label=label_, color=color_)
                self.axEQE_2.semilogy(wave, eqe, linewidth=3, linestyle=ls_, label=label_, color=color_)
            elif number == 1:
                self.axEQE_1.plot(energy, eqe, linewidth=3, linestyle=ls_, label=label_, color=color_)
                self.axEQE_2.semilogy(energy, eqe, linewidth=3, linestyle=ls_, label=label_, color=color_)

            return True

        else:
            return False

    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------

    # Page 3 - Fit EQE (Marcus Theory)
    # Page 5 - Fit EQE (MLJ Theory)

    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------

    # Function to select EQE to fit

    def pre_fit_EQE(self,
                    eqe_df,
                    startE,
                    stopE,
                    startFit,
                    stopFit,
                    startPlotFit,
                    stopPlotFit,
                    filename_Box,
                    label_Box,
                    color_Box,
                    file_no
                    ):
        """
        Wrapper function to fit EQE
        :param eqe_df: EQE data [dataFrame]
        :param startE: Start energy to plot EQE [float]
        :param stopE: Stop energy to plot EQE [float]
        :param startFit: Start energy to fit [float]
        :param stopFit: Stop energy to fit [float]
        :param startPlotFit: Start energy to plot fit [float]
        :param stopPlotFit: Stop energy to plot fit [float]
        :param filename_Box: GUI text box with filename to use for plot labeling [ui object]
        :param label_Box: GUI text box with plot label [ui object]
        :param color_Box: GUI text box with plot color [ui object]
        :param file_no: Number of EQE file to plot [int]
        :return: None
        """

        # Change this if you want to plot on the same graph 
        if self.fit_plot:
            self.axFit_1, self.axFit_2 = set_up_EQE_plot()  # Sets up plot of energy vs EQE / log(EQE)
            self.fit_plot = False

        ok_plot_Fit = self.plot_fit_EQE(eqe_df,
                                        startE,
                                        stopE,
                                        startFit,
                                        stopFit,
                                        startPlotFit,
                                        stopPlotFit,
                                        filename_Box,
                                        label_Box,
                                        color_Box,
                                        file_no
                                        )

        if ok_plot_Fit:
            self.axFit_1.legend()
            self.axFit_2.legend()
            plt.show()

    # -----------------------------------------------------------------------------------------------------------

    # Function to fit and plot EQE

    def plot_fit_EQE(self,
                     eqe_df,
                     startE,
                     stopE,
                     startFit,
                     stopFit,
                     startPlotFit,
                     stopPlotFit,
                     filename_Box,
                     label_Box,
                     color_Box,
                     file_no
                     ):
        """
        Function to fit and plot EQE
        :param eqe_df: EQE data [dataFrame]
        :param startE: Start energy to plot EQE [float]
        :param stopE: Stop energy to plot EQE [float]
        :param startFit: Start energy to fit [float]
        :param stopFit: Stop energy to fit [float]
        :param startPlotFit: Start energy to plot fit [float]
        :param stopPlotFit: Stop energy to plot fit [float]
        :param filename_Box: GUI text box with filename to use for plot labeling [ui object]
        :param label_Box: GUI text box with plot label [ui object]
        :param color_Box: GUI text box with plot color [ui object]
        :param file_no: Number of EQE file to plot [int]
        :return: True / False
        """

        include_Disorder = False
        fit_opticalPeak = False
        fit_ok = False

        startE = startE.value()  # Pick start energy
        stopE = stopE.value()  # Pick stop energy

        startFit = startFit.value()  # Pick start fit energy
        stopFit = stopFit.value()  # Pick stop fit energy

        startPlotFit = startPlotFit.value()
        stopPlotFit = stopPlotFit.value()

        if Fit_is_valid(eqe_df, startE, stopE, startFit, stopFit, file_no):  # Check if files are not empty and in range
            wave, energy, eqe, log_eqe = compile_EQE(eqe_df, startE, stopE, 1)  # Compile EQE file
            wave_fit, energy_fit, eqe_fit, log_eqe_fit = compile_EQE(eqe_df, startFit, stopFit, 1)  # Compile fit range

            label_ = pick_EQE_Label(label_Box, filename_Box)
            color_ = pick_EQE_Color(color_Box, file_no)

            # Fit EQE (Marcus Theory)
            if (str(file_no)).isnumeric():

                if file_no == 1:  # Data in first row

                    self.T_CT = self.ui.Temperature_1.value()
                    guessStart = self.ui.guessStart_1.value()
                    guessStop = self.ui.guessStop_1.value()
                    guessStart_sig = self.ui.guessStartSig_1.value()
                    guessStop_sig = self.ui.guessStopSig_1.value()
                    if self.ui.disorder_1.isChecked():
                        include_Disorder = True
                    if self.ui.OptButton_1.isChecked() and not self.ui.CTButton_1.isChecked():  # Check peak to fit
                        fit_opticalPeak = True
                    elif self.ui.OptButton_1.isChecked() and self.ui.CTButton_1.isChecked():
                        self.logger.error('Please select a valid peak to fit.')
                    elif not self.ui.OptButton_1.isChecked() and not self.ui.CTButton_1.isChecked():
                        self.logger.error('Please select a valid peak to fit.')

                elif file_no == 2:  # Data in second row
                    self.T_CT = self.ui.Temperature_2.value()
                    guessStart = self.ui.guessStart_2.value()
                    guessStop = self.ui.guessStop_2.value()
                    guessStart_sig = self.ui.guessStartSig_2.value()
                    guessStop_sig = self.ui.guessStopSig_2.value()
                    if self.ui.disorder_2.isChecked():
                        include_Disorder = True
                    if self.ui.OptButton_2.isChecked() and not self.ui.CTButton_2.isChecked():
                        fit_opticalPeak = True
                    elif self.ui.OptButton_2.isChecked() and self.ui.CTButton_2.isChecked():
                        self.logger.error('Please select a valid peak to fit.')
                    elif not self.ui.OptButton_2.isChecked() and not self.ui.CTButton_2.isChecked():
                        self.logger.error('Please select a valid peak to fit.')

                # Attempt peak fit:
                x_gaussian = linspace(startPlotFit, stopPlotFit, 50)

                ECT_guess = np.arange(guessStart, guessStop + 0.1, 0.05)  # Extract peak guess range
                Sig_guess = np.arange(guessStart_sig, guessStop_sig + 0.01, 0.01)  # Extract sigma guess range

                # Initialize parameters
                p0 = None

                if include_Disorder:
                    best_guess_df = pd.DataFrame()
                    p0_list = []
                    R2_list = []

                    # ECT_guess_list = []
                    # sig_guess_list = []
                    # f_list = []
                    # l_list = []
                    # ECT_list = []
                    # sig_list = []

                    # TODO: Replace with "guess_fit" function?
                    for ECT in ECT_guess:
                        for sig in Sig_guess:
                            try:
                                best_vals, covar, y_fit, r_squared = fit_model(self.gaussian_disorder,
                                                                               energy_fit,
                                                                               eqe_fit,
                                                                               p0=p0,
                                                                               include_disorder=True
                                                                               )
                                if r_squared > 0:
                                    p0_list.append(p0)
                                    R2_list.append(r_squared)

                                    # ECT_guess_list.append(ECT)
                                    # sig_guess_list.append(sig)
                                    # f_list.append(best_vals[0])
                                    # l_list.append(best_vals[1])
                                    # ECT_list.append(best_vals[2])
                                    # sig_list.append(best_vals[3])

                                else:
                                    raise Exception('Wrong fit determined.')
                                p0 = [0.001, 0.1, round(ECT, 3), round(sig, 3)] # NOTE: Modify guesses if fit unsuccessful
                            except:
                                p0 = [0.001, 0.1, round(ECT, 3), round(sig, 3)] # NOTE: Modify guesses if fit unsuccessful

                    best_guess_df['p0'] = p0_list
                    best_guess_df['R2'] = R2_list

                    # best_guess_df['ECT Guess'] = ECT_guess_list
                    # best_guess_df['Sigma Guess'] = sig_guess_list
                    # best_guess_df['f'] = f_list
                    # best_guess_df['l'] = l_list
                    # best_guess_df['ECT'] = ECT_list
                    # best_guess_df['Sigma'] = sig_list
                    # best_guess_df.to_csv('~/Desktop/Best_Guesses.csv')

                    best_R2 = max(best_guess_df['R2'])
                    best_p0 = best_guess_df['p0'][best_guess_df['R2'] == best_R2].values[0]  # Find best initial guess

                    # Determine fit values of fit with best intial guess
                    best_vals, covar, y_fit, r_squared = fit_model(self.gaussian_disorder,
                                                                   energy_fit,
                                                                   eqe_fit,
                                                                   p0=best_p0,
                                                                   include_disorder=True
                                                                   )
                    y_gaussian = [self.gaussian_disorder(value,
                                                         best_vals[0],
                                                         best_vals[1],
                                                         best_vals[2],
                                                         best_vals[3]
                                                         ) for value in x_gaussian]
                    fit_ok = True  # Accept fit

                else:  # if disorder is not to be included
                    # TODO: Replace with "guess_fit" function?
                    for ECT in ECT_guess:
                        y_gaussian = []
                        try:
                            best_vals, covar, y_fit, r_squared = fit_function(self.gaussian, energy_fit, eqe_fit, p0=p0)
                            if r_squared > 0:
                                y_gaussian = [self.gaussian(value,
                                                                best_vals[0],
                                                                best_vals[1],
                                                                best_vals[2]
                                                                ) for value in x_gaussian]
                                fit_ok = True
                                best_p0 = p0
                                break  # This breaks the for loop
                            else:
                                raise Exception('Wrong fit determined.')
                        except:
                            p0 = [0.001, 0.1, round(ECT, 3)] # NOTE: Modify guesses if fit unsuccessful

                if fit_ok:

                    self.logger.info('Fit Results: ')
                    print("")
                    print('Initial Guess (eV) : ', best_p0)

                    print('-' * 80)
                    print('Temperature [T] (K) : ', self.T_CT)
                    print('Oscillator Strength [f] (eV**2) : ', format(best_vals[0], '.6f'),
                          '+/-', format(math.sqrt(covar[0, 0]), '.6f'))
                    print('Reorganization Energy [l] (eV) : ', format(best_vals[1], '.6f'),
                          '+/-', format(math.sqrt(covar[1, 1]), '.6f'))

                    if fit_opticalPeak:
                        print('Optical Peak Energy [E_Opt] (eV) : ', format(best_vals[2], '.6f'),
                              '+/-', format(math.sqrt(covar[2, 2]), '.6f'))
                    else:
                        print('CT State Energy [ECT] (eV) : ', format(best_vals[2], '.6f'),
                              '+/-', format(math.sqrt(covar[2, 2]), '.6f'))

                    if include_Disorder:
                        print('Disorder (eV) : ', format(best_vals[3], '.6f'),
                              '+/-', format(math.sqrt(covar[3, 3]), '.6f'))

                    print('R_Squared : ', format(r_squared, '.6f'))
                    print('-' * 80)
                    print("")

                    # Plot EQE data and CT fit
                    if include_Disorder:
                        fit_label = 'Gaussian Fit + Disorder'
                    else:
                        fit_label = 'Gaussian Fit'

                    if fit_opticalPeak:
                        fit_linestyle = 'dotted'
                    else:
                        fit_linestyle = '--'

                    self.axFit_1.plot(energy,
                                      eqe,
                                      linewidth=3,
                                      label=label_,
                                      color=color_
                                      )
                    plt.draw()
                    self.axFit_1.plot(x_gaussian,
                                      y_gaussian,
                                      linewidth=2,
                                      label=fit_label,
                                      color='#000000',
                                      linestyle=fit_linestyle
                                      )

                    self.axFit_2.semilogy(energy,
                                          eqe,
                                          linewidth=3,
                                          label=label_,
                                          color=color_
                                          )
                    plt.draw()
                    self.axFit_2.plot(x_gaussian,
                                      y_gaussian,
                                      linewidth=2,
                                      label=fit_label,
                                      color='#000000',
                                      linestyle=fit_linestyle
                                      )

                    # Save fit data
                    if self.ui.save_gaussianFit.isChecked():
                        fit_file = pd.DataFrame()
                        fit_file['Energy'] = x_gaussian
                        fit_file['Signal'] = y_gaussian
                        fit_file['Temperature'] = self.T_CT
                        fit_file['Oscillator Strength (eV**2)'] = best_vals[0]
                        fit_file['Reorganization Energy (eV)'] = best_vals[1]

                        if fit_opticalPeak:
                            fit_file['Optical Peak Energy (eV)'] = best_vals[2]
                        else:
                            fit_file['CT State Energy (eV)'] = best_vals[2]

                        if include_Disorder:
                            fit_file['Sigma (eV)'] = best_vals[3]

                        save_fit_file = filedialog.asksaveasfilename()  # User to pick folder & name to save to
                        save_fit_path, save_fit_filename = os.path.split(save_fit_file)
                        if len(save_fit_path) != 0:  # Check if the user actually selected a path
                            os.chdir(save_fit_path)  # Change the working directory
                            fit_file.to_csv(save_fit_filename)  # Save data to csv
                            self.logger.info('Saving fit data to: %s' % str(save_fit_file))
                            os.chdir(self.data_dir)  # Change the directory back
                    return True

                else:
                    self.logger.info('Optimal parameters not found.')
                    return False

            # Fit EQE (MLJ Theory)
            elif file_no == 'x1':

                self.S_i = self.ui.Huang_Rhys.value()
                self.hbarw_i = self.ui.vib_Energy.value()
                self.T_x = self.ui.extra_Temperature.value()

                guessStart_CT = self.ui.extraGuessStart.value()
                guessStop_CT = self.ui.extraGuessStop.value()
                guessStart_sig = self.ui.extraGuessStart_sig.value()
                guessStop_sig = self.ui.extraGuessStop_sig.value()

                if self.ui.extra_static_Disorder.isChecked():
                    include_Disorder = True

                # Attempt peak fit:
                x_MLJ_theory = linspace(startPlotFit, stopPlotFit, 50)

                ECT_guess = np.arange(guessStart_CT, guessStop_CT + 0.1, 0.05)
                sig_guess = np.arange(guessStart_sig, guessStop_sig + 0.1, 0.2)

                # Initialize parameters
                r_squared = 0
                p0 = None
                bounds = None  # First try without bounds and initial guesses

                if include_Disorder:
                    # TODO: Replace with "guess_fit" function?
                    for ECT in ECT_guess:
                        for sig in sig_guess:
                            try:
                                best_vals, covar, y_fit, r_squared = fit_function(self.MLJ_gaussian_disorder,
                                                                                  energy_fit,
                                                                                  eqe_fit,
                                                                                  p0=p0,
                                                                                  bounds=bounds,
                                                                                  include_disorder=include_Disorder
                                                                                  )
                                y_MLJ_theory = [self.MLJ_gaussian_disorder(value,
                                                                           best_vals[0],
                                                                           best_vals[1],
                                                                           best_vals[2],
                                                                           best_vals[3]
                                                                           ) for value in x_MLJ_theory]
                            except:
                                p0 = [0.001, 0.1, round(ECT, 3), round(sig, 3)]  # NOTE: Modify guesses if fit unsuccessful
                                bounds = (0, [0.1, 0.4, 1.6, 0.2])  # NOTE: Modify bounds if fit unsuccessful
                else:
                    # TODO: Replace with "guess_fit" function?
                    for ECT in ECT_guess:
                        try:
                            best_vals, covar, y_fit, r_squared = fit_function(self.MLJ_gaussian,
                                                                              energy_fit,
                                                                              eqe_fit,
                                                                              p0=p0
                                                                              )
                            y_MLJ_theory = [self.MLJ_gaussian(value,
                                                              best_vals[0],
                                                              best_vals[1],
                                                              best_vals[2]
                                                              ) for value in x_MLJ_theory]
                        except:
                            p0 = [0.001, 0.1, round(ECT, 3)]  # NOTE: Modify guesses if fit unsuccessful

                if r_squared > 0:
                    self.logger.info('Fit Results: ')
                    print("")
                    print('Initial Guess (eV) : ', p0)

                    print('-' * 80)
                    print('Temperature [T] (K) : ', self.T_x)
                    print('Oscillator Strength [f] (eV**2) : ', format(best_vals[0], '.6f'), '+/-',
                          format(math.sqrt(covar[0, 0]), '.6f'))
                    print('Reorganization Energy [l] (eV) : ', format(best_vals[1], '.6f'), '+/-',
                          format(math.sqrt(covar[1, 1]), '.6f'))
                    print('CT State Energy [ECT] (eV) : ', format(best_vals[2], '.6f'), '+/-',
                          format(math.sqrt(covar[2, 2]), '.6f'))
                    if include_Disorder:
                        print('Disorder (eV) : ', format(best_vals[3], '.6f'),
                              '+/-', format(math.sqrt(covar[3, 3]), '.6f'))

                    print('R_Squared : ', format(r_squared, '.6f'))
                    print('-' * 80)
                    print("")

                    # Plot EQE data and CT fit
                    self.axFit_1.plot(energy,
                                      eqe,
                                      linewidth=3,
                                      label=label_,
                                      color=color_
                                      )
                    plt.draw()
                    if include_Disorder:
                        self.axFit_1.plot(x_MLJ_theory,
                                          y_MLJ_theory,
                                          linewidth=2,
                                          label='MLJ Fit + Disorder',
                                          color='#000000',
                                          linestyle='--'
                                          )
                    else:
                        self.axFit_1.plot(x_MLJ_theory,
                                          y_MLJ_theory,
                                          linewidth=2,
                                          label='MLJ Fit',
                                          color='#000000',
                                          linestyle='--'
                                          )
                    plt.draw()

                    self.axFit_2.semilogy(energy,
                                          eqe,
                                          linewidth=3,
                                          label=label_,
                                          color=color_
                                          )
                    plt.draw()
                    if include_Disorder:
                        self.axFit_2.plot(x_MLJ_theory,
                                          y_MLJ_theory,
                                          linewidth=2,
                                          label='MLJ Fit + Disorder',
                                          color='#000000',
                                          linestyle='--'
                                          )
                    else:
                        self.axFit_2.plot(x_MLJ_theory,
                                          y_MLJ_theory,
                                          linewidth=2,
                                          label='MLJ Fit',
                                          color='#000000',
                                          linestyle='--'
                                          )
                    plt.draw()

                    return True

                else:
                    self.logger.info('Optimal parameters not found.')
                    return False
        else:
            return False

    # -----------------------------------------------------------------------------------------------------------

    # Function to generate heat map of fitting values

    def heatMap(self,
                eqe_df,
                startStartE,
                startStopE,
                stopStartE,
                stopStopE,
                filename_Box,
                label_Box,
                color_Box,
                file_no
                ):
        """
        Function to generate heatmap for standard Marcus theory fitting
        :param eqe_df: EQE data [dataFrame]
        :param startStartE: Start of the start fit energy range [float]
        :param startStopE: Start of the stop fit energy range [float]
        :param stopStartE: Stop of the start fit energy range [float]
        :param stopStopE: Stop of the stop fit energy range [float]
        :param filename_Box: GUI text box with filename to use for plot labeling [ui object]
        :param label_Box: GUI text box with plot label [ui object]
        :param color_Box: GUI text box with plot color [ui object]
        :param file_no: Number of EQE file to perform heatmap fits for [int]
        :return: None
        """

        include_Disorder = False
        fit_opticalPeak = False

        startStartE = startStartE.value()
        startStopE = startStopE.value()
        stopStartE = stopStartE.value()
        stopStopE = stopStopE.value()

        startStart_ok = StartStop_is_valid(startStartE, startStopE)  # Check that start energy is lower than stop energy
        # startStop_ok = StartStop_is_valid(startStopE, stopStartE)
        stopStop_ok = StartStop_is_valid(stopStartE, stopStopE)

        if startStart_ok and stopStop_ok:  # If all operations are valid, proceed with heat map calculations

            startEnergies = []  # Create empty lists for start and stop energies
            stopEnergies = []

            start_df = []  # Create empty lists for fit data collection
            stop_df = []
            f_df = []
            l_df = []
            Ect_df = []
            R_df = []
            sig_df = []

            x = float(startStartE)
            y = float(stopStartE)

            while x <= float(startStopE):  # As long as start range start value is smaller than start range stop value
                startEnergies.append(round(x, 3))  # Round to three digits
                x += 0.005  # Adjust this number to change resolution

            while y <= float(stopStopE):
                stopEnergies.append(round(y, 3))
                y += 0.005

            for start in tqdm(startEnergies):  # Iterate through start energies
                for stop in tqdm(stopEnergies):  # Iterate through stop energies

                    wave_fit, energy_fit, eqe_fit, log_eqe_fit = compile_EQE(eqe_df, start, stop, 1)

                    # Fit EQE (Marcus Theory)
                    if (str(file_no)).isnumeric():

                        if file_no == 1:  # Range 1
                            self.T_CT = self.ui.Temperature_1.value()
                            guessStart = self.ui.guessStart_1.value()
                            guessStop = self.ui.guessStop_1.value()
                            guessStart_sig = self.ui.guessStartSig_1.value()
                            guessStop_sig = self.ui.guessStopSig_1.value()
                            if self.ui.disorder_1.isChecked():
                                include_Disorder = True
                            if self.ui.OptButton_1.isChecked() and not self.ui.CTButton_1.isChecked():  # Check peak to fit
                                fit_opticalPeak = True
                            elif self.ui.OptButton_1.isChecked() and self.ui.CTButton_1.isChecked():
                                self.logger.error('Please select a valid peak to fit.')
                            elif not self.ui.OptButton_1.isChecked() and not self.ui.CTButton_1.isChecked():
                                self.logger.error('Please select a valid peak to fit.')

                        elif file_no == 2:  # Range 2
                            self.T_CT = self.ui.Temperature_2.value()
                            guessStart = self.ui.guessStart_2.value()
                            guessStop = self.ui.guessStop_2.value()
                            guessStart_sig = self.ui.guessStartSig_2.value()
                            guessStop_sig = self.ui.guessStopSig_2.value()
                            if self.ui.disorder_2.isChecked():
                                include_Disorder = True
                            if self.ui.OptButton_2.isChecked() and not self.ui.CTButton_2.isChecked():
                                fit_opticalPeak = True
                            elif self.ui.OptButton_2.isChecked() and self.ui.CTButton_2.isChecked():
                                self.logger.error('Please select a valid peak to fit.')
                            elif not self.ui.OptButton_2.isChecked() and not self.ui.CTButton_2.isChecked():
                                self.logger.error('Please select a valid peak to fit.')

                        # Attempt peak fit:
                        ECT_guess = np.arange(guessStart, guessStop + 0.2, 0.2)  # Increased step to accelerate process
                        Sig_guess = np.arange(guessStart_sig, guessStop_sig + 0.01, 0.01)  # Extract sigma guess range

                        p0 = None

                        if include_Disorder:
                            best_guess_df = pd.DataFrame()
                            p0_list = []
                            R2_list = []
                            # TODO: Replace with "guess_fit" function?
                            for ECT in ECT_guess:
                                for sig in Sig_guess:
                                    try:
                                        best_vals, covar, y_fit, r_squared = fit_model(self.gaussian_disorder,
                                                                                       energy_fit,
                                                                                       eqe_fit,
                                                                                       p0=p0,
                                                                                       include_disorder=True
                                                                                       )
                                        if r_squared > 0:
                                            p0_list.append(p0)
                                            R2_list.append(r_squared)
                                        else:
                                            raise Exception('Wrong fit determined.')
                                        p0 = [0.001, 0.1, round(ECT, 3), round(sig, 3)]  # NOTE: Modify guesses if fit unsuccessful
                                    except:
                                        p0 = [0.001, 0.1, round(ECT, 3), round(sig, 3)]  # NOTE: Modify guesses if fit unsuccessful

                            if len(p0_list) != 0:  # Check that the list is not empty

                                best_guess_df['p0'] = p0_list
                                best_guess_df['R2'] = R2_list

                                best_R2 = max(best_guess_df['R2'])
                                best_p0 = best_guess_df['p0'][best_guess_df['R2'] == best_R2].values[
                                    0]  # Find best guess

                                # Determine fit values of fit with best intial guess
                                best_vals, covar, y_fit, r_squared = fit_model(self.gaussian_disorder,
                                                                               energy_fit,
                                                                               eqe_fit,
                                                                               p0=best_p0,
                                                                               include_disorder=True
                                                                               )
                                start_df.append(start)
                                stop_df.append(stop)
                                f_df.append(best_vals[0])
                                l_df.append(best_vals[1])
                                Ect_df.append(best_vals[2])
                                sig_df.append(best_vals[3])
                                R_df.append(r_squared)

                            # else:
                            #     self.logger.info('Optimal parameters not found.')

                        else:  # If disorder is not to be included
                            # TODO: Replace with "guess_fit" function?
                            for ECT in ECT_guess:
                                try:
                                    best_vals, covar, y_fit, r_squared = fit_function(self.gaussian, energy_fit,
                                                                                      eqe_fit, p0=p0)

                                    if r_squared > 0:
                                        start_df.append(start)
                                        stop_df.append(stop)
                                        f_df.append(best_vals[0])
                                        l_df.append(best_vals[1])
                                        Ect_df.append(best_vals[2])
                                        R_df.append(r_squared)

                                    break
                                except:
                                    p0 = [0.001, 0.1, round(ECT, 3)] # NOTE: Modify guesses if fit unsuccessful
                                    if ECT == ECT_guess[-1]:
                                        self.logger.info('Optimal parameters not found.')

                    # Fit EQE (MLJ Theory)
                    elif file_no == 'x1':

                        self.S_i = self.ui.Huang_Rhys.value()
                        self.hbarw_i = self.ui.vib_Energy.value()
                        self.T_x = self.ui.extra_Temperature.value()

                        guessStart = self.ui.extraGuessStart.value()
                        guessStop = self.ui.extraGuessStop.value()
                        guessStart_sig = self.ui.extraGuessStart_sig.value()
                        guessStop_sig = self.ui.extraGuessStop_sig.value()

                        if self.ui.extra_static_Disorder.isChecked():
                            include_Disorder = True

                        # Attempt peak fit:
                        ECT_guess = np.arange(guessStart, guessStop + 0.1, 0.05)
                        sig_guess = np.arange(guessStart_sig, guessStop_sig + 0.1, 0.2)

                        p0 = None
                        bounds = None
                        r_squared = 0

                        if include_Disorder:
                            # TODO: Replace with "guess_fit" function?
                            for ECT in ECT_guess:
                                for sig in sig_guess:
                                    try:
                                        best_vals, covar, y_fit, r_squared = fit_function(self.MLJ_gaussian_disorder,
                                                                                          energy_fit,
                                                                                          eqe_fit,
                                                                                          p0=p0,
                                                                                          bounds=bounds,
                                                                                          include_disorder=include_Disorder
                                                                                          )
                                        if r_squared > 0:
                                            start_df.append(start)
                                            stop_df.append(stop)
                                            f_df.append(best_vals[0])
                                            l_df.append(best_vals[1])
                                            Ect_df.append(best_vals[2])
                                            R_df.append(r_squared)
                                            sig_df.append(best_vals[3])
                                            break
                                    except:
                                        p0 = [0.001, 0.1, round(ECT, 3), round(sig, 3)]  # NOTE: Modify guesses if fit unsuccessful
                                        bounds = (0, [0.1, 0.4, 1.6, 0.2])  # NOTE: Modify bounds if fit unsuccessful
                                        if ECT == ECT_guess[-1]:
                                            self.logger.info('Optimal parameters not found.')
                                if r_squared > 0: # To break the second loop
                                    break
                        else:
                            # TODO: Replace with "guess_fit" function?
                            for ECT in ECT_guess:
                                try:
                                    best_vals, covar, y_fit, r_squared = fit_function(self.MLJ_gaussian,
                                                                                      energy_fit,
                                                                                      eqe_fit,
                                                                                      p0=p0
                                                                                      )
                                    if r_squared > 0:
                                        start_df.append(start)
                                        stop_df.append(stop)
                                        f_df.append(best_vals[0])
                                        l_df.append(best_vals[1])
                                        Ect_df.append(best_vals[2])
                                        R_df.append(r_squared)
                                        break
                                except:
                                    p0 = [0.001, 0.1, round(ECT, 3)]  # NOTE: Modify guesses if fit unsuccessful
                                    if ECT == ECT_guess[-1]:
                                        self.logger.info('Optimal parameters not found.')

            if len(R_df) != 0:  # Check that there are results to plot

                if include_Disorder:
                    parameter_df = pd.DataFrame(
                        {'Start': start_df, 'Stop': stop_df, 'f': f_df, 'l': l_df, 'Ect': Ect_df,
                         'R_Squared': R_df, 'Sig': sig_df})  # Create a dataFrame with all results
                else:
                    parameter_df = pd.DataFrame(
                        {'Start': start_df, 'Stop': stop_df, 'f': f_df, 'l': l_df, 'Ect': Ect_df,
                         'R_Squared': R_df})  # Create a dataFrame with all results
                max_index = parameter_df[parameter_df['R_Squared'] == max(parameter_df['R_Squared'])].index.values[0]

                self.logger.info('Fit Results: ')
                print("")
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

                if include_Disorder:
                    print('Average Gaussian Disorder [Sig] (eV) : ', format(parameter_df['Sig'].mean(), '.6f'), '+/-',
                          format(parameter_df['Sig'].std(), '.6f'))

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

                    if include_Disorder:
                        print('Average Gaussian Disorder [Sig] (eV) : ',
                              format(parameter_df['Sig'][parameter_df['R_Squared'] == 1.0].mean(), '.6f'), '+/-',
                              format(parameter_df['Sig'][parameter_df['R_Squared'] == 1.0].std(), '.6f'))

                else:
                    print('Max R_squared : ', format(max(parameter_df['R_Squared']), '.6f'))
                    print('Start Energy (eV) : ', parameter_df['Start'][max_index])
                    print('Stop Energy (eV) : ', parameter_df['Stop'][max_index])
                print('-' * 80)
                print("")

                f_df = parameter_df.pivot('Stop', 'Start', 'f')  # Pivot the dataFrame: x-value = Stop, y-value = Start, value = f
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

                if include_Disorder:
                    Sig_df = parameter_df.pivot('Stop', 'Start', 'Sig')
                    plt.ion()
                    plt.figure(figsize=(11, 9))
                    self.heatmap_5 = seaborn.heatmap(Sig_df, xticklabels=3, yticklabels=3)
                    plt.xlabel('Initial Energy Value (eV)', fontsize=17, fontweight='medium')
                    plt.ylabel('Final Energy Value (eV)', fontsize=17, fontweight='medium')
                    plt.title('Gaussian Disorder (eV)', fontsize=17, fontweight='medium')
                    cbar = self.heatmap_5.collections[0].colorbar
                    cbar.ax.tick_params(labelsize=15)
                    plt.yticks(rotation=360)
                    plt.xticks(rotation=90)
                    plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
                    plt.show()

            else:
                self.logger.info('No fits determined.')

    # -----------------------------------------------------------------------------------------------------------

    # Gaussian function

    def gaussian(self, E, f, l, Ect):
        """
        Standard gaussian function
        :param E: List of energy values [list of floats]
        :param f: Oscillator strength [float]
        :param l: Reorganization Energy [float]
        :param Ect: Charge Transfer State Energy [float]
        :return: EQE value [float]
        """
        return (f / (E * math.sqrt(4 * math.pi * l * self.T_CT * self.k))) * exp(
            -(Ect + l - E) ** 2 / (4 * l * self.k * self.T_CT))

    # Gaussian function including disorder

    # # Old function that I believe is incorrect
    # def gaussian_disorder(self, E, f, l, Ect):
    #     """
    #     Standard gaussian function including disorder
    #     :param E: List of energy values [list of floats]
    #     :param f: Oscillator strength [float]
    #     :param l: Reorganization Energy [float]
    #     :param Ect: Charge Transfer State Energy [float]
    #     :return: EQE value [float]
    #     """
    #     return (f / (E * math.sqrt(4 * math.pi * l * self.T_CT * self.k + 2 * self.sig ** 2))) * exp(
    #         -(Ect + l - E) ** 2 / (4 * l * self.k * self.T_CT + 2 * self.sig ** 2))

    def gaussian_disorder(self, E, f, l, Ect, sig):
        """
        New gaussian function including disorder
        :param E: List of energy values [list of floats]
        :param f: Oscillator strength [float]
        :param l: Reorganization Energy [float]
        :param Ect: Charge Transfer State Energy [float]
        :param sig: Gaussian disorder [float]
        :return: EQE value [float]
        """
        return (f / (E * math.sqrt(2 * math.pi * (2 * l * self.T_CT * self.k + sig ** 2)))) * exp(
            -((Ect - (sig ** 2 / (2 * self.k * self.T_CT)) + l + (sig ** 2 / (2 * self.k * self.T_CT)) - E) ** 2 / (
                    4 * l * self.k * self.T_CT + 2 * sig ** 2)))

    # -----------------------------------------------------------------------------------------------------------

    # MLJ function

    def MLJ_gaussian(self, E, f, l_o, Ect):  # Double check if this equation is correct
        """
        MLJ function
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

    def MLJ_gaussian_disorder(self, E, f, l, Ect, sig):  # Double check if this equation is correct
        """
        MLJ function including disorder
        :param E: List of energy values
        :param f: Oscillator strength
        :param l: Low frequency reorganization energy
        :param Ect: Charge Transfer State Energy
        :return: EQE value
        """
        EQE = 0
        for n in range(0, 6):
            EQE_n = (f / (E * math.sqrt(2 * math.pi * (2 * l * self.T_x * self.k + sig ** 2))) \
                     * (math.exp(-self.S_i) * self.S_i ** n / math.factorial(n)) \
                     * exp(-(Ect + l - E + n * self.hbarw_i) ** 2 \
                           / (4 * l * self.k * self.T_x + 2 * sig ** 2)))
            EQE += EQE_n
        return EQE

    # -----------------------------------------------------------------------------------------------------------

    # Page 6 - Fit EL and EQE

    # -----------------------------------------------------------------------------------------------------------

    # Function to scale and reduce EL and EQE

    def pre_plot_EL_EQE(self,
                        data_df,
                        startE,
                        stopE,
                        data_no,
                        fit=False
                        ):
        """
        Wrapper function to plot and fit EQE and EL data
        :param data_df: EQE or EL data [dataFrame]
        :param startE: Start energy [float]
        :param stopE: Stop energy [float]
        :param data_no: Number of data file to plot and fit [int]
        :param fit: Boolean value to specify whether to fit [bool]
        :return: None
        """

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
                    red_EL = [EL_signal[x] / (EL_energy[x])
                              for x in range(len(EL_signal))]  # Divide by energy to reduce (checked in Benduhn thesis)
                    red_EL_scaled = [red_EL[x] / scaleFactor
                                     for x in range(len(red_EL))]

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

                        # TODO: Old code to be removed?
                        # EL_abs = [EL_signal[x] / bb_dict[EL_energy[x]] for x in range(len(EL_energy))]
                        # red_EL_abs_scaled = [EL_abs[x] / (scaleFactor * EL_energy[x]) for x in range(len(EL_abs))]

                        scaleFactor_calc = self.ui.scalePlot_calc.value()

                        EL_abs_scaled = [EL_signal[x] / (scaleFactor * bb_dict[EL_energy[x]])
                                         for x in range(len(EL_energy))]
                        EL_abs_scaled_test = [EL_signal[x] / (scaleFactor_calc * bb_dict[EL_energy[x]])
                                              for x in range(len(EL_energy))]  # test which scale factor is correct

                        red_EL_abs_scaled = [EL_energy[x] * EL_signal[x] / (scaleFactor * bb_dict[EL_energy[x]])
                                             for x in range(len(EL_energy))]
                        red_EL_abs_scaled_test = [
                            EL_energy[x] * EL_signal[x] / (scaleFactor_calc * bb_dict[EL_energy[x]])
                            for x in range(len(EL_energy))]  # test which scale factor is correct

                        if not fit:
                            # label_ = pick_EQE_Label(self.ui.textBox_EL2, self.ui.textBox_EL1)
                            # color_ = pick_EQE_Color(self.ui.textBox_EL3, 100) # not currently used
                            label_ = '$\mathregular{Red. EQE_{cal}}$'
                            color_ = '#ff7716'  # Orange

                            # Experiments
                            # plot(self.axEL_1, self.axEL_2, EL_energy, red_EL_abs_scaled, label_, color_)

                            # plot(self.axEL_1, self.axEL_2, EL_energy, EL_abs_scaled,
                            #      label_='EQE,cal with EL SF', color_='blue')
                            # plot(self.axEL_1, self.axEL_2, EL_energy, EL_abs_scaled_test,
                            #      label_='EQE,cal with EQE,calc SF', color_='green')
                            # plot(self.axEL_1, self.axEL_2, EL_energy, red_EL_abs_scaled,
                            #      label_='Red. EQE,cal with EL SF', color_='orange')
                            plot(self.axEL_1, self.axEL_2, EL_energy, red_EL_abs_scaled_test,
                                 label_=label_, color_=color_)

                        elif fit:
                            self.fit_EL_EQE(EL_energy, red_EL_abs_scaled_test, self.ui.startFit_EL2,
                                            self.ui.stopFit_EL2, 1)

            else:
                self.logger.error('Please select a valid EL file.')

        elif data_no == 2:  # EQE Data

            if Data_is_valid(data_df, startE, stopE) and StartStop_is_valid(startE, stopE):

                self.Red_EQE_meas = pd.DataFrame()  # For determining the intersect between abs and emission
                EQE_wave, EQE_energy, EQE, EQE_log = compile_EQE(data_df, startE, stopE, 1)
                red_EQE = [EQE[x] * EQE_energy[x] for x in
                           range(len(EQE))]  # Multiplication confirmed in Benduhn thesis

                if not fit:
                    label_ = pick_EQE_Label(self.ui.textBox_EL5, self.ui.textBox_EL4)
                    # color_ = pick_EQE_Color(self.ui.textBox_EL6, 100)
                    color_ = '#000000'  # Black
                    plot(self.axEL_1, self.axEL_2, EQE_energy, red_EQE, label_, color_)

                elif fit:
                    self.fit_EL_EQE(EQE_energy, red_EQE, self.ui.startFit_EQE, self.ui.stopFit_EQE, 1)

    # -----------------------------------------------------------------------------------------------------------

    # Function to fit reduced EL and EQE

    def fit_EL_EQE(self,
                   energy,
                   y,
                   startE,
                   stopE,
                   data_no
                   ):
        """
        Function to fit EL and EQE data
        :param energy: List of energy values to fit [list]
        :param y: List of data values to fit [list]
        :param startE: Start energy to fit [float]
        :param stopE: Stop energy to fit [float]
        :param data_no: Number of data file to fit [int]
        :return: None
        """

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
                        # TODO: Fix formatting here!
                        if include_Disorder:  # Include Disorder in Fit
                            # y_fit_smooth = savgol_filter(y_fit, 51, 3) # In case you need to smooth the data
                            # y_fit_smooth = [x for x in y_fit_smooth]
                            # log_y_fit = [math.log(x) for x in y_fit]
                            # plot(self.axEL_1, self.axEL_2, energy_fit, y_fit_smooth, 'Smoothed Data', '#330000')
                            best_vals, covar = curve_fit(self.gaussian_EL_disorder,
                                                         energy_fit,
                                                         y_fit,
                                                         p0=[0.001, 0.1, self.ui.EL_CT_State.value()]
                                                         )
                            for value in x_gaussian:
                                y_gaussian.append(
                                    self.gaussian_EL_disorder(value, best_vals[0], best_vals[1], best_vals[2]))
                        else:  # Without Disorder
                            best_vals, covar = curve_fit(self.gaussian_EL,
                                                         energy_fit,
                                                         y_fit,
                                                         p0=[0.001, 0.1, self.ui.EL_CT_State.value()]
                                                         )
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

                    self.logger.info('Fit Results: ')
                    print("")
                    print('-' * 80)
                    print('Temperature [T] (K): ', self.T_EL)
                    print('Oscillator Strength [f] (eV**2) : ', format(best_vals[0], '.6f'), '+/-',
                          format(math.sqrt(covar[0, 0]), '.6f'))
                    print('Reorganization Energy [l] (eV) : ', format(best_vals[1], '.6f'), '+/-',
                          format(math.sqrt(covar[1, 1]), '.6f'))
                    print('CT State Energy [ECT] (eV) : ', format(best_vals[2], '.6f'), '+/-',
                          format(math.sqrt(covar[2, 2]), '.6f'))
                    print('-' * 80)
                    print("")

                    if include_Disorder:
                        self.axEL_1.plot(x_gaussian,
                                         y_gaussian,
                                         linewidth=2,
                                         label='Gaussian Fit + Disorder',
                                         color='#000000',
                                         linestyle='--'
                                         )
                        self.axEL_2.plot(x_gaussian,
                                         y_gaussian,
                                         linewidth=2,
                                         label='Gaussian Fit + Disorder',
                                         color='#000000',
                                         linestyle='--'
                                         )
                    else:
                        self.axEL_1.plot(x_gaussian,
                                         y_gaussian,
                                         linewidth=2,
                                         label='Gaussian Fit',
                                         color='#000000',
                                         linestyle='--'
                                         )
                        self.axEL_2.plot(x_gaussian,
                                         y_gaussian,
                                         linewidth=2,
                                         label='Gaussian Fit',
                                         color='#000000',
                                         linestyle='--'
                                         )
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

                    self.logger.info('Fit Results: ')
                    print("")
                    print('-' * 80)
                    print('Temperature [T] (K): ', self.T_EL)
                    print('Oscillator Strength [f] (eV**2) : ', format(best_vals[0], '.6f'), '+/-',
                          format(math.sqrt(covar[0, 0]), '.6f'))
                    print('Reorganization Energy [l] (eV) : ', format(best_vals[1], '.6f'), '+/-',
                          format(math.sqrt(covar[1, 1]), '.6f'))
                    print('CT State Energy [ECT] (eV) : ', format(best_vals[2], '.6f'), '+/-',
                          format(math.sqrt(covar[2, 2]), '.6f'))
                    print('-' * 80)
                    print("")

                    if include_Disorder:
                        self.axEL_1.plot(x_gaussian,
                                         y_gaussian,
                                         linewidth=2,
                                         label='MLJ Fit + Disorder',
                                         color='#000000',
                                         linestyle='--'
                                         )
                        self.axEL_2.plot(x_gaussian,
                                         y_gaussian,
                                         linewidth=2,
                                         label='MLJ Fit + Disorder',
                                         color='#000000',
                                         linestyle='--'
                                         )
                    else:
                        self.axEL_1.plot(x_gaussian,
                                         y_gaussian,
                                         linewidth=2,
                                         label='MLJ Fit',
                                         color='#000000',
                                         linestyle='--'
                                         )
                        self.axEL_2.plot(x_gaussian,
                                         y_gaussian,
                                         linewidth=2,
                                         label='MLJ Fit',
                                         color='#000000',
                                         linestyle='--'
                                         )
                    plt.legend()
                    plt.draw()

            except:
                self.logger.info('Optimal parameters not found.')

    # -----------------------------------------------------------------------------------------------------------

    # Gaussian function for reduced EL

    def gaussian_EL(self, E, f, l, Ect):
        """
        Standard gaussian to fit EL
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
        Standard gaussian including disorder to fit EL
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
        MLJ function to fit EL
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
        MLJ function including disorder to fit EL
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
        Standard gaussian function to fit reduced EQE
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
        Standard gaussian function including disorder to fit reduced EQE
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
        MLJ function to fit reduced EQE
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
        MLJ function including disorder to fit reduced EQE
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

    # Page 7 - Subtract and Add Peak Fits

    # -----------------------------------------------------------------------------------------------------------

    # Function to subtract peak fit from EQE

    def subtract_Fit(self,
                     data_Fit,
                     data_EQE,
                     label_Fit,
                     label_EQE,
                     color_Fit,
                     color_EQE
                     ):
        """
        Function to subtract fit from EQE data
        :param data_Fit: Fit data to subtract [dataFrame]
        :param data_EQE: EQE data to subtract fit from [dataFrame]
        :param label_Fit: GUI text box with fit plot label [ui object]
        :param label_EQE: GUI text box with EQE plot label [ui object]
        :param color_Fit: GUI text box with fit plot color [ui object]
        :param color_EQE: GUI text box with EQE plot color [ui object]
        :return:
        """

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
                self.logger.error('Please import a valid fit file.')

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
                    self.logger.info('Saving fit data to: %s' % str(save_sub_file))
                    os.chdir(self.data_dir)  # Change the directory back

            self.axSub_1, self.axSub_2 = set_up_EQE_plot()

            self.axSub_1.plot(data_Fit['Energy'],
                              data_Fit['Signal'],
                              linewidth=2,
                              linestyle='--',
                              color=color_fit,
                              label=label_fit
                              )
            self.axSub_1.plot(data_EQE['Energy'],
                              data_EQE['EQE'],
                              linewidth=2,
                              linestyle='-',
                              color=color_eqe,
                              label=label_eqe
                              )
            self.axSub_1.plot(data_EQE['Energy'],
                              sub_EQE,
                              linewidth=2,
                              linestyle='-',
                              color='#1f77b4',
                              label='Subtracted EQE'
                              )
            self.axSub_1.legend()

            self.axSub_2.plot(data_Fit['Energy'],
                              data_Fit['Signal'],
                              linewidth=2,
                              linestyle='--',
                              color=color_fit,
                              label=label_fit
                              )
            self.axSub_2.plot(data_EQE['Energy'],
                              data_EQE['EQE'],
                              linewidth=2,
                              linestyle='-',
                              color=color_eqe,
                              label=label_eqe
                              )
            self.axSub_2.plot(data_EQE['Energy'],
                              sub_EQE,
                              linewidth=2,
                              linestyle='-',
                              color='#1f77b4',
                              label='Subtracted EQE'
                              )
            self.axSub_2.legend()

    # -----------------------------------------------------------------------------------------------------------

    # Function to add peak fits

    def add_Fits(self,
                 data_OptFit,
                 data_CTFit,
                 data_EQE
                 ):
        """
        Function to add fit peaks
        :param data_OptFit: Optical peak fit data [dataFrame]
        :param data_CTFit: CT peak fit data [dataFrame]
        :param data_EQE: EQE data [dataFrame]
        :return: None
        """

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

        try:  # Check if fit for an optical peak was imported
            T_OptFit = data_OptFit['Temperature'][0]
            f_OptFit = data_OptFit['Oscillator Strength (eV**2)'][0]
            l_OptFit = data_OptFit['Reorganization Energy (eV)'][0]
            E_OptFit = data_OptFit['Optical Peak Energy (eV)'][0]
        except:
            self.logger.error('No optical peak fit imported.')

        try:  # Check if fit for a CT state fit was importated
            T_CTFit = data_CTFit['Temperature'][0]
            f_CTFit = data_CTFit['Oscillator Strength (eV**2)'][0]
            l_CTFit = data_CTFit['Reorganization Energy (eV)'][0]
            E_CTFit = data_CTFit['CT State Energy (eV)'][0]
        except:
            self.logger.error('No CT state fit imported.')

        if E_OptFit != 0 and E_CTFit != 0:  # Only progress if a valid energy was imported

            if len(data_EQE) != 0:

                for x in range(len(data_EQE['Energy'])):
                    if data_EQE['Energy'][x] < max(data_OptFit['Energy']):
                        OptFit_value = calculate_gaussian_absorption(data_EQE['Energy'][x],
                                                                     f_OptFit,
                                                                     l_OptFit,
                                                                     E_OptFit,
                                                                     T_OptFit
                                                                     )
                        CTFit_value = calculate_gaussian_absorption(data_EQE['Energy'][x],
                                                                    f_CTFit,
                                                                    l_CTFit,
                                                                    E_CTFit,
                                                                    T_CTFit
                                                                    )
                        add_Energy.append(data_EQE['Energy'][x])
                        add_Fits.append(OptFit_value + CTFit_value)

                self.axAdd_1, self.axAdd_2 = set_up_EQE_plot()

                self.axAdd_1.plot(data_OptFit['Energy'],
                                  data_OptFit['Signal'],
                                  linewidth=2,
                                  linestyle='--',
                                  color=color_OptFit,
                                  label=label_OptFit
                                  )
                self.axAdd_1.plot(data_CTFit['Energy'],
                                  data_CTFit['Signal'],
                                  linewidth=2,
                                  linestyle='--',
                                  color=color_CTFit,
                                  label=label_CTFit
                                  )
                self.axAdd_1.plot(add_Energy,
                                  add_Fits,
                                  linewidth=2,
                                  linestyle='dotted',
                                  color='grey',
                                  label='$\mathrm{S_1}$ + CT Fit'
                                  )
                self.axAdd_1.plot(data_EQE['Energy'],
                                  data_EQE['EQE'],
                                  linewidth=2,
                                  linestyle='-',
                                  color=color_EQE,
                                  label=label_EQE
                                  )
                self.axAdd_1.legend()

                self.axAdd_2.plot(data_OptFit['Energy'],
                                  data_OptFit['Signal'],
                                  linewidth=2,
                                  linestyle='--',
                                  color=color_OptFit,
                                  label=label_OptFit
                                  )
                self.axAdd_2.plot(data_CTFit['Energy'],
                                  data_CTFit['Signal'],
                                  linewidth=2,
                                  linestyle='--',
                                  color=color_CTFit,
                                  label=label_CTFit
                                  )
                self.axAdd_2.plot(add_Energy,
                                  add_Fits,
                                  linewidth=2,
                                  linestyle='dotted',
                                  color='grey',
                                  label='$\mathrm{S_1}$ + CT Fit'
                                  )
                self.axAdd_2.plot(data_EQE['Energy'],
                                  data_EQE['EQE'],
                                  linewidth=2,
                                  linestyle='-',
                                  color=color_EQE,
                                  label=label_EQE
                                  )
                self.axAdd_2.legend()

                df_add = pd.DataFrame()
                df_add['Energy'] = np.array(add_Energy)
                df_add['EQE'] = np.array(add_Fits)
                df_add.to_csv('Fit_sum.csv')  # TODO: make more versatile

            else:
                self.logger.error('Please import a valid EQE file.')

    # -----------------------------------------------------------------------------------------------------------

    # Page 4 - Extended Fits (Marcus Theory)

    # -----------------------------------------------------------------------------------------------------------

    # Separate Double Peak Fit

    # Function to compile fits for separate double peak fitting

    def double_fit(self):
        """
        Function to perform separate double fit of S1 and CT peaks
        :return: None
        """

        increase_factor = 1.05  # Adding 5% to the selected data
        include_disorder = False
        guessRange_Sig = None

        # Import relevant parameters

        if self.ui.bias_DoubleFit.isChecked():
            self.bias = True
            self.tolerance = float(self.ui.tolerance.value()) / 100
            self.logger.info('Constraining fit below EQE data.')
        else:
            self.bias = False
            self.logger.info('Not constraining fit.')

        if self.ui.disorder_DoubleFit.isChecked():
            include_disorder = True
            startGuess_Sig = float(self.ui.guessStartSig_CT.value())
            stopGuess_Sig = float(self.ui.guessStopSig_CT.value())
            guessSig_ok = StartStop_is_valid(startGuess_Sig, stopGuess_Sig)
            if guessSig_ok:
                guessRange_Sig = np.round(np.arange(startGuess_Sig, stopGuess_Sig + 0.01, 0.01), 3).tolist()
            else:
                guessRange_Sig = np.round(np.arange(startGuess_Sig, startGuess_Sig + 0.1, 0.01), 3).tolist()

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

        startCT_ok = StartStop_is_valid(startStart_CT, startStop_CT)
        stopCT_ok = StartStop_is_valid(stopStart_CT, stopStop_CT)

        guessOpt_ok = StartStop_is_valid(startGuess_Opt, stopGuess_Opt)
        guessCT_ok = StartStop_is_valid(startGuess_CT, stopGuess_CT)

        # Compile all start / stop energies for Opt and CT fit

        if startOpt_ok and stopOpt_ok and guessOpt_ok and startCT_ok and stopCT_ok and guessCT_ok:

            startRange_Opt = np.round(np.arange(startStart_Opt, startStop_Opt + 0.005, 0.01),
                                      3).tolist()  # Change step to 0.05
            stopRange_Opt = np.round(np.arange(stopStart_Opt, stopStop_Opt + 0.005, 0.01), 3).tolist()

            startRange_CT = np.round(np.arange(startStart_CT, startStop_CT + 0.005, 0.01), 3).tolist()
            stopRange_CT = np.round(np.arange(stopStart_CT, stopStop_CT + 0.005, 0.01), 3).tolist()

            guessRange_Opt = np.round(np.arange(startGuess_Opt, stopGuess_Opt + 0.1, 0.05), 2).tolist()
            guessRange_CT = np.round(np.arange(startGuess_CT, stopGuess_CT + 0.1, 0.05), 2).tolist()

            # Compile a dataFrame with all combinations of start / stop values for Opt and CT fit

            self.logger.info('Compiling Fit Ranges ...')

            df_Opt = pd.DataFrame()
            df_CT = pd.DataFrame()

            start_Opt_list = []
            stop_Opt_list = []
            start_CT_list = []
            stop_CT_list = []

            for startOpt in startRange_Opt:
                for stopOpt in stopRange_Opt:
                    start_Opt_list.append(startOpt)
                    stop_Opt_list.append(stopOpt)

            for startCT in startRange_CT:
                for stopCT in stopRange_CT:
                    start_CT_list.append(startCT)
                    stop_CT_list.append(stopCT)

            df_Opt['Start'] = start_Opt_list
            df_Opt['Stop'] = stop_Opt_list

            df_CT['Start'] = start_CT_list
            df_CT['Stop'] = stop_CT_list

            # Calculate all optical peak fits

            self.logger.info('Calculating Optical Peak Fits ...')

            cal_vals_Opt = list(map(lambda x: calculate_guess_fit(x=x,
                                                                  df=df_Opt,
                                                                  eqe=eqe,
                                                                  function=self.gaussian_double,
                                                                  guessRange=guessRange_Opt
                                                                  ), tqdm(range(len(df_Opt)))))

            best_vals_Opt = list(map(lambda list_: sep_list(list_, 0), cal_vals_Opt))

            R2_Opt = list(map(lambda list_: sep_list(list_, 1), cal_vals_Opt))

            df_Opt['Fit'] = best_vals_Opt
            df_Opt['R2'] = R2_Opt

            # Calculate CT state fits

            start_Opt_list = []
            stop_Opt_list = []
            start_CT_list = []
            stop_CT_list = []

            best_vals_Opt = []
            best_vals_CT = []

            R2_Opt = []
            R2_CT = []

            combined_R2_list = []

            Opt_Fit_list = []
            CT_Fit_list = []
            combined_Fit_list = []
            Energy_list = []
            EQE_list = []

            df_results = pd.DataFrame()

            self.logger.info('Calculating CT State Fits ...')

            if include_disorder:
                self.logger.info('Including CT State Disorder ...')

            # If Optical peak to be subtracted before CT fit
            if self.ui.subtract_DoubleFit.isChecked() and not self.ui.bestSubtract_DoubleFit_2.isChecked():
                self.logger.info('Subtracting All Optical Peak Fits ...')
                for x in tqdm(range(len(df_Opt))):
                    for y in tqdm(range(len(df_CT))):
                        if df_Opt['R2'][x] > 0:  # Check that the optical peak fit was successful
                            new_eqe = subtract_Opt(eqe, df_Opt['Fit'][x], T=self.T_double)

                            if include_disorder:
                                best_vals, r_squared = guess_fit(eqe=new_eqe,
                                                                 startE=df_CT['Start'][y],
                                                                 stopE=df_CT['Stop'][y],
                                                                 function=self.gaussian_disorder_double,
                                                                 guessRange=guessRange_CT,
                                                                 guessRange_sig=guessRange_Sig,
                                                                 include_disorder=True
                                                                 )
                            else:
                                best_vals, r_squared = guess_fit(eqe=new_eqe,
                                                                 startE=df_CT['Start'][y],
                                                                 stopE=df_CT['Stop'][y],
                                                                 function=self.gaussian_double,
                                                                 guessRange=guessRange_CT
                                                                 )
                        else:
                            best_vals = [0, 0, 0]
                            r_squared = 0

                        start_Opt_list.append(df_Opt['Start'][x])
                        stop_Opt_list.append(df_Opt['Stop'][x])
                        start_CT_list.append(df_CT['Start'][y])
                        stop_CT_list.append(df_CT['Stop'][y])
                        best_vals_Opt.append(df_Opt['Fit'][x])
                        best_vals_CT.append(best_vals)
                        R2_Opt.append(df_Opt['R2'][x])
                        R2_CT.append(r_squared)

                        # Calculate combined fit here
                        parameter_dict = calculate_combined_fit(stopE=df_Opt['Stop'][x],
                                                                best_vals_Opt=df_Opt['Fit'][x],
                                                                best_vals_CT=best_vals,
                                                                R2_Opt=df_Opt['R2'][x],
                                                                R2_CT=r_squared,
                                                                eqe=eqe,
                                                                T=self.T_double,
                                                                bias=self.bias,
                                                                tolerance=self.tolerance,
                                                                range=increase_factor,
                                                                include_disorder=include_disorder)

                        combined_R2_list.append(parameter_dict['R2_Combined'])
                        combined_Fit_list.append(parameter_dict['Combined_Fit'])
                        Opt_Fit_list.append(parameter_dict['Opt_Fit'])
                        CT_Fit_list.append(parameter_dict['CT_Fit'])
                        Energy_list.append(parameter_dict['Energy'])
                        EQE_list.append(parameter_dict['EQE'])

            # If only best Optical peak is to be subtracted before CT fit
            elif self.ui.bestSubtract_DoubleFit_2.isChecked() and not self.ui.subtract_DoubleFit.isChecked():
                self.logger.info('Subtracting Only Best Optical Peak Fit ...')

                # best_fit_index = df_Opt['Fit'][df_Opt['R2']==max(df_Opt['R2'])].index[0]
                # print(best_fit_index)

                # To avoid picking a fit that has a high R2 but moves above the data
                advanced_R2_list = []
                for x in range(len(df_Opt)):
                    wave_fit, energy_fit, eqe_fit, log_eqe_fit = compile_EQE(eqe,
                                                                             df_Opt['Start'][x],
                                                                             df_Opt['Stop'][x] * increase_factor,
                                                                             1
                                                                             )
                    y_fit = [self.gaussian_double(e,
                                                  df_Opt['Fit'][x][0],
                                                  df_Opt['Fit'][x][1],
                                                  df_Opt['Fit'][x][2]
                                                  ) for e in energy_fit]
                    advanced_R2_list.append(R_squared(eqe_fit, y_fit))

                df_Opt['Advanced R2'] = advanced_R2_list

                best_fit_index = df_Opt['Fit'][df_Opt['Advanced R2'] == max(df_Opt['Advanced R2'])].index[0]
                # print(best_fit_index)

                new_eqe = subtract_Opt(eqe, df_Opt['Fit'][best_fit_index], T=self.T_double)

                for y in tqdm(range(len(df_CT))):

                    if include_disorder:
                        best_vals, r_squared = guess_fit(eqe=new_eqe,
                                                         startE=df_CT['Start'][y],
                                                         stopE=df_CT['Stop'][y],
                                                         function=self.gaussian_disorder_double,
                                                         guessRange=guessRange_CT,
                                                         guessRange_sig=guessRange_Sig,
                                                         include_disorder=True
                                                         )
                    else:
                        best_vals, r_squared = guess_fit(eqe=new_eqe,
                                                         startE=df_CT['Start'][y],
                                                         stopE=df_CT['Stop'][y],
                                                         function=self.gaussian_double,
                                                         guessRange=guessRange_CT
                                                         )

                    start_Opt_list.append(df_Opt['Start'][best_fit_index])
                    stop_Opt_list.append(df_Opt['Stop'][best_fit_index])
                    start_CT_list.append(df_CT['Start'][y])
                    stop_CT_list.append(df_CT['Stop'][y])
                    best_vals_Opt.append(df_Opt['Fit'][best_fit_index])
                    best_vals_CT.append(best_vals)
                    R2_Opt.append(df_Opt['R2'][best_fit_index])
                    R2_CT.append(r_squared)

                    # Calculate combined fit here
                    parameter_dict = calculate_combined_fit(stopE=df_Opt['Stop'][x],
                                                            best_vals_Opt=df_Opt['Fit'][x],
                                                            best_vals_CT=best_vals,
                                                            R2_Opt=df_Opt['R2'][x],
                                                            R2_CT=r_squared,
                                                            eqe=eqe,
                                                            T=self.T_double,
                                                            bias=self.bias,
                                                            tolerance=self.tolerance,
                                                            range=increase_factor,
                                                            include_disorder=include_disorder)

                    combined_R2_list.append(parameter_dict['R2_Combined'])
                    combined_Fit_list.append(parameter_dict['Combined_Fit'])
                    Opt_Fit_list.append(parameter_dict['Opt_Fit'])
                    CT_Fit_list.append(parameter_dict['CT_Fit'])
                    Energy_list.append(parameter_dict['Energy'])
                    EQE_list.append(parameter_dict['EQE'])

            # If Optical peak not to be subtracted before CT fit
            elif not self.ui.subtract_DoubleFit.isChecked() and not self.ui.bestSubtract_DoubleFit_2.isChecked():
                self.logger.info('Not Subtracting Optical Peak Fits.')
                for x in tqdm(range(len(df_Opt))):
                    for y in tqdm(range(len(df_CT))):

                        if include_disorder:
                            best_vals, r_squared = guess_fit(eqe=eqe,
                                                             startE=df_CT['Start'][y],
                                                             stopE=df_CT['Stop'][y],
                                                             function=self.gaussian_disorder_double,
                                                             guessRange=guessRange_CT,
                                                             guessRange_sig=guessRange_Sig,
                                                             include_disorder=True
                                                             )
                        else:
                            best_vals, r_squared = guess_fit(eqe=eqe,
                                                             startE=df_CT['Start'][y],
                                                             stopE=df_CT['Stop'][y],
                                                             function=self.gaussian_double,
                                                             guessRange=guessRange_CT
                                                             )

                    start_Opt_list.append(df_Opt['Start'][x])
                    stop_Opt_list.append(df_Opt['Stop'][x])
                    start_CT_list.append(df_CT['Start'][y])
                    stop_CT_list.append(df_CT['Stop'][y])
                    best_vals_Opt.append(df_Opt['Fit'][x])
                    best_vals_CT.append(best_vals)
                    R2_Opt.append(df_Opt['R2'][x])
                    R2_CT.append(r_squared)

                    # Calculate combined fit here
                    parameter_dict = calculate_combined_fit(stopE=df_Opt['Stop'][x],
                                                            best_vals_Opt=df_Opt['Fit'][x],
                                                            best_vals_CT=best_vals,
                                                            R2_Opt=df_Opt['R2'][x],
                                                            R2_CT=r_squared,
                                                            eqe=eqe,
                                                            T=self.T_double,
                                                            bias=self.bias,
                                                            tolerance=self.tolerance,
                                                            range=increase_factor,
                                                            include_disorder=include_disorder)

                    combined_R2_list.append(parameter_dict['R2_Combined'])
                    combined_Fit_list.append(parameter_dict['Combined_Fit'])
                    Opt_Fit_list.append(parameter_dict['Opt_Fit'])
                    CT_Fit_list.append(parameter_dict['CT_Fit'])
                    Energy_list.append(parameter_dict['Energy'])
                    EQE_list.append(parameter_dict['EQE'])

            else:
                self.logger.info('Please select valid fit settings.')

            # The same code but using map functions
            ### NOTE: Disorder not yet included

            # # If Optical peak to be subtracted before CT fit
            #
            # if self.ui.subtract_DoubleFit.isChecked():
            #
            #     parameter_list = list(map(lambda x: map_fit(x=x,
            #                                                 df_Opt=df_Opt,
            #                                                 df_CT=df_CT,
            #                                                 eqe=eqe,
            #                                                 guessRange_CT=guessRange_CT,
            #                                                 function=self.gaussian_double,
            #                                                 T=self.T_double,
            #                                                 bias=self.bias,
            #                                                 tolerance=self.tolerance,
            #                                                 sub_fit=1),
            #                               tqdm(range(len(df_Opt)))))
            #
            #     best_vals_Opt = sep_list_list(list(map(lambda list_: sep_list(list_, 0), parameter_list)))
            #     R2_Opt = sep_list_list(list(map(lambda list_: sep_list(list_, 1), parameter_list)))
            #
            #     best_vals_CT = sep_list_list(list(map(lambda list_: sep_list(list_, 2), parameter_list)))
            #     R2_CT = sep_list_list(list(map(lambda list_: sep_list(list_, 3), parameter_list)))
            #
            #     start_Opt_list = sep_list_list(list(map(lambda list_: sep_list(list_, 4), parameter_list)))
            #     stop_Opt_list = sep_list_list(list(map(lambda list_: sep_list(list_, 5), parameter_list)))
            #
            #     start_CT_list = sep_list_list(list(map(lambda list_: sep_list(list_, 6), parameter_list)))
            #     stop_CT_list = sep_list_list(list(map(lambda list_: sep_list(list_, 7), parameter_list)))
            #
            #     combined_R2_list = sep_list_list(list(map(lambda list_: sep_list(list_, 8), parameter_list)))
            #     combined_Fit_list = sep_list_list(list(map(lambda list_: sep_list(list_, 9), parameter_list)))
            #     Opt_Fit_list = sep_list_list(list(map(lambda list_: sep_list(list_, 10), parameter_list)))
            #     CT_Fit_list = sep_list_list(list(map(lambda list_: sep_list(list_, 11), parameter_list)))
            #     Energy_list = sep_list_list(list(map(lambda list_: sep_list(list_, 12), parameter_list)))
            #     EQE_list = sep_list_list(list(map(lambda list_: sep_list(list_, 13), parameter_list)))
            #
            # else:
            #     parameter_list = list(map(lambda x: map_fit(x=x,
            #                                                 df_Opt=df_Opt,
            #                                                 df_CT=df_CT,
            #                                                 eqe=eqe,
            #                                                 guessRange_CT=guessRange_CT,
            #                                                 function=self.gaussian_double,
            #                                                 T=self.T_double,
            #                                                 bias=self.bias,
            #                                                 tolerance=self.tolerance,
            #                                                 sub_fit=0),
            #                               tqdm(range(len(df_Opt)))))
            #
            #     best_vals_Opt = sep_list_list(list(map(lambda list_: sep_list(list_, 0), parameter_list)))
            #     R2_Opt = sep_list_list(list(map(lambda list_: sep_list(list_, 1), parameter_list)))
            #
            #     best_vals_CT = sep_list_list(list(map(lambda list_: sep_list(list_, 2), parameter_list)))
            #     R2_CT = sep_list_list(list(map(lambda list_: sep_list(list_, 3), parameter_list)))
            #
            #     start_Opt_list = sep_list_list(list(map(lambda list_: sep_list(list_, 4), parameter_list)))
            #     stop_Opt_list = sep_list_list(list(map(lambda list_: sep_list(list_, 5), parameter_list)))
            #
            #     start_CT_list = sep_list_list(list(map(lambda list_: sep_list(list_, 6), parameter_list)))
            #     stop_CT_list = sep_list_list(list(map(lambda list_: sep_list(list_, 7), parameter_list)))
            #
            #     combined_R2_list = sep_list_list(list(map(lambda list_: sep_list(list_, 8), parameter_list)))
            #     combined_Fit_list = sep_list_list(list(map(lambda list_: sep_list(list_, 9), parameter_list)))
            #     Opt_Fit_list = sep_list_list(list(map(lambda list_: sep_list(list_, 10), parameter_list)))
            #     CT_Fit_list = sep_list_list(list(map(lambda list_: sep_list(list_, 11), parameter_list)))
            #     Energy_list = sep_list_list(list(map(lambda list_: sep_list(list_, 12), parameter_list)))
            #     EQE_list = sep_list_list(list(map(lambda list_: sep_list(list_, 13), parameter_list)))

            if len(best_vals_Opt) == len(best_vals_CT) and len(best_vals_Opt) != 0:  # Confirm lists are acceptable

                df_results['Start_Opt'] = start_Opt_list
                df_results['Stop_Opt'] = stop_Opt_list
                df_results['Fit_Opt'] = best_vals_Opt
                df_results['R2_Opt'] = R2_Opt
                df_results['Start_CT'] = start_CT_list
                df_results['Stop_CT'] = stop_CT_list
                df_results['Fit_CT'] = best_vals_CT
                df_results['R2_CT'] = R2_CT

                # Add combined fit to dataFrame
                df_results['Total_R2'] = combined_R2_list
                df_results['Total_Fit'] = combined_Fit_list
                df_results['Opt_Fit'] = Opt_Fit_list
                df_results['CT_Fit'] = CT_Fit_list
                df_results['Energy'] = Energy_list
                df_results['EQE'] = EQE_list

                # Find best fit

                self.logger.info('Determining Best Fit ...')
                self.logger.info('Fit Results: ')
                print("")

                label = pick_EQE_Label(self.ui.textBox_dF2, self.ui.textBox_dF1)

                for x in np.arange(1, 6, 1):
                    print('-' * 80)
                    print(('Best Fit No. {} : ').format(x))
                    df_results = find_best_fit(df_both=df_results,
                                               eqe=eqe,
                                               T=self.T_double,
                                               label=label,
                                               n_fit=x,
                                               include_disorder=include_disorder
                                               )

                print('-' * 80)
                print("")

        self.bias = False

    # -----------------------------------------------------------------------------------------------------------

    # Gaussian fitting function for double fit

    def gaussian_double(self, E, f, l, Ect):
        """
        Standard gaussian functions to separately fit double peaks
        :param x: List of energy values [list]
        :param f: Oscillator strength [float]
        :param l: Reorganization Energy [float]
        :param Ect: Peak Energy [float]
        :return: EQE value [list]
        """

        return (f / (E * math.sqrt(4 * math.pi * l * self.T_double * self.k))) * exp(
            -(Ect + l - E) ** 2 / (4 * l * self.k * self.T_double))

    # Gaussian fitting function including disorder for double fit

    def gaussian_disorder_double(self, E, f, l, Ect, sig):
        """
        Standard gaussian functions including disorder to separately fit double peaks
        :param E: List of energy values [list of floats]
        :param f: Oscillator strength [float]
        :param l: Reorganization Energy [float]
        :param Ect: Charge Transfer State Energy [float]
        :param sig: Gaussian disorder [float]
        :return: list of EQE values [list of floats]
        """

        return [(f / (e * math.sqrt(2 * math.pi * (2 * l * self.T_double * self.k + sig ** 2)))) * exp(
            -(Ect - (sig ** 2 / (2 * self.k * self.T_double)) + l + (
                    sig ** 2 / (2 * self.k * self.T_double)) - e) ** 2 / (
                    4 * l * self.k * self.T_double + 2 * sig ** 2)) for e in E]

    # -----------------------------------------------------------------------------------------------------------

    # Simultaneous Double Peak Fit

    # Function to perform simultaneous double peak fitting once

    def sim_double_fit_single(self):
        """
        Function to perform simultaneous double peak fitting once
        :return: None
        """
        # Import relevant parameters

        if self.ui.disorder_Sim.isChecked():
            include_disorder = True
        else:
            include_disorder = False

        eqe = self.data_sim
        self.T_sim = self.ui.Temperature_Sim.value()

        bound_dict = self.load_sim_dict()

        # Compile EQE data
        energy_fit, eqe_fit = compile_Data(energy=eqe['Energy'],
                                           y=eqe['EQE'],
                                           startE=bound_dict['start_fit'],
                                           stopE=bound_dict['stop_fit']
                                           )
        # Set plot range
        x_plot = linspace(bound_dict['start_plot'], bound_dict['stop_plot'], 50)

        if include_disorder:
            # NOTE: Adjust this to change initial guesses
            p0 = [0.001, 0.15, 1.30, 0.01, 0.150, 1.5, 0.1]

            best_vals, covar, y_fit, r_squared = fit_model_double(function=self.gaussian_disorder_double_sim,
                                                                  energy_fit=energy_fit,
                                                                  eqe_fit=eqe_fit,
                                                                  bound_dict=bound_dict,
                                                                  p0=p0,
                                                                  include_disorder=include_disorder,
                                                                  print_report=False
                                                                  )
            y_CT = [calculate_gaussian_disorder_absorption(i,
                                                           f=best_vals[0],
                                                           l=best_vals[1],
                                                           E=best_vals[2],
                                                           sig=best_vals[6],
                                                           T=self.T_sim
                                                           ) for i in x_plot]
            y_sum = [self.gaussian_disorder_double_sim(i,
                                                       fCT=best_vals[0],
                                                       lCT=best_vals[1],
                                                       ECT=best_vals[2],
                                                       fopt=best_vals[3],
                                                       lopt=best_vals[4],
                                                       Eopt=best_vals[5],
                                                       sig=best_vals[6]
                                                       ) for i in x_plot]
        else:
            # NOTE: Adjust this to change initial guesses
            p0 = [0.001, 0.15, 1.30, 0.01, 0.150, 1.5]

            best_vals, covar, y_fit, r_squared = fit_model_double(function=self.gaussian_double_sim,
                                                                  energy_fit=energy_fit,
                                                                  eqe_fit=eqe_fit,
                                                                  bound_dict=bound_dict,
                                                                  p0=p0,
                                                                  include_disorder=include_disorder,
                                                                  print_report=False
                                                                  )
            y_CT = [calculate_gaussian_absorption(i,
                                                  f=best_vals[0],
                                                  l=best_vals[1],
                                                  E=best_vals[2],
                                                  T=self.T_sim
                                                  ) for i in x_plot]
            y_sum = [self.gaussian_double_sim(i,
                                              fCT=best_vals[0],
                                              lCT=best_vals[1],
                                              ECT=best_vals[2],
                                              fopt=best_vals[3],
                                              lopt=best_vals[4],
                                              Eopt=best_vals[5]
                                              ) for i in x_plot]

        y_opt = [calculate_gaussian_absorption(i,
                                               f=best_vals[3],
                                               l=best_vals[4],
                                               E=best_vals[5],
                                               T=self.T_sim
                                               ) for i in x_plot]

        # TODO: Check that all prints have errors
        print('-' * 35)
        print('R_Squared : ', format(r_squared, '.6f'))
        print('-' * 35)
        print('f_Opt (eV**2) : ', format(best_vals[3], '.6f'), '+/-', format(math.sqrt(covar[3, 3]), '.6f'))
        print('l_Opt (eV) : ', format(best_vals[4], '.6f'), '+/-', format(math.sqrt(covar[4, 4]), '.6f'))
        print('E_Opt (eV) : ', format(best_vals[5], '.6f'), '+/-', format(math.sqrt(covar[5, 5]), '.6f'))
        print('-' * 35)
        print('f_CT (eV**2) : ', format(best_vals[0], '.6f'), '+/-', format(math.sqrt(covar[0, 0]), '.6f'))
        print('l_CT (eV) : ', format(best_vals[1], '.6f'), '+/-', format(math.sqrt(covar[1, 1]), '.6f'))
        print('E_CT (eV) : ', format(best_vals[2], '.6f'), '+/-', format(math.sqrt(covar[2, 2]), '.6f'))

        if include_disorder:
            print('Sigma (eV) : ', format(best_vals[6], '.6f'), '+/-', format(math.sqrt(covar[6, 6]), '.6f'))

        # print('Temperature [T] (K) : ', T)
        # print('-' * 80)

        label = pick_EQE_Label(self.ui.textBox_simFit_label, self.ui.textBox_simFit)

        axDouble_1, axDouble_2 = set_up_plot(flag='Energy')

        axDouble_1.plot(eqe['Energy'],
                        eqe['EQE'],
                        linewidth=2,
                        linestyle='-',
                        label=label,
                        color='black'
                        )
        axDouble_1.plot(x_plot,
                        y_opt,
                        linewidth=2,
                        linestyle='dotted',
                        label='Optical Peak Fit'
                        )
        axDouble_1.plot(x_plot,
                        y_CT,
                        linewidth=2,
                        linestyle='--',
                        label='CT State Fit'
                        )
        axDouble_1.plot(x_plot,
                        y_sum,
                        linewidth=2,
                        linestyle='dashdot',
                        label='Total Fit'
                        )
        axDouble_1.legend()

        axDouble_2.plot(eqe['Energy'],
                        eqe['EQE'],
                        linewidth=2,
                        linestyle='-',
                        label=label,
                        color='black'
                        )
        axDouble_2.plot(x_plot,
                        y_opt,
                        linewidth=2,
                        linestyle='--',
                        label='Optical Peak Fit'
                        )
        axDouble_2.plot(x_plot,
                        y_CT,
                        linewidth=2, linestyle='--',
                        label='CT State Fit'
                        )
        axDouble_2.plot(x_plot,
                        y_sum,
                        linewidth=2,
                        linestyle='dashdot',
                        label='Total Fit'
                        )
        axDouble_2.set_ylim([10 ** (-7), max(eqe['EQE']) * 1.4])
        axDouble_2.legend()

    # -----------------------------------------------------------------------------------------------------------

    # Function to perform simultaneous double peak fitting multiple times

    def sim_double_fit(self):
        """
        Function to perform simultaneous double peak fitting multiple times to find best fit
        :return: None
        """
        # Import relevant parameters

        eqe = self.data_sim
        self.T_sim = self.ui.Temperature_Sim.value()

        bound_dict = self.load_sim_dict()

        if self.ui.disorder_Sim.isChecked():
            include_disorder = True
        else:
            include_disorder = False

        if self.ui.bias_SimFit.isChecked():
            self.bias_sim = True
            self.tolerance_sim = float(self.ui.tolerance_Sim.value()) / 100
            self.logger.info('Constraining fit below EQE data.')
        else:
            self.bias_sim = False
            self.logger.info('Not constraining fit.')

        # # Create selection of initial guesses
        # # NOTE: Change this to adjust initial guesses
        # startGuess_Opt = 1.4
        # stopGuess_Opt = 1.8
        # startGuess_CT = 1.2
        # stopGuess_CT = 1.5
        # startGuess_sig = 0.01
        # stopGuess_sig = 0.2

        # Check that all start and stop energies are valid
        start_ok = StartStop_is_valid(bound_dict['start_start'], bound_dict['start_stop'])
        stop_ok = StartStop_is_valid(bound_dict['stop_start'], bound_dict['stop_stop'])

        # Compile all start / stop energies for CT fit
        if start_ok and stop_ok:

            startRange = np.round(np.arange(bound_dict['start_start'],
                                            bound_dict['start_stop'] + 0.005,
                                            0.01), 3).tolist()
            stopRange = np.round(np.arange(bound_dict['stop_start'],
                                           bound_dict['stop_stop'] + 0.005,
                                           0.01), 3).tolist()

            # guessRange_Opt = np.round(np.arange(startGuess_Opt,
            #                                     stopGuess_Opt + 0.05,
            #                                     0.1), 2).tolist()  # Reduced step size from 0.05 to 0.1
            # guessRange_CT = np.round(np.arange(startGuess_CT,
            #                                    stopGuess_CT + 0.05,
            #                                    0.1), 2).tolist()
            # guessRange_sig = np.round(np.arange(startGuess_sig,
            #                                     stopGuess_sig + 0.1,
            #                                     0.2), 2).tolist()  # Reduced step size from 0.05 to 0.2

            # Compile a dataFrame with all combinations of start / stop values for the optical and CT fit
            self.logger.info('Compiling Fit Ranges ...')

            df = pd.DataFrame()

            start_list = []
            stop_list = []

            for start in startRange:
                for stop in stopRange:
                    start_list.append(start)
                    stop_list.append(stop)

            df['Start'] = start_list
            df['Stop'] = stop_list

            self.logger.info('Calculating Fits ...')

            best_vals_Opt = []
            best_vals_CT = []

            start_list = []
            stop_list = []

            R2_Opt_list = []
            R2_CT_list = []
            R2_sum_list = []
            R2_average_list = []

            Opt_Fit_list = []
            CT_Fit_list = []
            combined_Fit_list = []
            Energy_list = []
            EQE_list = []

            df_results = pd.DataFrame()

            for x in tqdm(range(len(df))):
                if df['Start'][x] < df['Stop'][x]:
                    wave_fit, energy_fit, eqe_fit, log_eqe_fit = compile_EQE(eqe, df['Start'][x], df['Stop'][x], 1)

                    if include_disorder:
                        # NOTE: Adjust this to change initial guesses
                        # CAVEAT: I tried the guess_fit function to test other guesses but didn't achieve great results
                        p0 = [0.001, 0.15, 1.30, 0.01, 0.150, 1.5, 0.1]

                        try:
                            best_vals, covar, y_fit, r_squared = fit_model_double(
                                function=self.gaussian_disorder_double_sim,
                                energy_fit=energy_fit,
                                eqe_fit=eqe_fit,
                                bound_dict=bound_dict,
                                p0=p0,
                                include_disorder=include_disorder,
                                print_report=False
                                )

                            best_CT = [
                                best_vals[0],
                                best_vals[1],
                                best_vals[2],
                                best_vals[6],
                            ]
                            best_Opt = [
                                best_vals[3],
                                best_vals[4],
                                best_vals[5]
                            ]
                        except:
                            best_vals = [0, 0, 0, 0, 0, 0, 0]
                    else:
                        # NOTE: Adjust this to change initial guesses
                        # CAVEAT: I tried the guess_fit function to test other guesses but didn't achieve great results
                        p0 = [0.001, 0.15, 1.30, 0.01, 0.150, 1.5]

                        try:
                            best_vals, covar, y_fit, r_squared = fit_model_double(function=self.gaussian_double_sim,
                                                                                  energy_fit=energy_fit,
                                                                                  eqe_fit=eqe_fit,
                                                                                  bound_dict=bound_dict,
                                                                                  p0=p0,
                                                                                  include_disorder=include_disorder,
                                                                                  print_report=False
                                                                                  )

                            # # Old Code to test guess_fit function
                            # best_vals, r_squared = guess_fit(eqe=eqe,
                            #                                  startE=df['Start'][x],
                            #                                  stopE=df['Stop'][x],
                            #                                  function=self.gaussian_double_sim,
                            #                                  guessRange=guessRange_CT,
                            #                                  guessRange_opt=guessRange_Opt,
                            #                                  bounds=bound_dict,
                            #                                  simultaneous_double=True,
                            #                                  )

                            best_CT = [
                                best_vals[0],
                                best_vals[1],
                                best_vals[2]
                            ]
                            best_Opt = [
                                best_vals[3],
                                best_vals[4],
                                best_vals[5]
                            ]
                        except:
                            best_vals = [0, 0, 0, 0, 0, 0]

                    if sum(best_vals) != 0:  # Check that the fit was successful
                        # Calculate combined fit here
                        parameter_dict = calculate_combined_fit(eqe=eqe,
                                                                stopE=df['Stop'][x],
                                                                best_vals_Opt=best_Opt,
                                                                best_vals_CT=best_CT,
                                                                T=self.T_sim,
                                                                bias=self.bias_sim,
                                                                tolerance=self.tolerance_sim,
                                                                include_disorder=include_disorder)

                        start_list.append(df['Start'][x])
                        stop_list.append(df['Stop'][x])
                        best_vals_CT.append(best_CT)
                        best_vals_Opt.append(best_Opt)

                        R2_sum_list.append(parameter_dict['R2_Combined'])
                        R2_Opt_list.append(parameter_dict['R2_Opt'])
                        R2_CT_list.append(parameter_dict['R2_CT'])
                        R2_average_list.append(parameter_dict['R2_Average'])

                        combined_Fit_list.append(parameter_dict['Combined_Fit'])
                        Opt_Fit_list.append(parameter_dict['Opt_Fit'])
                        CT_Fit_list.append(parameter_dict['CT_Fit'])
                        Energy_list.append(parameter_dict['Energy'])
                        EQE_list.append(parameter_dict['EQE'])

                    else:  # If sum was unsuccessful, skip and move on
                        pass

                else:
                    pass

            if len(best_vals_Opt) == len(best_vals_CT) and len(best_vals_Opt) != 0:

                df_results['Start'] = start_list
                df_results['Stop'] = stop_list
                df_results['Fit_Opt'] = best_vals_Opt
                df_results['R2_Opt'] = R2_Opt_list
                df_results['Fit_CT'] = best_vals_CT
                df_results['R2_CT'] = R2_CT_list

                # Add combined fit to dataFrame
                df_results['Total_R2'] = R2_sum_list
                df_results['Total_Fit'] = combined_Fit_list
                df_results['Opt_Fit'] = Opt_Fit_list
                df_results['CT_Fit'] = CT_Fit_list
                df_results['Energy'] = Energy_list
                df_results['EQE'] = EQE_list
                df_results['Comp_R2'] = R2_average_list

                # Find best fit

                self.logger.info('Determining Best Fit ...')
                self.logger.info('Fit Results: ')
                print("")

                label = pick_EQE_Label(self.ui.textBox_simFit_label, self.ui.textBox_simFit)

                for x in np.arange(1, 5, 1):
                    print('-' * 80)
                    print(('Best Fit No. {} : ').format(x))
                    df_results = find_best_fit(df_both=df_results,
                                               eqe=eqe,
                                               T=self.T_sim,
                                               label=label,
                                               n_fit=x,
                                               include_disorder=include_disorder,
                                               simultaneous_double=True
                                               )

                print('-' * 80)
                print("")

    # -----------------------------------------------------------------------------------------------------------

    # Function to load dictionary of fit bounds

    def load_sim_dict(self):
        """
        Function to load dictionary of fit bounds from GUI paramaters
        :return: bound_dict: Dictionary of boundary values [dict]
        """

        start_fit = self.ui.startFit_Sim.value()
        stop_fit = self.ui.stopFit_Sim.value()
        start_plot = self.ui.startPlot_Sim.value()
        stop_plot = self.ui.stopPlot_Sim.value()

        start_start = self.ui.startStart_Sim.value()
        start_stop = self.ui.startStop_Sim.value()
        stop_start = self.ui.stopStart_Sim.value()
        stop_stop = self.ui.stopStop_Sim.value()

        start_Eopt = self.ui.startBound_Eopt.value()
        stop_Eopt = self.ui.stopBound_Eopt.value()
        start_fopt = self.ui.startBound_fopt.value()
        stop_fopt = self.ui.stopBound_fopt.value()
        start_lopt = self.ui.startBound_lopt.value()
        stop_lopt = self.ui.stopBound_lopt.value()

        start_ECT = self.ui.startBound_ECT.value()
        stop_ECT = self.ui.stopBound_ECT.value()
        start_fCT = self.ui.startBound_fCT.value()
        stop_fCT = self.ui.stopBound_fCT.value()
        start_lCT = self.ui.startBound_lCT.value()
        stop_lCT = self.ui.stopBound_lCT.value()

        start_sig = self.ui.startBound_sig.value()
        stop_sig = self.ui.stopBound_sig.value()

        bound_dict = {
            'start_Eopt': start_Eopt,
            'stop_Eopt': stop_Eopt,
            'start_fopt': start_fopt,
            'stop_fopt': stop_fopt,
            'start_lopt': start_lopt,
            'stop_lopt': stop_lopt,
            'start_ECT': start_ECT,
            'stop_ECT': stop_ECT,
            'start_fCT': start_fCT,
            'stop_fCT': stop_fCT,
            'start_lCT': start_lCT,
            'stop_lCT': stop_lCT,
            'start_sig': start_sig,
            'stop_sig': stop_sig,
            'start_fit': start_fit,
            'stop_fit': stop_fit,
            'start_plot': start_plot,
            'stop_plot': stop_plot,
            'start_start': start_start,
            'start_stop': start_stop,
            'stop_start': stop_start,
            'stop_stop': stop_stop
        }

        return bound_dict

    # # -----------------------------------------------------------------------------------------------------------

    # Gaussian fitting function for simultaneous double peak fit

    def gaussian_double_sim(self, E, fCT, lCT, ECT, fopt, lopt, Eopt):
        """
        Standard gaussian function to perform simultaneous double peak fitting
        :param E: List of energy values [list]
        :param fCT: CT Oscillator strength [float]
        :param lCT: CT Reorganization Energy [float]
        :param ECT: CT Peak Energy [float]
        :param fopt: S1 Oscillator strength [float]
        :param lopt: S1 Reorganization Energy [float]
        :param Eopt: S1 Peak Energy [float]
        :return: EQE value [list]
        """

        val_CT = (fCT / (E * math.sqrt(4 * math.pi * lCT * self.T_sim * self.k))) * exp(
            -(ECT + lCT - E) ** 2 / (4 * lCT * self.k * self.T_sim))
        val_opt = (fopt / (E * math.sqrt(4 * math.pi * lopt * self.T_sim * self.k))) * exp(
            -(Eopt + lopt - E) ** 2 / (4 * lopt * self.k * self.T_sim))

        return val_CT + val_opt

    # Gaussian fitting function for simultaneous double peak fit including disorder

    def gaussian_disorder_double_sim(self, E, fCT, lCT, ECT, fopt, lopt, Eopt, sig):
        """
        Standard gaussian function including disorder to perform simultaneous double peak fitting
        :param E: List of energy values [list]
        :param fCT: CT Oscillator strength [float]
        :param lCT: CT Reorganization Energy [float]
        :param ECT: CT Peak Energy [float]
        :param fopt: S1 Oscillator strength [float]
        :param lopt: S1 Reorganization Energy [float]
        :param Eopt: S1 Peak Energy [float]
        :return: EQE value [list]
        """

        val_CT = (fCT / (E * math.sqrt(2 * math.pi * (2 * lCT * self.T_sim * self.k + sig ** 2)))) * exp(
            -(ECT + lCT - E) ** 2 / (4 * lCT * self.k * self.T_sim + 2 * sig ** 2))
        val_opt = (fopt / (E * math.sqrt(4 * math.pi * lopt * self.T_sim * self.k))) * exp(
            -(Eopt + lopt - E) ** 2 / (4 * lopt * self.k * self.T_sim))

        return val_CT + val_opt

    # -----------------------------------------------------------------------------------------------------------

    # Functions to clear plots

    # -----------------------------------------------------------------------------------------------------------

    # Function to clear "Calculate EQE" plot

    def clear_plot(self):
        """
        Function to clear plot
        :return: None
        """

        plt.close()  # Close the current plot
        self.ax1, self.ax2 = set_up_plot()  # Setting up new plot is preferred over plt.clf() in case window was closed

    # -----------------------------------------------------------------------------------------------------------

    # Function to clear EQE plot

    def clear_EQE_plot(self):
        """
        Function to clear EQE fit plot
        :return: None
        """

        plt.close()
        plt.close()
        self.axFit_1, self.axFit_2 = set_up_EQE_plot()

    # -----------------------------------------------------------------------------------------------------------

    # Function to clear EL plot

    def clear_EL_plot(self):
        """
        Function to clear EL plot
        :return: None
        """

        plt.close()
        plt.close()
        self.axEL_1, self.axEL_2 = set_up_EL_plot()


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
