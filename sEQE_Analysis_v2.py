#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 11:59:40 2018

@author: jungbluth
"""

import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
import serial
import itertools
import os, sys
import tkinter as tk
from tkinter import filedialog
import xlrd
import random
from colour import Color

# for the gui
from PyQt5 import QtCore, QtGui, QtWidgets
import sEQE_Analysis_template_v2


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        
        QtWidgets.QMainWindow.__init__(self)
        
        
        # Set up the user interface from Designer
        
        self.ui = sEQE_Analysis_template_v2.Ui_MainWindow()
        self.ui.setupUi(self)         


        # Tkinter
        
        root = tk.Tk()
        root.withdraw()
        
        
        # Create lists for reference / data / EQEQ files

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
        
        
        # Handle Calculate Buttons
        
        self.ui.calculateButton_1.clicked.connect(lambda: self.pre_EQE(self.ref_1, self.data_1, self.ui.startNM_1, self.ui.stopNM_1, 1))
        self.ui.calculateButton_2.clicked.connect(lambda: self.pre_EQE(self.ref_2, self.data_2, self.ui.startNM_2, self.ui.stopNM_2, 2))
        self.ui.calculateButton_3.clicked.connect(lambda: self.pre_EQE(self.ref_3, self.data_3, self.ui.startNM_3, self.ui.stopNM_3, 3))
        self.ui.calculateButton_4.clicked.connect(lambda: self.pre_EQE(self.ref_4, self.data_4, self.ui.startNM_4, self.ui.stopNM_4, 4))
        self.ui.calculateButton_5.clicked.connect(lambda: self.pre_EQE(self.ref_5, self.data_5, self.ui.startNM_5, self.ui.stopNM_5, 5))

        
        # Handle Export All Data Button
        
        self.ui.exportButton.clicked.connect(self.export_EQE)
        
        
        # Handle Plot EQE Buttons
        
        self.ui.plotEQEButton_Wavelength.clicked.connect(lambda: self.pre_plot_EQE(0))
        self.ui.plotEQEButton_Energy.clicked.connect(lambda: self.pre_plot_EQE(1))
        
        
        # Handle Clear Plot Button
        
        self.ui.clearButton.clicked.connect(self.clear_plot)        
        
        
        # Import photodiode calibration files

        Si_file = pd.ExcelFile("FDS100-CAL.xlsx") # The files are in the sEQE Analysis folder, just as this program
#        print(Si_file.sheet_names)
        self.Si_cal = Si_file.parse('Sheet1')
#        print(self.Si_cal)
        
        InGaAs_file = pd.ExcelFile("FGA21-CAL.xlsx")
        self.InGaAs_cal = InGaAs_file.parse('Sheet1')
        
        
        # Data directory
        
        self.data_dir = '/home/jungbluth/Desktop/EQE Control Software/sEQE Data'   ### CHANGE THIS IF NEEDED LATER


        # Define variables        
        
        self.h = 6.626 * math.pow(10,-34) # [m^2 kg/s]
        self.c = 2.998 * math.pow(10,8) # [m/s]
        self.q = 1.602 * math.pow(10,-19) # [C]
        
        
        # To export calculated EQE files
        
        self.export = False
        
        self.do_plot = True


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
                self.ref_1 = pd.DataFrame.from_csv(file_) # Turn file into dataFrame  
#                print(self.ref_1)
#                print(len(self.ref_1)) 
                
            elif textBox_no == 3:
                self.ref_2 = pd.DataFrame.from_csv(file_) 

            elif textBox_no == 5:
                self.ref_3 = pd.DataFrame.from_csv(file_)    
                               
            elif textBox_no == 7:
                self.ref_4 = pd.DataFrame.from_csv(file_) 
                
            elif textBox_no == 9:
                self.ref_5 = pd.DataFrame.from_csv(file_) 
                
            # For data files:
                
            elif textBox_no == 2:
                self.data_1 = pd.DataFrame.from_csv(file_)
                
            elif textBox_no == 4:
                self.data_2 = pd.DataFrame.from_csv(file_) 

            elif textBox_no == 6:
                self.data_3 = pd.DataFrame.from_csv(file_)    
                               
            elif textBox_no == 8:
                self.data_4 = pd.DataFrame.from_csv(file_) 
                
            elif textBox_no == 10:
                self.data_5 = pd.DataFrame.from_csv(file_)
                
            # For EQE files:
                
            elif textBox_no == 11:
                self.EQE_1 = pd.DataFrame.from_csv(file_)
                
            elif textBox_no == 14:
                self.EQE_2 = pd.DataFrame.from_csv(file_) 

            elif textBox_no == 17:
                self.EQE_3 = pd.DataFrame.from_csv(file_)    
                               
            elif textBox_no == 20:
                self.EQE_4 = pd.DataFrame.from_csv(file_) 
                
            elif textBox_no == 23:
                self.EQE_5 = pd.DataFrame.from_csv(file_)
                
            elif textBox_no == 26:
                self.EQE_6 = pd.DataFrame.from_csv(file_)
                
            elif textBox_no == 29:
                self.EQE_7 = pd.DataFrame.from_csv(file_) 

            elif textBox_no == 32:
                self.EQE_8 = pd.DataFrame.from_csv(file_)    
                               
            elif textBox_no == 35:
                self.EQE_9 = pd.DataFrame.from_csv(file_) 
                
            elif textBox_no == 38:
                self.EQE_10 = pd.DataFrame.from_csv(file_)
           
# -----------------------------------------------------------------------------------------------------------        

    #### Functions to calculate EQE

# -----------------------------------------------------------------------------------------------------------
  
    ### Function to select data and reference file   
  
    def pre_EQE(self, ref_df, data_df, start, stop, range_no):
        
        startNM = start.value()
        stopNM = stop.value()
        if self.is_valid(ref_df, data_df, startNM, stopNM, range_no):
            self.calculate_EQE(ref_df, data_df, startNM, stopNM, range_no)
            
# ----------------------------------------------------------------------------------------------------------- 

    ### Function to calculate EQE
  
    def calculate_EQE (self, ref_df, data_df, startNM, stopNM, range_no):
                   
        power_dict = {}        
        Wavelength = []  
        Energy = []                  
        EQE = []
        log_EQE = []
        
        
        for x in range(len(ref_df['Wavelength'])): # Iterate through columns of reference file
            power_dict[ref_df['Wavelength'][x]] = ref_df['Power'][x] # Add wavelength and corresponding power to dictionary
             
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

    #### Function to export EQE data

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
            if self.is_valid(self.ref_1, self.data_1, startNM1, stopNM1, 1): # Check that files are non-empty and within wavelength range            
                Wave_1, Energy_1, EQE_1, log_EQE_1 = self.calculate_EQE(self.ref_1, self.data_1, startNM1, stopNM1) # Extract data
                export_1 = pd.DataFrame({'Wavelength': Wave_1, 'Energy': Energy_1, 'EQE': EQE_1, 'Log_EQE': log_EQE_1}) # Create dataFrame with EQE data
                wave_inc['1'] = Wave_1[0] # Add the first wavelength value to the wave_inc list
            else:
                ok_1 = False # Set variable to False if calculation is invalid
                
        if self.ui.exportBox_2.isChecked():
            startNM2 = self.ui.startNM_2.value() 
            stopNM2 = self.ui.stopNM_2.value()
            if self.is_valid(self.ref_2, self.data_2, startNM2, stopNM2, 2): 
                Wave_2, Energy_2, EQE_2, log_EQE_2 = self.calculate_EQE(self.ref_2, self.data_2, startNM2, stopNM2)           
                export_2 = pd.DataFrame({'Wavelength': Wave_2, 'Energy': Energy_2, 'EQE': EQE_2, 'Log_EQE': log_EQE_2})
                wave_inc['2'] = Wave_2[0]
            else:
                ok_2 = False              
                
        if self.ui.exportBox_3.isChecked():
            startNM3 = self.ui.startNM_3.value()
            stopNM3 = self.ui.stopNM_3.value() 
            if self.is_valid(self.ref_3, self.data_3, startNM3, stopNM3, 3): 
                Wave_3, Energy_3, EQE_3, log_EQE_3 = self.calculate_EQE(self.ref_3, self.data_3, startNM3, stopNM3)
                export_3 = pd.DataFrame({'Wavelength': Wave_3, 'Energy': Energy_3, 'EQE': EQE_3, 'Log_EQE': log_EQE_3})  
                wave_inc['3'] = Wave_3[0]
            else:
                ok_3 = False                   
            
        if self.ui.exportBox_4.isChecked():
            startNM4 = self.ui.startNM_4.value() 
            stopNM4 = self.ui.stopNM_4.value()
            if self.is_valid(self.ref_4, self.data_4, startNM4, stopNM4, 4): 
                Wave_4, Energy_4, EQE_4, log_EQE_4 = self.calculate_EQE(self.ref_4, self.data_4, startNM4, stopNM4) 
                export_4 = pd.DataFrame({'Wavelength': Wave_4, 'Energy': Energy_4, 'EQE': EQE_4, 'Log_EQE': log_EQE_4})
                wave_inc['4'] = Wave_4[0]
            else:
                ok_4 = False                   
        
        if self.ui.exportBox_5.isChecked():
            startNM5 = self.ui.startNM_5.value() 
            stopNM5 = self.ui.stopNM_5.value() 
            if self.is_valid(self.ref_5, self.data_5, startNM5, stopNM5, 5):
                Wave_5, Energy_5, EQE_5, log_EQE_5 = self.calculate_EQE(self.ref_5, self.data_5, startNM5, stopNM5)
                export_5 = pd.DataFrame({'Wavelength': Wave_5, 'Energy': Energy_5, 'EQE': EQE_5, 'Log_EQE': log_EQE_5}) 
                wave_inc['5'] = Wave_5[0]
            else:
                ok_5 = False   

       
        if ok_1 and ok_2 and ok_3 and ok_4 and ok_5: # Check if all operations are empty or if fields are left empty
            
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
 
# -----------------------------------------------------------------------------------------------------------        

    #### Function to plot EQE data

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
        
        if number == 0:
            self.set_up_EQE_plot(0)
        elif number == 1:
            self.set_up_EQE_plot(1)
                       
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
            ok_EQE_9 = self.plot_EQE(self.EQE_9, self.ui.startEQE_9, self.ui.stopEQE_9, self.ui.textBox_34, self.ui.textBox_35, self.ui.textBox_36, 1, number)        
                
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

       startEQE = startNM.value() # Pick start wavelength
       stopEQE = stopNM.value() # Pick stop wavelength  
       
       if self.EQE_is_valid(eqe_df, startEQE, stopEQE, file_no): # Check that files are non-empty and within wavelength range
           wave, energy, eqe, log_eqe = self.compile_EQE(eqe_df, startEQE, stopEQE)
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
                            
    ### Function to select EQE data  
  
    def compile_EQE(self, eqe_df, startNM, stopNM):
        
        Wavelength = []
        Energy = []
        EQE = []
        log_EQE = []
        
        for y in range(len(eqe_df['Wavelength'])): # Iterate through columns of EQE file                        
            if startNM <= eqe_df['Wavelength'][y] <= stopNM: # Compile EQE only if start <= wavelength <= stop, otherwise ignore
                Wavelength.append(eqe_df['Wavelength'][y])
                Energy.append(eqe_df['Energy'][y])          
                EQE.append(eqe_df['EQE'][y])   
                log_EQE.append(eqe_df['Log_EQE'][y])                  
        
        if len(Wavelength) == len(EQE) and len(Energy) == len(log_EQE): # Check that the lengths are the same
            return Wavelength, Energy, EQE, log_EQE
            
        else:
            print('Error Code 1: Length mismatch.')  
       

# -----------------------------------------------------------------------------------------------------------        
                
    def pick_Label(self, range_no, startNM, stopNM):
        
        label = str("Range") + str(range_no) + "_" + str(int(startNM)) + "nm" + "-" + str(int(stopNM)) + "nm"
        return label
                
# -----------------------------------------------------------------------------------------------------------        
                
    def pick_EQE_Label(self, label_Box, filename_Box):
        
        label = label_Box.toPlainText()
        filename = filename_Box.toPlainText()
                
        if len(label) != 0:
            return label
        else: # We don't need to check that there is a filename, as the "pick_label" function is only called after checking "EQE_is_valid"
            return filename
                
# -----------------------------------------------------------------------------------------------------------        
 
    def pick_EQE_Color(self, colour_Box, file_no):
        
        colour = colour_Box.toPlainText() 
        colour = colour.replace(" ", "") # If the user inputs "Sky Blue" instead of "SkyBlue" etc.     
#        print("Colour: ", colour)
        
        random_colour = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
#        print(random_colour)
                
        if len(colour) != 0:
            if self.is_Colour(colour):
                return colour
            else:
                print('Please name a valid colour for EQE File %s.' % str(file_no))
                return random_colour
        else:
            return random_colour  
            
# -----------------------------------------------------------------------------------------------------------                  

    def is_Colour(self, colour):
        
        try:
            Color(colour)
            return True
        except:
            return False             
               
# -----------------------------------------------------------------------------------------------------------                  

    #### Helper Functions

# -----------------------------------------------------------------------------------------------------------
       
    ### Function to check if reference & data files are non-empty and within wavelength range              
            
    def is_valid(self, ref_df, data_df, startNM, stopNM, range_no):
        
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
                
                
    ### Function to check if EQE files are non-empty and within wavelength range              
            
    def EQE_is_valid(self, eqe_df, startNM, stopNM, EQE_no):
        
        if len(eqe_df): 
            
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
               
    ### Function to set up plot
               
    def set_up_plot(self):                

#        style.use('ggplot')
        fig1 = plt.figure()
                    
        self.ax1 = fig1.add_subplot(2,1,1)
        plt.ylabel('EQE', fontsize=17, fontweight='medium')              
        plt.grid(True)
#        plt.box()
#        plt.title('EQE vs. Wavelength', fontsize=17, fontweight='medium')
#        plt.title('EQE vs. Energy', fontsize=17, fontweight='medium')
        plt.tick_params(labelsize=14)
        plt.minorticks_on()
        plt.rcParams['figure.facecolor']='xkcd:white'
        plt.rcParams['figure.edgecolor']='xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)

        self.ax2 = fig1.add_subplot(2,1,2)
        self.ax2.set_yscale('log') # To generate log scale axis
        plt.xlabel('Wavelength (nm)', fontsize=17, fontweight='medium')
#        plt.xlabel('Energy (eV)', fontsize=17, fontweight='medium') 
        plt.ylabel('Log(EQE)', fontsize=17, fontweight='medium')              
        plt.grid(True)
#        plt.box()
        plt.tick_params(labelsize=14)
        plt.minorticks_on()
        plt.rcParams['figure.facecolor']='xkcd:white'
        plt.rcParams['figure.edgecolor']='xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
        
        plt.show()
        
    
    ### Function to set up plot
               
    def set_up_EQE_plot(self, number):                

#        style.use('ggplot')

        fig3, self.ax3 = plt.subplots()        
        plt.ylabel('EQE', fontsize=17, fontweight='medium')   
        if number == 0:
            plt.xlabel('Wavelength (nm)', fontsize=17, fontweight='medium')
        elif number == 1:
            plt.xlabel('Energy (eV)', fontsize=17, fontweight='medium')            
#        plt.grid(True)
#        plt.box()
        plt.tick_params(labelsize=14)
        plt.minorticks_on()
        plt.rcParams['figure.facecolor']='xkcd:white'
        plt.rcParams['figure.edgecolor']='xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
#        plt.show()

        fig4, self.ax4 = plt.subplots()        
        self.ax4.set_yscale('log') # To generate log scale axis       
        plt.ylabel('Log EQE', fontsize=17, fontweight='medium') 
        if number == 0:
            plt.xlabel('Wavelength (nm)', fontsize=17, fontweight='medium')
        elif number == 1:
            plt.xlabel('Energy (eV)', fontsize=17, fontweight='medium')              
#        plt.grid(True)
#        plt.box()
        plt.tick_params(labelsize=14)
        plt.minorticks_on()
        plt.rcParams['figure.facecolor']='xkcd:white'
        plt.rcParams['figure.edgecolor']='xkcd:white'
        plt.tick_params(labelsize=15, direction='in', axis='both', which='major', length=8, width=2)
        plt.tick_params(labelsize=15, direction='in', axis='both', which='minor', length=4, width=2)
#        plt.show()
        
# ----------------------------------------------------------------------------------------------------------- 
        
    ### Function to clear plot        
        
    def clear_plot(self):
        plt.close() # Close the current plot
        self.set_up_plot() # Set up a new plot, this is preferred over plt.clf() in case the plot window was closed
#        self.do_plot = True        
        
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
        
def main():

  app = QtWidgets.QApplication(sys.argv)
  monoUI = MainWindow()
  monoUI.show()
  sys.exit(app.exec_())

if __name__ == "__main__":
  main()



