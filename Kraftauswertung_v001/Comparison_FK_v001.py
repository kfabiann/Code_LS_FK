"""

Modified: 04.10.2021
@author: Fabian Kneubuehler

Theme: Comparison of process forces 

#-*- coding:utf-8 -*-
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import sys
import openpyxl
from openpyxl import load_workbook
import numpy_indexed as npi

##########################################################################
##########################################################################
# General settings

# https://matplotlib.org/users/usetex.html
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Float format constrained to two decimal places
float_formatter = "{:.2f}".format
np.set_printoptions(formatter={'float_kind': float_formatter})

##########################################################################
##########################################################################
# Plot selection
CuttingSpeed = 1
Feed = 1
MeanValues = 0
MedianValues = 1

##########################################################################
##########################################################################
# Load data

# Load Excel file with process parameters
filename_Excel = "Zerspanungsversuche_Hobeln_Protokoll_Gesamt_03102021.xlsx"
workbook = load_workbook(
    filename=filename_Excel, data_only=True)
MainSheet = workbook.worksheets[0]

# Load numpy array
Datafile_path = r"Evaluation_files\Evaluation_Repetions\ProcessForceEvaluationRepetions.npy"

Force_eval = np.load(Datafile_path)

##########################################################################
##########################################################################
# Divide array into subarrays

Fc1 = Force_eval[Force_eval[:, 0] == 5, :]
Fc2 = Force_eval[Force_eval[:, 0] == 10, :]
Fc3 = Force_eval[Force_eval[:, 0] == 20, :]
Fc4 = Force_eval[Force_eval[:, 0] == 30, :]
Ff1 = Force_eval[Force_eval[:, 1] == 0.01, :]
Ff2 = Force_eval[Force_eval[:, 1] == 0.05, :]
Ff3 = Force_eval[Force_eval[:, 1] == 0.1, :]
Ff4 = Force_eval[Force_eval[:, 1] == 0.2, :]
Ff5 = Force_eval[Force_eval[:, 1] == 0.3, :]


Fc = np.array([Fc1, Fc2, Fc3, Fc4])
Ff = np.array([Ff1, Ff2, Ff3, Ff4, Ff5])


##########################################################################
##########################################################################

if (CuttingSpeed != 0):

    title = "Process forces over feed"

    fig, ax = plt.subplots(constrained_layout=True)

    if (MeanValues != 0):
        # Cutting force
        for iii in range(len(Fc)):
            ax.errorbar(Fc[iii][:, 1], Fc[iii][:, 2], Fc[iii][:, 3], label=(
                'Mean Fc $v_c$=' + str(Fc[iii][0, 0]) + ' m/min'), linestyle='-', capsize=10)
        # Feed force
        for iii in range(len(Fc)):
            ax.errorbar(Fc[iii][:, 1], Fc[iii][:, 7], Fc[iii][:, 8], label=(
                'Mean Ff $v_c$=' + str(Fc[iii][0, 0]) + ' m/min'), linestyle='-', capsize=10)
        # Passive force
        # for iii in range(len(Fc)):
        #     ax.errorbar(Fc[iii][:, 1], Fc[iii][:,12], Fc[iii][:, 13], label=(
        #         'Mean Fp $v_c$=' + str(Fc[iii][0, 0]) + ' m/min'), linestyle='-', capsize=10)

    if (MedianValues != 0):

        # Cutting force
        for iii in range(len(Fc)):
            lower_q = Fc[iii][:, 4] - Fc[iii][:, 5]  # lower quantile
            upper_q = Fc[iii][:, 6] - Fc[iii][:, 4]  # upper quantile
            ax.errorbar(Fc[iii][:, 1], Fc[iii][:, 4], [lower_q, upper_q], label=(
                'Median Fc $v_c$=' + str(Fc[iii][0, 0]) + ' m/min'), linestyle='-', capsize=10)

        # Feed force
        for iii in range(len(Fc)):
            lower_q = Fc[iii][:, 9] - Fc[iii][:, 10]  # lower quantile
            upper_q = Fc[iii][:, 11] - Fc[iii][:, 9]  # upper quantile
            ax.errorbar(Fc[iii][:, 1], Fc[iii][:, 9], [lower_q, upper_q], label=(
                'Median Ff $v_c$=' + str(Fc[iii][0, 0]) + ' m/min'), linestyle='-', capsize=10)

        # Passive force
        # for iii in range(len(Fc)):
        #     lower_q = Fc[iii][:, 14] - Fc[iii][:, 15]  # lower quantile
        #     upper_q = Fc[iii][:, 16] - Fc[iii][:, 14]  # upper quantile
        #     ax.errorbar(Fc[iii][:, 1], Fc[iii][:, 14], [lower_q, upper_q], label=(
        #         'Median Fp $v_c$=' + str(Fc[iii][0, 0]) + ' m/min'), linestyle='-', capsize=10)

    # Plot settings
    ax.set_xlabel(r'$f$ [mm]')
    ax.set_ylabel(r'Force [N/mm]')
    ax.set_title(title)
    #ax.legend(loc="upper left")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # Save plot
    title_mod = title.replace(" ", "_")

    if (MeanValues != 0):
        title_mod = title_mod + "_mean"
    if (MedianValues != 0):
        title_mod = title_mod + "_median"

    plt.savefig(r"Figures\PNG\Comparisons\ " + title_mod + ".png", dpi=300)
    plt.savefig(r"Figures\EPS\Comparisons\ " + title_mod + ".pdf", dpi=300)

    plt.show()

if (Feed != 0):

    title = "Process forces over cutting speed"
    fig, ax = plt.subplots(constrained_layout=True)

    if (MeanValues != 0):
        # Cutting force
        for iii in range(len(Ff)):
            ax.errorbar(Ff[iii][:, 0], Ff[iii][:, 2], Ff[iii][:, 3], label=(
                'Mean Fc $v_c$=' + str(Ff[iii][0, 1]) + ' mm'), linestyle='-', capsize=10)
        # Feed force
        for iii in range(len(Ff)):
            ax.errorbar(Ff[iii][:, 0], Ff[iii][:, 7], Ff[iii][:, 8], label=(
                'Mean Ff $v_c$=' + str(Ff[iii][0, 1]) + ' mm'), linestyle='-', capsize=10)
        # Passive force
        # for iii in range(len(Ff)):
        #     ax.errorbar(Ff[iii][:, 0], Ff[iii][:,12], Ff[iii][:, 13], label=(
        #         'Mean Fp $v_c$=' + str(Ff[iii][0, 1]) + ' mm'), linestyle='-', capsize=10)

    if (MedianValues != 0):

        # Cutting force
        for iii in range(len(Ff)):
            lower_q = Ff[iii][:, 4] - Ff[iii][:, 5]  # lower quantile
            upper_q = Ff[iii][:, 6] - Ff[iii][:, 4]  # upper quantile
            ax.errorbar(Ff[iii][:, 0], Ff[iii][:, 4], [lower_q, upper_q], label=(
                'Median Fc $v_c$=' + str(Ff[iii][0, 1]) + ' mm'), linestyle='-', capsize=10)

        # Feed force
        for iii in range(len(Ff)):
            lower_q = Ff[iii][:, 9] - Ff[iii][:, 10]  # lower quantile
            upper_q = Ff[iii][:, 11] - Ff[iii][:, 9]  # upper quantile
            ax.errorbar(Ff[iii][:, 0], Ff[iii][:, 9], [lower_q, upper_q], label=(
                'Median Ff $v_c$=' + str(Ff[iii][0, 1]) + ' mm'), linestyle='-', capsize=10)

    # Plot settings
    ax.set_xlabel(r'$v_c$ [m/min]')
    ax.set_ylabel(r'Force [N/mm]')
    ax.set_title(title)
    #ax.legend(loc="upper left")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # Save plot
    title_mod = title.replace(" ", "_")
    if (MeanValues != 0):
        title_mod = title_mod + "_mean"
    if (MedianValues != 0):
        title_mod = title_mod + "_median"

    plt.savefig(r"Figures\PNG\Comparisons\ " + title_mod + ".png", dpi=300)
    plt.savefig(r"Figures\EPS\Comparisons\ " + title_mod + ".pdf", dpi=300)

    plt.show()
