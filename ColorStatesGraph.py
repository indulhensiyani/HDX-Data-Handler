from __future__ import division
import pandas as pd
import os
import sys
import csv
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.interpolate import make_interp_spline
from scipy.interpolate import interp1d
import numpy as np
import seaborn as sns

# define the input and output
filepathOutput = ''  # "C:\\path\\to\\output\\plot\\directory\\"
filepathInput_SD = ''  # "C:\\path\\to\\experimental\\state\\data\\"

# Import csv files
inputcsvfileSASA = ''  # "csv_filename_for_SASA_linear"
inputcsvfileState_data = ''  # "csv_filename_for_experimental_data"
inputcsvfileMaxSASA = ''  # "csv_filename_for_maximum_SASA"

# Read the file, filter by Sequence and provide the PeptideList
State_data = pd.read_csv(filepathInput_SD + inputcsvfileState_data + '.csv', header=0)
State_data["Start-End"] = State_data["Start"].astype(str) + "-" + State_data["End"].astype(
    str)  # or use >>> pandas.concat([df['a'], df['b']]).unique()array(['a', 'b', 'd', 'e', 'g', 'h'], dtype=object)
StartEnd_SD = State_data['Start-End'].unique()
uniquePep_SD = State_data['Sequence'].unique()
uniqueState_SD = State_data['State'].unique()

# Read the file, filter by Sequence and provide the PeptideList
SASAlist = pd.read_csv(filepathOutput + inputcsvfileSASA + '.csv', header=0)

SASAmax = pd.read_csv(filepathOutput + inputcsvfileMaxSASA + '.csv', header=0)

for i, p in enumerate(uniquePep_SD):
    Seqname = SASAlist.iloc[i, 0]
    StartEnd = SASAlist.iloc[i, 3]
    MaxSASA = "{:.2f}".format(SASAmax.iloc[i, 4] * 2)

    fig, (plt3) = plt.subplots(ncols=1)

    # GLyPb plot settings
    for s in uniqueState_SD:
        x2 = State_data.loc[(((State_data['Sequence']) == p) & ((State_data['State']) == s)), 'Exposure']
        y2 = State_data.loc[(((State_data['Sequence']) == p) & ((State_data['State']) == s)), 'Uptake']
        max_x2 = max(x2)
        min_x2 = min(x2)
        print(x2)
        MaxD = State_data.loc[(((State_data['Sequence']) == p) & ((State_data['State']) == s)), 'MaxUptake'].unique()

        plt3.grid(color='lightgrey')
        plt3.plot(x2, y2, label=s, linewidth=3)
        plt3.scatter(x2, y2, linewidth=3, s=75)
        plt3.tick_params(axis='both', which='major', labelsize=16)
        plt3.tick_params(axis='x', colors="black", pad=2, labelsize=20, length=10)
        plt3.tick_params(axis='y', colors="black", pad=2, labelsize=20, length=10)
        plt3.xaxis.set_major_locator(ticker.MultipleLocator(1))
        plt3.set_xscale("log")
        plt3.set_xlim(left=min_x2, right=max_x2)
        plt3.set_xlabel("Time (ms)", fontsize=20)
        plt3.set_ylabel("Uptake (Da)", fontsize=20)
        plt3.set_ylim(0, MaxD)
        plt3.annotate(" and ".join([" MaxSASA "]) + " and ".join([" = "]) + str(MaxSASA), xy=(0.05, 0.9),
                      xycoords='axes fraction', fontsize=16)

    plt3.legend(loc='upper right', fontsize=12)
    fig.savefig(filepathOutput + StartEnd + Seqname + '.png')
