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
inputcsvfileSASA1 = ''  # "csv_filename_for_SASA_rigimol"
inputcsvfileSASA2 = ''  # "csv_filename_for_SASA_climber"

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
SASAlist1 = pd.read_csv(filepathOutput + inputcsvfileSASA1 + '.csv', header=0)
SASAlist2 = pd.read_csv(filepathOutput + inputcsvfileSASA2 + '.csv', header=0)

SASAmax = pd.read_csv(filepathOutput + inputcsvfileMaxSASA + '.csv', header=0)

for i, p in enumerate(uniquePep_SD):
    Seqname = SASAlist.iloc[i, 0]
    StartEnd = SASAlist.iloc[i, 3]
    MaxSASA = "{:.2f}".format(SASAmax.iloc[i, 4] * 2)

    # SASA list
    XSlin = list(range(0, len(SASAlist.iloc[0, 4:]), 1))
    YSlin = SASAlist.iloc[i, 4:]
    XSrig = list(range(0, len(SASAlist1.iloc[0, 4:]), 1))
    YSrig = SASAlist1.iloc[i, 4:]
    XScli = list(range(0, len(SASAlist2.iloc[0, 4:]), 1))
    YScli = SASAlist2.iloc[i, 4:]

    # SASA cubic interpolation
    arraySASAlinx = np.array(XSlin)
    arraySASAliny = np.array(YSlin)
    arraySASArigx = np.array(XSrig)
    arraySASArigy = np.array(YSrig)
    arraySASAclix = np.array(XScli)
    arraySASAcliy = np.array(YScli)

    linearScubic = interp1d(arraySASAlinx, arraySASAliny, kind="cubic")
    rigimolScubic = interp1d(arraySASArigx, arraySASArigy, kind="cubic")
    climberScubic = interp1d(arraySASAclix, arraySASAcliy, kind="cubic")

    fig, (plt1) = plt.subplots(ncols=1)

    # SASA plot interpolation
    curveSLx = np.linspace(arraySASAlinx.min(), arraySASAlinx.max(), 500)
    curveSLy = linearScubic(curveSLx)
    curveSRx = np.linspace(arraySASArigx.min(), arraySASArigx.max(), 500)
    curveSRy = rigimolScubic(curveSRx)
    curveSCx = np.linspace(arraySASAclix.min(), arraySASAclix.max(), 500)
    curveSCy = climberScubic(curveSCx)

    # himo ta ba ug balangaw char
    plt1.grid(color='lightgrey')
    plt1.plot(curveSLx, curveSLy, linewidth=2, color='#6D4F87', label='Linear')
    plt1.plot(curveSRx, curveSRy, linewidth=2, color='#23CA9C', label='Rigimol')
    plt1.plot(curveSCx, curveSCy, linewidth=2, color='#6D83F1', label='Climber')
    plt1.xaxis.set_major_locator(ticker.MultipleLocator(5))
    plt1.set_xlim(left=0, right=30)

    # SASA plot settings
    plt1.set_title([StartEnd, Seqname], fontsize=24)
    plt1.set_xlabel('#N of state', fontdict=None, labelpad=None, loc=None, fontsize=20)
    plt1.set_ylabel("SASA per peptide", color="black", fontsize=20)
    plt1.annotate(" and ".join([" MaxSASA "]) + " and ".join([" = "]) + str(MaxSASA), xy=(0.05, 0.05),
                  xycoords='axes fraction', fontsize=16)
    plt1.tick_params(axis='both', which='major', labelsize=16)
    plt1.tick_params(axis='x', colors="black", pad=2, labelsize=20, length=10)
    plt1.tick_params(axis='y', colors="black", pad=2, labelsize=20, length=10)

    plt1.legend(loc='upper right', fontsize=12)
    fig.savefig(filepathOutput + StartEnd + Seqname + '.png')
