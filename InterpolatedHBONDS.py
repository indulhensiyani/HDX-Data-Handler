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
inputcsvfileHbond = ''  # "csv_filename_for_Hbond_number_linear"
inputcsvfileHbond1 = ''  # "csv_filename_for_Hbond_number_rigimol"
inputcsvfileHbond2 = ''  # "csv_filename_for_Hbond_number_climber"

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
Hbondnumber = pd.read_csv(filepathOutput + inputcsvfileHbond + '.csv', header=0)
Hbondnumber1 = pd.read_csv(filepathOutput + inputcsvfileHbond1 + '.csv', header=0)
Hbondnumber2 = pd.read_csv(filepathOutput + inputcsvfileHbond2 + '.csv', header=0)

SASAmax = pd.read_csv(filepathOutput + inputcsvfileMaxSASA + '.csv', header=0)

for i, p in enumerate(uniquePep_SD):
    Seqname = SASAlist.iloc[i, 0]
    StartEnd = SASAlist.iloc[i, 3]
    MaxSASA = "{:.2f}".format(SASAmax.iloc[i, 4] * 2)

    # Hbonds list
    XHlin = list(range(0, len(Hbondnumber.iloc[0, 4:]), 1))
    YHlin = Hbondnumber.iloc[i, 4:]
    XHrig = list(range(0, len(Hbondnumber1.iloc[0, 4:]), 1))
    YHrig = Hbondnumber1.iloc[i, 4:]
    XHcli = list(range(0, len(Hbondnumber2.iloc[0, 4:]), 1))
    YHcli = Hbondnumber2.iloc[i, 4:]

    # Hbonds cubic interpolation
    arrayHbondlinx = np.array(XHlin)
    arrayHbondliny = np.array(YHlin)
    arrayHbondrigx = np.array(XHrig)
    arrayHbondrigy = np.array(YHrig)
    arrayHbondclix = np.array(XHcli)
    arrayHbondcliy = np.array(YHcli)

    linearHcubic = interp1d(arrayHbondlinx, arrayHbondliny, kind="cubic")
    rigimolHcubic = interp1d(arrayHbondrigx, arrayHbondrigy, kind="cubic")
    climberHcubic = interp1d(arrayHbondclix, arrayHbondcliy, kind="cubic")

    fig, (plt2) = plt.subplots(ncols=1)

    # Hbond plot interpolation
    curveHLx = np.linspace(arrayHbondlinx.min(), arrayHbondlinx.max(), 50)
    curveHLy = linearHcubic(curveHLx)
    curveHRx = np.linspace(arrayHbondrigx.min(), arrayHbondrigx.max(), 50)
    curveHRy = rigimolHcubic(curveHRx)
    curveHCx = np.linspace(arrayHbondclix.min(), arrayHbondclix.max(), 50)
    curveHCy = climberHcubic(curveHCx)

    plt2.grid(color='lightgrey')
    plt2.plot(curveHLx, curveHLy, linewidth=2, color='#6D4F87', label='Linear')
    plt2.plot(curveHRx, curveHRy, linewidth=2, color='#23CA9C', label='Rigimol')
    plt2.plot(curveHCx, curveHCy, linewidth=2, color='#6D83F1', label='Climber')
    plt2.xaxis.set_major_locator(ticker.MultipleLocator(5))
    plt2.set_xlim(left=0, right=30)

    # Hbond plot settings
    plt2.set_title([StartEnd, Seqname], fontsize=24)
    plt2.set_xlabel('#N of state', fontdict=None, labelpad=None, loc=None, fontsize=20)
    plt2.set_ylabel("#N of H bonds per peptide", color="black", fontsize=20)
    plt2.annotate(" and ".join([" MaxSASA "]) + " and ".join([" = "]) + str(MaxSASA), xy=(0.05, 0.05),
                  xycoords='axes fraction', fontsize=16)
    plt2.tick_params(axis='both', which='major', labelsize=16)
    plt2.tick_params(axis='x', colors="black", pad=2, labelsize=20, length=10)
    plt2.tick_params(axis='y', colors="black", pad=2, labelsize=20, length=10)

    plt2.legend(loc='upper right', fontsize=12)
    fig.savefig(filepathOutput + StartEnd + Seqname + '.png')
