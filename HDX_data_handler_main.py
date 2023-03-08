#!/usr/bin/python3
from __future__ import division

import __main__
import glob
import os
import sys
import subprocess
import tkinter as tk
from tkinter import ttk, filedialog

import numpy as np
import pandas as pd
import pymol.cmd
from epymol import rigimol
from pymol import *
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.interpolate import interp1d

__main__.pymol_argv = ['pymol']

# input filepath for 'script directory' folder
script_directory = ''  # C://path//to//HDX_data_handler_repo//script_directory//'


class int_two_test:

    def __init__(self, master=None):
        super().__init__()
        # build ui
        self.IntMain = tk.Tk() if master is None else tk.Toplevel(master)
        self.IntMain.configure(
            borderwidth=7,
            height=200,
            padx=15,
            pady=15,
            relief="groove",
            background="#0F5961")
        self.IntMain.geometry("500x450")
        self.IntMain.minsize(600, 550)
        self.IntMain.resizable(True, True)
        self.IntMain.title("HDX Data Handler")
        self.IntTitle = ttk.Label(self.IntMain)
        self.IntTitle.configure(
            background="#4FB9B6",
            font="{Segoe UI} 16 {}",
            justify="center",
            padding=10,
            relief="ridge",
            text='HDX Data Handler',
            width=22)
        self.IntTitle.grid(column=0, columnspan=2, padx=10, row=0, sticky="w")

        self.fileSASA = tk.Entry(self.IntMain)
        self.fileSASA.grid_forget()

        self.fileSASAshort = tk.Entry(self.IntMain)
        self.fileSASAshort.configure(background='#AEB2B5', foreground='#404040', justify="center", width=45)
        _text_ = 'enter filepath for state data (.csv)'
        self.fileSASAshort.delete("0", "end")
        self.fileSASAshort.insert("0", _text_)
        self.fileSASAshort.grid(column=1, columnspan=2, ipadx=2, ipady=2, row=3)

        # browse for state data (csv) files
        def browseStateD():
            self.fileSASAshort.delete("0", "end")
            csvfiletypes = [('CSV files', "*.csv"), ('All', "*.*")]
            file = filedialog.askopenfile(mode='r', filetypes=csvfiletypes)
            if file:
                pathStateD = os.path.abspath(file.name)
                self.fileSASA.insert(0, str(pathStateD))
                self.fileSASA.grid_forget()
                shortpathStateD = os.path.basename(os.path.abspath(file.name))
                self.fileSASAshort.insert(0, str(shortpathStateD))
                self.fileSASAshort.configure(state='disabled')

            # to retrieve state filepath, use:
            # variable name = (str(self.fileSASA.get()))

            # to retrieve state filename only, use:
            # variable name = (str(self.fileSASAshort.get()))

        BTcolor = ttk.Style()
        BTcolor.theme_use('alt')
        BTcolor.configure('TButton', background='#7FC7C7', foreground='black',
                          borderwidth=2, focuscolor='none', relief='groove')
        BTcolor.configure('A.TButton', background='#4FB9B6', foreground='black',
                          borderwidth=2, focuscolor='none', relief='groove')
        BTcolor.map('TButton', background=[('active', '#4FB9B6')])

        self.browseStateD = ttk.Button(self.IntMain, command=browseStateD)
        self.browseStateD.configure(text='Browse State Data', width=17)
        self.browseStateD.grid(
            column=0,
            ipadx=2,
            ipady=2,
            padx=10,
            pady=20,
            row=3)

        # graph three plots
        # still working on using the jupyter ver instead of the python ver
        def plotSASA():
            subprocess.run(['python', script_directory + 'InterpolatedSASA.py'],
                           check=True)

        self.PlotSASA = ttk.Button(self.IntMain)
        self.PlotSASA.configure(text='Plot SASA', width=11, command=plotSASA)
        self.PlotSASA.grid(column=3,
                           ipadx=2,
                           ipady=2,
                           padx=5,
                           row=13,
                           sticky="w")

        # graph Hbonds per peptide
        def plotHbonds():
            subprocess.run(['python', script_directory + 'InterpolatedHBONDS.py'],
                           check=True)

        self.PlotHBONDS = ttk.Button(self.IntMain)
        self.PlotHBONDS.configure(text='Plot H Bonds', width=11, command=plotHbonds)
        self.PlotHBONDS.grid(
            column=3,
            ipadx=2,
            ipady=2,
            padx=5,
            row=14,
            sticky="w")

        # graph three plots
        # still working on using the jupyter ver instead of the python ver
        def plotStates():
            subprocess.run(['python', script_directory + 'ColorStatesGraph.py'],
                           check=True)

        self.PlotSTATES = ttk.Button(self.IntMain)
        self.PlotSTATES.configure(text='Plot States', width=11, command=plotStates)
        self.PlotSTATES.grid(column=3,
                             ipadx=2,
                             ipady=2,
                             padx=5,
                             row=15,
                             sticky="w")

        # launches PyMOL and runs linear morphing
        def morphLinear():
            linearPDBdir = str(self.LinearDirectoryFull.get())

            state1 = str(self.state1.get())
            state2 = str(self.state2.get())

            pymol.finish_launching()
            pymol.cmd.delete("all")
            pymol.cmd.cd(linearPDBdir)

            pymol.cmd.fetch(code=state1, async_=0)
            pymol.cmd.fetch(code=state2, async_=0)

            pymol.cmd.remove(selection=state1 + ' and not c. A')
            pymol.cmd.remove(selection=state2 + ' and not c. A')
            pymol.cmd.remove('solvent or het')

            pymol.cmd.align(state1, state2, cycles=0)

            pymol.cmd.morph(name='linearout', sele1=state1, sele2=state2, refinement=2, method='linear')

            pymol.cmd.hide(selection=state1)
            pymol.cmd.hide(selection=state2)

        self.MorphLin = ttk.Button(self.IntMain, command=morphLinear)
        self.MorphLin.configure(text='Morph Linear', width=13)
        self.MorphLin.grid(column=2, ipadx=2, ipady=2, pady=5, row=6)

        # launches PyMOL and runs rigimol morphing
        def morphRigimol():
            rigimolPDBdir = str(self.RigimolDirectoryFull.get())

            state1 = str(self.state1.get())
            state2 = str(self.state2.get())

            pymol.finish_launching()
            pymol.cmd.delete("all")
            pymol.cmd.cd(rigimolPDBdir)

            pymol.cmd.fetch(code=state1, async_=0)
            pymol.cmd.fetch(code=state2, async_=0)

            pymol.cmd.remove(selection=state1 + ' and not c. A')
            pymol.cmd.remove(selection=state2 + ' and not c. A')
            pymol.cmd.remove('solvent or het')

            pymol.cmd.align(state1, state2, object='aln', cycles=0)

            pymol.cmd.create(name='rin', selection=state1 + ' in aln', source_state=1, target_state=1)
            pymol.cmd.create(name='rin', selection=state2 + ' in aln', source_state=1, target_state=2)

            rigimol.morph(source='rin', target='rigmorph', refinement=2)

        self.MorphRig = ttk.Button(self.IntMain, command=morphRigimol)
        self.MorphRig.configure(text='Morph Rigimol', width=13)
        self.MorphRig.grid(column=2, ipadx=2, ipady=2, row=7)

        def morphClimber(climberCode):
            S1 = (str(self.state1.get()))
            S2 = (str(self.state2.get()))

            # cd /mnt/c/path/to/climber
            # export CLIMBERDIR=/mnt/c/path/to/climber
            # mkdir RUNxxx
            # cd RUNxxx
            # curl -o S2_orig.pdb http://www.rcsb.org/pdb/files/S2.pdb
            # curl -o S1_orig.pdb http://www.rcsb.org/pdb/files/S1.pdb
            # awk '/ATOM/ {print}' < S1_orig.pdb > S1.pdb
            # awk '/ATOM/ {print}' < S2_orig.pdb | awk '$5=="A" {print}' > S2.pdb
            # ${CLIMBERDIR}/sh/rms.sh S1 S2 > morphx.al1
            # ${CLIMBERDIR}/sh/morphx.sh S1 S2 50
            # awk '/^imorph=/  {print $2 $12}' S1_S2_050step.log > steps_vs_drms
            # awk '/^imorph=/  {print $12 $19}' S1_S2_050step.log > drms_vs_energy
            pass

        self.MorphCli = ttk.Button(self.IntMain, command=morphClimber)
        self.MorphCli.configure(text='Morph Climber', width=13)
        self.MorphCli.grid(column=2, ipadx=2, ipady=2, pady=5, row=8)
        self.MorphCli.grid_forget()

        # selects and reads state data file
        def selectStateD():
            pass

        self.selectStateD = ttk.Button(self.IntMain, command=selectStateD)
        self.selectStateD.configure(text='Select File')
        self.selectStateD.grid(column=3, ipadx=2, ipady=2, row=3)
        self.selectStateD.grid_forget()

        self.state1 = tk.Entry(self.IntMain)
        self.state1.configure(background='#aeb2b5', foreground='#404040', justify="center", takefocus=False, width=20)
        _text_ = 'Enter State1 (.pdb)'
        self.state1.delete("0", "end")
        self.state1.insert("0", _text_)
        self.state1.grid(column=0, ipadx=2, ipady=2, pady=5, row=1, sticky="e")

        self.state2 = tk.Entry(self.IntMain)
        self.state2.configure(background='#aeb2b5', foreground='#404040', justify="center", width=20)
        _text_ = 'Enter State2 (.pdb)'
        self.state2.delete("0", "end")
        self.state2.insert("0", _text_)
        self.state2.grid(
            column=1,
            ipadx=2,
            ipady=2,
            padx=10,
            pady=5,
            row=1,
            sticky="w")

        # confirm State 1 PDB
        def selectPDB_S1():
            if len(self.state1.get()) == 4:
                self.state1.configure(state='normal')
            else:
                self.state1.delete(0, 'end')
                self.state1.configure(state='normal')

        self.selectState1 = ttk.Button(self.IntMain, command=selectPDB_S1)
        self.selectState1.configure(text='Enter PDB1')
        self.selectState1.grid(column=0, padx=10, row=2, sticky="w")

        # confirm State 1 PDB
        def selectPDB_S2():
            if len(self.state1.get()) == 4:
                self.state2.configure(state='normal')
            else:
                self.state2.delete(0, 'end')
                self.state2.configure(state='normal')

        self.SelectState2 = ttk.Button(self.IntMain, command=selectPDB_S2)
        self.SelectState2.configure(text='Enter PDB2')
        self.SelectState2.grid(column=1,
                               ipadx=2,
                               ipady=2,
                               padx=10,
                               pady=5,
                               row=2,
                               sticky="w")

        # select Linear PDB directory
        def set_LinPDB_dir():
            self.LinearDirectory.delete("0", "end")
            directory = filedialog.askdirectory()
            if directory:
                pathLinPDB = os.path.basename(os.path.abspath(directory))
                self.LinearDirectory.insert(0, str(pathLinPDB))
                pathLinPDB_Full = os.path.abspath(directory)
                self.LinearDirectoryFull.insert(0, str(pathLinPDB_Full) + '\\')
                self.LinearDirectoryFull.grid_forget()

            # to retrieve Linear directory use:
            # linearPDBdir = str(self.LinearDirectoryFull.get())

        self.pdb_dir_Lin = ttk.Button(self.IntMain, command=set_LinPDB_dir)
        self.pdb_dir_Lin.configure(text='Select Folder', width=17)
        self.pdb_dir_Lin.grid(column=0, pady=5, row=6)

        # select Rigimol PDB directory
        def set_RigPDB_dir():
            self.RigimolDirectory.delete("0", "end")
            directory = filedialog.askdirectory()
            if directory:
                pathRigPDB = os.path.basename(os.path.abspath(directory))
                self.RigimolDirectory.insert(0, str(pathRigPDB))
                pathRigPDB_Full = os.path.abspath(directory)
                self.RigimolDirectoryFull.insert(0, str(pathRigPDB_Full))
                self.RigimolDirectoryFull.grid_forget()

            # to retrieve Linear directory use:
            # rigimolPDBdir = str(self.RigimolDirectoryFull.get())

        self.pdb_dir_Rig = ttk.Button(self.IntMain, command=set_RigPDB_dir)
        self.pdb_dir_Rig.configure(text='Select Folder', width=17)
        self.pdb_dir_Rig.grid(column=0, row=7)

        # select Climber PDB directory
        def set_CliPDB_dir():
            self.RigimolDirectory.delete("0", "end")
            directory = filedialog.askdirectory()
            if directory:
                pathCliPDB = os.path.basename(os.path.abspath(directory))
                self.ClimberDirectory.insert(0, str(pathCliPDB))
                pathCliPDB_Full = os.path.abspath(directory)
                self.ClimberDirectoryFull.insert(0, str(pathCliPDB_Full))
                self.ClimberDirectoryFull.grid_forget()

            climberPDBdir = str(self.ClimberDirectoryFull.cget('text'))
            print(climberPDBdir)

        self.pdb_dir_Cli = ttk.Button(self.IntMain, command=set_CliPDB_dir)
        self.pdb_dir_Cli.configure(text='Select Folder', width=17)
        self.pdb_dir_Cli.grid(column=0, pady=5, row=8)
        self.pdb_dir_Cli.grid_forget()

        self.LinearDirectory = tk.Entry(self.IntMain)
        self.LinearDirectory.configure(background='#AEB2B5', foreground='#404040',
                                       justify="left", width=27)
        _text_ = 'Linear Directory'
        self.LinearDirectory.delete("0", "end")
        self.LinearDirectory.insert("0", _text_)
        self.LinearDirectory.grid(column=1, row=6, sticky="w")

        self.LinearDirectoryFull = ttk.Entry(self.IntMain)
        self.LinearDirectoryFull.grid_forget()

        self.RigimolDirectory = tk.Entry(self.IntMain)
        self.RigimolDirectory.configure(background='#AEB2B5', foreground='#404040',
                                        justify="left", width=27)
        _text_ = 'Rigimol Directory'
        self.RigimolDirectory.delete("0", "end")
        self.RigimolDirectory.insert("0", _text_)
        self.RigimolDirectory.grid(column=1, row=7, sticky="w")

        self.RigimolDirectoryFull = ttk.Entry(self.IntMain)
        self.RigimolDirectoryFull.grid_forget()

        self.ClimberDirectory = tk.Entry(self.IntMain)
        self.ClimberDirectory.configure(background='#AEB2B5', foreground='#404040',
                                        justify="left", width=27)
        _text_ = 'Climber Directory'
        self.ClimberDirectory.delete("0", "end")
        self.ClimberDirectory.insert("0", _text_)
        self.ClimberDirectory.grid(column=1, row=8, sticky="w")
        self.ClimberDirectory.grid_forget()

        self.ClimberDirectoryFull = ttk.Entry(self.IntMain)
        self.ClimberDirectoryFull.grid_forget()

        # select CSV input directory
        def set_csv_indir():
            self.csv_inpath.delete("0", "end")
            directory = filedialog.askdirectory()
            if directory:
                csv_dirpath = os.path.basename(os.path.abspath(directory))
                self.csv_inpath.insert(0, (str(csv_dirpath)))
                csv_InPathFull = os.path.abspath(directory)
                self.csv_inpathFull.insert(0, str(csv_InPathFull))
                self.csv_inpathFull.grid_forget()

            # to retrieve Input SASA directory use:
            # variable name = str(self.csv_inpathFull.get())

        self.csvin_dir = ttk.Button(self.IntMain, command=set_csv_indir)
        self.csvin_dir.configure(text='Select Folder', width=17)
        self.csvin_dir.grid(column=0, row=9, pady=15)

        self.csv_inpath = tk.Entry(self.IntMain)
        self.csv_inpath.configure(background='#AEB2B5', foreground='#404040',
                                  justify="left", width=27)
        _text_ = 'CSV Input Directory'
        self.csv_inpath.delete("0", "end")
        self.csv_inpath.insert("0", _text_)
        self.csv_inpath.grid(column=1, row=9, sticky="w")

        self.csv_inpathFull = ttk.Entry(self.IntMain)
        self.csv_inpathFull.grid_forget()

        # select CSV output directory
        def set_csv_outdir():
            self.csv_outpath.delete("0", "end")
            directory = filedialog.askdirectory()
            if directory:
                csv_dirpath = os.path.basename(os.path.abspath(directory))
                self.csv_outpath.insert(0, (str(csv_dirpath)))
                csv_fullpath = os.path.abspath(directory)
                self.csv_outpathFull.insert(0, str(csv_fullpath))
                self.csv_outpathFull.grid_forget()

            # to retrieve Output SASA directory use:
            # variable name = str(self.csv_outpathFull.get())

        self.csvout_dir = ttk.Button(self.IntMain, command=set_csv_outdir)
        self.csvout_dir.configure(text='Select Folder', width=17)
        self.csvout_dir.grid(column=0, row=10)

        self.csv_outpath = tk.Entry(self.IntMain)
        self.csv_outpath.configure(background='#AEB2B5', foreground='#404040',
                                   justify="left", width=27)
        _text_ = 'CSV Output Directory'
        self.csv_outpath.delete("0", "end")
        self.csv_outpath.insert("0", _text_)
        self.csv_outpath.grid(column=1, row=10, sticky="w")

        self.csv_outpathFull = ttk.Entry(self.IntMain)
        self.csv_outpathFull.grid_forget()

        # prints out linear PDBs
        def exportLin():
            linearPDBdir = str(self.LinearDirectoryFull.get())

            pymol.cmd.cd(linearPDBdir)
            pymol.cmd.multifilesave(filename='{name}_{state}.pdb', selection='linearout', state=0, format='pdb')

            pymol.cmd.delete("all")

        self.ExportLin = ttk.Button(self.IntMain, command=exportLin)
        self.ExportLin.configure(text='Export Linear', width=13)
        self.ExportLin.grid(column=3, ipadx=2, ipady=2, row=6)

        # prints out rigimol PDBs
        def exportRig():
            rigimolPDBdir = str(self.RigimolDirectoryFull.get())

            pymol.cmd.cd(rigimolPDBdir)
            pymol.cmd.multifilesave(filename='{name}_{state}.pdb', selection='rigmorph', state=0, format='pdb')

            pymol.cmd.delete("all")

        self.ExportRig = ttk.Button(self.IntMain, command=exportRig)
        self.ExportRig.configure(text='Export Rigimol', width=13)
        self.ExportRig.grid(column=3, ipadx=2, ipady=2, row=7)

        # prints out climber PDBs
        def exportCli():
            pass

        self.ExportCli = ttk.Button(self.IntMain, command=exportCli)
        self.ExportCli.configure(text='Export Climber', width=13)
        self.ExportCli.grid(column=3, ipadx=2, ipady=2, row=8)
        self.ExportCli.grid_forget()

        # calculate SASA for Linear-generated states
        def calculateLinSASA():
            from pymol import cmd
            import pandas as pd
            import os
            import sys
            import csv
            import glob

            # define the input and output
            inputcsvfile = (str(self.fileSASA.get()))
            inputLinearPDB = str(self.LinearDirectoryFull.get())

            outputSASApath = str(self.csv_outpathFull.get())
            outputcsvfile = "Linear_SASA_Values"
            outputcsvfilemax = "Max_SASA_Values"

            PeptideListSASA = pd.read_csv(inputcsvfile, header=0,
                                          usecols=["Sequence", "Start", "End"])
            PeptideListSASA = PeptideListSASA.drop_duplicates('Sequence', keep='first')
            PeptideListSASA = pd.DataFrame(PeptideListSASA, columns=['Sequence', 'Start', 'End'])
            PeptideListSASA["Start-End"] = PeptideListSASA["Start"].astype(str) + "-" + PeptideListSASA["End"].astype(
                str)

            def SASAperPep_def(inputLinearPDB, PeptideListSASA):

                file_list = glob.glob(inputLinearPDB + '*.pdb')
                if file_list:
                    file_list.sort()
                    for count, p in enumerate(file_list):
                        pdb_file = (os.path.split(p)[-1])
                        cmd.load(p)

                        cmd.set("dot_density", 2)  # high sampling density can increase, but it will get slower max 4
                        cmd.set("dot_solvent", 1)  # solvent accessible rather than molecular surface
                        cmd.set("solvent_radius", 1.4)
                        cmd.h_add("all")

                        SASAperPeptide = []
                        for i in range(len(PeptideListSASA)):
                            Selection = " and ".join([" resi "]) + PeptideListSASA.iloc[i, 3]
                            SASAperPeptide.append(cmd.get_area(Selection))
                        PeptideListSASA[str(pdb_file)] = SASAperPeptide
                        cmd.delete("all")

                    PeptideListSASA.to_csv(outputSASApath + '\\' + outputcsvfile + '.csv', index=None, header=True)

            def UpperSASAlimit(PeptideListSASA):

                SASAlimitperPeptide = []

                for i in range(len(PeptideListSASA)):
                    Selection = PeptideListSASA.iloc[i, 0]
                    cmd.fab(Selection, str(i),
                            hydro=1)
                    SASAlimitperPeptide.append(cmd.get_area('all') * 2)
                    cmd.delete(i)
                PeptideListSASA['MaxSASA'] = SASAlimitperPeptide
                PeptideListSASA.to_csv(outputSASApath + '\\' + outputcsvfilemax + '.csv', index=None, header=True)

            cmd.extend("SASAperPeptide_def", SASAperPep_def)
            SASAperPep_def(inputLinearPDB, PeptideListSASA)

            cmd.extend("UpperSASAlimit", UpperSASAlimit)
            UpperSASAlimit(PeptideListSASA)
            self.CalcSASA_Lin.configure(state='disabled')

        self.CalcSASA_Lin = ttk.Button(self.IntMain, command=calculateLinSASA)
        self.CalcSASA_Lin.configure(text='Calculate SASA', width=13)
        self.CalcSASA_Lin.grid(column=0, columnspan=2, ipadx=2, ipady=2, row=13)

        # calculate SASA for Rigimol-generated states
        def calculateRigSASA():
            from pymol import cmd
            import pandas as pd
            import os
            import sys
            import csv
            import glob

            # define the input and output
            inputcsvfile = (str(self.fileSASA.get()))
            inputRigimolPDB = str(self.RigimolDirectoryFull.get())

            outputSASApath = str(self.csv_outpathFull.get())
            outputcsvfile = "Rigimol_SASA_Values"

            PeptideListSASA = pd.read_csv(inputcsvfile, header=0,
                                          usecols=["Sequence", "Start", "End"])
            PeptideListSASA = PeptideListSASA.drop_duplicates('Sequence', keep='first')
            PeptideListSASA = pd.DataFrame(PeptideListSASA, columns=['Sequence', 'Start', 'End'])
            PeptideListSASA["Start-End"] = PeptideListSASA["Start"].astype(str) + "-" + PeptideListSASA["End"].astype(
                str)

            def SASAperPep_def(inputLinearPDB, PeptideListSASA):

                file_list = glob.glob(inputLinearPDB + '*.pdb')
                if file_list:
                    file_list.sort()
                    for count, p in enumerate(file_list):
                        pdb_file = (os.path.split(p)[-1])
                        cmd.load(p)

                        cmd.set("dot_density", 2)  # high sampling density can increase, but it will get slower max 4
                        cmd.set("dot_solvent", 1)  # solvent accessible rather than molecular surface
                        cmd.set("solvent_radius", 1.4)
                        cmd.h_add("all")

                        SASAperPeptide = []
                        for i in range(len(PeptideListSASA)):
                            Selection = " and ".join([" resi "]) + PeptideListSASA.iloc[i, 3]
                            SASAperPeptide.append(cmd.get_area(Selection))
                        PeptideListSASA[str(pdb_file)] = SASAperPeptide
                        cmd.delete("all")

                    PeptideListSASA.to_csv(outputSASApath + '\\' + outputcsvfile + '.csv', index=None, header=True)

            cmd.extend("SASAperPeptide_def", SASAperPep_def)
            SASAperPep_def(inputRigimolPDB, PeptideListSASA)
            self.CalcSASA_Rig.configure(state='normal')

        self.CalcSASA_Rig = ttk.Button(self.IntMain, command=calculateRigSASA)
        self.CalcSASA_Rig.configure(text='Calculate SASA', width=13)
        self.CalcSASA_Rig.grid(column=0, columnspan=2, ipadx=2, ipady=2, pady=5, row=14)

        # calculate SASA for Rigimol-generated states
        def calculateCliSASA():
            from pymol import cmd
            import pandas as pd
            import os
            import sys
            import csv
            import glob

            # define the input and output
            inputcsvfile = (str(self.fileSASA.get()))
            inputClimberPDB = str(self.ClimberDirectoryFull.get())

            outputSASApath = str(self.csv_outpathFull.get())
            outputcsvfile = "Climber_SASA_Values"

            PeptideListSASA = pd.read_csv(inputcsvfile, header=0,
                                          usecols=["Sequence", "Start", "End"])
            PeptideListSASA = PeptideListSASA.drop_duplicates('Sequence', keep='first')
            PeptideListSASA = pd.DataFrame(PeptideListSASA, columns=['Sequence', 'Start', 'End'])
            PeptideListSASA["Start-End"] = PeptideListSASA["Start"].astype(str) + "-" + PeptideListSASA["End"].astype(
                str)

            def SASAperPep_def(inputClimberPDB, PeptideListSASA):

                file_list = glob.glob(inputClimberPDB + '*.pdb')
                if file_list:
                    file_list.sort()
                    for count, p in enumerate(file_list):
                        pdb_file = (os.path.split(p)[-1])
                        cmd.load(p)

                        cmd.set("dot_density", 2)  # high sampling density can increase, but it will get slower max 4
                        cmd.set("dot_solvent", 1)  # solvent accessible rather than molecular surface
                        cmd.set("solvent_radius", 1.4)
                        cmd.h_add("all")

                        SASAperPeptide = []
                        for i in range(len(PeptideListSASA)):
                            Selection = " and ".join([" resi "]) + PeptideListSASA.iloc[i, 3]
                            SASAperPeptide.append(cmd.get_area(Selection))
                        PeptideListSASA[str(pdb_file)] = SASAperPeptide
                        cmd.delete("all")

                    PeptideListSASA.to_csv(outputSASApath + '\\' + outputcsvfile + '.csv', index=None, header=True)

            cmd.extend("SASAperPeptide_def", SASAperPep_def)
            SASAperPep_def(inputClimberPDB, PeptideListSASA)
            self.CalcSASA_Cli.configure(state='normal')

        self.CalcSASA_Cli = ttk.Button(self.IntMain, command=calculateCliSASA)
        self.CalcSASA_Cli.configure(text='Calculate SASA', width=13)
        self.CalcSASA_Cli.grid(column=0, columnspan=2, ipadx=2, ipady=2, row=15)

        # calculate H Bonds for Linear-generated states
        def calculateLinHBONDS():
            import pymol
            from pymol import cmd
            import pandas as pd
            import os
            import sys
            import csv
            import glob

            # define the input and output
            inputcsvfile = (str(self.fileSASA.get()))
            inputLinearPDB = str(self.LinearDirectoryFull.get())

            outputHBONDSpath = str(self.csv_outpathFull.get())
            outputList = "Linear_HBONDS_List"
            outputNum = "Linear_HBONDS_Number"

            # Read the file, filter by Sequence and provide the PeptideList
            PeptideListHbond = pd.read_csv(inputcsvfile, header=0,
                                           usecols=["Sequence", "Start", "End"])
            PeptideListHbond = PeptideListHbond.drop_duplicates('Sequence', keep='first')
            PeptideListHbond = pd.DataFrame(PeptideListHbond, columns=['Sequence', 'Start', 'End'])
            PeptideListHbond["Start:End"] = PeptideListHbond["Start"].astype(str) + "-" + PeptideListHbond[
                "End"].astype(str)

            # Another data frame for H-bond list if needed
            PeptideNumberHbond = pd.DataFrame(PeptideListHbond, columns=['Sequence', 'Start', 'End'])
            PeptideNumberHbond["Start:End"] = PeptideListHbond["Start"].astype(str) + "-" + PeptideListHbond[
                "End"].astype(str)

            def HbondPerPeptide_def(inputLinearPDB, PeptideListHbond):
                # Import pdb files
                file_list = glob.glob(inputLinearPDB + '*.pdb')
                if file_list:
                    file_list.sort()
                for count, p in enumerate(file_list):
                    pdb_file = (os.path.split(p)[-1])
                    pdb_file = os.path.splitext(pdb_file)[0]
                    cmd.load(p)

                    # Second loop for calculating number of H bonds
                    # Function will generate H-bonds per peptide stretch for a given selection
                    cutoff = 3.2
                    angle = 45
                    hb_list_name = 'hbonds'
                    Selection2 = pdb_file + "& e. o"
                    mode = 1

                    NumberHbond = []
                    ListHbond = []
                    for i in range(len(PeptideListHbond)):
                        selection = PeptideListHbond.iloc[i, 3]
                        Selection = " and ".join([" resi "]) + selection + " & name n"

                        hb = cmd.find_pairs(Selection, Selection2, mode=mode, cutoff=cutoff, angle=angle)
                        # convert hb list to set to remove duplicates
                        hb_set = set()
                        for atoms in hb:
                            a = [atoms[0], atoms[1]]
                            a.sort()
                            hb_set.add(tuple(a))
                        # convert set back to list and sort for easier reading
                        hb = list(hb_set)
                        hb.sort(key=lambda x: x[0][1])
                        pymol.stored.listA = []
                        pymol.stored.listB = []
                        pymol.stored.listC = []
                        C = []

                        for pairs in hb:
                            cmd.iterate("%s and index %s" % (pairs[0][0], pairs[0][1]),
                                        'pymol.stored.listA.append( "%1s/%3s`%s/%s/%i " % (chain, resn, resi, name, ID),)')

                            cmd.iterate("%s and index %s" % (pairs[1][0], pairs[1][1]),
                                        'pymol.stored.listB.append( "%1s/%3s`%s/%s/%i " % (chain, resn, resi, name, ID),)')
                            pymol.stored.listC.append(
                                cmd.distance(hb_list_name, "%s and index %s" % (pairs[0][0], pairs[0][1]),
                                             "%s and index %s" % (pairs[1][0], pairs[1][1])))
                            for line in enumerate(pymol.stored.listA):  # there is a problem here ask Joshy :)
                                C.append("%s   %s   %.2f" % (
                                    pymol.stored.listA[line[0]], pymol.stored.listB[line[0]], pymol.stored.listC[line[0]]))

                        ListHbond.append(';'.join(C))
                        NumberHbond.append(len(pymol.stored.listA))

                    PeptideListHbond[str(pdb_file)] = ListHbond
                    PeptideNumberHbond[str(pdb_file)] = NumberHbond
                    cmd.delete("all")

                PeptideListHbond.to_csv(outputHBONDSpath + '\\' + outputList + '.csv', index=None, header=True)
                PeptideNumberHbond.to_csv(outputHBONDSpath + '\\' + outputNum + '.csv', index=None, header=True)
                cmd.extend("HbondPerPeptide_def", HbondPerPeptide_def)

            HbondPerPeptide_def(inputLinearPDB, PeptideListHbond)
            self.CalcHBONDS_Lin.configure(state='normal')

        self.CalcHBONDS_Lin = ttk.Button(self.IntMain, command=calculateLinHBONDS)
        self.CalcHBONDS_Lin.configure(text='Calculate Hbonds', width=15, state='normal')
        self.CalcHBONDS_Lin.grid(column=1, ipadx=2, ipady=2, row=13, sticky="e")

        # calculate H Bonds for Rigimol-generated states
        def calculateRigHBONDS():
            import pymol
            from pymol import cmd
            import pandas as pd
            import os
            import sys
            import csv
            import glob

            # define the input and output
            inputcsvfile = (str(self.fileSASA.get()))
            inputRigimolPDB = str(self.RigimolDirectoryFull.get())

            outputHBONDSpath = str(self.csv_outpathFull.get())
            outputList = "Rigimol_HBONDS_List"
            outputNum = "Rigimol_HBONDS_Number"

            # Read the file, filter by Sequence and provide the PeptideList
            PeptideListHbond = pd.read_csv(inputcsvfile, header=0,
                                           usecols=["Sequence", "Start", "End"])
            PeptideListHbond = PeptideListHbond.drop_duplicates('Sequence', keep='first')
            PeptideListHbond = pd.DataFrame(PeptideListHbond, columns=['Sequence', 'Start', 'End'])
            PeptideListHbond["Start:End"] = PeptideListHbond["Start"].astype(str) + "-" + PeptideListHbond[
                "End"].astype(str)

            # Another data frame for H-bond list if needed
            PeptideNumberHbond = pd.DataFrame(PeptideListHbond, columns=['Sequence', 'Start', 'End'])
            PeptideNumberHbond["Start:End"] = PeptideListHbond["Start"].astype(str) + "-" + PeptideListHbond[
                "End"].astype(str)

            def HbondPerPeptide_def(inputRigimolPDB, PeptideListHbond):
                # Import pdb files
                file_list = glob.glob(inputRigimolPDB + '*.pdb')
                if file_list:
                    file_list.sort()
                for count, p in enumerate(file_list):
                    pdb_file = (os.path.split(p)[-1])
                    pdb_file = os.path.splitext(pdb_file)[0]
                    cmd.load(p)

                    # Second loop for calculating number of H bonds
                    # Function will generate H-bonds per peptide stretch for a given selection
                    cutoff = 3.2
                    angle = 45
                    hb_list_name = 'hbonds'
                    Selection2 = pdb_file + "& e. o"
                    mode = 1

                    NumberHbond = []
                    ListHbond = []
                    for i in range(len(PeptideListHbond)):
                        selection = PeptideListHbond.iloc[i, 3]
                        Selection = " and ".join([" resi "]) + selection + " & name n"

                        hb = cmd.find_pairs(Selection, Selection2, mode=mode, cutoff=cutoff, angle=angle)
                        # convert hb list to set to remove duplicates
                        hb_set = set()
                        for atoms in hb:
                            a = [atoms[0], atoms[1]]
                            a.sort()
                            hb_set.add(tuple(a))
                        # convert set back to list and sort for easier reading
                        hb = list(hb_set)
                        hb.sort(key=lambda x: x[0][1])
                        pymol.stored.listA = []
                        pymol.stored.listB = []
                        pymol.stored.listC = []
                        C = []

                        for pairs in hb:
                            cmd.iterate("%s and index %s" % (pairs[0][0], pairs[0][1]),
                                        'pymol.stored.listA.append( "%1s/%3s`%s/%s/%i " % (chain, resn, resi, name, ID),)')

                            cmd.iterate("%s and index %s" % (pairs[1][0], pairs[1][1]),
                                        'pymol.stored.listB.append( "%1s/%3s`%s/%s/%i " % (chain, resn, resi, name, ID),)')
                            pymol.stored.listC.append(
                                cmd.distance(hb_list_name, "%s and index %s" % (pairs[0][0], pairs[0][1]),
                                             "%s and index %s" % (pairs[1][0], pairs[1][1])))
                            for line in enumerate(pymol.stored.listA):  # there is a problem here ask Joshy :)
                                C.append("%s   %s   %.2f" % (
                                    pymol.stored.listA[line[0]], pymol.stored.listB[line[0]],
                                    pymol.stored.listC[line[0]]))

                        ListHbond.append(';'.join(C))
                        NumberHbond.append(len(pymol.stored.listA))

                    PeptideListHbond[str(pdb_file)] = ListHbond
                    PeptideNumberHbond[str(pdb_file)] = NumberHbond
                    cmd.delete("all")

                PeptideListHbond.to_csv(outputHBONDSpath + '\\' + outputList + '.csv', index=None, header=True)
                PeptideNumberHbond.to_csv(outputHBONDSpath + '\\' + outputNum + '.csv', index=None, header=True)
                cmd.extend("HbondPerPeptide_def", HbondPerPeptide_def)

            HbondPerPeptide_def(inputRigimolPDB, PeptideListHbond)
            self.CalcHBONDS_Rig.configure(state='normal')

        self.CalcHBONDS_Rig = ttk.Button(self.IntMain, command=calculateRigHBONDS)
        self.CalcHBONDS_Rig.configure(text='Calculate Hbonds', width=15)
        self.CalcHBONDS_Rig.grid(column=1, ipadx=2, ipady=2, row=14, sticky="e")

        # calculate H Bonds for Climber-generated states
        def calculateCliHBONDS():
            import pymol
            from pymol import cmd
            import pandas as pd
            import os
            import sys
            import csv
            import glob

            # define the input and output
            inputcsvfile = (str(self.fileSASA.get()))
            inputClimberPDB = str(self.ClimberDirectoryFull.get())

            outputHBONDSpath = str(self.csv_outpathFull.get())
            outputList = "Climber_HBONDS_List"
            outputNum = "Climber_HBONDS_Number"

            # Read the file, filter by Sequence and provide the PeptideList
            PeptideListHbond = pd.read_csv(inputcsvfile, header=0,
                                           usecols=["Sequence", "Start", "End"])
            PeptideListHbond = PeptideListHbond.drop_duplicates('Sequence', keep='first')
            PeptideListHbond = pd.DataFrame(PeptideListHbond, columns=['Sequence', 'Start', 'End'])
            PeptideListHbond["Start:End"] = PeptideListHbond["Start"].astype(str) + "-" + PeptideListHbond[
                "End"].astype(str)

            # Another data frame for H-bond list if needed
            PeptideNumberHbond = pd.DataFrame(PeptideListHbond, columns=['Sequence', 'Start', 'End'])
            PeptideNumberHbond["Start:End"] = PeptideListHbond["Start"].astype(str) + "-" + PeptideListHbond[
                "End"].astype(str)

            def HbondPerPeptide_def(inputClimberPDB, PeptideListHbond):
                # Import pdb files
                file_list = glob.glob(inputClimberPDB + '*.pdb')
                if file_list:
                    file_list.sort()
                for count, p in enumerate(file_list):
                    pdb_file = (os.path.split(p)[-1])
                    pdb_file = os.path.splitext(pdb_file)[0]
                    cmd.load(p)

                    # Second loop for calculating number of H bonds
                    # Function will generate H-bonds per peptide stretch for a given selection
                    cutoff = 3.2
                    angle = 45
                    hb_list_name = 'hbonds'
                    Selection2 = pdb_file + "& e. o"
                    mode = 1

                    NumberHbond = []
                    ListHbond = []
                    for i in range(len(PeptideListHbond)):
                        selection = PeptideListHbond.iloc[i, 3]
                        Selection = " and ".join([" resi "]) + selection + " & name n"

                        hb = cmd.find_pairs(Selection, Selection2, mode=mode, cutoff=cutoff, angle=angle)
                        # convert hb list to set to remove duplicates
                        hb_set = set()
                        for atoms in hb:
                            a = [atoms[0], atoms[1]]
                            a.sort()
                            hb_set.add(tuple(a))
                        # convert set back to list and sort for easier reading
                        hb = list(hb_set)
                        hb.sort(key=lambda x: x[0][1])
                        pymol.stored.listA = []
                        pymol.stored.listB = []
                        pymol.stored.listC = []
                        C = []

                        for pairs in hb:
                            cmd.iterate("%s and index %s" % (pairs[0][0], pairs[0][1]),
                                        'pymol.stored.listA.append( "%1s/%3s`%s/%s/%i " % (chain, resn, resi, name, ID),)')

                            cmd.iterate("%s and index %s" % (pairs[1][0], pairs[1][1]),
                                        'pymol.stored.listB.append( "%1s/%3s`%s/%s/%i " % (chain, resn, resi, name, ID),)')
                            pymol.stored.listC.append(
                                cmd.distance(hb_list_name, "%s and index %s" % (pairs[0][0], pairs[0][1]),
                                             "%s and index %s" % (pairs[1][0], pairs[1][1])))
                            for line in enumerate(pymol.stored.listA):  # there is a problem here ask Joshy :)
                                C.append("%s   %s   %.2f" % (
                                    pymol.stored.listA[line[0]], pymol.stored.listB[line[0]],
                                    pymol.stored.listC[line[0]]))

                        ListHbond.append(';'.join(C))
                        NumberHbond.append(len(pymol.stored.listA))

                    PeptideListHbond[str(pdb_file)] = ListHbond
                    PeptideNumberHbond[str(pdb_file)] = NumberHbond
                    cmd.delete("all")

                PeptideListHbond.to_csv(outputHBONDSpath + '\\' + outputList + '.csv', index=None, header=True)
                PeptideNumberHbond.to_csv(outputHBONDSpath + '\\' + outputNum + '.csv', index=None, header=True)
                cmd.extend("HbondPerPeptide_def", HbondPerPeptide_def)

            HbondPerPeptide_def(inputClimberPDB, PeptideListHbond)
            self.CalcHBONDS_Cli.configure(state='normal')

        self.CalcHBONDS_Cli = ttk.Button(self.IntMain, command=calculateCliHBONDS)
        self.CalcHBONDS_Cli.configure(text='Calculate Hbonds', width=15)
        self.CalcHBONDS_Cli.grid(column=1, ipadx=2, ipady=2, row=15, sticky="e")

        self.CalculateSubheading = ttk.Label(self.IntMain)
        self.CalculateSubheading.configure(
            background="#c0c0c0",
            font="{Calibri Light} 11 {}",
            justify="center",
            padding=10,
            relief="ridge",
            state="normal",
            text='Calculate SASA and H Bonds',
            width=20)
        self.CalculateSubheading.grid(
            column=0,
            columnspan=2,
            padx=10,
            pady=13,
            row=12,
            sticky="w")

        self.GraphsSubheading = ttk.Label(self.IntMain)
        self.GraphsSubheading.configure(
            background="#c0c0c0",
            font="{Calibri Light} 11 {}",
            justify="center",
            padding=10,
            relief="ridge",
            state="normal",
            text='Plot Graphs',
            width=22)
        self.GraphsSubheading.grid(
            column=2,
            columnspan=2,
            padx=20,
            pady=15,
            row=12,
            sticky="w")

        self.SASAgraph = ttk.Label(self.IntMain)
        self.SASAgraph.configure(
            justify="center",
            padding=3,
            state="normal",
            text='SASA',
            width=13)
        self.SASAgraph.grid(column=2, padx=20, row=13, sticky="w")

        self.HBONDgraph = ttk.Label(self.IntMain)
        self.HBONDgraph.configure(
            justify="left",
            padding=3,
            text='H Bonds',
            width=13)
        self.HBONDgraph.grid(column=2, padx=20, row=14, sticky="w")

        self.STATESgraph = ttk.Label(self.IntMain)
        self.STATESgraph.configure(
            justify="left",
            padding=3,
            text='HDX Activity',
            width=13)
        self.STATESgraph.grid(column=2, padx=20, row=15, sticky="w")

        self.LinearCalc = ttk.Label(self.IntMain)
        self.LinearCalc.configure(
            justify="left",
            padding=3,
            text='Linear',
            width=10)
        self.LinearCalc.grid(column=0, padx=10, row=13, sticky="w")

        self.RigimolCalc = ttk.Label(self.IntMain)
        self.RigimolCalc.configure(
            justify="left",
            padding=3,
            text='Rigimol',
            width=10)
        self.RigimolCalc.grid(column=0, padx=10, row=14, sticky="w")

        self.ClimberCalc = ttk.Label(self.IntMain)
        self.ClimberCalc.configure(
            justify="left",
            padding=3,
            text='Climber',
            width=10)
        self.ClimberCalc.grid(column=0, padx=10, row=15, sticky="w")

        # Main widget
        self.mainwindow = self.IntMain

    def run(self):
        self.mainwindow.mainloop()


if __name__ == "__main__":
    app = int_two_test()
    app.run()
