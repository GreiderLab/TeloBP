import tkinter as tk
from tkinter import filedialog
from Bio import SeqIO
import sys
import os
import gzip
import numpy as np
import pandas as pd



sys.path.insert(0, '../TeloBP')
from TeloBP import *
from constants import errorReturns
import subprocess

fileMode = False
dataDir = None
# Create a variable to store the output directory path, default to current directory
outputDir = os.getcwd()

def update_button_state():
    print(f"Data directory: {dataDir}")
    if fileModeVar.get() == True:
        button1.config(state=tk.NORMAL)
        button2.config(state=tk.DISABLED)
    else:
        button1.config(state=tk.DISABLED)
        button2.config(state=tk.NORMAL)
    if dataDir is None:
        print("Data directory is not selected")
        buttonRunAnalysis.config(state=tk.DISABLED)
    else:
        print(f"Data directory selected: {dataDir}")
        buttonRunAnalysis.config(state=tk.NORMAL)

def open_file_explorer():
    filepath = filedialog.askopenfilename()
    global fileMode
    fileMode = True
    #print(f'File selected: {filepath}')
    global dataDir
    dataDir = filepath
    update_button_state()
    # Process the selected file

def open_folder_explorer():
    folderpath = filedialog.askdirectory()
    global fileMode
    fileMode = False
    # print(f'Folder selected: {folderpath}')
    global dataDir
    dataDir = folderpath
    update_button_state()
    # Process the selected folder



root = tk.Tk()
root.geometry("1000x500")  # Set the UI size to 1000 x 500

# Create a frame to hold the buttons on the left side
button_frame = tk.Frame(root)
button_frame.pack(side=tk.LEFT, pady=20)

# Define a Tkinter variable
teloNP = tk.BooleanVar(value = True)
fileModeVar = tk.BooleanVar(value = True)

# Create a Checkbutton
toggle = tk.Checkbutton(button_frame, text="TeloNP Mode", variable=teloNP)
toggle.pack()

toggle1 = tk.Checkbutton(button_frame, text="Input Single File Mode", variable=fileModeVar, command=update_button_state)
toggle1.pack()

# Create a label for the file selection
label = tk.Label(button_frame, text="Select input fastq(s)")
label.pack()


button1 = tk.Button(button_frame, text='Open File Explorer', command=open_file_explorer)
button1.pack()

button2 = tk.Button(button_frame, text='Open Folder Explorer', command=open_folder_explorer)
button2.pack()


labelDataType = tk.Label(button_frame, text=f"Data type: {'File mode' if fileMode else 'Folder mode'}")
labelDataType.pack()


def select_output_directory():
    global outputDir
    outputDir = filedialog.askdirectory()
    print(f"Output directory selected: {outputDir}")
# Create a button for output directory selectionE
buttonOutputDir = tk.Button(button_frame, text='Select Output Directory', command=select_output_directory)
buttonOutputDir.pack()

def run_anal():
    teloNPBool = teloNP.get()
    # run_analysis(dataDir, fileMode, teloNPBool, outputDir, progressLabel, output_frame)
    print(f"Data directory: {dataDir}")
    print(f"Output directory: {outputDir}")
    print(f"File mode: {fileMode}")
    print(f"TeloNP mode: {teloNPBool}")
    print(f"python teloBPCmd.py '{dataDir}' '{outputDir}' {'--fileMode' if fileMode else ''} {'--teloNPBool' if teloNPBool else ''}")
    subprocess.run(['python', 'teloBPCmd.py', dataDir, outputDir, '--fileMode' if fileMode else '', '--teloNPBool' if teloNPBool else ''])
    # print(f"python teloBPCmd.py {dataDir} {outputDir} {'--fileMode' if fileMode else ''} {'--teloNPBool' if teloNPBool else ''}")
    # subprocess.run(['python', 'teloBPCmd.py', dataDir, outputDir, '--fileMode' if fileMode else '', '--teloNPBool' if teloNPBool else ''])


buttonRunAnalysis = tk.Button(button_frame, text='Run Analysis', command=run_anal, state=tk.DISABLED)
buttonRunAnalysis.pack()

update_button_state()
# Create a frame to hold the text outputs on the right side
output_frame = tk.Frame(root)
output_frame.pack(side=tk.RIGHT)

progressLabel = tk.Label(output_frame, text=f"Progress: ")
progressLabel.pack()



'''
text2 = tk.Text(output_frame, height=10, width=50)
text2.pack()
'''


root.mainloop()

root.mainloop()
