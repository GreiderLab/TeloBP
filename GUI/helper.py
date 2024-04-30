import os
import tkinter as tk
from tkinter import filedialog
from Bio import SeqIO
import gzip
import numpy as np
import pandas as pd
from pandarallel import pandarallel
import multiprocessing as mp


from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

import sys

sys.path.insert(0, '../TeloBP')
from TeloBP import *
from constants import errorReturns

def graph(labels, sizes, colors, outputText, output_frame):

    # Create figure and axis
    fig = Figure(figsize=(4, 4), dpi=100)
    ax = fig.add_subplot(111)

    # Plot pie chart
    ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140)

    # Create canvas and draw pie chart onto canvas
    canvas = FigureCanvasTkAgg(fig, master=output_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()

    # Create text widgets in the output frame
    outputTextLabel = tk.Label(output_frame, text=outputText)
    outputTextLabel.pack()

