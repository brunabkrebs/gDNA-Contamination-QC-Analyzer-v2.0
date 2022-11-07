'''
File name: gDNAcont-analyzer.py
Author: Bruna Krebs Kutche
Date created: 10/30/2022
Python version: 3.9
'''

'''
###############
##  READ ME  ##
###############

Instructions:

> Pre-process data on Biorad's CFX Maestro software: Subtract baseline, Determine Cq via Single Threshold
> Export all files as .csv to a *** local directory ***
> Generate Layout file according to each experiment using only the allowed labels
> Result files will be saved in the same directory as the Quantification Amplification Results file, which MUST be a local directory


'''

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import pandas as pd
import tkinter
from tkinter import filedialog
import math
import numpy as np
import os
import shutil
from tkinter import *
from tkinter import messagebox
from PIL import Image, ImageTk
import matplotlib.ticker as tck


### Get input files ###

gui = tkinter.Tk()
gui.wm_iconbitmap('BKK_Lapps_ico2.ico')
gui.title( 'E.coli Genomic DNA Contamination Assay Analyzer' )
gui.geometry( "750x300" )
gui.grid_columnconfigure( 0, minsize=30 )
gui.grid_columnconfigure( 2, minsize=10 )

for rows in range( 0, 11 ):
    gui.grid_rowconfigure( rows, minsize=20 )
title = Label( gui, text='Instructions:', font="Arial 14", width='15', height='1', justify='left',
                  anchor='w' ).grid( row=1, column=1, sticky=tkinter.W )
body = Label( gui, text= '> Pre-process data on Biorad\'s CFX Maestro software:\n        * Subtract baseline\n        * Determine Cq via Single Threshold'
                         '\n\n> Export all files as .csv to a local directory.'
                         '\n\n> Generate Layout file according to each experiment using only the allowed labels.'
                         '\n\n> Result files will be saved in the same directory as the .csv files.'
                         '\n\n\n                    Close this window to start analyzing',
              font="Arial 11", justify='left', anchor='w' ).grid( row=2, column=1, sticky=tkinter.W )


# Add logo
image = Image.open('BKK_Lapps_logo.png')
image = image.resize((150, 200))
image2 = ImageTk.PhotoImage(image)
panel = Label(gui, image=image2)
panel.grid(row=1, column=2, rowspan=2, sticky=tkinter.E)
taglabel = Label(gui,text='by Bruna Krebs Kutche', font='Arial 9', fg='#AAB7B8',justify='right', anchor='w').grid(row=3, column=2,sticky=tkinter.W, pady=5)

gui.wait_window()

# Read Biorad's output .csv file
filename1 = filedialog.askopenfilename(initialdir='C:/Users/BrunaK/Desktop', filetypes=[('CSV files', '*.csv')],
                                       title="Open the Quantification Amplification Results file")

file1 = pd.read_csv(filename1)
data = pd.DataFrame(file1)

# Read Biorad's Cq values file
filename3 = filedialog.askopenfilename(filetypes=[('CSV files', '*.csv')],
                                       title="Open the Quantification Cq Results file")
file3 = pd.read_csv(filename3)
cqtable = pd.DataFrame(file3)
cqtable = cqtable[['Well', 'Cq']]


# Read Biorad's Tmelt file - if it exists
tmelt = ()
filename4 = filedialog.askopenfilename(filetypes=[('CSV files', '*.csv')],
                                       title="Open the Melt Curve Derivative Results file")
file4 = pd.read_csv(filename4)
tmelttable = pd.DataFrame(file4)


# Read layout file
filename2 = filedialog.askopenfilename(filetypes=[('CSV files', '*.csv')],
                                       title='Open the Layout file')
file2 = pd.read_csv(filename2)
layout = pd.DataFrame(file2)
layout = layout.drop(layout.columns[0], axis=1)

# Determine directory where results will be saved
savedir = os.path.dirname(filename1)

### Process data ###

# Delete unused columns from CSV files
data = data.drop(data.columns[0], axis=1)
data = data.drop(data.columns[0], axis=1)
cqtable = cqtable.sort_values(by=['Well'])
tmelttable = tmelttable.drop(tmelttable.columns[0], axis=1)


# Rename columns according to layout

# Create list with legend names
legend = []
i = 0
while i < 8:
    legend.append(list(layout.iloc[i]))
    i = i + 1

legend = [x for y in legend for x in y]

# Rename cells
data.columns = legend
cqtable.index = legend
tmelttable = tmelttable.set_index('Temperature')
tmelttable.columns = legend

# Keep only used columns according to layout
data = data.loc[:, data.columns.notnull()]
cqtable = cqtable.loc[cqtable.index.notnull(), :]
cqtable = cqtable.drop(cqtable.columns[0], axis=1)
cqtable = cqtable.sort_index()


# Calculate Cq Mean and StDev
cqmean = cqtable.groupby(level=0).mean()
cqmean = cqmean.rename(columns={"Cq": "Cq Mean"})

cqstdev = cqtable.groupby(level=0).std()
cqstdev = cqstdev.rename(columns={"Cq": "Cq St.Dev."})

cqsummary = cqmean
cqsummary['Cq St.Dev.'] = cqstdev['Cq St.Dev.'].values
cqsummary = cqsummary.round(decimals = 2)
cqsummary = cqsummary.reindex(['NTC','10pg gDNA','1pg gDNA','0.1pg gDNA','Reference + 1pg gDNA','Test + 1pg gDNA','Reference only','Test only'])

# Export files to .csv
data.to_csv(os.path.join(savedir,r'Raw_Data.csv'))
cqtable.to_csv(os.path.join(savedir,r'Cq_Table.csv'))
cqsummary.to_csv(os.path.join(savedir,r'Cq_Summary.csv'))

# Make table for each criteria
criteria1tb = cqsummary.loc[['NTC', '1pg gDNA'], ['Cq Mean']]
criteria2tb = cqsummary.loc[['Reference + 1pg gDNA', 'Test + 1pg gDNA', '1pg gDNA'], ['Cq Mean']]
criteria3tb = cqsummary.loc[['Reference only', 'Test only', '1pg gDNA'], ['Cq Mean']]

criteria1 = ''
if cqsummary.loc['NTC']['Cq Mean'] > cqsummary.loc['1pg gDNA']['Cq Mean']:
    criteria1 = 'Pass'
else:
    criteria1 = 'Fail'

criteria2 = ''
if cqsummary.loc['1pg gDNA']['Cq Mean'] - 1 <= cqsummary.loc['Reference + 1pg gDNA']['Cq Mean'] <= cqsummary.loc['1pg gDNA']['Cq Mean'] + 1:
    criteria2 = 'Reference Enzyme: Pass'
else:
    criteria2 = 'Reference Enzyme: Fail'

criteria2b = ''
if cqsummary.loc['1pg gDNA']['Cq Mean'] - 1 <= cqsummary.loc['Test + 1pg gDNA']['Cq Mean'] <= cqsummary.loc['1pg gDNA']['Cq Mean'] + 1:
    criteria2b = 'Test Enzyme: Pass'
else:
    criteria2b = 'Test Enzyme: Fail'

criteria3 = ''
if cqsummary.loc['Reference only']['Cq Mean'] > cqsummary.loc['1pg gDNA']['Cq Mean']:
    criteria3 = 'Reference Enzyme: Pass'
else:
    criteria3 = 'Reference Enzyme: Fail'

criteria3b = ''
if cqsummary.loc['Test only']['Cq Mean'] > cqsummary.loc['1pg gDNA']['Cq Mean']:
    criteria3b = 'Test Enzyme: Pass'
else:
    criteria3b = 'Test Enzyme: Fail'

### Plot Results ###

plt.rcParams['font.size'] = '16'
datamax = data.max().max()
ylimit = int(math.ceil(datamax / 1000)) * 1000

# Plot raw traces for standard QC experiments
plt.figure(1)
plt.plot(data['NTC'], color='black')
plt.plot(data['10pg gDNA'], color='red')
plt.plot(data['1pg gDNA'], color='blue')
plt.plot(data['0.1pg gDNA'], color='green')
plt.plot(data['Reference + 1pg gDNA'], color='#DA54A3')
plt.plot(data['Test + 1pg gDNA'], color='#00CC99')
plt.plot(data['Reference only'], color='#0071BC')
plt.plot(data['Test only'], color='orange')

plt.title('gDNA contamination assay')
plt.xlabel('Cycle')
plt.ylabel('RFU')
plt.axis([0, 40, 0, ylimit])
plt.grid(color='#C8CCD0', linewidth=0.5)
plt.tight_layout()
plt.savefig(os.path.join(savedir,r'Raw traces.png'))
plt.savefig('Raw traces.png')
plt.close()

# Legend figure
plt.figure(2)
black_square = mpatches.Patch(facecolor='black', edgecolor='black', label='NTC')
red_square = mpatches.Patch(facecolor='red', edgecolor='black', label='10pg gDNA')
blue_square = mpatches.Patch(facecolor='blue', edgecolor='black', label='1pg gDNA')
green_square = mpatches.Patch(facecolor='green', edgecolor='black', label='0.1pg gDNA')
square1 = mpatches.Patch(facecolor='#DA54A3', edgecolor='black', label='Reference + 1pg gDNA')
square2 = mpatches.Patch(facecolor='#00CC99', edgecolor='black', label='Test + 1pg gDNA')
square3 = mpatches.Patch(facecolor='#0071BC', edgecolor='black', label='Reference only')
square4 = mpatches.Patch(facecolor='orange', edgecolor='black', label='Test only')
plt.legend(handles=[black_square, red_square, blue_square, green_square, square1, square2, square3, square4], loc='center left')
plt.grid(False)
plt.axis('off')
plt.savefig(os.path.join(savedir,r'Legend.png'))
plt.savefig('Legend.png')
plt.close()


# # Plot Melt Peak
plt.figure(3)
plt.plot(tmelttable['NTC'], color='black')
plt.plot(tmelttable['10pg gDNA'], color='red')
plt.plot(tmelttable['1pg gDNA'], color='blue')
plt.plot(tmelttable['0.1pg gDNA'], color='green')
plt.plot(tmelttable['Reference + 1pg gDNA'], color='#DA54A3')
plt.plot(tmelttable['Test + 1pg gDNA'], color='#00CC99')
plt.plot(tmelttable['Reference only'], color='#0071BC')
plt.plot(tmelttable['Test only'], color='orange')

plt.title('Melt Peak')
plt.xlabel('Temperature (C)')
plt.ylabel('-d(RFU)/dT')
plt.grid(color='#C8CCD0', linewidth=0.5)
plt.minorticks_on()
plt.tight_layout()
plt.savefig(os.path.join(savedir,r'Melt peak.png'))
plt.savefig('Melt peak.png')
plt.close()


# ## Get Threshold value on BioRad software ##
# th = float(input("What is the Threshold Value determined by BioRad?"))
#
# ## Determine Cq of samples ##
#
#
# ncols = data.shape[1] #Determines how many columns in the Data matrix
#
# data[data.iloc[:, 1].gt(th)].index[0] #prints Cq of column 1
#
# #prints Cq of all columns
# cq=[]
# i=0
# while i < data.shape[1]:
#     cq.append((data[data.iloc[:, i].gt(th)].index[0]))
#     i = i+1
#
# data.loc[len(data)]=cq #Adds Cq values to bottom of each column
#
# #Table with Cq values of each sample
# cqtable = data.loc[40]
#
# cqavg = cqtable.groupby(level=0).mean()
# cqstdev = cqtable.groupby(level=0).std()
# summary = pd.DataFrame({'Cq mean':cqavg, 'Cq StDev':cqstdev})
#
# summary.to_csv('summary.csv')
#
# print(summary)


# Set up the HTML page

html = f'''
    <html>
        <head>
            <title>Report</title>
        </head>
        <body>
            <h1>E.coli gDNA Contamination Assay QC - Report</h1>
            <h2>Raw Traces</h2>
            <img src='Raw Traces.png')' width="600">
            <img src='Legend.png' width="400">
            <h2>Cq Values</h2>
            {cqsummary.to_html()}
            <br>
            <h2>Criteria to Pass QC:</h2>
            <p>1. 'NTC' should have higher average Cq values than '1pg gDNA'</p>
            {criteria1tb.to_html()}
            <p>{criteria1}</p>
            <br>
            <p>2. Average Cq values for 'Enzyme + 1pg gDNA' should be within 1 Cq of '1pg gDNA'</p>
            {criteria2tb.to_html()}
            <p>{criteria2}</p>
            <p>{criteria2b}</p>
            <br>
            <p>3. 'Enzyme only' should have higher average Cq than '1pg gDNA'</p>
            {criteria3tb.to_html()}
            <p>{criteria3}</p>
            <p>{criteria3b}</p>
            <br>
            <h2>Melt Peak</h2>
            <img src='Melt peak.png')' width="600">
            <br>
            <br>
            <br>
        </body>
    </html>
    '''

# Write the html string as an HTML file
with open('html_report.html', 'w') as f:
    f.write(html)

# Open html file on browser
os.system("start html_report.html")
shutil.copy('html_report.html', savedir)