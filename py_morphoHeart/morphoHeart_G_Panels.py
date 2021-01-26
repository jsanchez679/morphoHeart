# -*- coding: utf-8 -*-
"""
Script that renames ISH images according to genotype.
Images to be renamed need to be inside RAW_ISHXX-XX folder, which should be inside ISHXX-XX folder
Creates figures with All ISH Images and Grouped by Genotype

@author: Juliana Sanchez-Posada
Date: August 04, 2020
"""

#%% Importing python packages

import os
import pandas as pd
import matplotlib.pyplot as plt

from os.path import join, isdir
from pathlib import Path
from skimage import io
from pylab import *

rc('axes', linewidth=1, edgecolor = 'w')

#%% Verify working dir
root_path = os.getcwd()
wd = r"D:\Documents JSP\Dropbox\Dropbox_Juliana\PhD_Thesis\Data_ongoing\ISH_ongoing\A_ISH_Analysis\py_ISHAnalysis"
if root_path != wd:
    os.chdir(wd)
    root_path = os.getcwd()
    
# Directory of ISH_ongoing  
dir_ISH_ongoing = Path(root_path).parent.parent
print("ISH_ongoing:", dir_ISH_ongoing)

#%% Get information abour ISH to process
# Ask for ISH number
ISH_No = input("Enter ISH Number (ISHXX): ")
# Ask for ISH number
ISH_Ref = input("Enter ISH Number (ISH00-XX): ")
# Ask for ISH tube number
ISH_tNo = input("Enter ISH Tube Number (tXXX): ")

#%% Get ISH-Folder Directories
# Directory of ISHXX
ISHXX = "ISH"+ISH_No
dir_ISHXX = join(dir_ISH_ongoing,ISHXX)
# ISHXX-XX
ISHXX_XX = ISHXX+"-"+ISH_Ref
# Directory of ISHXX-XX-tXXX
ISHXX_t = ISHXX_XX+"_t"+ISH_tNo
dir_ISHXX_t = join(dir_ISHXX,ISHXX_t)
print("ISH_tube:", dir_ISHXX_t)

# Directory of folder with rotated and cropped images (RotCrop_ISHXX-XX)
dir_RotCrop_ISHXXt = os.path.join(dir_ISHXX_t,"RotCrop_"+ISHXX_XX)
print("RotCrop_ISHXX-XX dir: ",dir_RotCrop_ISHXXt)

# Directory of images Rotated and Cropped
if isdir(dir_RotCrop_ISHXXt) == False:
    print("Directory with images Rotated and Cropped not found!!!")

# Directory of ISHXX_Genotype.txt
file_ISH_Gen = os.path.join(dir_ISHXX,"ISH"+ISH_No+"_Genotype.txt")
print("ISHXX_Genotype.txt dir: ",file_ISH_Gen)

# Directory of INFO_ISHXX.txt
file_ISH_INFO = os.path.join(dir_ISHXX,"INFO_ISH"+ISH_No+".txt")
print("INFO_ISHXX.txt dir: ",file_ISH_INFO)

print("Dir to save images: ", dir_ISHXX_t)

#%% Get dataframe with genotype
#Select the columns to get data from
usecols_gen = ["Pos", "ISH"+ISH_No+"-"+ISH_Ref]
# Create dataframe with genotype
df_ISH_Gen = pd.read_csv(file_ISH_Gen, sep="\t", header=0, usecols=usecols_gen, encoding = "ISO-8859-1", index_col ="Pos")

#%% Get dataframe with ISH info
#Select the columns to get data from
#usecols_info = ["ISH No", "Probe used", "Stage ", "Strain"]

# Create dataframe with INFO of ISH
df_ISH_INFO = pd.read_csv(file_ISH_INFO, sep="\t", header=0, index_col ="ISH No", encoding = "ISO-8859-1")

Strain = df_ISH_INFO.loc[ISHXX_XX,"Strain"][0:23]
Stage = df_ISH_INFO.loc[ISHXX_XX,"Stage "]
Probe = df_ISH_INFO.loc[ISHXX_XX,"Probe used"]
label_plt = Probe + ", " + Strain + ", " + Stage
plt_title = ISHXX_t

#%% Rename files and count images per genotype
sep = os.sep
tif_in_dir = Path(dir_RotCrop_ISHXXt)
dirs2Im = []

ht_count = 0
dirs2Im_ht = []
mt_count= 0
dirs2Im_mt = []
wt_count = 0
dirs2Im_wt = []

for num, filename in enumerate(tif_in_dir.glob('*.TIF*')):
    filename_str = str(filename)
    #print('filename_str: ', filename_str)
    filename_split = filename_str.split(sep)
    filename_tt = filename_split[-1]
    imCode = filename_tt[9:12]
    print('>> Image code/ref: ', imCode)   
    genotype = df_ISH_Gen.loc[imCode, ISHXX_XX]
    print(">> Genotype: ", genotype)
    dir2file = sep.join(filename_split[:-1])
    dir2file = os.path.join(dir2file,filename_tt[0:12])
    if genotype == 'ht' or genotype == 'mt' or genotype == 'wt':
        new_name = dir2file+"_"+genotype+".tif"
       # print('New_name: ',new_name)
        os.rename(filename,new_name)
        if genotype == 'ht':
            ht_count += 1
            dirs2Im_ht.append(new_name)
        elif genotype == 'mt':
            mt_count += 1
            dirs2Im_mt.append(new_name)
        elif genotype == 'wt':
            wt_count += 1
            dirs2Im_wt.append(new_name)
    else: 
        new_name = dir2file+"_ck.tif"
        #print('New_name (ck): ',new_name)
        os.rename(filename,new_name)
    dirs2Im.append(new_name)
    
print("All images have been renamed!")

# Total number of images in dir
TotIm = num

print("total:", TotIm ,", ht:",ht_count, ", mt:", mt_count, ", wt:", wt_count)

#%% Function definition to create panel
def panelCreate(Im_size, n_im, n_cols, n_rows, dir2Plot, plt_title, fig_title, label_plots, saveFig):
    #Plot
    h_fig11 = n_rows*Im_size
    w_fig11 = n_cols*Im_size
        # 1     2    3    4     5    6
    y = [1.05,1.00,0.95,0.93,0.92,0.91]
    
    fig11 = plt.figure(figsize=(w_fig11, h_fig11), constrained_layout=False)
    
    # gridspec inside gridspec
    grid = fig11.add_gridspec(n_rows, n_cols, wspace=0.0, hspace=0.0)
    
    for i in range(n_im):
        #Get Image and Label
        dir2Image = Path(dir2Plot[i])
        myIm = io.imread(dir2Image)
        txt = str(dir2Image.name)[9:12]
        # Plot
        ax = fig11.add_subplot(grid[i])
        ax.imshow(myIm)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.95, 0.05, txt,
                verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes,
                color='white', fontsize=10, fontweight='bold')
        
    plt.suptitle(plt_title+" - "+label_plots, y=y[n_rows-1], style='italic')
    if saveFig :
        plt.savefig(fig_title, bbox_inches='tight', pad_inches=0, dpi = 300)
    plt.show()
    
#%% Definitions for all panels
Im_size = 1.5

# All ISH images
exact = TotIm % 12
rows = TotIm // 12
cols = 12

if exact != 0:
    rows = rows + 1

# Heterozygous
exact_ht = ht_count % 6
rows_ht = ht_count // 6
cols_ht = 6
    
if exact_ht != 0:
    rows_ht = rows_ht + 1
    
if ht_count < cols_ht:
    cols_ht = ht_count
    
# Mutants type
exact_mt = mt_count % 6
rows_mt = mt_count // 6
cols_mt = 6
    
if exact_mt != 0:
    rows_mt = rows_mt + 1
    
if mt_count < cols_mt:
    cols_mt = mt_count
    
# Wild type
exact_wt = wt_count % 6
rows_wt = wt_count // 6
cols_wt = 6
    
if exact_wt != 0:
    rows_wt = rows_wt + 1
    
if wt_count < cols_wt:
    cols_wt = wt_count
    
saveFig = True

all_title = os.path.join(dir_ISHXX_t, plt_title+"_All.tif")
ht_title = os.path.join(dir_ISHXX_t, plt_title+"_ht.tif")
mt_title = os.path.join(dir_ISHXX_t, plt_title+"_mt.tif")
wt_title = os.path.join(dir_ISHXX_t, plt_title+"_wt.tif")

panelCreate(Im_size, TotIm, cols, rows, dirs2Im, ISHXX_t+" [All]", all_title, label_plt, saveFig)
panelCreate(Im_size, ht_count, cols_ht, rows_ht, dirs2Im_ht, ISHXX_t+" [ht]", ht_title, label_plt, saveFig)
panelCreate(Im_size, mt_count, cols_mt, rows_mt, dirs2Im_mt, ISHXX_t+" [mt]", mt_title, label_plt, saveFig)
panelCreate(Im_size, wt_count, cols_wt, rows_wt, dirs2Im_wt, ISHXX_t+" [wt]", wt_title, label_plt, saveFig)

   
