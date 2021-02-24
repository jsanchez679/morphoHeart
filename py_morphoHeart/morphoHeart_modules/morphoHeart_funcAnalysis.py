# -*- coding: utf-8 -*-
"""
morphoHeart_funcAnalysis

Version: Nov, 2020
@author: Juliana Sanchez-Posada

"""
#%% Importing python packages
import os
import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from itertools import count
from progress.bar import Bar
suffix = '%(index)d/%(max)d - %(elapsed)ds'

#%% Importing morphoHeart packages
from .morphoHeart_funcBasics import alert

#%% func - getVariables
def getVariables (list_vars, name):
    """


    Parameters
    ----------
    list_vars : TYPE
        DESCRIPTION.
    name : TYPE
        DESCRIPTION.

    Returns
    -------
    variables2loop : TYPE
        DESCRIPTION.

    """

    print('\nVariables:')
    for c, value in enumerate(list_vars, 1):
        print(c-1, value)
    input_var = input('Select the '+ name +' you would like to process: ')

    if input_var == 'All':
        var_num = list(range(0,len(list_vars),1))

    else:
        var_num = []
        comma_split = input_var.split(',')

        for string in comma_split:
            if '-' in string:
                minus_split = string.split('-')
                #print(minus_split)
                for n in list(range(int(minus_split[0]),int(minus_split[1])+1,1)):
                    #print(n)
                    var_num.append(n)
            else:
                var_num.append(int(string))

    variables2loop = []
    for i, num in enumerate(var_num):
        variables2loop.append(list_vars[num])

    return variables2loop

#%% func - def_variables
def def_variables(module):
    """
    

    Parameters
    ----------
    module : TYPE
        DESCRIPTION.

    Returns
    -------
    variables : TYPE
        DESCRIPTION.
    ylabels : TYPE
        DESCRIPTION.

    """
    if module == 'morphoHeart_D_AnalyseAllData':
        variables = ["SurfArea_Myoc","SurfArea_Int.Myoc","SurfArea_Ext.Myoc",
                     "SurfArea_Endo","SurfArea_Int.Endo","SurfArea_Ext.Endo",
                     "SurfArea_CJ","SurfArea_Int.CJ","SurfArea_Ext.CJ",
                     "Vol_Int.Myoc","Vol_Ext.Myoc",
                     "Vol_Int.Endo","Vol_Ext.Endo",
                     "linLine_Int.Myoc(Cut)","linLine_Ext.Endo(Cut)",
                     "Length_CL_Int.Myoc(Cut)","Length_CL_Ext.Endo(Cut)",
                     "Vol_Myoc","Vol_Atr.Myoc","Vol_Vent.Myoc",
                     "Vol_Endo","Vol_Atr.Endo","Vol_Vent.Endo",
                     "Vol_CJ","Vol_Atr.CJ","Vol_Vent.CJ",
                     "ang_Heart","ang_Atr","ang_Vent","ang_BtwChambers",
                     'Looping_Ratio_Myoc','Looping_Ratio_Endo',
                     'Vol_Atr.ExtMyoc','Vol_Vent.ExtMyoc','Vol_Atr.IntEndo','Vol_Vent.IntEndo']

        ylabels  = ["Surface Area\nMyocardium [um$^2$]","Surface Area\nInt.Myocardium [um$^2$]","Surface Area\nExt.Myocardium [um$^2$]",
                     "Surface Area\nEndocardium [um$^2$]","Surface Area\nInt.Endocardium [um$^2$]","Surface Area\nExt. Endocardium [um$^2$]",
                     "Surface Area\nCardiac Jelly [um$^2$]","Surface Area\nInt.Cardiac Jelly [um$^2$]","Surface Area\nExt.Cardiac Jelly [um$^2$]",
                     "Volume Int.Myocardium [um$^3$]","Heart Volume [um$^3$]",
                     "Heart Lumen Volume [um$^3$]","Volume Ext.Endocardium [um$^3$]",
                     "Linear Heart Length (Int.Myoc) [um]","Linear Heart Length (Ext.Endo) [um]",
                     "Looped Heart Length (Int.Myoc) [um]","Linear Heart Length (Ext.Endo) [um]",
                     "Volume Myocardium [um$^3$]","Atrial Volume Myocardium [um$^3$]","Ventricular Volume Myocardium [um$^3$]",
                     "Volume Endocardium [um$^3$]","Atrial Volume Endocardium [um$^3$]","Ventricular Volume Endocardium [um$^3$]",
                     "Volume CJ [um$^3$]","Atrial Volume CJ [um$^3$]","Ventricular Volume CJ [um$^3$]",
                     "Angle Heart wrt Sample (\N{DEGREE SIGN})","Atrial Angle (\N{DEGREE SIGN})","Ventricular Angle (\N{DEGREE SIGN})","Angle between Chambers (\N{DEGREE SIGN})",
                     'Looping Ratio (Int.Myoc)','Looping Ratio (Ext.Endo)',
                     'Atrial Volume [um$^3$]','Ventricular Volume [um$^3$]','Atrial Lumen Volume [um$^3$]','Ventricular Lumen Volume [um$^3$]']
    
    return variables, ylabels

#%% func - getVarsANDLabels_UserInput
def getVarsANDLabels_UserInput (variables, labels):
    """
    

    Parameters
    ----------
    variables : TYPE
        DESCRIPTION.
    labels : TYPE
        DESCRIPTION.

    Returns
    -------
    vars2loop : TYPE
        DESCRIPTION.
    labels2loop : TYPE
        DESCRIPTION.

    """

    print('\nVariables:')
    for c, value in enumerate(variables, 1):
        print(c-1, value)
    input_var = input('Select the variables you would like to process: ')

    if input_var == 'All':
        var_num = list(range(0,len(variables),1))

    else:
        var_num = []
        comma_split = input_var.split(',')

        for string in comma_split:
            if '-' in string:
                minus_split = string.split('-')
                #print(minus_split)
                for n in list(range(int(minus_split[0]),int(minus_split[1])+1,1)):
                    #print(n)
                    var_num.append(n)
            else:
                var_num.append(int(string))

    vars2loop = []
    labels2loop = []
    for i, num in enumerate(var_num):
        vars2loop.append(variables[num])
        labels2loop.append(labels[num])

    return vars2loop, labels2loop

#%% func - getVarsANDLabels_Autom
def getVarsANDLabels_Autom (variables, labels, input_var):
    """
    

    Parameters
    ----------
    variables : TYPE
        DESCRIPTION.
    labels : TYPE
        DESCRIPTION.

    Returns
    -------
    vars2loop : TYPE
        DESCRIPTION.
    labels2loop : TYPE
        DESCRIPTION.

    """

    var_num = []
    comma_split = input_var.split(',')

    for string in comma_split:
        if '-' in string:
            minus_split = string.split('-')
            #print(minus_split)
            for n in list(range(int(minus_split[0]),int(minus_split[1])+1,1)):
                #print(n)
                var_num.append(n)
        else:
            var_num.append(int(string))

    vars2loop = []
    labels2loop = []
    for i, num in enumerate(var_num):
        vars2loop.append(variables[num])
        labels2loop.append(labels[num])

    return vars2loop, labels2loop

#%% func - kde_sklearn
def kde_sklearn(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scikit-learn"""
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x[:, np.newaxis])
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return np.exp(log_pdf)

#%% func - kdeThPlots
def kdeThPlots(filename, df_file, file_num, variable, thData, dir2save, save = True):
    """
    Function that creates Kernel Density Estimation Plots for each region of the heart (e.g, atrium, ventricle, 
    left, right, dorsal and ventral and then

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    df_file : TYPE
        DESCRIPTION.
    file_num : TYPE
        DESCRIPTION.
    variable : TYPE
        DESCRIPTION.
    thData : TYPE
        DESCRIPTION.
    dir2save : TYPE
        DESCRIPTION.
    save : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    df_pdfs : TYPE
        DESCRIPTION.
        
    Some links for reference:
        #https://stackabuse.com/kernel-density-estimation-in-python-using-scikit-learn/
        #https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/
        #https://jakevdp.github.io/PythonDataScienceHandbook/05.13-kernel-density-estimation.html

    """
    if variable == 'cj_thickness':
        title = filename + ' - Cardiac Jelly Thickness [um]'
        xlabel = 'Cardiac Jelly Thickness [um]'
        step = 0.05
    elif variable == 'myoc_intBall' :
        title = filename + ' - Myoc.Int Ballooning [um]'
        xlabel = 'Myoc.Int Ballooning [um]'
        step = 0.25
    
    regions_div = [['AtrVent','atrium', 'ventricle'], ['LeftRight','left','right'],['DorsVent','dorsal','ventral']]
    color = [['','tomato','gold'],['','deepskyblue', 'darkblue'],['','greenyellow', 'darkgreen']]
    num_linspace = int((max(round(thData[variable]))/step)+2)
    x_grid = np.linspace(0, max(round(thData[variable]))+step, num_linspace)

    genotype = df_file.loc[file_num,'Gene_A']+':'+df_file.loc[file_num,'Genotype_A']
    if df_file.loc[file_num,'Gene_B'] != '-':
        genotype = str(genotype+'/'+df_file.loc[file_num,'Gene_B']+':'+df_file.loc[file_num,'Genotype_B'])
        
    df_pdfs = pd.DataFrame(columns = ['Filename','Strain','Stage','Genotype', 'x_grid','AtrVent','atrium', 'ventricle', 'LeftRight','left','right','DorsVent','dorsal','ventral'])
    df_pdfs['Filename'] = [filename for j in range(len(x_grid))]
    df_pdfs['Strain'] = [df_file.loc[file_num,'Strain'] for j in range(len(x_grid))]
    df_pdfs['Stage'] = [df_file.loc[file_num,'Stage'] for j in range(len(x_grid))]
    df_pdfs['Genotype'] = [genotype for j in range(len(x_grid))]
    df_pdfs['x_grid'] = x_grid
    plot_dir = os.path.join(dir2save, filename+"_")
    
    bar = Bar('Creating density plots', max=3, suffix = suffix, check_tty=False, hide_cursor=False)
    for i, reg, col in zip(count(), regions_div, color):
        df_one = thData[thData[reg[0]] == reg[1]]
        th_one = df_one[variable]
        bw_one = len(th_one)**(-1./(1+4))*np.std(th_one)
        pdf_one = kde_sklearn(th_one, x_grid, bandwidth=bw_one)
        
        df_two = thData[thData[reg[0]] == reg[2]]
        th_two = df_two[variable]
        bw_two = len(th_two)**(-1./(1+4))*np.std(th_two)
        pdf_two = kde_sklearn(th_two, x_grid, bandwidth=bw_two)
        
        print('\n - bandwidths:', format(bw_one,'.2f'), '-',format(bw_two, '.2f'))
        pdf_comb = np.add(pdf_one, pdf_two)
        
        pct_one = pdf_one/pdf_comb
        # pct_two = pdf_two/pdf_comb
    
        ones = np.ones((len(x_grid),))
        
        fig, ax = plt.subplots(figsize=(8,5))
        ax.fill_between(x_grid, pct_one, alpha=0.5, label= reg[1], color = col[1])
        ax.fill_between(x_grid, pct_one, ones, alpha=0.5, label= reg[2], color= col[2])
        # ax.plot(x_grid, pct_one, linewidth=1, alpha=0.5, label= reg[1])
        # ax.plot(x_grid, pct_two, linewidth=1, alpha=0.5, label= reg[2])
        ax.set_xlim(0, x_grid[-1])
        ax.set_ylim(0, 1)        
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.set_title(title, fontsize = 10)
        ax.set_xlabel(xlabel, fontsize=10)
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        if save: 
            plt.savefig(plot_dir+"kde"+reg[0]+".png", dpi=300, bbox_inches='tight', transparent=True)
        plt.show()

        # fig, ax = plt.subplots()
        # ax.plot(x_grid, pdf_one, linewidth=1, alpha=0.5, label= 'pdf_'+reg[1])
        # ax.plot(x_grid, pdf_two, linewidth=1, alpha=0.5, label= 'pdf_'+reg[2])
        # ax.legend(loc='upper right')
        
        df_pdfs[reg[0]] = pct_one
        df_pdfs[reg[1]] = pdf_one
        df_pdfs[reg[2]] = pdf_two
        
        bar.next()
        
    bar.finish()
    alert('wohoo',1)
    
    return df_pdfs

#%% func - compareKDE
# import seaborn as sns
# sns.set_theme(style="darkgrid")

# # Load an example dataset with long-form data
# fmri = sns.load_dataset("fmri")

# # Plot the responses for different events and regions
# sns.lineplot(x="timepoint", y="signal",
#              hue="region", style="event",
#              data=fmri)

# flights = sns.load_dataset("flights")
# flights.head()

# flights_wide = flights.pivot("year", "month", "passengers")
# flights_wide.head()
# sns.lineplot(data=flights_wide)

# sns.lineplot(data=flights, x="year", y="passengers")
# sns.lineplot(data=flights, x="year", y="passengers", hue="month")

#%% - ALERT WHEN IMPORTED
print ("IMPORTED: morphoHeart_funcAnalysis")
alert('jump',1)
