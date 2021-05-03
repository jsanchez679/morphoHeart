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
# from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
# import random

import pandas as pd
import seaborn as sns
from itertools import count
from progress.bar import Bar
import math

suffix = '%(index)d/%(max)d - %(elapsed)ds'

#%% Importing morphoHeart packages
from .morphoHeart_funcBasics import alert, ask4input, getInputNumbers, loadDF

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

    print('\n- Variables:')
    for c, value in enumerate(list_vars, 1):
        print(c-1, value)
    input_var = ask4input('Select the '+ name +' you would like to process: ', str)

    if input_var == 'all':
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

#%% func - getGenotypeAll
def getGenotypeAll(df_input):
    df_input['GenotA'] = df_input['Gene_A']+':'+df_input['Genotype_A']
    df_input['GenotB'] = df_input['Gene_B']+':'+df_input['Genotype_B']
    df_input['GenotypeAll'] = df_input['GenotA']+'/'+df_input['GenotB']
    df_input['Ref'] = df_input['LS_Session']+'_'+df_input['Fish_ref']
    
    # Create new column with complete genotype
    genotypeAll = [] 
    for i, genotA, genotB in zip(count(), df_input["GenotA"], df_input["GenotB"]): 
        if genotB == "-:-": 
            genotypeAll.append(genotA) 
        else: 
            genotypeAll.append(genotA+'/'+genotB) 
    df_input["GenotypeAll"] = genotypeAll 
    
    return df_input
    
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
    if module == 'morphoHeart_D':
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
                     "Volume Myocardium [um$^3$]","Atrial Volume\nMyocardium [um$^3$]","Ventricular Volume\nMyocardium [um$^3$]",
                     "Volume Endocardium [um$^3$]","Atrial Volume\nEndocardium [um$^3$]","Ventricular Volume\nEndocardium [um$^3$]",
                     "Volume Cardiac Jelly [um$^3$]","Atrial Volume\nCardiac Jelly [um$^3$]","Ventricular Volume\nCardiac Jelly [um$^3$]",
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

    print('\- nVariables:')
    for c, value in enumerate(variables, 1):
        print(c-1, value)
    input_var = ask4input('Select the variables you would like to process: ', str)

    if input_var == 'all':
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
    
    print('\n- Creating density plots... this process takes a while, about 8-12 min/plot out of 3 plots, be patient :)' )
    bar = Bar('- Creating density plots', max=3, suffix = suffix, check_tty=False, hide_cursor=False)
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

#%% func - filterDF2Plot
def filterDF2Plot (df_input, df_type = 'meas'):
    """
    Define which df to plot

    Parameters
    ----------
    df_input : TYPE
        DESCRIPTION.

    Returns
    -------
    df_filt : TYPE
        DESCRIPTION.
    genots : TYPE
        DESCRIPTION.
    strains : TYPE
        DESCRIPTION.
    stages : TYPE
        DESCRIPTION.

    """
    
    if df_type != 'kde':
        genots_or = sorted(df_input.GenotypeAll.unique(), reverse = True)
    else: 
        genots_or = sorted(df_input.Genotype.unique(), reverse = True)
    strains_or = sorted(df_input.Strain.unique())
    stages_or = sorted(df_input.Stage.unique())
    
    print('- Information found in the imported dataframe:')
    print(' - Genotypes: ', genots_or)
    print(' - Strains: ', strains_or)      
    print(' - Stages: ', stages_or)
    
    filtDF = ask4input('Do you want to filter dataframe to create plots containing just one group of data? \n\t(e.g. only wild-types, only one strain or stage) [0]:no/[1]:yes: ', bool)
    
    df_filt = df_input
    if df_type != 'kde':
        genots = sorted(df_filt.GenotypeAll.unique(), reverse = True)
    else: 
        genots = sorted(df_filt.Genotype.unique(), reverse = True)
        
    strains = sorted(df_filt.Strain.unique())
    stages = sorted(df_filt.Stage.unique())
    
    while filtDF:

        filter_type = ask4input('Select the category by which you want to filter the data \n\t[0]: Genotype\n\t[1]: Strain\n\t[2]: Stage\t>> :', int)
        
        if filter_type == 0:
            str_print = 'Select genotypes to include '
            for i, gen in enumerate(genots):
                str_print = str_print+'\n\t['+str(i)+']: '+ gen
            str_print = str_print+'\n\t[all]: All >> : '
            gen_filt = ask4input(str_print, str)
            filt_genot = getInputNumbers(gen_filt, genots)
            sel_genots = [genots[j] for j in filt_genot]
            if df_type != 'kde':
                df_filt = df_filt[df_filt['GenotypeAll'].isin(sel_genots)]
            else: 
                df_filt = df_filt[df_filt['Genotype'].isin(sel_genots)]
                
        if filter_type == 1:
            str_print = 'Select strains to include '
            for i, strain in enumerate(strains):
                str_print = str_print+'\n\t['+str(i)+']: '+ strain
            str_print = str_print+'\n\t[all]: All >> : '
            strain_filt = ask4input(str_print, str)
            filt_strain = getInputNumbers(strain_filt, strains)
            sel_strains = [strains[j] for j in filt_strain]
            df_filt = df_filt[df_filt['Strain'].isin(sel_strains)]
            
        if filter_type == 2:
            str_print = 'Select stages to include '
            for i, stg in enumerate(stages):
                str_print = str_print+'\n\t['+str(i)+']: '+ stg
            str_print = str_print+'\n\t[all]: All >> : '
            stage_filt = ask4input(str_print, str)
            filt_stage = getInputNumbers(stage_filt, stages)
            sel_stages = [stages[j] for j in filt_stage]
            df_filt = df_filt[df_filt['Stage'].isin(sel_stages)]
        
        if df_type != 'kde':
            print('\n- Filtered dataframe: ')
            print(df_filt[['LS_Session','Fish_ref','Strain','Stage','GenotypeAll']])
        else: 
            print('\n- Sample Filtered dataframe: ')
            print(df_filt[['Strain','Stage','Genotype']].sample(10))
            
        if df_type != 'kde':
            genots = sorted(df_filt.GenotypeAll.unique(), reverse = True)
        else: 
            genots = sorted(df_filt.Genotype.unique(), reverse = True)
            
        strains = sorted(df_filt.Strain.unique())
        stages = sorted(df_filt.Stage.unique())
        
        filtDF = ask4input('Do you want to filter the dataframe using another category [0]:no/[1]:yes? : ', bool)
    
    return df_filt, genots, strains, stages

#%% func - filterR_Autom 
def filterR_Autom (df_input, filters, group, col_out = ''):
        
    df_pivot = df_input.groupby(filters)
    df_out = df_pivot.get_group(group)
    if col_out != '':
        df_out = df_out[col_out]
    
    return df_out
    
#%% func - getLegends
def getLegends(genots, strains, stages):
    """
    

    Parameters
    ----------
    genots : TYPE
        DESCRIPTION.
    strains : TYPE
        DESCRIPTION.
    stages : TYPE
        DESCRIPTION.

    Returns
    -------
    out_genots : TYPE
        DESCRIPTION.
    out_strains : TYPE
        DESCRIPTION.
    out_stages : TYPE
        DESCRIPTION.

    """
    
    all_genots = ['hapln1a:wt', 'hapln1a:ht', 'hapln1a:mt', 
                  'hapln1a:wt/spaw:wt', 'hapln1a:wt/spaw:ht', 'hapln1a:wt/spaw:mt', 
                  'hapln1a:ht/spaw:wt', 'hapln1a:ht/spaw:ht', 'hapln1a:ht/spaw:mt',
                  'hapln1a:mt/spaw:wt', 'hapln1a:mt/spaw:ht', 'hapln1a:mt/spaw:mt',
                  'hapln1a_OE:gu', 'hapln1a_OE:g', 'hapln1a_OE:u', 'hapln1a_OE:wt',
                  'vcana:wt', 'vcana:ht', 'vcana:mt',]
    
    leg_genots = ['$hapln1a^{+/+}$', '$hapln1a^{+/-}$', '$hapln1a^{-/-}$', 
                  '$hapln1a^{+/+}/spaw^{+/+}$', '$hapln1a^{+/+}/spaw^{+/-}$', '$hapln1a^{+/+}/spaw^{-/-}$', 
                  '$hapln1a^{+/-}/spaw^{+/+}$', '$hapln1a^{+/-}/spaw^{+/-}$', '$hapln1a^{+/-}/spaw^{-/-}$',
                  '$hapln1a^{-/-}/spaw^{+/+}$', '$hapln1a^{-/-}/spaw^{+/-}$', '$hapln1a^{-/-}/spaw^{-/-}$',
                  'hapln1a_OE:gu', 'hapln1a_OE:g', 'hapln1a_OE:u', 'hapln1a_OE:wt',
                  '$vcana^{+/+}$', '$vcana^{+/-}$', '$vcana^{-/-}$',]
    
    out_genots = []
    for gen in genots:
        out_genots.append(leg_genots[all_genots.index(gen)])
    
    all_strains = ['hapln1a prom187/+ (F2s) InX','hapln1a prom187/+ (F3s) InX', 
                   'hapln1a prom241/+ (F2s) InX', 'hapln1a prom241/+ (F3s) InX',
                   'spaw+/-; hapln1a prom241/+ InX, vcana prom365/+ (F2s) InX']
    
    leg_strains = ["$hapln1a^{\Delta 187} (F2s)$", "$hapln1a^{\Delta 187} (F3s)$", 
                   "$hapln1a^{\Delta 241} (F2s)$", "$hapln1a^{\Delta 241} (F3s)$",
                   "$spaw^{+/-}; hapln1a^{\Delta 241}$", "$vcana^{\Delta 365}$"]
    
    out_strains = []
    for strain in strains:
        out_strains.append(leg_strains[all_strains.index(strain)])
    
    all_stages = ['32-34', '48-50', '58-60', '72-74']
    
    leg_stages = ['32-34hpf', '48-50hpf', '58-60hpf', '72-74hpf']
    
    out_stages = []
    for stg in stages:
        try: 
            out_stages.append(leg_stages[all_stages.index(stg)])
        except :
            out_stages.append('')

    return out_genots, out_strains, out_stages
        
#%% func - plotInGroups
def plotInGroups(script, input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
                     h_plot, w_plot, save, dir2save, info, dpi = 300, sharey = False, h_add = 5, w_add = 1, ext = 'png'):
    
    styles = ['o', '^', 's', 'v', 'D', '<', 'p', '>'] # https://matplotlib.org/stable/api/markers_api.html
    variables, ylabels = def_variables(script)
    
    for i, input_var, title in zip(count(), input_vars, titles):
        vars2plot, labels2plot = getVarsANDLabels_Autom(variables, ylabels, input_var)
        sns.set_context('poster') # notebook, talk, poster, paper
        # Set up the matplotlib figure
        num_vars = len(vars2plot)
        plots_per_col = 3
        plots_per_row = math.ceil(num_vars/plots_per_col)
        
        index_no_plot = list(range(3,(plots_per_col+1)*plots_per_row,4))
        for index in index_no_plot:
            vars2plot.insert(index, '')
            labels2plot.insert(index, '')
            
        # Set up the matplotlib figure
        size_col = (plots_per_col+1)*h_plot+h_add
        size_row = plots_per_row*w_plot+w_add
        
        # Genotypes and Strains being plotted 
        genots = sorted(df2plot.GenotypeAll.unique(), reverse=True)
        strains = sorted(df2plot.Strain.unique())
        stages = sorted(df2plot.Stage.unique())
        
        if i == 0: 
            print('- Genotypes: ', genots)
            print('- Strains: ', strains)
            print('- Stages: ', stages)
        
        palettes = ['mediumturquoise', 'darkmagenta']
       
        #  Create figure  - plt.clf()
        gridkw = dict(width_ratios=[1,1,1,0.2])
        fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col+1, figsize=(size_col, size_row), sharex=False, sharey=sharey, gridspec_kw=gridkw)
        fig.subplots_adjust(hspace=0.5, wspace=0.5)
        
        sns.set_style("ticks")
        sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
                                                        "ytick.labelsize":8, "lines.linewidth": 2.5,
                                                        "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
        # Define legends for strain and genotype
        legend_elem_gen = []
        for gen, gen_lab in enumerate(gen_legend):
            legend_elem_gen.append(Line2D([0], [0], marker='h', color='w', label=gen_lab,
                                    markerfacecolor=palettes[gen], markersize=20))
        space = [Line2D([0], [0], marker='o', color='w', label='',
                                    markerfacecolor='w', markersize=20)]
        legend_elem_strain = []
        for n_str, str_lab, mark in zip(count(), strain_legend, styles):
            legend_elem_strain.append(Line2D([0], [0], marker=mark, color='w', label=str_lab,
                                    markerfacecolor='k', markersize=15))
            
        handle_new = legend_elem_gen+space+legend_elem_strain
        legend_new = gen_legend+['']+strain_legend
        
        marker_size = 12
        dodge = True
        jitter = 0.2
        for n, ax, var, ylabel in zip(count(), axes.flatten(), vars2plot, labels2plot):
            if n in index_no_plot:
                if n == 3:
                    ax.set_axis_off()
                    ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(-1.8, 1), frameon = False)
                else: 
                    ax.remove()
            else: 
                for j, strain, style in zip(count(), strains, styles):
                    df_plot = df2plot[df2plot['Strain'] == strain]
                    m = sns.stripplot(x="Stage", y=var, hue="GenotypeAll", hue_order = genots, data=df_plot, ax = ax, order=['32-34','48-50','72-74'],
                                  marker = style, palette = palettes, jitter=jitter, dodge = dodge, size = marker_size)
                    box = ax.get_position()
                    ax.set(xlabel="Stage [hpf]", ylabel=ylabel);
                    ax.set_position([box.x0, box.y0, box.width*1, box.height])
                    ax.get_legend().remove()
                    
                    # ax.legend(handles0[:],labels0[:],loc='center left', bbox_to_anchor=(1, 0.5))
                    sns.despine()
                    
                    if n == 0:
                        handles, labels = m.get_legend_handles_labels()
                        #print(handles)
                        print(labels)

        fig.suptitle(title, fontsize = 30, y=1)
        dir2savef = os.path.join(dir2save, 'meas_all', 'R_')
        if info != '':
            fig_title = dir2savef+info+"_"+title+"."+ext
        else: 
            fig_title = dir2savef+title+"."+ext
        
        if save: 
            plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
            
#%% func - plotPerVariable
def plotPerVariable(script, input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
                     h_plot, w_plot, save, dir2save, info, dpi = 300, h_add = 5, w_add = 1, ext = 'png'):
    
    styles = ['o', '^', 's', 'v', 'D', '<', 'p', '>'] # https://matplotlib.org/stable/api/markers_api.html
    variables, ylabels = def_variables(script)
    
    for i, input_var, title in zip(count(), input_vars, titles):
        vars2plot, labels2plot = getVarsANDLabels_Autom(variables, ylabels, input_var)
        
        # Set up the matplotlib figure
        stages = sorted(df2plot.Stage.unique())
        plots_per_col = len(stages)
        plots_per_row = 1
        stages.append('')
        
        # Set up the matplotlib figure
        size_col = (plots_per_col+1)*h_plot-12
        size_row = plots_per_row*w_plot+w_add
        
        # Genotypes and Strains being plotted 
        genots = sorted(df2plot.GenotypeAll.unique(), reverse=True)
        strains = sorted(df2plot.Strain.unique())

        index_no_plot = list(range(3,(plots_per_col+1)*plots_per_row,4))

        if i == 0: 
            print('- Genotypes: ', genots)
            print('- Strains: ', strains)
            print('- Stages: ', stages)
        
        palettes = ['mediumturquoise', 'darkmagenta']
        
        for nn, var, ylabel in zip(count(), vars2plot, labels2plot):
            #  Create figure  - plt.clf()
            gridkw = dict(width_ratios=[1,1,1,0.2])
            fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col+1, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
            fig.subplots_adjust(hspace=1.5, wspace=0.05)
            
            sns.set_style("ticks")
            sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
                                                            "ytick.labelsize":'small', "lines.linewidth": 2.5,
                                                            "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
            # Define legends for strain and genotype
            legend_elem_gen = []
            for gen, gen_lab in enumerate(gen_legend):
                legend_elem_gen.append(Line2D([0], [0], marker='h', color='w', label=gen_lab,
                                        markerfacecolor=palettes[gen], markersize=10))
            handle_new = legend_elem_gen
            legend_new = gen_legend
            
            marker_size = 10
            dodge = False
            jitter = 0.2
            for n, ax, stg in zip(count(), axes.flatten(), stages):
                # print(n)
                if n in index_no_plot:
                    if n == 3:
                        ax.set_axis_off()
                        ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(0, 1), frameon = False)
                    else: 
                        ax.remove()
                else: 
                    df_plot = df2plot[df2plot['Stage'] == stg]
                    m = sns.stripplot(x="Strain", y=var, hue="GenotypeAll", hue_order = genots, data=df_plot, ax = ax, order=strains,
                                  marker = styles[0], palette = palettes, jitter=jitter, dodge = dodge, size = marker_size)
                    box = ax.get_position()
                    m.set_xticklabels(strain_legend, rotation=50)
                    # m.set_axis_labels(strains)
                    ax.set(xlabel="\nStage: "+stg+'hpf')
                    if n == 0:
                        ax.set(ylabel='\n'+ylabel+'\n')
                    else: 
                        ax.set(ylabel = '')
                        ax.tick_params(axis = 'y', labelcolor='w', width=0.1)
                        ax.spines['left'].set_color('gray')
                        ax.spines['left'].set_linestyle((0,(4,4)))#"dashed")
                    ax.set_position([box.x0, box.y0, box.width*1, box.height])
                    ax.get_legend().remove()
                    
                    # ax.legend(handles0[:],labels0[:],loc='center left', bbox_to_anchor=(1, 0.5))
                    sns.despine()
                    
                    if n == 0:
                        handles, labels = m.get_legend_handles_labels()
                        # print(handles)
                        print(var, '- Genotypes: ',labels)
    
            fig.suptitle(title, fontsize = 30, y=1)
            dir2savef = os.path.join(dir2save, 'meas_Ind', 'R_')
            if info != '':
                fig_title = dir2savef+info+"_Ind_"+var+"."+ext
            else: 
                fig_title = dir2savef+"Ind_"+var+"."+ext
            
            if save: 
                plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)
                
#%% func - plotPerVariableLabels
def plotPerVariableLabels(script, input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
                     h_plot, w_plot, save, dir2save, info, dpi = 300, h_add = 5, w_add = 1):
    
    styles = ['o', '^', 's', 'v', 'D', '<', 'p', '>'] # https://matplotlib.org/stable/api/markers_api.html
    variables, ylabels = def_variables(script)
    
    for i, input_var, title in zip(count(), input_vars, titles):
        vars2plot, labels2plot = getVarsANDLabels_Autom(variables, ylabels, input_var)
        
        # Set up the matplotlib figure
        stages = sorted(df2plot.Stage.unique())
        plots_per_col = len(stages)
        plots_per_row = 1
        stages.append('')
        
        # Set up the matplotlib figure
        size_col = (plots_per_col+1)*h_plot-12
        size_row = plots_per_row*w_plot+w_add
        
        # Genotypes and Strains being plotted 
        genots = sorted(df2plot.GenotypeAll.unique(), reverse=True)
        strains = sorted(df2plot.Strain.unique())

        index_no_plot = list(range(3,(plots_per_col+1)*plots_per_row,4))

        if i == 0: 
            print('- Genotypes: ', genots)
            print('- Strains: ', strains)
            print('- Stages: ', stages)
        
        palettes = ['mediumturquoise', 'darkmagenta']
        
        for nn, var, ylabel in zip(count(), vars2plot, labels2plot):
            #  Create figure  - plt.clf()
            gridkw = dict(width_ratios=[1,1,1,0.2])
            fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col+1, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
            fig.subplots_adjust(hspace=1.5, wspace=0.05)
            
            sns.set_style("ticks")
            sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
                                                            "ytick.labelsize":'small', "lines.linewidth": 2.5,
                                                            "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
            # Define legends for strain and genotype
            legend_elem_gen = []
            for gen, gen_lab in enumerate(gen_legend):
                legend_elem_gen.append(Line2D([0], [0], marker='h', color='w', label=gen_lab,
                                        markerfacecolor=palettes[gen], markersize=10))
            handle_new = legend_elem_gen
            legend_new = gen_legend
            
            marker_size = 10
            dodge = False
            jitter = 0.2
            for n, ax, stg in zip(count(), axes.flatten(), stages):
                # print(n)
                if n in index_no_plot:
                    if n == 3:
                        ax.set_axis_off()
                        ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(0, 1), frameon = False)
                    else: 
                        ax.remove()
                else: 
                    df_plot = df2plot[df2plot['Stage'] == stg]
                    m = sns.stripplot(x="Strain", y=var, hue="GenotypeAll", hue_order = genots, data=df_plot, ax = ax, order=strains,
                                  marker = styles[0], palette = palettes, jitter=jitter, dodge = dodge, size = marker_size)
                    box = ax.get_position()
                    m.set_xticklabels(strain_legend, rotation=50)
                    ax.set(xlabel="\nStage: "+stg+'hpf')
                    ymin, ymax = ax.get_ylim()
                    up = (ymax-ymin)//20
                    # print(up)
                    df_plot = df_plot.sort_values(by=['Stage','Strain', var])
                    fish_refs = df_plot.Ref.unique()
                    ha = ['right', 'left', 'center']*50
                    for j, ref in enumerate(fish_refs):
                        strain_pos =strains.index(df_plot['Strain'].values[j])
                        gen_pos = genots.index(df_plot['GenotypeAll'].values[j])
                        # print(ref, ha[j])
                        # rand_x = random.uniform(0.95, 1.05)
                        # rand_y = random.uniform(0.95, 1.05)
                        # ax.text(x=strain_pos*rand_x, y=df_plot[var].values[j]*rand_y, s=ref, horizontalalignment=ha[j], size=10, color=palettes[gen_pos])
                        ax.text(x=strain_pos, y=df_plot[var].values[j]+up, s=ref, horizontalalignment=ha[j], size=10, color=palettes[gen_pos])
    
                    if n == 0:
                        ax.set(ylabel='\n'+ylabel+'\n')
                    else: 
                        ax.set(ylabel = '')
                        ax.tick_params(axis = 'y', labelcolor='w', width=0.1)
                        ax.spines['left'].set_color('gray')
                        ax.spines['left'].set_linestyle((0,(4,4)))#"dashed")
                    ax.set_position([box.x0, box.y0, box.width*1, box.height])
                    ax.get_legend().remove()
                    
                    # ax.legend(handles0[:],labels0[:],loc='center left', bbox_to_anchor=(1, 0.5))
                    sns.despine()
                    
                    if n == 0:
                        handles, labels = m.get_legend_handles_labels()
                        # print(handles)
                        print(var, '- Genotypes: ',labels)
    
            fig.suptitle(title, fontsize = 30, y=1)
            dir2savef = os.path.join(dir2save, 'meas_Lab', 'R_')
            if info != '':
                fig_title = dir2savef+"Lab_"+info+"_"+var+".png"
            else: 
                fig_title = dir2savef+"Lab_"+var+".png"
            
            if save: 
                plt.savefig(fig_title, dpi=dpi, bbox_inches='tight', transparent=True)

#%% func - plotKDEs
def plotKDEs(classif, classif_lab, df_PDF, save, dir2save, info, ext, dpi = 300):
    
    stages = sorted(df_PDF.Stage.unique())
    n_stages = len(stages)
    strains = sorted(df_PDF.Strain.unique())
    genots = sorted(df_PDF.Genotype.unique(), reverse = True)
    gen_legend, strain_legend, stage_legend = getLegends(genots, strains, stages)
    stages.append('')

    pal = [['tomato','gold'],['deepskyblue', 'darkblue'],['greenyellow','darkgreen']]
    err_pal = [['lightcoral','goldenrod','lightcoral','goldenrod'],['cyan', 'steelblue','cyan', 'steelblue'],['yellowgreen','limegreen','yellowgreen','limegreen']]
    
    width_ratios=[1 for n in range(n_stages)]+[0.1]
    gridkw = dict(width_ratios=width_ratios)
    for j, cl, cl_lab in zip(count(), classif, classif_lab): 
        xgrid_max = df_PDF['x_grid'].max()
        # print('xgrid_max:', xgrid_max)
        fig, axes = plt.subplots(nrows=1, ncols=n_stages+1, figsize=(13*n_stages, 10), gridspec_kw=gridkw, sharex = False, sharey = True)
        fig.subplots_adjust(wspace=0.1)
        for n, ax, stg in zip(count(), axes.flatten(), stages):
            if n < n_stages:
                df_plot = df_PDF[df_PDF['Stage'] == stg]
                m = sns.lineplot(data=df_plot, x="x_grid", y=cl, hue="Genotype", hue_order = genots, 
                                 estimator='mean', ci=95, ax = ax, palette = pal[j], err_kws = dict(alpha = 0.3))
                fill_under_lines(color = err_pal[j], ax = ax, alpha = 0.2)
                ax.set_xlim(0, xgrid_max)
                ax.set_ylim(0, 1)        
                ax.get_legend().remove()
                ax.set(xlabel='Cardiac jelly thickness [um]', ylabel='\nThickness distribution\n'+ cl_lab + '\n')
                ax.title.set_text("Stage: "+stg+'hpf\n')
                sns.despine()
                if n == 0:
                    handles, labels = m.get_legend_handles_labels()
                    print(cl_lab,'- Genotypes: ',labels)
                # if n == 2: 
                #     ax.legend(handles, gen_legend, loc='lower right', frameon = False)
            else: 
                ax.set_axis_off()
                ax.legend(handles, gen_legend, loc='upper left', frameon = False)
        
        dir2savef = os.path.join(dir2save, 'kde', 'R_')
        if save: 
            plt.savefig(dir2savef+info+"kdeAll_"+cl+"."+ext, dpi=300, bbox_inches='tight', transparent=True)

#%% func - plotKDEIndiv
def plotKDEIndiv(classif, classif_lab, df_PDF, save, dir2save, info, ext, dpi = 300):
    
    stages_o = sorted(df_PDF.Stage.unique())
    n_stages = len(stages_o)
    stages = []
    for n in range(n_stages):
        stages.append(stages_o[n])
        stages.append('')
        
    strains = sorted(df_PDF.Strain.unique())
    genots = sorted(df_PDF.Genotype.unique(), reverse = True)
    gen_legend, strain_legend, stage_legend = getLegends(genots, strains, stages)
    gridkw = dict(height_ratios=[1 for n in range(n_stages)], width_ratios=[1,0.3])
    
    for j, cl, cl_lab in zip(count(), classif, classif_lab): 
        xgrid_max = df_PDF['x_grid'].max()
        # print('xgrid_max:', xgrid_max)
        fig, axes = plt.subplots(nrows=n_stages, ncols=2, figsize=(12, 10*n_stages), gridspec_kw=gridkw, sharex = False, sharey = True)
        fig.subplots_adjust(hspace=0.5)
        for n, ax, stg in zip(count(), axes.flatten(), stages):
            # print(n)
            if n in [0,2,4,6]:
                df_plot = df_PDF[df_PDF['Stage'] == stg]
                m = sns.lineplot(data=df_plot, x="x_grid", y=cl, hue='Filename', style = "Genotype",# hue_order = genots, 
                                 estimator='mean', ci=95, ax = ax, #palette = pal[j], 
                                 err_kws = dict(alpha = 0.3))
                # fill_under_lines(color = err_pal[j], ax = ax, alpha = 0.2)
                ax.set_xlim(0, xgrid_max)
                ax.set_ylim(0, 1)        
                ax.get_legend().remove()
                ax.set(xlabel='Cardiac jelly thickness [um]', ylabel='\nThickness distribution\n'+ cl_lab + '\n')
                ax.title.set_text("Stage: "+stg+'hpf\n')
                sns.despine()
                handles, labels = m.get_legend_handles_labels()
            else: 
                # print('aja')
                ax.set_axis_off()
                ax.legend(handles, labels, loc='upper left', fontsize = 'xx-small', frameon = False)
                
        dir2savef = os.path.join(dir2save, 'kde', 'R_')
        if save: 
            plt.savefig(dir2savef+info+"kdeIndivAll_"+cl+"."+ext, dpi=300, bbox_inches='tight', transparent=True)
            
#%% func - fill_under_lines
def fill_under_lines(color, ax=None, alpha=.2, **kwargs):
        if ax is None:
            ax = plt.gca()
        # print('ja:',len(ax.lines))
        for i, line in enumerate(ax.lines):
            x, y = line.get_xydata().T
            # print(color[i])
            ax.fill_between(x, 0, y, color=color[i], alpha=alpha, **kwargs)

#%% func - getHeatmaps2Unify
def getHeatmaps2Unify(folders, chamber, thickness, dir_R_hmf, operation = 'mean'):
    
    dfs_hmf = []
    num = 0
    for n, file in enumerate(folders):
        # print(file)
        hmf_file = 'hmf_unloop'+chamber+thickness
        try: 
            dfs_hmf.append(loadDF(file[2:], hmf_file, dir_R_hmf))
            num += 1
        except: 
            print('-No hmf heatmap dataframe found for '+thickness+' within '+file+' results folder!')
            continue
        
    if operation == 'std':
        df_hmf = pd.concat(dfs_hmf).groupby(level=0).std()
    elif operation == 'mean':
        df_hmf = pd.concat(dfs_hmf).groupby(level=0).mean()
    elif operation == 'sem':
        df_hmf = pd.concat(dfs_hmf).groupby(level=0).sem()
    
    
    return df_hmf, num

#%% func - unifyHeatmap
def unifyHeatmap(df, chamber, stage, genotype, gen_info, thickness, vmin, vmax, n_val, dir2save, savePlot, cmap = 'jet'):
    
    stage = stage+'hpf'
    if thickness == 'CjTh':
        title = 'Cardiac jelly thickness [um] - ('+chamber+', '+stage+', '+genotype+' - n='+str(n_val)+')'
    elif thickness == 'myocIntBall': 
        title = 'Myocardium ballooning [um] - ('+chamber+', '+stage+', '+genotype+' - n='+str(n_val)+')'
    
    # Make figure 
    df_max = format(df.max().max(), '.2f')
    fig, ax = plt.subplots(figsize=(16, 10))
    if 'std' in gen_info:
        title = '(STD - max:'+df_max+'um) '+ title 
    elif 'sem' in gen_info:
        title = '(SEM - max:'+df_max+'um) '+ title 
    ax = sns.heatmap(df, cmap=cmap, vmin = vmin, vmax = vmax)
    
    max_val = df.index.max()
    if max_val > 1.1: 
        y_labels = [1,2]
        y_text = '[Inflow tract >> Valve]'
    else: # Ventricle
        y_labels = [0,1]
        y_text = '[Valve >> Outflow tract]'
        
    x_pos = ax.get_xticks()
    # x_lab = ax.get_xticklabels()
    x_pos_new = np.linspace(x_pos[0], x_pos[-1], 19)
    x_lab_new = np.arange(-180,200,20)
    ax.set_xticks(x_pos_new) 
    ax.set_xticklabels(x_lab_new, rotation=30)
    
    y_pos = ax.get_yticks()
    # y_lab = ax.get_yticklabels()
    y_pos_new = np.linspace(y_pos[0], y_pos[-1], 11)
    y_lab_new = np.linspace(y_labels[0],y_labels[1],11)
    y_lab_new = [format(y,'.2f') for y in y_lab_new]
    ax.set_yticks(y_pos_new) 
    ax.set_yticklabels(y_lab_new, rotation=0)
    
    plt.ylabel('Centreline position '+y_text+'\n', fontsize=10)
    plt.xlabel('Angle (\N{DEGREE SIGN}) [Dorsal >> Right >> Ventral >> Left >> Dorsal]', fontsize=10)
    plt.title(title, fontsize = 15)
    
    dir4heatmap = os.path.join(dir2save,'hmf', 'hmfAll_'+thickness+'_'+chamber+'_'+stage+'_'+gen_info+'.png')

    if savePlot: 
        plt.savefig(dir4heatmap, dpi=300, bbox_inches='tight', transparent=True)
    

#%% Palette stuff
    # Save a palette to a variable:
    # palette = sns.color_palette("husl", 8*n_gen*n_strain)
    # sns.palplot(palette)
    # for gen in range(n_gen):
    #     palettes.append(palette[gen*len(palette)//n_gen+2])
    
#%% Others -  Plots of all groups of variables x3 (stages)

    # for i, input_var, title in zip(count(), input_vars, titles):
    #     vars2plot, labels2plot = fcAn.getVarsANDLabels_Autom(variables, ylabels, input_var)
        
    #     # Set up the matplotlib figure
    #     num_vars = len(vars2plot)
    #     plots_per_col = len(stages)
    #     plots_per_row = num_vars
    #     stages = sorted(df2plot.Stage.unique())
    #     n_stages = len(stages)
        
    #     vars2plotX3 = []
    #     labels2plotX3 = []
    #     for ii, var, lab in zip(count(), vars2plot, labels2plot): 
    #         for n in range(n_stages):
    #             vars2plotX3.append(var)
    #             labels2plotX3.append(lab)
        
    #     index_first = list(range(0,(plots_per_col+1)*plots_per_row,4))
    #     stagesX3 = stages*num_vars
    #     index_no_plot = list(range(3,(plots_per_col+1)*plots_per_row,4))
    #     for index in index_no_plot:
    #         vars2plotX3.insert(index, '')
    #         labels2plotX3.insert(index, '')
    #         stagesX3.insert(index,'')
            
    #     # Set up the matplotlib figure
    #     size_col = (plots_per_col+1)*h_plot+h_add-5
    #     size_row = plots_per_row*w_plot+w_add+20
        
    #     # Number of genotypes: 
    #     genots = sorted(df2plot.GenotypeAll.unique(), reverse=True)
    #     n_gen = len(genots)
    #     strains = sorted(df2plot.Strain.unique())
    #     n_strain = len(strains)

    #     if i == 0: 
    #         print('Genotypes: ', genots)
    #         print('Strains: ', strains)
    #         print('Stages: ', stages)
        
    #     palettes = ['mediumturquoise', 'darkmagenta']
        
    #     #  Create figure  - plt.clf()
    #     gridkw = dict(width_ratios=[1,1,1,0.2])
    #     fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col+1, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
    #     fig.subplots_adjust(hspace=1.5, wspace=0)
        
    #     sns.set_style("ticks")
    #     sns.set_context('poster', font_scale = 1, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
    #                                                     "ytick.labelsize":'small', "lines.linewidth": 2.5,
    #                                                     "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
    #     # Define legends for strain and genotype
    #     legend_elem_gen = []
    #     for gen, gen_lab in enumerate(gen_legend):
    #         legend_elem_gen.append(Line2D([0], [0], marker='h', color='w', label=gen_lab,
    #                                 markerfacecolor=palettes[gen], markersize=20))
    #     handle_new = legend_elem_gen
    #     legend_new = gen_legend
        
    #     marker_size = 12
    #     dodge = False
    #     jitter = 0.2
    #     for n, ax, stg, var, ylabel in zip(count(), axes.flatten(), stagesX3, vars2plotX3, labels2plotX3):
    #         if n in index_no_plot:
    #             if n == 3:
    #                 ax.set_axis_off()
    #                 ax.legend(handle_new, legend_new, loc='upper left', bbox_to_anchor=(0, 1), frameon = False)
    #             else: 
    #                 ax.remove()
    #         else: 
    #             df_plot = df2plot[df2plot['Stage'] == stg]
    #             m = sns.stripplot(x="Strain", y=var, hue="GenotypeAll", hue_order = genots, data=df_plot, ax = ax, order=strains,
    #                           marker = style, palette = palettes, jitter=jitter, dodge = dodge, size = marker_size)
    #             box = ax.get_position()
                
    #             m.set_xticklabels(strain_legend, rotation=45)
    #             # m.set_axis_labels(strains)
    #             if n in index_first:
    #                 ax.set(ylabel=ylabel)
    #             else: 
    #                 ax.set(ylabel = '')
    #                 ax.tick_params(axis = 'y', labelcolor='w', width=0.1)
    #                 ax.spines['left'].set_color('gray')
    #                 ax.spines['left'].set_linestyle("dashed")
    #             ax.set(xlabel="\nStage: "+stg+'hpf')
    #             ax.set_position([box.x0, box.y0, box.width*1, box.height])
    #             ax.get_legend().remove()
                
    #             # ax.legend(handles0[:],labels0[:],loc='center left', bbox_to_anchor=(1, 0.5))
    #             sns.despine()
                
    #             if n == 0:
    #                 handles, labels = m.get_legend_handles_labels()
    #                 print(handles)
    #                 print(labels)

    #     fig.suptitle(title, fontsize = 30, y=0.9)
        
        # if save: 
        #     plt.savefig(dir_pl_meas+title+"StgStrainSplit.png", dpi=300, bbox_inches='tight', transparent=True)
        
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
