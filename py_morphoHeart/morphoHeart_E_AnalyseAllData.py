# -*- coding: utf-8 -*-
"""
morphoHeart - E. ANALYSE ALL DATA

@author: Juliana Sanchez-Posada
"""

#%% Importing python packages
import os
# import numpy as np
# from skimage import measure
import pandas as pd
import glob
from itertools import count
import seaborn as sns

# Verify working dir
def setWorkingDir (root_path, init = False):
    if not init:
        wd = os.path.dirname(os.path.abspath(__file__))
        if root_path != wd:
            os.chdir(wd)
            root_path = os.getcwd()
    # init = True
    print("Current working directory: {0}".format(os.getcwd()))

    return root_path, init

root_path, init = setWorkingDir(os.getcwd())

#%% Start D_AnalyseData
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcAnalysis as fcAn
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics

    #%% Get main directories 
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    dir_R_meas = os.path.join(dir_data2Analyse,'R_All', 'df_meas')
    dir_pl_meas = os.path.join(dir_data2Analyse,'R_All', 'pl_meas')
    dir_R_cjPdfs = os.path.join(dir_data2Analyse,'R_All', 'df_cjPDFs')
    dir_R_hmf = os.path.join(dir_data2Analyse,'R_All', 'df_hmf')
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = 'R', out_type = 'xlsx')

    #%% Create plots for df_measurements 
    # Get directories
    all_CSVs = glob.glob(dir_R_meas + "/*.csv")
    df_meas = pd.concat((pd.read_csv(f) for f in all_CSVs))
    df_meas['Looping_Ratio_Myoc'] = df_meas['Length_CL_Int.Myoc(Cut)']/ df_meas['linLine_Int.Myoc(Cut)']
    df_meas['Looping_Ratio_Endo'] = df_meas['Length_CL_Ext.Endo(Cut)']/ df_meas['linLine_Ext.Endo(Cut)']
    df_meas = fcAn.getGenotypeAll(df_meas)
    
    # Create new column with complete genotype
    genotypeAll = [] 
    for i, genotA, genotB in zip(count(), df_meas["GenotA"], df_meas["GenotB"]): 
        if genotB == "-:-": 
            genotypeAll.append(genotA) 
        else: 
            genotypeAll.append(genotA+'/'+genotB) 
    df_meas["GenotypeAll"] = genotypeAll 
    
    #% Filter if needed
    df2plot, genots, strains, stages = fcAn.filterDF2Plot(df_meas)
    gen_legend, strain_legend , stage_legend = fcAn.getLegends(genots, strains, stages)
    variables, ylabels = fcAn.def_variables('morphoHeart_D')
    
    #% Plot results
    save = fcBasics.ask4input('Do you want to save the created figures? [0]:no/[1]:yes: ', bool)
    if save: 
        info = fcBasics.ask4input('Add info for plot name: ', str)
        svgpng = ['png','svg']
        svgORpng = fcBasics.ask4input('Do you want to save the created figures as [0]:png, [1]:svg, [2]:both? : ', int)
        if svgORpng == 2:
            svgORpng = svgpng
        else:
            svgORpng = [svgpng[svgORpng]]
    else: 
        info = ''
        svgORpng = ''
    
    titles = ['Surface Areas','Heart Size','Lumen Size', 'Heart Looping', 'Tissue Layer Volumes - Myocardium',
              'Tissue Layer Volumes - Endocardium','Tissue Layer Volumes - Cardiac Jelly','Angles']
    input_vars = ['0-8','10,32,33','11,34,35','13,15,30','17-19','20-22','23-25','27-29']
    fcAn.plotInGroups('morphoHeart_D', input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
                     h_plot = 8, w_plot = 6, save = save, dir2save = dir_pl_meas, info = info, dpi = 300, ext=svgORpng)#, sharey = True)
    
#%%
    titles = ['Heart and Lumen Size','Tissue Layer Volumes']
    input_vars = ['10,32,33,11,34,35','17-19,20-22,23-25']
    fcAn.plotPerVariable('morphoHeart_D', input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
                     h_plot = 8, w_plot = 6, save = save, dir2save = dir_pl_meas, info = info, dpi = 300, ext=svgORpng)
    
    #%% Plots with embryo ref labels
    titles = ['Heart and Lumen Size','Tissue Layer Volumes','Surface Areas','Heart Looping','Angles']
    input_vars = ['10,32,33,11,34,35','17-19,20-22,23-25','0-8','13,15,30','27-29']
    fcAn.plotPerVariableLabels('morphoHeart_D', input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
                     h_plot = 8, w_plot = 6, save = save, dir2save = dir_pl_meas, info = info, dpi = 300)
    
    #%% MODIFIED WHEN CREATING CJ PLOTS - CHECK!
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set_context('notebook')
    genots = sorted(df2plot.GenotypeAll.unique(), reverse=True)
    strains = sorted(df2plot.Strain.unique())
    stages = sorted(df2plot.Stage.unique())
    info = 'all'
    plots_per_col = len(stages)
    plots_per_row = 1
    stages.append('')
    h_plot = 8; w_plot = 6
    w_add = 1
        
    # Set up the matplotlib figure
    size_col = (plots_per_col+1)*h_plot-12
    size_row = plots_per_row*w_plot+w_add
    
    # fcAn.getVarsANDLabels_UserInput (variables, ylabels)
    input_vars = ['10,32,33,30,13,23-25']
    for j, input_var in enumerate(input_vars):
        vars2plot, labels2plot = fcAn.getVarsANDLabels_Autom(variables, ylabels, input_var)
        for nn, var, ylabel in zip(count(), vars2plot, labels2plot):
            #  Create figure  - plt.clf()
            gridkw = dict(width_ratios=[1,1,1,0.2])
            fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col+1, figsize=(size_col, size_row), sharex=False, sharey=True, gridspec_kw=gridkw)
            fig.subplots_adjust(hspace=1.5, wspace=0.05)
            # m = sns.pointplot(x='Stage', y=var, hue="GenotypeAll", data=df2plot, dodge=True, join=False, ci=95, palette="Set2", 
            #                   order=stages, ax=axes)
            m2 = sns.stripplot(x='Stage', y=var, hue="GenotypeAll", data=df2plot, dodge=True, palette="Set2", ax=axes,
                               order=stages)
            
            dir2savef = os.path.join(dir_pl_meas, 'meas_svg', 'R_'+info)
            fig_title = dir2savef+var+".svg"
            # plt.savefig(fig_title, dpi=300, bbox_inches='tight', transparent=True)
            
    #%%
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set_context('notebook')
    genots = sorted(df2plot.GenotypeAll.unique(), reverse=True)
    strains = sorted(df2plot.Strain.unique())
    stages = sorted(df2plot.Stage.unique())
    info = 'all'
    
    # fcAn.getVarsANDLabels_UserInput (variables, ylabels)
    input_vars = ['10,32,33,30,13,23-25']
    for j, input_var in enumerate(input_vars):
        vars2plot, labels2plot = fcAn.getVarsANDLabels_Autom(variables, ylabels, input_var)
        for i, var in enumerate(vars2plot):
            print(var)
            fig, axes = plt.subplots(1,1)
            # m = sns.pointplot(x='Stage', y=var, hue="GenotypeAll", data=df2plot, dodge=True, join=False, ci=95, palette="Set2", 
            #                   order=stages, ax=axes)
            m2 = sns.stripplot(x='Stage', y=var, hue="GenotypeAll", data=df2plot, dodge=True, palette="Set2", ax=axes,
                               order=stages)
            
            dir2savef = os.path.join(dir_pl_meas, 'meas_svg', 'R_'+info)
            fig_title = dir2savef+var+".svg"
            # plt.savefig(fig_title, dpi=300, bbox_inches='tight', transparent=True)
    
    #%% Create plots for df_cjPDFs 
    # Get directories
    all_CSVs = glob.glob(dir_R_cjPdfs + "/*.csv")
    df_cjPDF = pd.concat((pd.read_csv(f) for f in all_CSVs))
    df_cjPDF2plot, genots, strains, stages = fcAn.filterDF2Plot(df_cjPDF, 'kde')
    
    #% Plot results
    save = fcBasics.ask4input('Do you want to save the created figures? [0]:no/[1]:yes: ', bool)
    if save: 
        info = fcBasics.ask4input('Add info for plot name: ', str)
        svgORpng = fcBasics.ask4input('Do you want to save the created figures as png or svg? : ', str)
    else: 
        info = ''
        svgORpng = ''
        
    classif = ['AtrVent', 'LeftRight', 'DorsVent']
    classif_lab = ['Atrium-Ventricle', 'Left-Right', 'Dorsal-Ventral']
    fcAn.plotKDEs(classif, classif_lab, df_PDF = df_cjPDF2plot, save = save, dir2save = dir_pl_meas, 
                  info = info+'_cj', ext=svgORpng, dpi = 300)
    fcAn.plotKDEIndiv(classif, classif_lab, df_PDF = df_cjPDF2plot, save = save, dir2save = dir_pl_meas, 
                      info = info+'_cj', ext=svgORpng, dpi = 300)
    
    #%% Just for wt
    #% Filter if needed
    df2plot, genots, strains, stages = fcAn.filterDF2Plot(df_meas)
    gen_legend, strain_legend , stage_legend = fcAn.getLegends(genots, strains, stages)
    variables, ylabels = fcAn.def_variables('morphoHeart_D')
    
    #% Plot results
    save = fcBasics.ask4input('Do you want to save the created figures? [0]:no/[1]:yes: ', bool)
    info = fcBasics.ask4input('Add info for plot name: ', str)
    
    titles = ['Surface Areas','Heart Size','Lumen Size', 'Heart Looping', 'Tissue Layer Volumes - Myocardium',
              'Tissue Layer Volumes - Endocardium','Tissue Layer Volumes - Cardiac Jelly','Angles']
    input_vars = ['0-8','10,32,33','11,34,35','13,15,30','17-19','20-22','23-25','27-29']
    fcAn.plotInGroups('morphoHeart_D', input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
                     h_plot = 8, w_plot = 6, save = save, dir2save = dir_pl_meas, info ='_'+info, dpi = 300)#, sharey = True)
    
    
    #%% Heatmaps
    df_dataset = fcAn.getGenotypeAll(df_dataset)
    # dir_R_hmf = os.path.join(dir_data2Analyse,'R_All', 'df_hmN')
    
    filters = ['Stage', 'GenotypeAll']
    # groups = [('32-34','hapln1a:wt'), ('32-34','hapln1a:mt'), 
    #           ('48-50','hapln1a:wt'), ('48-50','hapln1a:mt'), 
    #           ('72-74','hapln1a:wt'), ('72-74','hapln1a:mt')]
    groups = [('32-34','hapln1a:wt'), ('32-34','hapln1a:mt')] 
    # groups = [('48-50','hapln1a:wt'), ('48-50','hapln1a:mt')]
    # groups = [('72-74','hapln1a:wt'), ('72-74','hapln1a:mt')]
    
    for group in groups:
        # print(group)
        folders = fcAn.filterR_Autom (df_dataset, filters, group, col_out = 'Folder')
        for chamber in ['Atr', 'Vent']:
            print('Unifying group: ', group, '- chamber: ', chamber)
            thickness = 'CjTh'
            df_hmf, num = fcAn.getHeatmaps2Unify(folders, chamber, thickness, dir_R_hmf)
            # df_hmf_std, num = fcAn.getHeatmaps2Unify(folders, chamber, thickness, dir_R_hmf, operation = 'std')
            # df_hmf_sem, num = fcAn.getHeatmaps2Unify(folders, chamber, thickness, dir_R_hmf, operation = 'sem')
            gen_info = list(group[1])
            if '/' in gen_info:
                ind_sep = gen_info.index('/') 
                gen_info = list(gen_info); gen_info[ind_sep] = '_'; 
            ind_gen = list(filter(lambda x: gen_info[x] == ':', range(len(gen_info)))) 
            for ind in ind_gen:
                gen_info[ind+1] = gen_info[ind+1].upper()
            if len(ind_gen) == 2:
                ind_gen[1] -= 1
            [gen_info.pop(ind) for ind in ind_gen]
            gen_info  = ''.join(gen_info)
            vmax = 0; vmin = 25
            fcAn.unifyHeatmap(df_hmf, chamber, genotype=group[1], gen_info = gen_info+'_241', stage=group[0], thickness= thickness, 
                      vmin=vmin, vmax=vmax, n_val = num, dir2save = dir_pl_meas, savePlot = True, cmap = 'jet')
            # fcAn.unifyHeatmap(df_hmf_std, chamber, genotype=group[1], gen_info = gen_info+'_241stdN', stage=group[0], thickness= thickness, 
            #           vmin=vmin, vmax=vmax, n_val = num, dir2save = dir_pl_meas, savePlot = True, cmap = 'jet')
            # fcAn.unifyHeatmap(df_hmf_sem, chamber, genotype=group[1], gen_info = gen_info+'_241sem', stage=group[0], thickness= thickness, 
            #           vmin=vmin, vmax=vmax, n_val = num, dir2save = dir_pl_meas, savePlot = True, cmap = 'jet')
            
            
    for group in groups:
        folders = fcAn.filterR_Autom (df_dataset, filters, group, col_out = 'Folder')
        for chamber in ['Atr', 'Vent']:
            print('Unifying group: ', group, '- chamber: ', chamber)
            thickness = 'myocIntBall'
            df_hmf, num = fcAn.getHeatmaps2Unify(folders, chamber, thickness, dir_R_hmf)
            # df_hmf_std, num = fcAn.getHeatmaps2Unify(folders, chamber, thickness, dir_R_hmf, operation = 'std')
            # df_hmf_sem, num = fcAn.getHeatmaps2Unify(folders, chamber, thickness, dir_R_hmf, operation = 'sem')
            gen_info = list(group[1])
            if '/' in gen_info:
                ind_sep = gen_info.index('/') 
                gen_info = list(gen_info); gen_info[ind_sep] = '_'; 
            ind_gen = list(filter(lambda x: gen_info[x] == ':', range(len(gen_info)))) 
            for ind in ind_gen:
                gen_info[ind+1] = gen_info[ind+1].upper()
            if len(ind_gen) == 2:
                ind_gen[1] -= 1
            [gen_info.pop(ind) for ind in ind_gen]
            gen_info  = ''.join(gen_info)
            vmax = 0; vmin = 100
            fcAn.unifyHeatmap(df_hmf, chamber, genotype=group[1], gen_info = gen_info+'_241', stage=group[0], thickness= thickness, 
                      vmin=vmin, vmax=vmax, n_val = num, dir2save = dir_pl_meas, savePlot = True, cmap = 'jet')
            # fcAn.unifyHeatmap(df_hmf_std, chamber, genotype=group[1], gen_info = gen_info+'_241std', stage=group[0], thickness= thickness, 
            #           vmin=vmin, vmax=vmax, n_val = num, dir2save = dir_pl_meas, savePlot = True, cmap = 'jet')
            # fcAn.unifyHeatmap(df_hmf_sem, chamber, genotype=group[1], gen_info = gen_info+'_241sem', stage=group[0], thickness= thickness, 
            #           vmin=vmin, vmax=vmax, n_val = num, dir2save = dir_pl_meas, savePlot = True, cmap = 'jet')
            
    #%% Cardiac jelly in a beating heart
    #% Create plots for df_measurements 
    # Get directories
    all_CSVs = glob.glob(dir_R_meas + "/*.csv")
    df_meas = pd.concat((pd.read_csv(f) for f in all_CSVs))
    df_meas['Looping_Ratio_Myoc'] = df_meas['Length_CL_Int.Myoc(Cut)']/ df_meas['linLine_Int.Myoc(Cut)']
    df_meas['Looping_Ratio_Endo'] = df_meas['Length_CL_Ext.Endo(Cut)']/ df_meas['linLine_Ext.Endo(Cut)']
    df_meas = fcAn.getGenotypeAll(df_meas)
    
    # Create new column with complete genotype
    genotypeAll = [] 
    for i, genotA, genotB in zip(count(), df_meas["GenotA"], df_meas["GenotB"]): 
        if genotB == "-:-": 
            genotypeAll.append(genotA) 
        else: 
            genotypeAll.append(genotA+'/'+genotB) 
    df_meas["GenotypeAll"] = genotypeAll 
    
    #% Filter if needed
    df2plot, genots, strains, stages = fcAn.filterDF2Plot(df_meas)
    gen_legend, strain_legend , stage_legend = fcAn.getLegends(genots, strains, stages)
    variables, ylabels = fcAn.def_variables('morphoHeart_D')
    
    df2plot['time_point'] = df2plot.Folder.str[14:18]
    titles = ['Heart Size','Lumen Size', 'Heart Looping', 'Tissue Layer Volumes - Myocardium',
              'Tissue Layer Volumes - Endocardium','Tissue Layer Volumes - Cardiac Jelly','Angles']
    input_vars = ['10,32,33','11,34,35','13,15,30','17-19','20-22','23-25','27-29']
    
    info = 'CJ_DATA'
    fcAn.plotRelPlotInGroups('morphoHeart_D', input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
                     h_plot = 8, w_plot = 6, save = save, dir2save = dir_pl_meas, info = info, dpi = 300, sharey = False, ext=svgORpng)#, sharey = True)
    
   
#%%
init = True
    
                    

