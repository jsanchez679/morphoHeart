# -*- coding: utf-8 -*-
"""
morphoHeart - E. ANALYSE ALL DATA

@author: Juliana Sanchez-Posada
"""

#%% Importing python packages
import os
import pandas as pd
import glob
from itertools import count
# import seaborn as sns
# import math
# import matplotlib.pyplot as plt

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

    #%% Get main directories and measurements dataframe
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    dir_R_meas = os.path.join(dir_data2Analyse,'R_All','df_all', 'df_meas')
    dir_pl_meas = os.path.join(dir_data2Analyse,'R_All','plots_all')
    dir_R_cjPdfs = os.path.join(dir_data2Analyse,'R_All','df_all','df_cjPDFs')
    dir_R_hmf = os.path.join(dir_data2Analyse,'R_All','df_all','df_hmf')
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = 'R', out_type = 'xlsx')
    all_CSVs = glob.glob(dir_R_meas + "/*.csv")
    df_meas = pd.concat((pd.read_csv(f) for f in all_CSVs))
    df_meas = df_meas.loc[:, ~df_meas.columns.str.contains('^Unnamed')]
    df2plot = []
    
    #%% Add all other measurements to dataframe
    df_meas = fcAn.sortDFCols(fcAn.getVarRatios(fcAn.getMainStrain(fcAn.getGenotypeAll(df_meas))))
    fcAn.printDFINfo(df_meas)
    # fcBasics.saveDF('All', df_meas, 'df_meas', os.path.join(dir_data2Analyse,'R_All', 'df_meas','R_temp'))
    
    #%% Filter dataset
    df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
    fcAn.printDFINfo(df2plot)
    save, info, ext = fcAn.q_savePlot()
    
    #%% Plot in groups or per variable
    groups = ['Surface Area', 'Heart Size','Lumen Size', 'Heart Looping', 'Tissue Myocardium',
              'Tissue Endocardium','Tissue Cardiac Jelly','Angles','Volume Percentages','Tissue Volume Percentages'] 
    
    pl_groups = fcAn.plot_groups()
    vars_dict = fcAn.def_variables(plot_type = 'strip_plots')
    # x_var = 'Stage'; hue_var = 'GenotypeAll'; shape_var = 'Strain_o'
    x_var = 'GenotypeAll'; hue_var = 'Stage'; shape_var = 'Strain_o' #***
    for group in groups: 
        vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group])
        fcAn.plotInGroups (df2plot, vars2plot = vars2plot, x_var =  x_var, hue_var =  hue_var, shape_var = shape_var, 
                           title = pl_groups[group]['title'], labels2plot = labels2plot, ips = (6,6), dir2save = dir_pl_meas,
                           n_cols = pl_groups[group]['n_cols'], h_add = 5, w_add = 1, sharey = False, 
                           yticks_lab = pl_groups[group]['yticks_lab'], info =info, save = save, dpi = 300, ext = ext)
        
    #%% Tissue composition (whole, atrium and ventricle)
    print('\n=> TISSUE COMPOSITION ANALYSIS')
    fcAn.printDFINfo(df2plot)

    colours = ['lightseagreen','darkorange','darkmagenta' ]
    vars2plot_all = [['Vol_Myoc','Vol_CJ','Vol_Endo'],
                     ['Vol_Atr.Myoc','Vol_Atr.CJ','Vol_Atr.Endo'],
                     ['Vol_Vent.Myoc','Vol_Vent.CJ','Vol_Vent.Endo']]
    add_vars2plot_all = ['Vol_Int.Endo','Vol_Atr.IntEndo','Vol_Vent.IntEndo']
    title_sp_all = ['Whole Heart ','Atrial ', 'Ventricular ']
    
    save, info, ext = fcAn.q_savePlot()
    # x = 'GenotypeAll'; col = 'Stage'; per = 'per Stage and Genotype'
    x = 'Stage'; col = 'GenotypeAll'; per = 'per Genotype and Stage'
    group_vars = ['Stage', 'GenotypeAll']
    txt_title = fcAn.get_txt_title(['Strain_o'], df2plot)
    for n, vars2plot, add_var, title_sp in zip(count(), vars2plot_all, add_vars2plot_all, title_sp_all):
        for m, stack100, unit in zip(count(), [True, False], [' (%)', ' [um$^3$]']):
            # => Only tissues 
            fcAn.barPlots(df2plot = df2plot, vars2plot = vars2plot, group_vars = group_vars, 
                          x_var = x, col_var = col, colours = colours, 
                          title = title_sp +'Tissue Composition '+ per, txt_title = txt_title,
                          ylabel = title_sp +'Tissue Composition'+unit, dir2save = dir_pl_meas,
                          info = info, stack100 = stack100, sub_bar_lab = True, save = save, ext = ext)
            # => Including lumen 
            fcAn.barPlots(df2plot = df2plot, vars2plot = vars2plot + [add_var], group_vars = group_vars, 
                          x_var = x, col_var = col, colours = colours + ['tomato'], 
                          title = title_sp +'Composition '+ per, txt_title = txt_title,
                          ylabel = title_sp +'Composition'+unit, dir2save = dir_pl_meas,
                          info = info, stack100 = stack100, sub_bar_lab = True, save = save, ext = ext)
    
    #%% Tissue Expansion or Shrinkage
    print('\n=> TISSUE EXPANSION OR SHRINKAGE ANALYSIS')
    fcAn.printDFINfo(df2plot)
    
    group_vars = ['GenotypeAll', 'Stage'] # ['GenotypeAll', 'Strain_o', 'Stage']
    txt_title = fcAn.get_txt_title(['Strain_o'], df2plot)
    for  n, vars2plot, add_var, title_sp in zip(count(), vars2plot_all, add_vars2plot_all, title_sp_all):
        _, df4plot_pct_change = fcAn.get_dfPctChange(df2plot, vars2plot, group_vars)
        # => x = 'Stage', col = 'GenotypeAll'
        fcAn.pctChange_barPlots(df2plot = df4plot_pct_change, vars2plot = [var+'_Change' for var in vars2plot], 
                                group_vars = group_vars, x_var = 'Stage', col_var = 'GenotypeAll', colours = colours,
                                title = title_sp + 'Percentage Changes\nin Tissue Composition per Genotype and Stage', 
                                txt_title = txt_title, ylabel = 'Percentage changes in tissue composition \nwith respect to previous stage(%)', 
                                dir2save = dir_pl_meas, info=info, sub_bar_lab = True, save = save, ext = ext)
        
        # => x = 'GenotypeAll', col = 'Stage'
        # fcAn.pctChange_barPlots(df2plot = df4plot_pct_change, vars2plot = [var+'_Change' for var in vars2plot], 
        #                         group_vars = group_vars, x = 'GenotypeAll', col = 'Stage', colours = colours,
        #                         title = title_sp + 'Percentage Changes\nin Tissue Composition per Genotype and Stage', 
        #                         txt_title = txt_title, ylabel = 'Percentage changes in tissue composition \nwith respect to previous stage(%)', 
        #                         dir2save = dir_pl_meas, info=info, sub_bar_lab = True, save = save)
        
        # => Including lumen 
        # vars2plot = vars2plot + [add_var]
        # _, df4plot_pct_change = fcAn.get_df_pctChange(df2plot, vars2plot, group_vars)
        # fcAn.pctChange_barPlots(df2plot = df4plot_pct_change, vars2plot = [var+'_Change' for var in vars2plot], 
        #                         group_vars = group_vars, x = 'Stage', col = 'GenotypeAll',colours = colours+['tomato'], 
        #                         title = title_sp + 'Percentage Changes\nin Composition per Genotype and Stage', 
        #                         txt_title = txt_title, ylabel = 'Percentage changes in tissue composition \nwith respect to previous stage(%)', 
        #                         dir2save = dir_pl_meas, info=info, sub_bar_lab = True, save = save)

    #%% Heatmaps 
    # - turbo colormap: https://ai.googleblog.com/2019/08/turbo-improved-rainbow-colormap-for.html
    # - https://github.com/matplotlib/matplotlib/issues/7081/
    # - https://matplotlib.org/stable/tutorials/colors/colormaps.html
    
    print('\n=> THICKNESS AND BALLOONING HEATMAP ANALYSIS')
    df_dataset = fcAn.getGenotypeAll(df_dataset)
    df_dataset = fcAn.getMainStrain(df_dataset)
    fcAn.printDFINfo(df_dataset)
    df_dataset_hm, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_dataset, [])
    fcAn.printDFINfo(df_dataset_hm)

    save, info, ext = fcAn.q_savePlot()
    filters = ['Stage', 'GenotypeAll']
    groups = [('32-34','hapln1a:wt'), ('32-34','hapln1a:mt'), 
              ('48-50','hapln1a:wt'), ('48-50','hapln1a:mt'), 
              ('72-74','hapln1a:wt'), ('72-74','hapln1a:mt')]
    # groups = [('32-34','hapln1a:wt'), ('32-34','hapln1a:mt')] 
    # groups = [('48-50','hapln1a:wt'), ('48-50','hapln1a:mt')]
    # groups = [('72-74','hapln1a:wt'), ('72-74','hapln1a:mt')]
    
    normalise = True; perChamber = False; norm_type = 'opt_div'; opt_norm = [normalise, perChamber, norm_type]
    for variable in ['CjTh', 'myocIntBall']:
        for chamber in ['Atr', 'Vent']:
            fcAn.meanHM(df_dataset_hm = df_dataset_hm, filters = filters, groups = groups, chamber = chamber, 
                        variable = variable, opt_norm = opt_norm,  dir2load_df = dir_R_hmf, dir2save_hmf = dir_pl_meas, 
                        dir_data2Analyse = dir_data2Analyse, save = save)
                
    #%% Create plots for df_cjPDFs 
    print('\n=> KERNEL DENSITY ESTIMATE (KDE) PLOT ANALYSIS')
    all_CSVs = glob.glob(dir_R_cjPdfs + "/*.csv")
    df_cjPDF = pd.concat((pd.read_csv(f) for f in all_CSVs))
    df_cjPDF = df_cjPDF.loc[:, ~df_cjPDF.columns.str.contains('^Unnamed')]
    df_cjPDF = fcAn.getMainStrain(df_cjPDF)
    fcAn.printDFINfo(df_cjPDF,  df_type = 'kde')
    df_cjPDF2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_cjPDF, df2plot = [], df_type = 'kde')
    
    #% Plot results
    save, info, ext = fcAn.q_savePlot()
    classif = ['AtrVent', 'LeftRight', 'DorsVent']
    classif_lab = ['Atrium-Ventricle', 'Left-Right', 'Dorsal-Ventral']
    fcAn.plotKDEs(classif, classif_lab, df_PDF = df_cjPDF2plot, save = save, dir2save = dir_pl_meas, 
                  info = info+'_cj', ext=ext, dpi = 300)
    fcAn.plotKDEIndiv(classif, classif_lab, df_PDF = df_cjPDF2plot, save = save, dir2save = dir_pl_meas, 
                      info = info+'_cj', ext=ext, dpi = 300)
  
    #%% Cardiac jelly in a beating heart
    print('\n=> CARDIAC JELLY IN THE BEATING HEART ANALYSIS')
    df2relPlot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot = [])
    fcAn.printDFINfo(df2relPlot)
    dict_legends = fcAn.def_legends(df2relPlot)
    df2relPlot['time_point'] = df2relPlot.Folder.str[14:18]
    save, info, ext = fcAn.q_savePlot()
    
    groups = ['Surface Area', 'Heart Size','Lumen Size', 'Heart Looping', 'Tissue Myocardium',
             'Tissue Endocardium','Tissue Cardiac Jelly','Angles','Volume Percentages','Tissue Volume Percentages'] 
    
    pl_groups = fcAn.plot_groups()
    vars_dict = fcAn.def_variables(plot_type = 'strip_plots')
    x_var = 'time_point'; hue_var = 'Stage'; shape_var = 'Strain_o' #***
    for group in groups: 
        vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group])
        fcAn.relPlotInGroups(df2relPlot, vars2plot = vars2plot, x_var =  x_var, hue_var =  hue_var, shape_var = shape_var, 
                           title = pl_groups[group]['title'], labels2plot = labels2plot, ips = (6,6), dir2save = dir_pl_meas,
                           n_cols = pl_groups[group]['n_cols'], h_add = 5, w_add = 1, sharey = False, 
                           yticks_lab = pl_groups[group]['yticks_lab'], info =info, save = save, dpi = 300, ext = ext)
    
    #%% - STATISTICS!!!
    print("Select genotypes to include in analysis: ")
    genot2filter = input ("[A] = 'wt','ht','mt' / [B] = 'wt','mt': ")

    if genot2filter == 'A':
        genotypes2filter = ['wt','ht','mt']
        box_pairs = [("wt", "ht"), ("ht", "mt"), ("wt", "mt")]
    else:
        genotypes2filter = ['wt','mt']
        box_pairs = [("wt", "mt")] 
    
    df2analyse_o = df2analyse_o[df2analyse_o['Genotype'].isin(genotypes2filter)]
    order2plot = genotypes2filter
    gen2txt = ''.join(genotypes2filter)
    
    alpha = 0.05
    x ='Genotype'
    vars2loop = stAn.getVariables(variables)
    for var in vars2loop:
        print(">>>> Running normality tests for variable:", var, "/ allele: ", alleleName, '\n')
        pvalues_norm, txt_normtest, all_data_normal = fcAn.normalityTests (alpha, x, var, genotypes2filter, df2plot, 'ShapiroWilk')
    
        test2use, multcomp_txt, txt_testSelected = fcAn.selectStatisticalTest (genotypes2filter, all_data_normal, var)
        
        x_label= 'Genotype'
        y_label= titles[variables.index(var)]+"\n"
        
        if test2use == 'One-way-ANOVA':
            pval_OWA, pval_Tukey_o, pval_Tukey = stAn.stTest_OWANOVA (x, var, df2plot)
            test_res = stAn.plotStatistics (df2plot, x, var, x_label, y_label, box_pairs, figsize_St, genot2filter,
                                 colors_bgd, colors_pts, order2plot,
                                 alleleName, txt_normtest, multcomp_txt, txt_testSelected, txt_note,
                                 test2use, pval_OWA, pval_Tukey,
                                 plots_dir, saveFig_St, ext='svg')
        elif test2use == 'Kruskal': 
            pval_Kruskal, pval_Dunn_o, pval_Dunn = stAn.stTest_Kruskal (x, var, df2plot)
            test_res = stAn.plotStatistics (df2plot, x, var, x_label, y_label, box_pairs, figsize_St, genot2filter,
                                 colors_bgd, colors_pts, order2plot,
                                 alleleName, txt_normtest, multcomp_txt, txt_testSelected, txt_note,
                                 test2use, pval_Kruskal, pval_Dunn,
                                 plots_dir, saveFig_St, ext='svg')
        else:
            test_res = stAn.plotStatistics (df2plot, x, var, x_label, y_label, box_pairs, figsize_St, genot2filter,
                                 colors_bgd, colors_pts, order2plot,
                                 alleleName, txt_normtest, multcomp_txt, txt_testSelected, txt_note,
                                 test2use, [], [],
                                 plots_dir, saveFig_St, ext='svg')
        
    #%% OTHERS TO CHECK LATER!!!
    #%% Plots with embryo ref labels
    # save = True
    # titles = ['Heart and Lumen Size','Tissue Layer Volumes','Surface Areas','Heart Looping','Angles', 'Volume Percentages']
    # input_vars = ['10,32,33,11,34,35','17-19,20-22,23-25','0-8','13,15,30','27-29', '36-38']
    # fcAn.plotPerVariableLabels('morphoHeart_D', input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
    #                  h_plot = 8, w_plot = 6, save = save, dir2save = dir_pl_meas, info = info, dpi = 300)
    
    
#%%
init = True
