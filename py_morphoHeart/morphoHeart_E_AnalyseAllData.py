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

# Verify working dir
def setWorkingDir (root_path, init = False):
    if not init:
        wd = os.path.dirname(os.path.abspath(__file__))
        if root_path != wd:
            os.chdir(wd)
            root_path = os.getcwd()
    print("- Current working directory: {0}".format(os.getcwd()))
    if not init: 
        print('\nIMPORTANT NOTE:\n- Remember to start running from cell %Start D_AnalyseData.\n- NEVER run as an individual cell the cell called %Importing python packages')
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
    dir_R_hm_reg = os.path.join(dir_data2Analyse,'R_All','df_all','df_hm_reg')
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = 'R', out_type = 'xlsx')
    all_CSVs = glob.glob(dir_R_meas + "/*.csv")
    #% df_meas - df containing all .csv files in the R_all folder within im_morphoHeart
    df_meas = pd.concat((pd.read_csv(f) for f in all_CSVs))
    df_meas = df_meas.loc[:, ~df_meas.columns.str.contains('^Unnamed')]
    df2plot = []; df_cjPDF2plot = []
    # Add all other measurements to dataframe and apply format
    df_meas = fcAn.spAnalysis(fcAn.getMainStrain(fcAn.getGenotypeAll(df_meas)))
    df_meas = fcAn.cleanDfCols(df_meas)
    df_meas = fcAn.sortDFCols(fcAn.getVarRatios(df_meas))

    # fcAn.printDFINfo(df_meas)
    fcBasics.saveDF('All', df_meas, 'df_meas', os.path.join(dir_data2Analyse,'R_All', 'df_all','df_meas','R_temp'))
    # fcBasics.saveDF('spaw-mts_hapln1a-wts', df2plot, 'df2plot', os.path.join(dir_data2Analyse,'R_All', 'df_all','df_meas','R_temp'))
    
    #%% Set context 
    thesis_contxt_input = True#fcAn.ask4input('Define context: \n\t- [0]: presentation or poster \n\t- [1]: document:', bool)
    thesis_use = fcAn.Thesis_contxt(thesis_contxt_input)
    ips = [(0.8,6), (0.2,1.5)]
    plotInGroups = False

    pl_indiv = fcAn.plot_indiv()
    for n, var in enumerate(pl_indiv.keys()):
        print(str(n)+'. -'+var)
        
    #%% Plot Individual Plots for Group of Filtered Dataset - working for document style!
    newGroup = True
    # groupsA = ['h1a241vs187_mxWts_50hpf','h1a241vs187_NotMxWts_50hpf',#]#,
    #              'h1a241_mtVwt_time', 
    #              'h1aOE_mxCtrls_50hpf', 'h1aOE_NotMxCtrls_50hpf', 
    #              'h1a_oeVmt_mxCtrls_50hpf','h1a_oeVmt_NotMxCtrls_50hpf', 
    #              'vcana_mtVwt_50hpf',
    #              'spaw_mtVsibs_time_shape', 
    #              'spaw-h1a_mtVsibs_time_shape']
    groupsA = ['h1a241_mtVwt_time']
    groupsA = ['h1a241vs187_mxWts_50hpf']
    groupsA = ['h1a241vs187_NotMxWts_50hpf']
    
    groupsA = ['vcana_mtVwt_50hpf']
    
    groupsA = ['h1a241vs187_NotMxWts_50hpf']
    
    groupsA = ['h1a241vs187_NotMxWts_50hpf']
    
    
    groupsA =['spaw_mtVsibs_time_shape']
    groupsA =['spaw-h1a_mtVsibs_time_shape']
    groupsA = ['all_NotMx']
    #%%
    # while newGroup: 
    for info in groupsA:
        yes_plot = fcBasics.ask4input('>> Do you want to plot: '+ info+'? [0]:no, go to the next group / [1]: yes! >> :', bool)
        if yes_plot:
            #% Filter dataset
            print('>> FILTER FOR: ', info)
            df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
            fcAn.printDFINfo(df2plot)
            df2plot = fcAn.modifySibsGenot(df2plot)
            df2plot = fcAn.modifyOEControlGenot(df2plot)
            df2plot = fcAn.normalAsDextLoopers(df2plot)
            
            # save, info, ext = fcAn.q_savePlot()
            save = True; 
            # info = input('> Figures title: ')
            if 'NotMx' in info:
                GenotVar = 'GenotypeAll'
            else: 
                GenotVar = 'GenotypeF'
            
            # Settings for statistical analysis
            run_stats = fcBasics.ask4input('Do you want to run statistical analysis? [0]:no / [1]: yes! >> : ', bool)
            filters = ['Stage', GenotVar]; alpha = 0.05
            btw_x = False; btw_hue = True
            #%%
            # Individual Plots
            x_var = GenotVar
            hue_var = 'Stage'; 
            if 'shape' in info:
                shape_var = 'spAnalysis' 
            else: 
                shape_var = 'Strain_o' 
            # Settings for plots
            pl_indiv = fcAn.plot_indiv()
            vars_dict = fcAn.def_variables(plot_type = 'strip_plots')
            groups = list(pl_indiv.keys())[63:]#+list(pl_indiv.keys())[63:66]
            # groups = ['Vol_Ext.Myoc']#,'Vol_Atr.ExtMyoc','Vol_Vent.ExtMyoc']#,'Looping_Ratio_Myoc']
        
            for group in groups: 
                print('\n>>> ',group)
                vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group], 'indiv')
                dict_legends = fcAn.def_legends(df2plot)
                if 'shape' in info:
                    fcAn.plotIndivperXShape(# General Plot Settings
                           df2plot, vars2plot = vars2plot, x_var = x_var, hue_var = hue_var, shape_var = shape_var, 
                           # Size, title, labels and legends
                           title = pl_indiv[group]['title'], labels2plot = labels2plot, 
                           dict_legends = dict_legends, ips = ips[thesis_contxt_input], 
                           # Other plot settings
                           suptitle = False, right_legend = False, ctx_style = thesis_use.get_context(),
                           yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
                           yset = pl_indiv[group]['yset'],
                           box_plot = True, show_fliers = False,
                           # Saving settings
                           save = save, dpi = 300, ext = ['png', 'svg'], info = info+'-'+pl_indiv[group]['graph_no'], 
                           dir2save = os.path.join(dir_pl_meas,'pl_indiv'))
                else: 
                    fcAn.plotIndivperX(# General Plot Settings
                            df2plot, vars2plot = vars2plot, x_var = x_var, hue_var = hue_var, shape_var = shape_var, 
                            # Size, title, labels and legends
                            title = pl_indiv[group]['title'], labels2plot = labels2plot, 
                            dict_legends = dict_legends, ips = ips[thesis_contxt_input], 
                            # Other plot settings
                            suptitle = False, right_legend = False, ctx_style = thesis_use.get_context(),
                            yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
                            yset = pl_indiv[group]['yset'],
                            box_plot = True, show_fliers = False,
                            # Saving settings
                            save = save, dpi = 300, ext = ['png', 'svg'], info = info+'-'+pl_indiv[group]['graph_no'], 
                            dir2save = os.path.join(dir_pl_meas,'pl_indiv'))

                if run_stats: 
                    _, hue_values, _ = fcAn.props_ordered(df2plot, x_var, hue_var, shape_var)
                    dicts_stats = []
                    for hue_value in hue_values:
                        df_xfilt = df2plot[df2plot[hue_var] == hue_value]
                        box_pairs_all, box_pairs_f = fcAn.def_box_pairs(df_xfilt, x_var, hue_var, btw_x, btw_hue)
                        dicts_stats.append(fcAn.runStatisticalTests(data = df_xfilt, filters = filters, norm_test = 'Shapiro-Wilk', 
                                                              box_pairs_all = box_pairs_all, box_pairs_f = box_pairs_f, 
                                                              vars_group = vars2plot, alpha = alpha))
                    stats = True; stats_set = [stats, dicts_stats, alpha]
                    # Plot using the defined statistics in the dictionary
                    r_plStats = fcAn.plotIndivStatsTxtAllX(# General Plot Settings
                                        df2plot, vars2plot = vars2plot, x_var =  x_var, hue_var = hue_var, shape_var = shape_var, 
                                        # Size, title, labels and legends
                                        title = pl_indiv[group]['title'], 
                                        dict_legends = dict_legends,
                                        # Statistic Settings
                                        stats_set = stats_set, 
                                        # Other plot settings
                                        suptitle = True, ctx_style = thesis_use.get_context(),
                                        # Saving settings
                                        save = save, dpi = 300, ext = ['png'], info = info+'-'+pl_indiv[group]['graph_no'], 
                                        dir2save = os.path.join(dir_pl_meas,'pl_indiv'))
                if 'time' in info: 
                    fcAn.plotMultTimeCourse (# General Plot Settings
                            df2plot = df2plot, vars2plot = vars2plot, x_var=hue_var, hue_var=x_var,
                            # Size, title, labels and legends
                            title = pl_indiv[group]['title'], labels2plot = labels2plot, 
                            dict_legends = dict_legends, ips = (0.2,1.7), 
                            # Other plot settings
                            suptitle = False, right_legend = True, ctx_style = thesis_use.get_context(),
                            yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
                            yset = pl_indiv[group]['yset'],
                            # Saving settings
                            save = save, dpi = 300, ext = ['png', 'svg'], info = info+'-'+pl_indiv[group]['graph_no'], 
                            dir2save = os.path.join(dir_pl_meas,'pl_timecourse'))
        
        # newGroup = fcBasics.ask4input('Plot another filtered group? [0]: no, [1]: yes! >:', bool)
    
    #%% Looping Dir changing Shape
    #% Filter dataset
    # df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
    # fcAn.printDFINfo(df2plot)
    # df2plot = fcAn.modifySibsGenot(df2plot)
    # df2plot = fcAn.modifyOEControlGenot(df2plot)
    # df2plot = fcAn.normalAsDextLoopers(df2plot)
    
    # # save, info, ext = fcAn.q_savePlot()
    # save = True; 
    # info = input('> Figures title: ')
    # GenotVar = 'GenotypeF'
    # # GenotVar = 'GenotypeAll'
    
    # # Individual Plots
    # if not plotInGroups:  
    #     x_var = GenotVar
    #     hue_var = 'Stage'; shape_var = 'spAnalysis' 
    #     # Settings for plots
    #     pl_indiv = fcAn.plot_indiv()
    #     vars_dict = fcAn.def_variables(plot_type = 'strip_plots')
    #     # groups = list(pl_indiv.keys())[22:23]#[36:]
    #     groups = list(pl_indiv.keys())[63:65]#+list(pl_indiv.keys())[63:]#[36:]
    #     # groups = ['Vol_Ext.Myoc']#,'Vol_Atr.ExtMyoc','Vol_Vent.ExtMyoc']#,'Looping_Ratio_Myoc']
    
    #     for group in groups: 
    #         print('\n>>> ',group)
    #         vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group], 'indiv')
    #         dict_legends = fcAn.def_legends(df2plot)
            
    #         fcAn.plotIndivperXShape(# General Plot Settings
    #                            df2plot, vars2plot = vars2plot, x_var = x_var, hue_var = hue_var, shape_var = shape_var, 
    #                            # Size, title, labels and legends
    #                            title = pl_indiv[group]['title'], labels2plot = labels2plot, 
    #                            dict_legends = dict_legends, ips = ips[thesis_contxt_input], 
    #                            # Other plot settings
    #                            suptitle = False, right_legend = False, ctx_style = thesis_use.get_context(),
    #                            yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
    #                            yset = pl_indiv[group]['yset'],
    #                            box_plot = True, show_fliers = False,
    #                            # Saving settings
    #                            save = save, dpi = 300, ext = ['png', 'svg'], info = info+'-'+pl_indiv[group]['graph_no'], 
    #                            dir2save = os.path.join(dir_pl_meas,'pl_indiv'))


    #%% Looping Dir as x_var
    # xvar_loopingDir = False
    # if xvar_loopingDir: 
    #     #% Filter dataset
    #     df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
    #     fcAn.printDFINfo(df2plot)
    #     df2plot = fcAn.modifySibsGenot(df2plot)
    #     df2plot = fcAn.modifyOEControlGenot(df2plot)
        
    #     # save, info, ext = fcAn.q_savePlot()
    #     save = True; 
    #     info = input('> Figures title: ')
    #     GenotVar = 'GenotypeF'
    #     # GenotVar = 'GenotypeAll'
        
    #     # Individual Plots
    #     if not plotInGroups:  
    #         x_var = 'spAnalysis'
    #         hue_var = GenotVar; shape_var = 'Strain_o' 
    #         # Settings for plots
    #         pl_indiv = fcAn.plot_indiv()
    #         vars_dict = fcAn.def_variables(plot_type = 'strip_plots')
    #         groups = list(pl_indiv.keys())[15:25]+list(pl_indiv.keys())[63:]#[23:25]#[36:]
    #         # groups = ['Vol_Ext.Myoc']#,'Vol_Atr.ExtMyoc','Vol_Vent.ExtMyoc']#,'Looping_Ratio_Myoc']
        
    #         for group in groups: 
    #             print('\n>>> ',group)
    #             vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group], 'indiv')
    #             dict_legends = fcAn.def_legends(df2plot)
                
    #             fcAn.plotIndivperX(# General Plot Settings
    #                                df2plot, vars2plot = vars2plot, x_var = x_var, hue_var = hue_var, shape_var = shape_var, 
    #                                # Size, title, labels and legends
    #                                title = pl_indiv[group]['title'], labels2plot = labels2plot, 
    #                                dict_legends = dict_legends, ips = ips[thesis_contxt_input], 
    #                                # Other plot settings
    #                                suptitle = False, right_legend = False, ctx_style = thesis_use.get_context(),
    #                                yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
    #                                yset = pl_indiv[group]['yset'],
    #                                box_plot = True, show_fliers = False,
    #                                # Saving settings
    #                                save = save, dpi = 300, ext = ['png', 'svg'], info = info+'-'+pl_indiv[group]['graph_no'], 
    #                                dir2save = os.path.join(dir_pl_meas,'pl_indiv'))
               
#%%
    pl_indiv = fcAn.plot_indiv()
    for n, var in enumerate(pl_indiv.keys()):
        print(str(n)+'. -'+var)
        
    #%%  Plots through time for just one genotype!
    #%Plot for filtered wild-types (all strains) - working for document style!
    groupsB = [
                # 'all_mixedWt_time'#,
                # 'h1a241_mts_time', 
                'spaw_mts_time_shape', 
                # 'spaw-h1a_mts_time_shape']
                ]
    # groupsB = ['spaw_mts_time_shape', 
    #             'spaw-h1a_mts_time_shape']
    # groupsB = ['h1a241_mts_time']
    
    for info in groupsB:
        print(info)
        if info == 'all_mixedWt_time':
            # % Filter dataset (Option 1) - all wild-types combined
            filters = ['GenotypeF']; groups = [('wt:wt')]
            df2plot = fcAn.filterR_Autom(df_meas, filters, groups[0])
            df2plot = df2plot[df2plot['GenotypeAll'] != 'wt_tc:wt']
            uni_color = False
            fcAn.printDFINfo(df2plot)
        elif info == 'h1a241_mts_time':
            # % Filter dataset (Option 2) - only hapln1a241 mts
            filters = ['GenotypeAll']; groups = [('hapln1a241:mt')]
            df2plot = fcAn.filterR_Autom(df_meas, filters, groups[0])
            uni_color = True
            fcAn.printDFINfo(df2plot)
        elif info == 'spaw_mts_time_shape':
            # % Filter dataset (Option 2) - only spaw mts
            filters = ['GenotypeAll']; groups = [('hapln1a241:wt/spaw:mt')]
            df2plot = fcAn.filterR_Autom(df_meas, filters, groups[0])
            uni_color = True
            fcAn.printDFINfo(df2plot)
        elif info == 'spaw-h1a_mts_time_shape':
            # % Filter dataset (Option 2) - only spaw mts
            filters = ['GenotypeAll']; groups = [('hapln1a241:mt/spaw:mt')]
            df2plot = fcAn.filterR_Autom(df_meas, filters, groups[0])
            uni_color = True
            fcAn.printDFINfo(df2plot)
        else: 
            df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
            fcAn.printDFINfo(df2plot)
            uni_color = True
            
        save = True; ext = ['png', 'svg']; 
        # info = input('> Figures title: ')
        
        # Settings for statistical analysis
        run_stats = True#fcBasics.ask4input('Do you want to run statistical analysis? [0]:no / [1]: yes! >> : ', bool)
        filters = ['Stage', 'GenotypeF']; alpha = 0.05
        btw_x = False; btw_hue = True
        x_var = 'Stage'; 
        hue_var = 'GenotypeF'; 
        if 'shape' in info:
            shape_var = 'spAnalysis' 
        else: 
            shape_var = 'Strain_o' 
        # Settings for plots
        pl_indiv = fcAn.plot_indiv()
        vars_dict = fcAn.def_variables(plot_type = 'strip_plots')
        groups = list(pl_indiv.keys())[63:]
        # groups = ['Vol_Ext.Myoc']#,'Vol_Atr.ExtMyoc','Vol_Vent.ExtMyoc']#,'Looping_Ratio_Myoc']
        for group in groups: 
            print(group)
            vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group], 'indiv')
            dict_legends = fcAn.def_legends(df2plot)
            if 'shape' in info:
                fcAn.plotIndivperXShape(# General Plot Settings
                                df2plot, vars2plot = vars2plot, x_var = x_var, hue_var = hue_var, shape_var = shape_var, 
                                # Size, title, labels and legends
                                title = pl_indiv[group]['title'], labels2plot = labels2plot, 
                                dict_legends = dict_legends, ips = ips[thesis_contxt_input], 
                                # Other plot settings
                                suptitle = False, right_legend = False, 
                                ctx_style = thesis_use.get_context(), uni_color = uni_color, 
                                yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
                                yset = pl_indiv[group]['yset'],
                                box_plot = True, show_fliers = False,
                                # Saving settings
                                save = save, dpi = 300, ext = ['png', 'svg'], info = info+'-'+pl_indiv[group]['graph_no'], 
                                dir2save = os.path.join(dir_pl_meas,'pl_indiv'))
            else: 
                fcAn.plotIndivperX(# General Plot Settings
                                df2plot, vars2plot = vars2plot, x_var = x_var, hue_var = hue_var, shape_var = shape_var, 
                                # Size, title, labels and legends
                                title = pl_indiv[group]['title'], labels2plot = labels2plot, 
                                dict_legends = dict_legends, ips = ips[thesis_contxt_input], 
                                # Other plot settings
                                suptitle = False, right_legend = False, 
                                ctx_style = thesis_use.get_context(), uni_color = uni_color, 
                                yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
                                yset = pl_indiv[group]['yset'],
                                box_plot = True, show_fliers = False,
                                # Saving settings
                                save = save, dpi = 300, ext = ['png', 'svg'], info = info+'-'+pl_indiv[group]['graph_no'], 
                                dir2save = os.path.join(dir_pl_meas,'pl_indiv'))
            if run_stats: 
                _, hue_values, _ = fcAn.props_ordered(df2plot, x_var, hue_var, shape_var)
                dicts_stats = []
                for hue_value in hue_values:
                    df_xfilt = df2plot[df2plot[hue_var] == hue_value]
                    box_pairs_all, box_pairs_f = fcAn.def_box_pairs(df_xfilt, x_var, hue_var, btw_x, btw_hue)
                    dicts_stats.append(fcAn.runStatisticalTests(data = df_xfilt, filters = filters, norm_test = 'Shapiro-Wilk', 
                                                          box_pairs_all = box_pairs_all, box_pairs_f = box_pairs_f, 
                                                          vars_group = vars2plot, alpha = alpha))
                stats = True; stats_set = [stats, dicts_stats, alpha]
                # Plot using the defined statistics in the dictionary
                r_plStats = fcAn.plotIndivStatsTxtAllX(# General Plot Settings
                                    df2plot, vars2plot = vars2plot, x_var =  x_var, hue_var = hue_var, shape_var = shape_var, 
                                    # Size, title, labels and legends
                                    title = pl_indiv[group]['title'], 
                                    dict_legends = dict_legends,
                                    # Statistic Settings
                                    stats_set = stats_set, 
                                    # Other plot settings
                                    suptitle = True, ctx_style = thesis_use.get_context(),
                                    # Saving settings
                                    save = save, dpi = 300, ext = ['png'], info = info+'-'+pl_indiv[group]['graph_no'], 
                                    dir2save = os.path.join(dir_pl_meas,'pl_indiv'))
            if 'time' in info: 
                fcAn.plotMultTimeCourse (# General Plot Settings
                        df2plot = df2plot, vars2plot = vars2plot, x_var=x_var, hue_var=hue_var,
                        # Size, title, labels and legends
                        title = pl_indiv[group]['title'], labels2plot = labels2plot, 
                        dict_legends = dict_legends, ips = (0.2,1.7), 
                        # Other plot settings
                        suptitle = False, right_legend = True, ctx_style = thesis_use.get_context(),
                        yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
                        yset = pl_indiv[group]['yset'],
                        # Saving settings
                        save = save, dpi = 300, ext = ['png', 'svg'], info = info+'-'+pl_indiv[group]['graph_no'], 
                        dir2save = os.path.join(dir_pl_meas,'pl_timecourse'))
    
    #%% Pointplots of mts and wts through time
#     df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
#     # Setting spaw sibs (spaw:ht, h1a:wt, as wts to mix with h1a wts)
#     df2plot = fcAn.modifySibsGenot(df2plot)
#     fcAn.printDFINfo(df2plot)
#     save = True; ext = ['png','svg']
#     info = input('> Figures title: ')
#     x_var = 'Stage'; 
#     hue_var = 'GenotypeF'; shape_var = 'Strain_o' 
    
#     pl_indiv = fcAn.plot_indiv()
#     vars_dict = fcAn.def_variables(plot_type = 'strip_plots')
# #Add function to select variable(s) to plot
#     groups = list(pl_indiv.keys())[0:1]
    
#     for group in groups: 
#         print(group)
#         vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group], 'indiv')
#         dict_legends = fcAn.def_legends(df2plot)
#         fcAn.plotMultTimeCourse (# General Plot Settings
#                     df2plot = df2plot, vars2plot = vars2plot, x_var=x_var, hue_var=hue_var,
#                     # Size, title, labels and legends
#                     title = pl_indiv[group]['title'], labels2plot = labels2plot, 
#                     dict_legends = dict_legends, ips = (0.2,1.7), 
#                     # Other plot settings
#                     suptitle = False, right_legend = True, ctx_style = thesis_use.get_context(),
#                     yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
#                     yset = pl_indiv[group]['yset'],
#                     # Saving settings
#                     save = save, dpi = 300, ext = ext, info = info+'-'+pl_indiv[group]['graph_no'], 
#                     dir2save = os.path.join(dir_pl_meas,'pl_timecourse'))
    
    #%% Time course (agarose experiment) - working for document style!
    df2plot_tc = df_meas[df_meas['GenotypeAll'] == 'wt_tc:wt']
    fcAn.printDFINfo(df2plot_tc)
    save = True; info = 'wt_timecourse'; ext = ['png', 'svg']
    # save, info, ext = fcAn.q_savePlot()

    x_var = 'Stage'; hue_var = 'Fish_ref'
    pl_indiv = fcAn.plot_indiv()
    vars_dict = fcAn.def_variables(plot_type = 'strip_plots')
#Add function to select variable(s) to plot
    groups = list(pl_indiv.keys())[0:39]
    
    for group in groups: 
        print(group)
        vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group], 'indiv')
        dict_legends = fcAn.def_legends(df2plot_tc)
        fcAn.plotIndivTimeCourse (# General Plot Settings
                    df2plot = df2plot_tc, vars2plot = vars2plot, x_var=x_var, hue_var=hue_var,
                    # Size, title, labels and legends
                    title = pl_indiv[group]['title'], labels2plot = labels2plot, 
                    dict_legends = dict_legends, ips = (0.2,1.7), 
                    # Other plot settings
                    suptitle = False, right_legend = True, ctx_style = thesis_use.get_context(),
                    yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
                    # Saving settings
                    save = save, dpi = 300, ext = ext, info = info+'-'+pl_indiv[group]['graph_no'], 
                    dir2save = os.path.join(dir_pl_meas,'pl_timecourse'))
    
    #%% Tissue composition (whole, atrium and ventricle)
    print('\n=> TISSUE COMPOSITION ANALYSIS')
    df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
    fcAn.printDFINfo(df2plot)
    # save, info, ext = fcAn.q_savePlot()
    save = True; info = fcBasics.ask4input(' Title/info:', str) 
    ext = ['png','svg']
    
    # all_h1a_wtVmt
    # all_h1a187_wtVmt
    # all_h1aOE
    # all_spaw,h1a_sibsVmt
    # all_spaw_sibsVmt
    # all_wts
    
    colours = ['lightseagreen','darkorange','darkmagenta' ]
    vars2plot_all = [['Vol_Myoc','Vol_CJ','Vol_Endo'],
                     ['Vol_Atr.Myoc','Vol_Atr.CJ','Vol_Atr.Endo'],
                     ['Vol_Vent.Myoc','Vol_Vent.CJ','Vol_Vent.Endo']]
    add_vars2plot_all = ['Vol_Int.Endo','Vol_Atr.IntEndo','Vol_Vent.IntEndo']
    title_sp_all = ['Whole Heart ','Atrial ', 'Ventricular ']
    ylim_all = [[(0,1500e3),(0,2500e3)],[(0,1000e3),(0,1500e3)],[(0,600e3),(0,1000e3)]]
# Add function to select GenotypeAll or GenotypeF (joined wts) and x_var/hue_var
    # x = 'GenotypeAll'; col = 'Stage'; per = 'per Stage and Genotype'
    genotVar = 'GenotypeF'# 'GenotypeAll'
    x = 'Stage'; col = genotVar; per = 'per Genotype and Stage'#**
    group_vars = ['Stage', genotVar]
    txt_title = fcAn.get_txt_title(['Strain_o'], df2plot)
    dict_legends = fcAn.def_legends(df2plot)
    for n, vars2plot, add_var, title_sp, ylim in zip(count(), vars2plot_all, add_vars2plot_all, title_sp_all, ylim_all):
        for m, stack100, unit in zip(count(), [True, False], [' (%)', ' [$\mu$m$^3$]']):

            # => Only tissues 
            fcAn.barPlots(# General Plot Settings
                          df2plot = df2plot, vars2plot = vars2plot, group_vars = group_vars, 
                          x_var = x, col_var = col, 
                          # General Plot Settings
                          title = title_sp +'Tissue Composition '+ per, txt_title = txt_title,
                          ylabel = title_sp +'Tissue \nComposition '+unit, colours = colours, 
                          dict_legends = dict_legends, 
                          # Other plot settings
                          yticks_lab = '1e3 - d.', ylim = ylim[0], stack100 = stack100, sub_bar_lab = True,
                          ctx_style = thesis_use.get_context(), bot_legend = False, 
                          # Saving settings
                          save = save, ext = ext, info = info, dir2save = dir_pl_meas)
            # => Including lumen 
            fcAn.barPlots(# General Plot Settings
                          df2plot = df2plot, vars2plot = vars2plot + [add_var], group_vars = group_vars, 
                          x_var = x, col_var = col, 
                          # General Plot Settings
                          title = title_sp +'Composition '+ per, txt_title = txt_title,
                          ylabel = title_sp +'\nComposition'+unit, colours = colours + ['tomato'], 
                          dict_legends = dict_legends,
                          # Other plot settings
                          yticks_lab = '1e3 - d.', ylim = ylim[1], stack100 = stack100, sub_bar_lab = True,
                          ctx_style = thesis_use.get_context(), bot_legend = False, 
                          # Saving settings
                          save = save, ext = ext, info = info, dir2save = dir_pl_meas)
    
    #%% Tissue Expansion or Shrinkage (not corrected after dictLegends changed!)
    print('\n=> TISSUE EXPANSION OR SHRINKAGE ANALYSIS')
    df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
    fcAn.printDFINfo(df2plot)
    
    save = True; info = fcBasics.ask4input(' Title/info:', str) 
    ext = ['png','svg']
    # all_h1a_wtVmt
    # all_spaw,h1a_sibsVmt
    # all_spaw_sibsVmt
    # all_wts
    
    genot_var = 'GenotypeAll'
    group_vars = [genot_var, 'Stage'] # ['GenotypeAll', 'Strain_o', 'Stage']
    txt_title = fcAn.get_txt_title(['Strain_o'], df2plot)
    for  n, vars2plot, add_var, title_sp in zip(count(), vars2plot_all, add_vars2plot_all, title_sp_all):
        # print(n,vars2plot, add_var, title_sp)
        _, df4plot_pct_change = fcAn.get_dfPctChange(df2plot, vars2plot, group_vars)
        dict_legends = fcAn.def_legends(df4plot_pct_change.reset_index(), df_type = 'changes')
        # print(dict_legends)
        # => x = 'Stage', col = 'GenotypeAll'
        fcAn.pctChange_barPlots(# General Plot Settings
                                df2plot = df4plot_pct_change, vars2plot = [var+'_Change' for var in vars2plot], 
                                group_vars = group_vars, x_var = 'Stage', col_var = genot_var, 
                                # General Plot Settings
                                title = title_sp + 'Percentage Changes\nin Tissue Composition per Genotype and Stage', 
                                txt_title = txt_title, 
                                ylabel = 'Changes in tissue composition \nwith respect to previous stage (%)', 
                                colours = colours, dict_legends = dict_legends,
                                # Other plot settings
                                sub_bar_lab = True, 
                                ctx_style = thesis_use.get_context(), bot_legend = False, 
                                # Saving settings
                                save = save, ext = ext, info=info, dir2save = dir_pl_meas)
        
        # => x = 'GenotypeAll', col = 'Stage'
        # fcAn.pctChange_barPlots(df2plot = df4plot_pct_change, vars2plot = [var+'_Change' for var in vars2plot], 
        #                         group_vars = group_vars, x = 'GenotypeAll', col = 'Stage', colours = colours,
        #                         title = title_sp + 'Percentage Changes\nin Tissue Composition per Genotype and Stage', 
        #                         txt_title = txt_title, ylabel = 'Percentage changes in tissue composition \nwith respect to previous stage(%)', 
        #                         dir2save = dir_pl_meas, info=info, sub_bar_lab = True, save = save)
        
        # => Including lumen 
        vars2plot = vars2plot + [add_var]
        _, df4plot_pct_change = fcAn.get_dfPctChange(df2plot, vars2plot, group_vars)
        fcAn.pctChange_barPlots(# General Plot Settings
                                df2plot = df4plot_pct_change, vars2plot = [var+'_Change' for var in vars2plot], 
                                group_vars = group_vars, x_var = 'Stage', col_var = genot_var,
                                # General Plot Settings
                                title = title_sp + 'Percentage Changes\nin Composition per Genotype and Stage', 
                                txt_title = txt_title, 
                                ylabel = 'Changes in tissue composition \nwith respect to previous stage (%)', 
                                colours = colours+['tomato'], dict_legends = dict_legends,
                                # Other plot settings
                                sub_bar_lab = True, 
                                ctx_style = thesis_use.get_context(), bot_legend = False, 
                                # Saving settings
                                save = save, ext = ext, info=info, dir2save = dir_pl_meas)


    #%% Heatmaps 
    # - turbo colormap: https://ai.googleblog.com/2019/08/turbo-improved-rainbow-colormap-for.html
    # - https://github.com/matplotlib/matplotlib/issues/7081/
    # - https://matplotlib.org/stable/tutorials/colors/colormaps.html
    
    print('\n=> THICKNESS AND BALLOONING HEATMAP ANALYSIS')
    df_dataset = fcBasics.exportDatasetCSV(dir_data2Analyse, end_name = 'R', out_type = 'xlsx')
    df_dataset = fcAn.getMainStrain(fcAn.getGenotypeAll(df_dataset))
    # fcAn.printDFINfo(df_dataset)
    df_dataset_hm, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_dataset, [])
    fcAn.printDFINfo(df_dataset_hm)

    save, info, ext = fcAn.q_savePlot()
    filters = ['Stage', 'GenotypeAll', 'Strain_o']
    
    # groups = [('32-34','hapln1a241:wt','hapln1a prom241/+ InX'),('32-34','hapln1a241:mt','hapln1a prom241/+ InX')]
    # groups = [('48-50','hapln1a241:wt','hapln1a prom241/+ InX'),('48-50','hapln1a241:mt','hapln1a prom241/+ InX')]
    # groups = [('72-74','hapln1a241:wt','hapln1a prom241/+ InX'),('72-74','hapln1a241:mt','hapln1a prom241/+ InX')]
    
    # groups = [('32-34','hapln1a241:wt/spaw:ht','spaw+/-; hapln1a prom241/+ InX')]
    # groups = [('48-50','hapln1a241:wt/spaw:ht','spaw+/-; hapln1a prom241/+ InX')]
    groups = [('72-74','hapln1a241:wt/spaw:ht','spaw+/-; hapln1a prom241/+ InX')]
    # groups = [('48-50','vcana365:wt','vcana prom365/+ InX'),('48-50','vcana365:mt','vcana prom365/+ InX')]
    
    # groups = [('48-50','hapln1a187:wt','hapln1a prom187/+ InX'),('48-50','hapln1a187:mt','hapln1a prom187/+ InX')]
    
    
    # groups = [('48-50','galff:+/uas:+','myl7:galFF; UAS:TFP x UAS:hapln1a, cryaa:CFP'),
    #           ('48-50','galff:-/uas:+','myl7:galFF; UAS:TFP x UAS:hapln1a, cryaa:CFP'),
    #           ('48-50','galff:+/uas:-','myl7:galFF; UAS:TFP x UAS:hapln1a, cryaa:CFP')]
    
    # All wts
    df_dataset_hm = df_dataset_hm[df_dataset_hm['GenotypeAll'] != 'wt_tc:wt']
    filters = ['Stage', 'GenotypeF']
    groups = [('48-50','wt:wt')]
    # groups = [('32-34','wt:wt')]

    normalise = False; perChamber = False; norm_type = 'opt_div'; opt_norm = [normalise, perChamber, norm_type]
    for variable in ['CjTh', 'myocIntBall','MyocTh', 'EndoTh']:#'CjTh']:#, 
        for chamber in ['Atr', 'Vent']:
            fcAn.meanHM(df_dataset_hm = df_dataset_hm, filters = filters, groups = groups, chamber = chamber, 
                        variable = variable, opt_norm = opt_norm,  dir2load_df = dir_R_hmf, dir2save_hmf = dir_pl_meas, 
                        dir_data2Analyse = dir_data2Analyse, save = save, info = info)
            

    #%% Multiple linear regression for heatmaps
    # https://medium.com/swlh/interpreting-linear-regression-through-statsmodels-summary-4796d359035a
    angle_div = 90; position_div = 0.5
    for variable in ['CjTh']:#, 'myocIntBall','MyocTh', 'EndoTh']:
        df_mean, reglab = fcAn.stackRegHM(angle_div = angle_div, position_div = position_div, 
                        df_dataset_hm = df_dataset, df_meas = df_meas, var2an = variable, 
                        dir2load_df = dir_R_hmf, save = True, dir2save = dir_R_hm_reg)
        
    from statsmodels.formula.api import ols
    import statsmodels.api as sm
    # https://stackoverflow.com/questions/50733014/linear-regression-with-dummy-categorical-variables
    
    from sklearn import linear_model
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import r2_score
    from sklearn.metrics import mean_squared_error
    import seaborn as sns
            
    filters = ['Stage', 'Strain_o','Chamber']
    chambers = ['Atr']#, 'Vent']
    groups = [('32-34','hapln1a prom241/+ InX', i) for i in chambers]
    genotype_var = 'GenotypeAll'
              
    for variable in ['CjTh']:#, 'myocIntBall','MyocTh', 'EndoTh']:
        for group in groups:
            df_filt = fcAn.filterR_Autom (df_mean, filters, group)
            coVar =df_filt['coVar2add'].unique()[0]
            print('\n >> Variable: ', variable,' - Group: ', group, ' - Co-Variable: ', coVar)
            coVar_val = df_filt['coVar_value']
            df_reg = df_filt[[genotype_var, 'coVar_value']+reglab]
            df_reg["id"] = df_reg.index
            df_long = pd.wide_to_long(df_reg, stubnames='R-', i=['id',genotype_var], j='Region')
            df_long = df_long.rename(columns={'R-': variable+'_value'})
            df_long = df_long.reset_index().drop(columns = ['id'])
            df_long['Region'] = ['No'+str(num) for num in df_long['Region']]
            df_long_dum = pd.get_dummies(data = df_long, drop_first=True)
            formula = variable+'_value ~ ' 
            for var in df_long_dum.columns:
                print (var)
                if genotype_var in var:
                    df_long_dum = df_long_dum.rename(columns={var: var.replace(':', '')})
                if var != variable+'_value': 
                    formula = formula +' + ' + var

            fit = ols(formula, data=df_long_dum).fit() 
            fit.summary()
            
            x_df = df_long[[genotype_var,'Region','coVar_value']]
            y_df = df_long[[variable+'_value']]
            xd_df = pd.get_dummies(data = x_df, drop_first=True)

            x_train, x_test, y_train, y_test = train_test_split(xd_df, y_df, test_size = .20, random_state = 40)
    
            regr = linear_model.LinearRegression() # Do not use fit_intercept = False if you have removed 1 column after dummy encoding
            regr.fit(x_train, y_train)
            y_prediction = regr.predict(x_test)

            sns.regplot(y_test, y_prediction)
            
            # The coefficients
            print("Coefficients: \n", regr.coef_)
            # The mean squared error
            print("Mean squared error: %.2f" % mean_squared_error(y_test, y_prediction))
            # The coefficient of determination: 1 is perfect prediction
            print("Coefficient of determination: %.2f" % r2_score(y_test, y_prediction))
            # print the intercept
            print('Intercept: %.2f' % regr.intercept_)
            
            coeff_parameter = pd.DataFrame(list(regr.coef_[0]),list(xd_df.columns),columns=['Coefficient'])
            print(coeff_parameter)
            
            x_train_sm= sm.add_constant(x_train)
            x_train_sm= sm.add_constant(x_train)
            df_long_dum_train = x_train
            df_long_dum_train[variable+'_value'] = y_train
            formula = variable+'_value ~ ' 
            for var in df_long_dum_train.columns:
                print (var)
                if genotype_var in var:
                    df_long_dum_train = df_long_dum_train.rename(columns={var: var.replace(':', '')})
                if var != variable+'_value': 
                    formula = formula +' + ' + var
            
            ls = ols(formula, data=df_long_dum_train).fit()
            
            print(ls.summary())

    #%% Plotting different user's data
    user_analysis = False
    if user_analysis:
        df_meas = pd.concat((pd.read_csv(f) for f in all_CSVs))
        df_meas = df_meas.loc[:, ~df_meas.columns.str.contains('^Unnamed')]
        df_meas = fcAn.spAnalysis(fcAn.getMainStrain(fcAn.getGenotypeAll(df_meas)))
        
        #% Import Emma's wt data
        df_meas['User'] = pd.Series(['Juliana' for x in range(len(df_meas.index))])
        filters = ['GenotypeF']; groups = [('wt:wt')]
        df_meas = fcAn.filterR_Autom(df_meas, filters, groups[0])
        df_meas = df_meas[df_meas['GenotypeAll'] != 'wt_tc:wt']
        df_measEmma = pd.read_csv(os.path.join(dir_data2Analyse,'R_All', 'df_all','df_meas','R_temp','All_df_meas_EmmaWT.csv'))
        df_measEmma['User'] = pd.Series(['Emma' for x in range(len(df_measEmma.index))])
        df_meas = pd.concat([df_meas,df_measEmma])
        
        df_meas = df_meas.loc[:, ~df_meas.columns.str.contains('^Unnamed')]
        df2plot = []; df_cjPDF2plot = []
        # Add all other measurements to dataframe and apply format
        df_meas = fcAn.sortDFCols(fcAn.getVarRatios(fcAn.spAnalysis(fcAn.getMainStrain(fcAn.getGenotypeAll(df_meas)))))
        fcAn.printDFINfo(df_meas)
        
        #% Filter dataset
        # df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
        df2plot = df_meas[df_meas['Stage'] == '72-74']
        fcAn.printDFINfo(df2plot)
        save = True; run_stats = True
        info = 'UserComparison'
        GenotVar = 'GenotypeF'
        filters = ['User', GenotVar]; alpha = 0.05
        btw_x = True; btw_hue = False
            
        x_var = 'User'
        hue_var = GenotVar; shape_var = 'Strain_o' 
        # Settings for plots
        pl_indiv = fcAn.plot_indiv()
        vars_dict = fcAn.def_variables(plot_type = 'strip_plots')
        groups = list(pl_indiv.keys())[0:3]+list(pl_indiv.keys())[11:18]
    
        for group in groups: 
            print('\n>>> ',group)
            vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group], 'indiv')
            dict_legends = fcAn.def_legends(df2plot)
            
            fcAn.plotIndivperX(# General Plot Settings
                               df2plot, vars2plot = vars2plot, x_var = x_var, hue_var = hue_var, shape_var = shape_var, 
                               # Size, title, labels and legends
                               title = pl_indiv[group]['title'], labels2plot = labels2plot, 
                               dict_legends = dict_legends, ips = ips[thesis_contxt_input], 
                               # Other plot settings
                               suptitle = False, right_legend = False, ctx_style = thesis_use.get_context(),
                               yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
                               yset = pl_indiv[group]['yset'],
                               box_plot = True, show_fliers = False,
                               # Saving settings
                               save = save, dpi = 300, ext = ['png', 'svg'], info = info+'-'+pl_indiv[group]['graph_no'], 
                               dir2save = os.path.join(dir_pl_meas,'pl_indiv'))
        
            if run_stats: 
                _, hue_values, _ = fcAn.props_ordered(df2plot, x_var, hue_var, shape_var)
                dicts_stats = []
                for hue_value in hue_values:
                    print(hue_value)
                    df_xfilt = df2plot[df2plot[hue_var] == hue_value]
                    box_pairs_all, box_pairs_f = fcAn.def_box_pairs(df_xfilt, x_var, hue_var, btw_x, btw_hue)
                    box_pairs_all = [(('Emma', 'wt:wt'), ('Juliana', 'wt:wt'))]
                    box_pairs_f = [[(('Emma', 'wt:wt'), ('Juliana', 'wt:wt'))]]
                    dicts_stats.append(fcAn.runStatisticalTests(data = df_xfilt, filters = filters, norm_test = 'Shapiro-Wilk', 
                                                          box_pairs_all = box_pairs_all, box_pairs_f = box_pairs_f, 
                                                          vars_group = vars2plot, alpha = alpha))
                stats = True; stats_set = [stats, dicts_stats, alpha]
                # Plot using the defined statistics in the dictionary
                r_plStats = fcAn.plotIndivStatsTxtAllX(# General Plot Settings
                                    df2plot, vars2plot = vars2plot, x_var =  x_var, hue_var = hue_var, shape_var = shape_var, 
                                    # Size, title, labels and legends
                                    title = pl_indiv[group]['title'], 
                                    dict_legends = dict_legends,
                                    # Statistic Settings
                                    stats_set = stats_set, 
                                    # Other plot settings
                                    suptitle = True, ctx_style = thesis_use.get_context(),
                                    # Saving settings
                                    save = save, dpi = 300, ext = ['png'], info = info+'-'+pl_indiv[group]['graph_no'], 
                                    dir2save = os.path.join(dir_pl_meas,'pl_indiv'))
    
    #%% Plotting different tg lines data
    tg_analysis = False
    if tg_analysis:
        
        df_meas = pd.concat((pd.read_csv(f) for f in all_CSVs))
        df_meas = df_meas.loc[:, ~df_meas.columns.str.contains('^Unnamed')]
        df_meas = fcAn.spAnalysis(fcAn.getMainStrain(fcAn.getGenotypeAll(df_meas)))
        df_meas['User'] = pd.Series(['Juliana' for x in range(len(df_meas.index))])
        filters = ['GenotypeF']; groups = [('wt:wt')]
        df_meas = fcAn.filterR_Autom(df_meas, filters, groups[0])
        df_meas = df_meas[df_meas['GenotypeAll'] != 'wt_tc:wt']
        
        #% Import CJ heart data
        df_measTg = pd.read_csv(os.path.join(dir_data2Analyse,'R_All', 'df_all','df_meas','R_temp','LSCJ_F01_ResultsDF.csv'))
        df_measTg['User'] = pd.Series(['Anjalie' for x in range(len(df_measTg.index))])
        df_measTg = df_measTg.loc[:, ~df_measTg.columns.str.contains('^Unnamed')]
        df_meas = pd.concat([df_meas,df_measTg])
        
        df_meas = fcAn.getMainStrain(fcAn.getGenotypeAll(df_meas))
        df_meas = df_meas[df_meas['Stage'] == '48-50']
    
        df2plot = []; df_cjPDF2plot = []
        # Add all other measurements to dataframe and apply format
        df2plot = fcAn.sortDFCols((fcAn.getMainStrain(fcAn.getGenotypeAll(df_meas))))
        fcAn.printDFINfo(df_meas)
        
        #% Filter dataset
        # df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
        fcAn.printDFINfo(df2plot)
        save = True;
        info = 'TgComparison'
        GenotVar = 'GenotypeF'
            
        x_var = 'User'
        hue_var = GenotVar; shape_var = 'Strain_o'
        # Settings for plots
        pl_indiv = fcAn.plot_indiv()
        vars_dict = fcAn.def_variables(plot_type = 'strip_plots')
        groups = list(pl_indiv.keys())[0:3]+list(pl_indiv.keys())[11:18]
    
        for group in groups: 
            print('\n>>> ',group)
            vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group], 'indiv')
            dict_legends = fcAn.def_legends(df2plot)
            
            fcAn.plotIndivperX(# General Plot Settings
                               df2plot, vars2plot = vars2plot, x_var = x_var, hue_var = hue_var, shape_var = shape_var, 
                               # Size, title, labels and legends
                               title = pl_indiv[group]['title'], labels2plot = labels2plot, 
                               dict_legends = dict_legends, ips = ips[thesis_contxt_input], 
                               # Other plot settings
                               suptitle = False, right_legend = False, ctx_style = thesis_use.get_context(),
                               yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
                               yset = pl_indiv[group]['yset'],
                               box_plot = True, show_fliers = False,
                               # Saving settings
                               save = save, dpi = 300, ext = ['png', 'svg'], info = info+'-'+pl_indiv[group]['graph_no'], 
                               dir2save = os.path.join(dir_pl_meas,'pl_indiv'))

            
#%% OTHERS
others = False

if others:
    #%% Create plots for df_cjPDFs 
    print('\n=> KERNEL DENSITY ESTIMATE (KDE) PLOT ANALYSIS')
    all_CSVs = glob.glob(dir_R_cjPdfs + "/*.csv")
    df_cjPDF = pd.concat((pd.read_csv(f) for f in all_CSVs))
    df_cjPDF = df_cjPDF.loc[:, ~df_cjPDF.columns.str.contains('^Unnamed')]
    df_cjPDF = fcAn.getMainStrain(df_cjPDF)
    fcAn.printDFINfo(df_cjPDF,  df_type = 'kde')
    df_cjPDF2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_cjPDF, df_cjPDF2plot, df_type = 'kde')
    
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
    
    #%% - STATISTICS EASY FOR FUTURE DEVELOPMENT!!!
    import matplotlib.pyplot as plt
    import seaborn as sns
    from statannot import add_stat_annotation
    
    df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
    fcAn.printDFINfo(df2plot)
    
    groups = ['Surface Area', 'Heart Size','Lumen Size', 'Heart Looping', 'Tissue Myocardium',
             'Tissue Endocardium','Tissue Cardiac Jelly','Angles','Volume Percentages','Tissue Volume Percentages'] 
    group = groups[fcBasics.ask4input('Select variable from groups [int]:', int)]
    test2use = fcBasics.ask4input("test2use (Two vars: 't-test_ind'/'Mann-Whitney' - Three vars: 'Kruskal'): ", str, True)
    vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group])
    for var in vars2plot: 
        print(var)
        x = 'Stage' # factor
        y = var
        hue = 'GenotypeAll'
        hue_order = ['hapln1a:wt', 'hapln1a:mt']
        order = ['32-34', '48-50', '72-74']
        box_pairs = [(('32-34','hapln1a:wt'), ('32-34','hapln1a:mt')),
                      (('48-50','hapln1a:wt'), ('48-50','hapln1a:mt')),
                      (('72-74','hapln1a:wt'), ('72-74','hapln1a:mt'))]#, 
                      # (('32-34','hapln1a:wt'), ('48-50','hapln1a:wt')),
                      # (('48-50','hapln1a:wt'), ('72-74','hapln1a:wt')),
                      # (('32-34','hapln1a:wt'), ('72-74','hapln1a:wt'))]
        
        fig = plt.subplots()
        ax = sns.swarmplot(data = df2plot, x=x, y=y, hue=hue, order = order, dodge = True)
        ax1, test_results = add_stat_annotation(ax, data=df2plot, x=x, y=y, hue=hue, order = order,
                                                box_pairs=box_pairs,perform_stat_test=True, test = test2use,
                                                comparisons_correction=None,
                                                line_offset_to_box=0.1, line_offset=0.01,
                                                line_height=0.015, 
                                                text_format='star', loc='inside', verbose=2);
    
    #%% OTHERS TO CHECK LATER!!!
    #%% Plots with embryo ref labels
    # save = True
    # titles = ['Heart and Lumen Size','Tissue Layer Volumes','Surface Areas','Heart Looping','Angles', 'Volume Percentages']
    # input_vars = ['10,32,33,11,34,35','17-19,20-22,23-25','0-8','13,15,30','27-29', '36-38']
    # fcAn.plotPerVariableLabels('morphoHeart_D', input_vars, titles, df2plot, gen_legend, strain_legend , stage_legend,
    #                  h_plot = 8, w_plot = 6, save = save, dir2save = dir_pl_meas, info = info, dpi = 300)
    
    
#%%
init = True

#%% OLD
#     df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
#     fcAn.printDFINfo(df2plot)
#     save, info, ext = fcAn.q_savePlot()
    
        
#         #%%
    
#     df2plot, genots, strains, strains_o, stages = fcAn.filterDF2Plot(df_meas, df2plot)
#     fcAn.printDFINfo(df2plot)
#     save, info, ext = fcAn.q_savePlot()
#     # %%
# # Add function to select GenotypeAll or GenotypeF (joined wts) and x_var/hue_var
#     # x_var = 'Stage'; hue_var = 'GenotypeAll'; shape_var = 'Strain_o'
#     # x_var = 'GenotypeAll'; hue_var = 'Stage'; shape_var = 'Strain_o' #*** WTS ONLY THIS for all plots
    
#     # x_var = 'spAnalysis'; hue_var = 'Stage'; shape_var = 'Strain_o' # wt, lf, rt, ct
    
#     # Settings for statistical analysis
#     run_stats = fcBasics.ask4input('Do you want to run statistical analysis? [0]:no / [1]: yes! >> : ', bool)
#     filters = ['Stage', 'GenotypeAll']; alpha = 0.05
#     btw_x = False; btw_hue = True
    
#     # Select for indiv or group plots
#     plotInGroups = fcBasics.ask4input('Graph [0] individual or [1] group plots? >> : ', bool)
    
#     #%%
#     # Individual Plots
#     if not plotInGroups:  
#         x_var = 'GenotypeAll'; 
#         hue_var = 'Stage'; shape_var = 'Strain_o' #*** WTS ONLY THIS for all plots
#         # Settings for plots
#         pl_indiv = fcAn.plot_indiv()
#         vars_dict = fcAn.def_variables(plot_type = 'strip_plots')
# #Add function to select variable(s) to plot
#         groups = list(pl_indiv.keys())[:]
#         groups = ['Vol_Ext.Myoc']#,'Vol_Atr.ExtMyoc','Vol_Vent.ExtMyoc']#,'Looping_Ratio_Myoc']
    
#         for group in groups: 
#             print(group)
#             vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group], 'indiv')
#             dict_legends = fcAn.def_legends(df2plot)
#             if not run_stats: 
#                 fcAn.plotIndivperX(# General Plot Settings
#                                    df2plot, vars2plot = vars2plot, x_var = x_var, hue_var = hue_var, shape_var = shape_var, 
#                                    # Size, title, labels and legends
#                                    title = pl_indiv[group]['title'], labels2plot = labels2plot, 
#                                    dict_legends = dict_legends, ips = (0.8,6), 
#                                    # Other plot settings
#                                    suptitle = False, right_legend = False,
#                                    yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
#                                    yset = pl_indiv[group]['yset'],
#                                    box_plot = True, show_fliers = False,
#                                    # Saving settings
#                                    save = save, dpi = 300, ext = ext, info = info, 
#                                    dir2save = os.path.join(dir_pl_meas,'pl_indiv'))
#             else: 
#                 _, hue_values, _ = fcAn.props_ordered(df2plot, x_var, hue_var, shape_var)
#                 dicts_stats = []
#                 for hue_value in hue_values:
#                     df_xfilt = df2plot[df2plot[hue_var] == hue_value]
#                     box_pairs_all, box_pairs_f = fcAn.def_box_pairs(df_xfilt, x_var, hue_var, btw_x, btw_hue)
#                     dicts_stats.append(fcAn.runStatisticalTests(data = df_xfilt, filters = filters, norm_test = 'Shapiro-Wilk', 
#                                                           box_pairs_all = box_pairs_all, box_pairs_f = box_pairs_f, 
#                                                           vars_group = vars2plot, alpha = alpha))
#                 stats = True; stats_set = [stats, dicts_stats, alpha]
#                 # Plot using the defined statistics in the dictionary
#                 r_plStats = fcAn.plotIndivStatsperX(# General Plot Settings
#                                     df2plot, vars2plot = vars2plot, x_var =  x_var, hue_var = hue_var, shape_var = shape_var, 
#                                     # Size, title, labels and legends
#                                     title = pl_indiv[group]['title'], labels2plot = labels2plot, 
#                                     dict_legends = dict_legends, ips = (1,6), 
#                                     # Statistic Settings
#                                     stats_set = stats_set, statsTxt = True, 
#                                     # Other plot settings
#                                     suptitle = False, 
#                                     yticks_lab = pl_indiv[group]['yticks_lab'], ylim = pl_indiv[group]['ylim'], 
#                                     yset = pl_indiv[group]['yset'],
#                                     box_plot = True, show_fliers = False,
#                                     # Saving settings
#                                     save = save, dpi = 300, ext = ext, info = info, 
#                                     dir2save = os.path.join(dir_pl_meas,'pl_indiv_stats'))
                    
#                 # test_res, box_pairs_all = r_plStats
#             # input()
#     #%%
#     # Group Plots # NOT WORKING!!!!!
#     else:
#         x_var = 'Stage'; hue_var = 'GenotypeAll'; shape_var = 'Strain_o'
#         # Settings for plots
#         pl_groups = fcAn.plot_groups()
#         vars_dict = fcAn.def_variables(plot_type = 'strip_plots')
# # Add function to select group
#         groups = list(pl_groups.keys())[:]
#         # ['Surface Area', 'Heart Size','Lumen Size', 'Heart Looping', 'Tissue Myocardium',
#         #           'Tissue Endocardium','Tissue Cardiac Jelly','Ratios Cardiac Jelly', 
#         #           'Sagittal Angles','Ventral Angles', 'Volume Percentages','Tissue Volume Percentages'] 
#         groups = ['Heart Size']#['Surface Area CJ']#,'Lumen Size', 'Heart Looping','Tissue Cardiac Jelly (AtrVent)']
    
#         for group in groups: 
#             vars2plot, labels2plot = fcAn.selectVariables_auto(vars_dict, [group], 'group')
#             dict_legends = fcAn.def_legends(df2plot)
#             if not run_stats: 
#                 fcAn.plotInGroups(#General Plot Settings
#                                   df2plot, vars2plot = vars2plot, x_var =  x_var, hue_var =  hue_var, shape_var = shape_var, 
#                                   # Size, title, labels and legends
#                                   title = pl_groups[group]['title'], labels2plot = labels2plot, 
#                                   dict_legends = dict_legends, ips = (6,6), 
#                                   # Other plot settings
#                                   n_cols = pl_groups[group]['n_cols'], h_add = 5, w_add = 1, sharey = False, 
#                                   yticks_lab = pl_groups[group]['yticks_lab'], ylim = pl_groups[group]['ylim'], 
#                                   # Saving settings
#                                   save = save, dpi = 300, ext = ext, info =info, 
#                                   dir2save = os.path.join(dir_pl_meas,'pl_groups'))
    
#             else: 
#                 # Define all the multiple comparisons to calculate defining box_pairs 
#                 box_pairs_all, box_pairs_f = fcAn.def_box_pairs(df2plot, x_var, hue_var, btw_x, btw_hue)
#                 dict_stats = fcAn.runStatisticalTests(data = df2plot, filters = filters, norm_test = 'Shapiro-Wilk', 
#                                                       box_pairs_all = box_pairs_all, box_pairs_f = box_pairs_f, 
#                                                       vars_group = vars2plot, alpha = alpha)
#                 stats = True; stats_set = [stats, dict_stats, alpha]
#                 # Plot using the defined statistics in the dictionary
#                 r_plStats= fcAn.plotInGroupsStats(# General Plot Settings
#                                     df2plot, vars2plot = vars2plot, x_var =  x_var, hue_var = hue_var, shape_var = shape_var, 
#                                     # Size, title, labels and legends
#                                     title = pl_groups[group]['title'], labels2plot = labels2plot, 
#                                     dict_legends = dict_legends, ips = (6,6), 
#                                     # Statistic Settings
#                                     stats_set = stats_set,
#                                     # Other plot settings
#                                     n_cols = pl_groups[group]['n_cols'], h_add = 5, w_add = 2, sharey = False, 
#                                     yticks_lab = pl_groups[group]['yticks_lab'], ylim = '', 
#                                     # Saving settings
#                                     save = save, dpi = 300, ext = ext, info = info, 
#                                     dir2save = os.path.join(dir_pl_meas,'pl_groups_stats'))
                    
#                 test_res, box_pairs_all = r_plStats

