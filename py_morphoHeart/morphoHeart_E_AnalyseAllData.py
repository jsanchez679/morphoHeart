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
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import count
import math

init = False
# Verify working dir
def setWorkingDir (root_path, init):
    if not init:
        wd = os.path.dirname(os.path.abspath(__file__))
        if root_path != wd:
            os.chdir(wd)
            root_path = os.getcwd()
    # init = True
    print("Current working directory: {0}".format(os.getcwd()))

    return root_path, init

root_path, init = setWorkingDir(os.getcwd(),init)

#%% Start D_AnalyseData
if init:
    # Importing morphoHeart packages
    from morphoHeart_modules import morphoHeart_funcAnalysis as fcAn
    from morphoHeart_modules import morphoHeart_funcBasics as fcBasics
    from morphoHeart_modules import morphoHeart_funcPlot as fcPlot
    # import morphoHeart_funcContours as fcCont
    # import morphoHeart_funcMeshes as fcMeshes

    #%% Get main directories (check which ones are actually used)
    _, _, dir_data2Analyse = fcBasics.getMainDirectories(root_path)
    dir_R_meas = os.path.join(dir_data2Analyse,'R_All', 'df_meas')
    dir_pl_meas = os.path.join(dir_data2Analyse,'R_All', 'pl_meas','R_')
    dir_R_cjPdfs = os.path.join(dir_data2Analyse,'R_All', 'df_cjPDFs')

    #%% Create plots for df_measurements 
    
    # Get directories
    all_CSVs = glob.glob(dir_R_meas + "/*.csv")
    df_meas = pd.concat((pd.read_csv(f) for f in all_CSVs))
    df_meas['Looping_Ratio_Myoc'] = df_meas['Length_CL_Int.Myoc(Cut)']/ df_meas['linLine_Int.Myoc(Cut)']
    df_meas['Looping_Ratio_Endo'] = df_meas['Length_CL_Ext.Endo(Cut)']/ df_meas['linLine_Ext.Endo(Cut)']
    
    df_meas['GenotA'] = df_meas['Gene_A']+':'+df_meas['Genotype_A']
    df_meas['GenotB'] = df_meas['Gene_B']+':'+df_meas['Genotype_B']
    df_meas['GenotypeAll'] = df_meas['GenotA']+'/'+df_meas['GenotB']
    
    # Generate result using pandas 
    genotypeAll = [] 
    for i, genotA, genotB in zip(count(), df_meas["GenotA"], df_meas["GenotB"]): 
        if genotB == "-:-": 
            genotypeAll.append(genotA) 
        else: 
            genotypeAll.append(genotA+'/'+genotB) 
    df_meas["GenotypeAll"] = genotypeAll 
    
    #%% Plot results
    variables, ylabels = fcAn.def_variables('morphoHeart_D_AnalyseAllData')
    save = True
    
    #%%
    titles = ['Surface Areas','Heart and Lumen Size','Heart Looping','Tissue Layer Volumes','Angles']
    input_vars = ['0-8','10,32,33,11,34,35','13,15,30','17-19,20-22,23-25','27-29']
    
    #ask4context = # notebook, talk, poster, paper
    for i, input_var, title in zip(count(), input_vars, titles):
        vars2plot, labels2plot = fcAn.getVarsANDLabels_Autom(variables, ylabels, input_var)
        
        sns.set_style("ticks")
        # Set up the matplotlib figure
        num_vars = len(vars2plot)
        plots_per_col = 3
        plots_per_row = math.ceil(num_vars/plots_per_col)
        size_col = plots_per_col*8+5
        size_row = plots_per_row*6+1
        
        # Number of genotypes: 
        genots = sorted(df_all.GenotypeAll.unique())
        n_gen = len(genots)
        strains = sorted(df_all.Strain.unique())
        n_strain = len(strains)
        
        # Save a palette to a variable:
        palette = sns.color_palette("husl", 8*n_gen*n_strain)
        # sns.palplot(palette)
        
        palettes = []
        for strain in range(n_strain):
            palettes.append(palette[strain*len(palette)//n_strain+2])
        
        # plt.clf()
        fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col, figsize=(size_col, size_row), sharex=False, sharey=False)
        fig.subplots_adjust(hspace=0.5, wspace=0.8)
        
        sns.set_context('poster', font_scale = 0.8, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
                                                        "ytick.labelsize":'small', "lines.linewidth": 2.5,
                                                        "xtick.major.size": 10, "ytick.major.size": 10, "figure.titlesize" :"large"})
        
        marker_size = 10
        dodge = False
        jitter = 0.2
        if n_strain == 2:
            df_plot_0 = df_meas[df_meas['Strain'] == strains[0]]
            df_plot_1 = df_meas[df_meas['Strain'] == strains[1]]
            for ax, var, ylabel in zip(axes.flatten(), vars2plot, labels2plot):
    
                m_0 = sns.stripplot(x="Stage", y=var, hue="GenotypeAll", data=df_plot_0, ax = ax, order=['32-34','48-50','72-74'],
                              marker = 'o', palette = palettes, jitter=jitter, dodge = dodge, size = marker_size)
    
                m_1 = sns.stripplot(x="Stage", y=var, hue="GenotypeAll", data=df_plot_1, ax = ax, order=['32-34','48-50','72-74'],
                              marker = '^', palette = palettes, jitter=jitter, dodge = dodge, size = marker_size)
                box = ax.get_position()
                ax.set(xlabel="Stage [hpf]", ylabel=ylabel);
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                
                handles0, labels0 = m_0.get_legend_handles_labels()
                handles1, labels1 = m_1.get_legend_handles_labels()
            
                # ax.legend(handles0, labels0, loc='center left', bbox_to_anchor=(1, 0.5))
                ax.legend(handles0[1:],labels[1:],loc='center left', bbox_to_anchor=(1, 0.5))
                sns.despine()
                
        fig.suptitle(title, fontsize = 30, y=1)
        
        if save: 
            plt.savefig(dir_pl_meas+title+".png", dpi=300, bbox_inches='tight', transparent=True)

#%% 
    # perhaps use hue for strain and the different symbpl for the genotype nad use swarmplots instead? and in this case dodge = True
    #https://seaborn.pydata.org/generated/seaborn.swarmplot.html 
    # ax = sns.swarmplot(x="day", y="total_bill", hue="smoker",
    #                data=tips, palette="Set2", dodge=True)
    
    #%%
    vars2plot, labels2plot = fcAn.getVarsANDLabels_Autom(variables, ylabels, input_var)
    var = vars2plot[0]
    fig, axes = plt.subplots(figsize=(12, 5), sharex=False, sharey=False)
    g = sns.catplot(x="Strain", y=var,
                    hue="GenotypeAll", col="Stage", hue_order = sorted(df_all.GenotypeAll.unique()),# col_order = sorted(df_all.Strain.unique()),
                    data=df_meas, kind="strip", aspect = 0.5 )#,
    # locs, labels = plt.xticks()
    # plt.setp(labels, rotation=45)
    g.set_xticklabels(rotation=30)
    g.set_axis_labels("Survived","All")
    

                   # height=4, aspect=.7);
    #https://dev.to/thalesbruno/subplotting-with-matplotlib-and-seaborn-5ei8
    #https://seaborn.pydata.org/tutorial/axis_grids.html
    #https://stackoverflow.com/questions/56788245/is-there-a-restriction-on-catplot-with-subplot
    
#%%
    plt.clf()
    thu_fri_sat = tips[tips['time'] == 'Lunch']#tips[(tips['day']=='Thur') | (tips['day']=='Fri') | (tips['day']=='Sat')]
    colors = ['blue','yellow','green','red']
    m = sns.stripplot('size','total_bill',hue='day',
                      marker='o',data=thu_fri_sat, jitter=0.1, 
                      palette=sns.xkcd_palette(colors),
                      dodge=True,linewidth=2,edgecolor="gray")
    
    sun = tips[tips['time'] == 'Dinner']#tips[tips['day']=='Sun']
    n = sns.stripplot('size','total_bill',color='red',hue='day',alpha=0.5,
                      marker='^',data=sun, jitter=0.1, 
                      dodge=True,linewidth=0)
    handles, labels = n.get_legend_handles_labels()
    n.legend(handles[:4], labels[:4])
    
#%% https://stackoverflow.com/questions/38650895/how-do-i-add-multiple-markers-to-a-stripplot-in-seaborn
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    tips = sns.load_dataset("tips")
    
    plt.clf()
    thu_fri_sat = tips[(tips['day']=='Thur') | (tips['day']=='Fri') | (tips['day']=='Sat')]
    colors = ['blue','yellow','green','red']
    m = sns.stripplot('size','total_bill',hue='day',
                      marker='o',data=thu_fri_sat, jitter=0.1, 
                      palette=sns.xkcd_palette(colors),
                      dodge=True,linewidth=2,edgecolor="gray")
    
    sun = tips[tips['day']=='Sun']
    n = sns.stripplot('size','total_bill',color='red',hue='day',alpha=0.5,
                      marker='^',data=sun, jitter=0.1, 
                      dodge=True,linewidth=0)
    handles, labels = n.get_legend_handles_labels()
    n.legend(handles[:4], labels[:4])


#%% # colors_pts = sns.color_palette("husl", 8)# sns.color_palette("tab10")#sns.hls_palette(2, h=.5, l=.4)
    # #Define number of variations/categories to plot
    # colors_pts = [colors_pts[0], colors_pts[4]]
    
    # sns.set_style("ticks", {"xtick.major.size": 4, "ytick.major.size": 4})
    # sns.set_style("darkgrid")
    # sns.set(rc={"xtick.bottom" : True, "ytick.left" : True, "xtick.minor.width":0.3, "xtick.major.size": 2,"ytick.labelsize":5, "lines.linewidth": 1.5})
    sns.set_context('poster', font_scale = 0.8, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
                                                    "ytick.labelsize":'small', "lines.linewidth": 2.5,
                                                    "xtick.major.size": 10, "ytick.major.size": 10})

    fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col, figsize=(size_col, size_row), sharex=False, sharey=False)
    fig.subplots_adjust(hspace=0.5, wspace=0.8)
    fig.suptitle(title)#, fontsize = 20)

    for ax, var, ylabel in zip(axes.flatten(), vars2plot, labels2plot):
        # sns.barplot(x="Stage", y=var, hue="Genotype_A", data=df_all, ax = ax, order=['32-34','48-50','72-74'],
        #               palette = colors_pts)
        sns.stripplot(x="Stage", y=var, hue="Genotype_A", data=df_all, ax = ax, order=['32-34','48-50','72-74'],
                      palette = colors_pts, jitter=0.1, dodge = True, size = 8)
        ax.legend(bbox_to_anchor=(1, 1), loc='best')
        # Add x-axis and y-axis labels
        ax.set(xlabel="Stage [hpf]",
                  ylabel=ylabel);
        sns.despine()
        # ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#%% for looping
f, ax = plt.subplots()
sns.violinplot(data=data)
sns.despine(offset=10, trim=True);


ax =  sns.stripplot(x="day", y="total_bill", hue="smoker",

                   data=tips, palette="Set2", size=20, marker="D",

                   edgecolor="gray", alpha=.25)

# Save a palette to a variable:
palette = sns.color_palette("Set2", 16)
 
# Use palplot and pass in the variable:
sns.palplot(palette)
    #%%
    # Set up the matplotlib figure
    num_vars = len(vars2plot)
    plots_per_col = 2
    plots_per_row = math.ceil(num_vars/plots_per_col)
    size_col = plots_per_col*3.5
    size_row = plots_per_row*3.5

    colors_pts = sns.hls_palette(2, h=.5, l=.4)

    # sns.set_style("ticks", {"xtick.major.size": 4, "ytick.major.size": 4})
    sns.set_style("darkgrid")
    # sns.set(rc={"xtick.bottom" : True, "ytick.left" : True, "xtick.minor.width":0.3, "xtick.major.size": 2,"ytick.labelsize":5, "lines.linewidth": 1.5})
    sns.set_context('poster', font_scale = 0.5, rc={"grid.linewidth": 0.7,"xtick.bottom" : True, "ytick.left" : True,
                                                    "ytick.labelsize":'small', "lines.linewidth": 1.5})

    fig, axes = plt.subplots(nrows=plots_per_row, ncols=plots_per_col, figsize=(size_col, size_row), sharex=False, sharey=True)
    fig.subplots_adjust(hspace=0.5)
    fig.suptitle(title)

    for ax, var, ylabel in zip(axes.flatten(), vars2plot, labels2plot):
        sns.stripplot(x="Stage", y=var, hue="Genotype_A", data=df_all, ax = ax, order=['32-34','48-50','72-74'],
                      palette = colors_pts)
        ax.legend(bbox_to_anchor=(1, 1), loc='best')
        # Add x-axis and y-axis labels
        ax.set(xlabel="Stage [hpf]",
                  ylabel=ylabel);
        sns.despine()

  



    #%%
    sns.set_style("white")
    A_vL = df_filtAtr.loc[df_filtAtr.Class=='A_vL', 'Thickness']
    A_vR = df_filtAtr.loc[df_filtAtr.Class=='A_vR', 'Thickness']
    A_dL = df_filtAtr.loc[df_filtAtr.Class=='A_dL', 'Thickness']
    A_dR = df_filtAtr.loc[df_filtAtr.Class=='A_dR', 'Thickness']

    # Plot
    kwargs = dict(hist_kws={'alpha':.6}, kde_kws={'linewidth':2})

    plt.figure(figsize=(10,7), dpi= 80)
    sns.distplot(A_vL, color="dodgerblue", label="A_vL", **kwargs)
    sns.distplot(A_vR, color="orange", label="A_vR", **kwargs)
    sns.distplot(A_dL, color="deeppink", label="A_dL", **kwargs)
    sns.distplot(A_dR, color="lime", label="A_dR", **kwargs)
    plt.xlim(-5,25)
    plt.ylim(0,0.6)
    plt.suptitle(title, y=0.95, size=12)
    plt.legend();

init = True
    #%%
#     import pandas as pd
#     import matplotlib.pyplot as plt
#
#     # Import Data
#     df = pd.read_csv('https://raw.githubusercontent.com/selva86/datasets/master/diamonds.csv')
#
#     # Plot
#     fig, axes = plt.subplots(1, 5, figsize=(10,2.5), dpi=100, sharex=True, sharey=True)
#     colors = ['tab:red', 'tab:blue', 'tab:green', 'tab:pink', 'tab:olive']
#
#     for i, (ax, cut) in enumerate(zip(axes.flatten(), df.cut.unique())):
#         print(cut)
#         x = df.loc[df.cut==cut, 'depth']
#         ax.hist(x, alpha=0.5, bins=100, density=True, stacked=True, label=str(cut), color=colors[i])
#         ax.set_title(cut)
#
#     plt.suptitle('Probability Histogram of Diamond Depths', y=1.05, size=16)
#     ax.set_xlim(50, 70); ax.set_ylim(0, 1);
#     plt.tight_layout();
#
#     import seaborn as sns
#     sns.set_style("white")
#
#     # Import data
#     #df = pd.read_csv('https://raw.githubusercontent.com/selva86/datasets/master/diamonds.csv')
#     x1 = df.loc[df.cut=='Ideal', 'depth']
#     x2 = df.loc[df.cut=='Fair', 'depth']
#     x3 = df.loc[df.cut=='Good', 'depth']
#
#     # Plot
#     kwargs = dict(hist_kws={'alpha':.6}, kde_kws={'linewidth':2})
#
#     plt.figure(figsize=(10,7), dpi= 80)
#     sns.distplot(x1, color="dodgerblue", label="Compact", **kwargs)
#     sns.distplot(x2, color="orange", label="SUV", **kwargs)
#     sns.distplot(x3, color="deeppink", label="minivan", **kwargs)
#     plt.xlim(50,75)
#     plt.legend();
#
# #  -------------------------- REVISAMEEEEE ------------------------- pallollito
#     # #%%
#     # g = sns.catplot(x="Stage", y="SurfArea_Myoc", hue="Genotype_A", data=df_all)
#     #
#     #                 col="diet",
#
#
#     fig, axs = plt.subplots(2,5, figsize=(15, 6), facecolor='w', edgecolor='k')
#     fig.subplots_adjust(hspace = .5, wspace=.001)
#
#     axs = axs.ravel()
#
#     for i in range(10):
#
#         axs[i].contourf(np.random.rand(10,10),5,cmap=plt.cm.Oranges)
#         axs[i].set_title(str(250+i))
#     #%%
#     fig = plt.figure()
#     fig.add_subplot(2,2,1)   #top left
#     fig.add_subplot(2,2,2)   #top right
#     fig.add_subplot(2,2,3)   #bottom left
#     fig.add_subplot(2,2,4)   #bottom right
#     plt.show()
#
#     import seaborn as sns
#
#     sns.set_theme(style="ticks")
#
#     exercise = sns.load_dataset("exercise")
#
#     g = sns.catplot(x="time", y="pulse", hue="kind", data=exercise)
#
#     #%%
#     import matplotlib.pyplot as plt
#     import seaborn as sns
#     import numpy as np
#     from sklearn.datasets import load_iris
#
#     data = load_iris()
#     fig, axes = plt.subplots(nrows=2, ncols=2)
#     fig.subplots_adjust(hspace=0.5)
#     fig.suptitle('Distributions of Iris Features')
#
#     for ax, feature, name in zip(axes.flatten(), data.data.T, data.feature_names):
#         sns.distplot(feature, ax=ax, bins=len(np.unique(data.data.T[0]))//2)
#         ax.set(title=name[:-4].upper(), xlabel='cm')
#
#
#     #%% Old codes
#     ax1 = sns.countplot(x="LeftRight""," data=df_ptsAtr)
#     #%%
#     g = sns.catplot(x="DorsVent", hue="LeftRight",# col="survived",
#                     data=df_ptsAtr, kind="count",
#                     height=4, aspect=.7);
#
#     # ax = df_ptsAtr.hist(column='Thickness', by='LeftRight', bins=10, grid=False,
#     #              figsize=(8,10), layout=(3,1), sharex=True, color='#86bf91',
#     #              zorder=2, rwidth=0.9)
#
#     # for i,x in enumerate(ax):
#
#     #     # Despine
#     #     x.spines['right'].set_visible(False)
#     #     x.spines['top'].set_visible(False)
#     #     x.spines['left'].set_visible(False)
#
#     #     # Switch off ticks
#     #     x.tick_params(axis="both", which="both", bottom="off", top="off",
#     #                   labelbottom="on", left="off", right="off", labelleft="on")
#
#     #     # Draw horizontal axis lines
#     #     vals = x.get_yticks()
#     #     for tick in vals:
#     #         x.axhline(y=tick, linestyle='dashed', alpha=0.4, color='#eeeeee', zorder=1)
#
#     #     # Set x-axis label
#     #     x.set_xlabel("Thickness", labelpad=20, weight='bold', size=12)
#
#     #     # Set y-axis label
#     #     if i == 1:
#     #         x.set_ylabel("Left-Right", labelpad=50, weight='bold', size=12)
#
#     #     # Format y-axis label
#     #     x.yaxis.set_major_formatter(StrMethodFormatter('{x:,g}'))
#     #     x.tick_params(axis='x', rotation=0)
#
#     #%% Filter dataframe #cardiac jelly thickness
#     df_filtThick = df_ptsClassAll[df_ptsClassAll['Thickness'] > 0]
#     df_filtThick.sample(20)
#
#     classif_list = []
#
#     for index, row in df_filtThick.iterrows():
#         classif = ''
#         AorV = row['Atr-Vent']
#         #Get Atrium
#         if AorV == 'atrium':
#             classif = classif +'A_'
#             #Dorsal-Ventral
#             DorV = row['Dors-Vent_Atr']
#             if DorV == 'dorsal':
#                 classif = classif +'d'
#             else:
#                 classif = classif +'v'
#         #Get Ventricle
#         if AorV == 'ventricle':
#             classif = classif +'V_'
#             #Dorsal-Ventral
#             DorV = row['Dors-Vent_Vent']
#             if DorV == 'dorsal':
#                 classif = classif +'d'
#             else:
#                 classif = classif +'v'
#
#         #Left-right
#         LorR = row['Left-Right']
#         if LorR == 'left':
#             classif = classif +'L'
#         else:
#             classif = classif +'R'
#
#         classif_list.append(classif)
#     cjf4.alert("wohoo",1)
#
#     df_filtThick['Class'] = classif_list
#     df_filtThick.drop(columns=['X_pts_cjOut', 'Y_pts_cjOut', 'Z_pts_cjOut'])
#
#     df_filtThick_title = filename[0]+'_df_filtThick.csv'
#     df_filtThick_dir = os.path.join(dir_txtNnpy, df_filtThick_title)
#     df_filtThick.to_csv(df_filtThick_dir)
#     cjf4.alert("countdown",1)
#
#     #%% Plot
#     df_filtAtr = df_filtThick[df_filtThick['Atr-Vent']=='atrium']
#
#     #%%
#     fig, axes = plt.subplots(1, 4, figsize=(10,2.5), dpi=100, sharex=True, sharey=True)
#     colors = ['tab:blue', 'tab:orange', 'tab:pink', 'tab:green', 'tab:olive']
#
#     for i, (ax, Class) in enumerate(zip(axes.flatten(), df_filtAtr.Class.unique())):
#         print(Class)
#         x = df_filtAtr.loc[df_filtAtr.Class==Class, 'Thickness']
#         ax.hist(x, alpha=0.7, bins=100, density=True, stacked=True, label=str(Class), color=colors[i])
#         ax.set_title(Class)
#
#     title = filename[0] + ' - Probability Histogram of Thickness'
#     plt.suptitle(title, y=1.05, size=12)
#     ax.set_xlim(-5, 25);
#     ax.set_ylim(0, 0.6);
#     plt.tight_layout();
#
#     #%%
#     sns.set_style("white")
#     A_vL = df_filtAtr.loc[df_filtAtr.Class=='A_vL', 'Thickness']
#     A_vR = df_filtAtr.loc[df_filtAtr.Class=='A_vR', 'Thickness']
#     A_dL = df_filtAtr.loc[df_filtAtr.Class=='A_dL', 'Thickness']
#     A_dR = df_filtAtr.loc[df_filtAtr.Class=='A_dR', 'Thickness']
#
#     # Plot
#     kwargs = dict(hist_kws={'alpha':.6}, kde_kws={'linewidth':2})
#
#     plt.figure(figsize=(10,7), dpi= 80)
#     sns.distplot(A_vL, color="dodgerblue", label="A_vL", **kwargs)
#     sns.distplot(A_vR, color="orange", label="A_vR", **kwargs)
#     sns.distplot(A_dL, color="deeppink", label="A_dL", **kwargs)
#     sns.distplot(A_dR, color="lime", label="A_dR", **kwargs)
#     plt.xlim(-5,25)
#     plt.ylim(0,0.6)
#     plt.suptitle(title, y=0.95, size=12)
#     plt.legend();
#
#     #%%
#     # import pandas as pd
#     # import matplotlib.pyplot as plt
#
#     # # Import Data
#     # df = pd.read_csv('https://raw.githubusercontent.com/selva86/datasets/master/diamonds.csv')
#
#     # # Plot
#     # fig, axes = plt.subplots(1, 5, figsize=(10,2.5), dpi=100, sharex=True, sharey=True)
#     # colors = ['tab:red', 'tab:blue', 'tab:green', 'tab:pink', 'tab:olive']
#
#     # for i, (ax, cut) in enumerate(zip(axes.flatten(), df.cut.unique())):
#     #     print(cut)
#     #     x = df.loc[df.cut==cut, 'depth']
#     #     ax.hist(x, alpha=0.5, bins=100, density=True, stacked=True, label=str(cut), color=colors[i])
#     #     ax.set_title(cut)
#
#     # plt.suptitle('Probability Histogram of Diamond Depths', y=1.05, size=16)
#     # ax.set_xlim(50, 70); ax.set_ylim(0, 1);
#     # plt.tight_layout();
#
#     # import seaborn as sns
#     # sns.set_style("white")
#
#     # # Import data
#     # #df = pd.read_csv('https://raw.githubusercontent.com/selva86/datasets/master/diamonds.csv')
#     # x1 = df.loc[df.cut=='Ideal', 'depth']
#     # x2 = df.loc[df.cut=='Fair', 'depth']
#     # x3 = df.loc[df.cut=='Good', 'depth']
#
#     # # Plot
#     # kwargs = dict(hist_kws={'alpha':.6}, kde_kws={'linewidth':2})
#
#     # plt.figure(figsize=(10,7), dpi= 80)
#     # sns.distplot(x1, color="dodgerblue", label="Compact", **kwargs)
#     # sns.distplot(x2, color="orange", label="SUV", **kwargs)
#     # sns.distplot(x3, color="deeppink", label="minivan", **kwargs)
#     # plt.xlim(50,75)
#     # plt.legend();
#
#     #%%
#     # #Hi @jsanchez679 ,
#     # have you already tried using vmin and vmax ? As in:
#     # https://github.com/marcomusy/vedo/blob/master/vedo/examples/basic/mesh_sharemap.py
#
#     # Another way could be using:
#     # mesh.mapper().SetScalarRange(x0,x1)
