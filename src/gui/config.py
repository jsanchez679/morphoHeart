# -*- coding: utf-8 -*-
'''
morphoHeart_config

@author: Juliana Sanchez-Posada
'''
#%% Imports - ########################################################
import os
from pathlib import Path
from pandas import read_csv
from collections import defaultdict


#%% ##### - Authorship - #####################################################
__author__     = 'Juliana Sanchez-Posada'
__license__    = 'MIT'
__maintainer__ = 'J. Sanchez-Posada'
__email__      = 'jsanchezposadam@gmail.com'
__website__    = 'https://github.com/jsanchez679/morphoHeart'

#%% config class
class mH_Config():
    def __init__(self):
        self.version = '3.0.0'
        self.gui_sound = (True, 'All')
        self.theme = 'Light'
        self.heart_default = False

        #Config for vedo plots
        self.txt_font = 'Dalim'
        self.leg_font = 'LogoType' # 'Quikhand' 'LogoType'  'Dalim'
        self.leg_width = 0.18
        self.leg_height = 0.2
        self.txt_size = 0.8
        self.txt_color = '#696969'
        self.txt_slider_size = 0.8

        self.path_o = os.path.abspath(__file__)
        self.path_mHImages = Path(self.path_o).parent.parent.parent / 'images'
        self.path_ui = Path(self.path_o).parent.parent / 'gui' / 'ui'
        self.path_themes = Path(self.path_o).parent.parent / 'gui' / 'themes'
        self.path_templates = Path(self.path_o).parent.parent.parent / 'db' / 'templates'
        self.dev = False
        self.dev_plots = False
        self.dev_hm3d2d = False
        self.path_gui = Path(self.path_o).parent.parent / 'gui' 

        df_links = read_csv(self.path_gui / 'mH_links.csv', header=0, sep=';', usecols=["key0", "key1", "link"])
        d = defaultdict(dict)
        for a, b, c in df_links.itertuples(index=False):
            d[a][b] = c

        self.dict_links = dict(d)

        self.link2docs = self.dict_links['General']['user_manual'] #'https://drive.google.com/file/d/1-w9N3_SNzqNrpCAmTTvwc5br08KA6IbV/view?usp=drive_link'
        self.link2github = self.dict_links['General']['github'] #'https://github.com/jsanchez679/morphoHeart'
        self.link2paper = self.dict_links['General']['pre_print'] #'https://www.biorxiv.org/content/10.1101/2024.02.19.580991v2'
        self.link2vedo_plotter = 'https://vedo.embl.es/docs/vedo/plotter.html#Plotter'

mH_config = mH_Config()

print('morphoHeart! - Loaded config')

