# -*- coding: utf-8 -*-
'''
morphoHeart_config

@author: Juliana Sanchez-Posada
'''
#%% Imports - ########################################################
import os
from pathlib import Path

#%% ##### - Authorship - #####################################################
__author__     = 'Juliana Sanchez-Posada'
__license__    = 'MIT'
__maintainer__ = 'J. Sanchez-Posada'
__email__      = 'jsanchezposadam@gmail.com'
__website__    = 'https://github.com/jsanchez679/morphoHeart'

#%% config class
class mH_Config():
    def __init__(self):
        self.version = '2.1.0'
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
        self.dev = False
        self.dev_plots = False
        self.dev_hm3d2d = False

        self.link2docs = 'https://www.dropbox.com/scl/fi/ft7i6t8d0hc859t7ii9a3/User-Manual-v1o0o0.pdf?rlkey=ebyp3xyjhau78kc3tqu5td8s0&dl=0'
        self.link2github = 'https://github.com/jsanchez679/morphoHeart'
        self.link2paper = 'https://www.biorxiv.org/content/10.1101/2024.02.19.580991v2'

mH_config = mH_Config()

print('morphoHeart! - Loaded config')

