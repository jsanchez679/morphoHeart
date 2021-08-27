# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 09:07:28 2021

@author: mdp18js
"""


from vedo import *

em0 = Volume("LS34_F06_V_SR_0956_ch0.tif").isosurface(80)
em1 = Volume("LS34_F06_V_SR_0956_ch0.tif").isosurface(threshold = True)
em2 = Volume("LS34_F06_V_SR_0956_ch0.tif").isosurface(500)
em3 = Volume("LS34_F06_V_SR_0956_ch0.tif").isosurface(700)

vp = Plotter(N=4, axes = 4)
vp.show(em0, at = 0)
vp.show(em1, at = 1)
vp.show(em2, at = 2)
vp.show(em3, at = 3, interactive = True)


em0 = Volume("LS34_F06_V_SR_0956_ch0.tif", spacing = (0.228,0.228,0.469)).isosurface(400)
em1 = Volume("LS34_F06_V_SR_0956_ch0.tif", spacing = (0.228,0.228,0.469)).isosurface(300)
em2 = Volume("LS34_F06_V_SR_0956_ch0.tif", spacing = (0.228,0.228,0.469)).isosurface(450)


vp = Plotter(N=3, axes = 7)
vp.show(em1, at = 0)
vp.show(em2, at = 1)
vp.show(em0, at=2, interactive = True)


