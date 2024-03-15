'''
morphoHeart

@author: Juliana Sanchez-Posada
'''

#%% Imports - ########################################################
import sys
from PyQt6 import QtWidgets

#%% morphoHeart Imports - ##################################################
print('Welcome to morphoHeart!')
from src.gui.config import mH_config 
from src.gui.gui_classes import *
from src.modules import mH_classes_new as mHC
from src import mH_api as mA

#%% API
class Controller: 
    def __init__(self):
        self.welcome_win = None
        self.new_proj_win = None
        self.new_proj_win_from_temp = None
        self.meas_param_win = None
        self.load_proj_win = None
        self.load_multi_proj_win = None
        self.new_organ_win = None
        self.main_win = None
        self.multip_analysis_win = None
        self.load_s3s = None
        self.proj_settings_win = None
        self.organ_settings_win = None
        self.about_window = None

        self.wins = ['new_proj_win','meas_param_win','load_proj_win','load_multi_proj_win', 
                     'new_organ_win','main_win', 'multip_analysis_win', 'load_s3s', 
                     'proj_settings_win', 'organ_settings_win', 'about_window']

    def show_welcome(self):
        #Close previous windows if existent
        if self.new_proj_win != None:
            self.new_proj_win.close()
            self.new_proj_win = None
        if self.new_proj_win_from_temp != None:
            self.new_proj_win_from_temp.close()
            self.new_proj_win_from_temp = None
        if self.load_proj_win != None:
            self.load_proj_win.close()
            self.load_proj_win = None
        if self.load_multi_proj_win != None:
            self.load_multi_proj_win.close()
            self.load_multi_proj_win = None
        if self.multip_analysis_win != None: 
            self.multip_analysis_win.close()
            self.multip_analysis_win = None
        #Create welcome window and show
        if self.new_proj_win == None: 
            self.welcome_win = WelcomeScreen()
        self.welcome_win.show()

        #Connect Buttons
        # -Create new project
        self.welcome_win.button_new_proj.clicked.connect(lambda: self.show_create_new_proj())
        # -Load project
        self.welcome_win.button_load_proj.clicked.connect(lambda: self.show_load_proj())
        # -Create new project from template
        self.welcome_win.button_new_proj_from_template.clicked.connect(lambda: self.show_proj_from_temp())
        # -Multi-project analysis
        self.welcome_win.button_multi_proj.clicked.connect(lambda: self.show_multi_proj())

    def show_create_new_proj(self):
        #Close welcome window
        self.sound = self.welcome_win.sound
        self.welcome_win.close()
        #Create new proj window and show
        if self.new_proj_win == None: 
            self.new_proj_win = CreateNewProj(controller=self)
            self.init_create_new_proj()
        self.new_proj_win.show()

    def show_proj_from_temp(self): 
        #Close welcome window
        self.welcome_win.close()
        #Create new proj window and show
        if self.new_proj_win_from_temp == None: 
            self.new_proj_win_from_temp = NewProjFromTemp(controller=self)
            self.init_new_proj_from_temp()
        self.new_proj_win_from_temp.show()

    def show_meas_param(self):
        #Create meas param window and show
        if self.new_proj_win.check_to_set_params(): 
            if self.meas_param_win == None: 
                self.meas_param_win = SetMeasParam(mH_settings = self.new_proj_win.mH_settings, 
                                                parent=self.new_proj_win)
            else:
                #If the window to set measurements parameters had already been created, re initialise it with the channels selected. 
                self.meas_param_win.ch_all = self.new_proj_win.mH_settings['name_chs']
                if isinstance(self.new_proj_win.mH_settings['chNS'], dict) and len(list(self.new_proj_win.mH_settings['chNS'].keys()))>0:
                    if self.new_proj_win.mH_settings['chNS']['layer_btw_chs']:
                        self.meas_param_win.ch_all['chNS'] = self.new_proj_win.mH_settings['chNS']['user_nsChName']
                print(self.meas_param_win.ch_all)
                self.meas_param_win.set_meas_param_table()
            if self.new_proj_win.template: 
                self.meas_param_win.init_template(self.proj_temp)
            self.meas_param_win.show()
            self.meas_param_win.button_set_params.clicked.connect(lambda: self.set_proj_meas_param())
            self.new_proj_win.set_meas_param_all.setChecked(True)

    def show_load_proj(self): 
        #Close welcome window
        self.welcome_win.close()
        #Create Load Project Window and show
        if self.load_proj_win == None:
            self.load_proj_win = LoadProj() 
            self.init_load_proj()
        self.load_proj_win.show()
    
    def show_multi_proj(self):
        #Close welcome window
        self.welcome_win.close()
        #Create Load Multi-Project Window and show
        if self.load_multi_proj_win == None:
            self.load_multi_proj_win = Load_MultiProj() 
            self.init_load_multi_proj()
        self.load_multi_proj_win.show()

    def show_new_organ(self, parent_win:str):
        #Identify parent and close it
        if parent_win == 'new_proj_win':
            if self.new_proj_win.button_new_proj.isChecked():
                self.new_proj_win.close()
            else: 
                error_txt = "*Please create the New Project first before adding an organ to it."
                self.new_proj_win.win_msg(error_txt)
                self.new_proj_win.button_new_proj.setChecked(False)
                return
        elif parent_win == 'load_proj_win':
            self.load_proj_win.close()
        else: 
            print('Controller: show_new_organ: Other parent window?')
            alert('bubble')
        print('Controller, show_new_organ > parent_win:', parent_win)

        #Create new organ window and show
        if self.new_organ_win == None:
            self.new_organ_win = NewOrgan(proj = self.proj)
            self.init_new_organ_win(parent_win=parent_win)
        self.new_organ_win.show()

    def show_parent(self, win:str, parent:str):
        try:
            parent_win = getattr(self, parent)
            parent_win.show()
            getattr(self, win).close()
            setattr(self, win, None)

        except: 
            self.clear_win_show_welcome(parent=None)

    def show_main_window(self, parent_win:str):
        #Close new organ or load organ window 
        win = getattr(self, parent_win)
        if parent_win == 'new_organ_win':
            if self.new_organ_win.button_create_new_organ.isChecked():

                pass# self.new_organ_win.close()
            else: 
                error_txt = '*You need to first create the organ to continue.'
                self.new_organ_win.win_msg(error_txt)
                return
            
        elif parent_win == 'load_proj_win':
            self.load_proj_win.check_unique_organ_selected(self.proj)
            if self.load_proj_win.organ_selected != None:
                self.organ_to_analyse = self.load_proj_win.organ_selected#.replace(' ', '_')
                self.load_proj_win.win_msg('Loading organ '+self.organ_to_analyse+'...')
                self.load_organ(proj = self.proj, organ_to_load = self.organ_to_analyse)
                # self.load_proj_win.close()
            else: 
                if len(self.proj.organs) == 0: 
                    error_txt = "!The project selected does not contain organs. Add a new organ to this project by selecting 'Create New Organ'."
                    self.load_proj_win.win_msg(error_txt)
                else: 
                    error_txt = '*Please select one organ to open.'
                    self.load_proj_win.win_msg(error_txt)
                    print('Controller, show_main_window: Error in loading window')
                return
        else: 
            print('Controller, show_main_window: Other parent window?')
            alert('bubble')
        print('Controller, show_main_window > parent_win:', parent_win)

        #Create Main Project Window and show
        if self.main_win == None:
            self.main_win = MainWindow(proj = self.proj, organ = self.organ, controller=self) 
            self.init_main_win()
        win.close()
        setattr(self, parent_win, None)
        self.main_win.show()

    def show_analysis_window(self, parent_win:str, single_proj:bool): 
        print('This takes to Analysis Window from '+parent_win)
        win = getattr(self, parent_win)
        print('win.multi_organs_added:', win.multi_organs_added, '\n', len(win.multi_organs_added))
        if len(win.multi_organs_added) > 0:
            win.win_msg('Loading Organs in Analysis Window...')
            df_pando, dict_projs, dict_organs = self.load_multip_proj_and_organs(proj_org = win.multi_organs_added, single_proj=single_proj)
            # win.close()
        else: 
            error_txt = '*Please add organs to "Organs Added to Analysis" table to include in the combinatorial analysis.'
            win.win_msg(error_txt, win.go_to_analysis_window)
            return
        
        #Create Main Project Window and show
        if self.multip_analysis_win == None:
            if single_proj: 
                self.multip_analysis_win = MultipAnalysisWindow(projs = dict_projs, organs = dict_organs, df_pando= df_pando, controller=self) 
                self.init_multip_analysis_win()
                win.close()
                self.multip_analysis_win.show()
            else:
                pass 
                # self.multip_analysis_win = MainWindow(proj = self.proj, organ = self.organ, controller=self) 
                # self.init_multip_analysis_win()
                # self.multip_analysis_win.show()
    
    def show_load_closed_stacks(self):
        if self.load_s3s == None:
            self.load_s3s = Load_S3s(proj = self.proj, organ = self.organ, parent_win=self.main_win) 
        else: 
            self.load_s3s.show()
        
    def show_proj_settings(self, parent_win):
        if self.proj_settings_win == None:
            self.proj_settings_win = ProjSettings(proj = self.proj, controller=self) 
        self.proj_settings_win.show()
        getattr(parent_win, 'button_see_proj_settings').setChecked(False)
    
    def show_organ_settings(self, parent_win):
        if self.organ_settings_win == None: 
            self.organ_settings_win = OrganSettings(proj = self.proj, organ = self.organ, controller=self)
        self.organ_settings_win.show()
        getattr(parent_win, 'button_see_organ_settings').setChecked(False)
    
    def show_about(self): 
        if self.about_window == None:
            self.about_window = AboutScreen() 
        self.about_window.show()

    def show_template(self): 
        selection = self.new_proj_win_from_temp.cB_temp_names.currentText()
        if  selection != '--select--': 
            jsonDict_name = selection+'.json'
            json2open_dir = mH_config.path_templates / jsonDict_name
            if json2open_dir.is_file():
                with open(json2open_dir, "r") as read_file:
                    print(">> "+jsonDict_name+": Opening JSON encoded data")
                    load_dict = json.load(read_file)
                
            load_dict = make_Paths(load_dict)
            self.proj_temp = mHC.TempProj(proj_dict = load_dict)

            if self.proj_settings_win == None:
                self.proj_settings_win = ProjSettings(proj = self.proj_temp, controller=self, template=selection) 
            self.proj_settings_win.show()
            getattr(self.new_proj_win_from_temp, 'btn_view_temp').setChecked(False)
        else: 
            self.new_proj_win_from_temp.win_msg('*Select a valid template to be able to view its settings...', self.new_proj_win_from_temp.btn_view_temp)
            return

    #Inititalise windows
    def init_load_proj(self): 
        #Connect buttons
        # -Go Back
        self.load_proj_win.button_go_back.clicked.connect(lambda: self.show_welcome())
        self.load_proj_win.button_go_back_comb.clicked.connect(lambda: self.show_welcome())
        # -Browse project
        self.load_proj_win.button_browse_proj.clicked.connect(lambda: self.load_proj(parent_win='load_proj_win'))
        # -Show New Organ Window 
        self.load_proj_win.button_add_organ.clicked.connect(lambda: self.show_new_organ(parent_win='load_proj_win'))
        # -Go to main_window
        self.load_proj_win.go_to_main_window.clicked.connect(lambda: self.show_main_window(parent_win='load_proj_win'))
        # -Go to analysis window
        self.load_proj_win.go_to_analysis_window.clicked.connect(lambda: self.show_analysis_window(parent_win='load_proj_win', single_proj=True))
        # -See proj settings
        self.load_proj_win.button_see_proj_settings.clicked.connect(lambda: self.show_proj_settings(parent_win=self.load_proj_win))

    def init_load_multi_proj(self): 
        #Connect buttons
        # -Go Back
        self.load_multi_proj_win.button_go_back.clicked.connect(lambda: self.show_welcome())
        # -Browse project
        self.load_multi_proj_win.button_browse_proj.clicked.connect(lambda: self.load_proj(parent_win='load_multi_proj_win'))
       # -Go to analysis window
        self.load_multi_proj_win.go_to_analysis_window.clicked.connect(lambda: self.show_analysis_window(parent_win='load_multi_proj_win', single_proj=False))
        # -See proj settings
        self.load_multi_proj_win.button_see_proj_settings.clicked.connect(lambda: self.show_proj_settings(parent_win=self.load_proj_win))

    def init_multi_proj(self): 
        #Connect buttons
        # -Go Back
        self.load_multi_proj_win.button_go_back.clicked.connect(lambda: self.show_welcome())
        # -Browse project
        self.load_multi_proj_win.button_browse_proj.clicked.connect(lambda: self.load_proj(parent_win='load_multi_proj_win'))
        # -See proj settings
        self.load_multi_proj_win.button_see_proj_settings.clicked.connect(lambda: self.show_proj_settings(parent_win=self.load_multi_proj_win))

    def init_create_new_proj(self): 
        #Connect Buttons
        # -Go Back 
        self.new_proj_win.button_go_back.clicked.connect(lambda: self.show_welcome())
        # -Set Measurement Parameters 
        self.new_proj_win.set_meas_param_all.clicked.connect(lambda: self.show_meas_param())
        # -Create New Project 
        self.new_proj_win.button_new_proj.clicked.connect(lambda: self.new_proj())
        # -Show New Organ Window 
        self.new_proj_win.button_add_organ.clicked.connect(lambda: self.show_new_organ(parent_win='new_proj_win'))

    def init_new_proj_from_temp(self): 
        # -Go Back
        self.new_proj_win_from_temp.button_go_back.clicked.connect(lambda: self.show_welcome())
        self.new_proj_win_from_temp.btn_view_temp.clicked.connect(lambda: self.show_template())
        self.new_proj_win_from_temp.button_select_temp.clicked.connect(lambda: self.select_template())

    def init_new_organ_win(self, parent_win=None): 
        #Connect Buttons
        # -Go Back 
        self.new_organ_win.button_go_back.clicked.connect(lambda: self.show_parent(win='new_organ_win', parent=parent_win))
        # - Create New Organ
        self.new_organ_win.button_create_new_organ.clicked.connect(lambda: self.new_organ())
        # -Go to main_window
        self.new_organ_win.go_to_main_window.clicked.connect(lambda: self.show_main_window(parent_win='new_organ_win'))
        # -See proj settings
        self.new_organ_win.button_see_proj_settings.clicked.connect(lambda: self.show_proj_settings(parent_win=self.new_organ_win))

    def init_main_win(self): 
        
        self.main_win.button_see_proj_settings.clicked.connect(lambda: self.show_proj_settings(parent_win=self.main_win))
        self.main_win.button_see_organ_settings.clicked.connect(lambda: self.show_organ_settings(parent_win=self.main_win))
        self.init_segmentation_tab()
        self.init_morphoHeart_tab()
        self.init_morphoCell_tab()

        #Action buttons
        self.main_win.actionOpen_a_new_Project_and_Organ.triggered.connect(self.open_new_organ_and_project)
        self.main_win.actionCreate_new_Project.triggered.connect(self.create_new_project)
        self.main_win.actionOpen_another_organ_from_current_project.triggered.connect(self.open_another_organ_same_project)
        self.main_win.actionCreate_new_organ_within_the_current_project.triggered.connect(self.create_new_organ_same_project)
        self.main_win.actionAbout_morphoHeart.triggered.connect(self.show_about)

    def init_segmentation_tab(self): 
        #Segmentation Tab
        self.main_win.all_closed.clicked.connect(lambda: self.show_load_closed_stacks())

        #MASKING
        self.main_win.mask_ch1_play.clicked.connect(lambda: self.mask_ch('ch1'))
        self.main_win.mask_ch2_play.clicked.connect(lambda: self.mask_ch('ch2'))
        self.main_win.mask_ch3_play.clicked.connect(lambda: self.mask_ch('ch3'))
        self.main_win.mask_ch4_play.clicked.connect(lambda: self.mask_ch('ch4'))

        #MASKING NPY
        self.main_win.npy_mask_ch1_play.clicked.connect(lambda: self.mask_npy('ch1'))
        self.main_win.npy_mask_ch2_play.clicked.connect(lambda: self.mask_npy('ch2'))
        self.main_win.npy_mask_ch3_play.clicked.connect(lambda: self.mask_npy('ch3'))
        self.main_win.npy_mask_ch4_play.clicked.connect(lambda: self.mask_npy('ch4'))

        #AUTOMATICALLY CLOSE CONTOURS
        self.main_win.autom_close_ch1_play.clicked.connect(lambda: self.autom_close_contours('ch1'))
        self.main_win.autom_close_ch2_play.clicked.connect(lambda: self.autom_close_contours('ch2'))
        self.main_win.autom_close_ch3_play.clicked.connect(lambda: self.autom_close_contours('ch3'))
        self.main_win.autom_close_ch4_play.clicked.connect(lambda: self.autom_close_contours('ch4'))

        #MANUALLY CLOSE CONTOURS
        self.main_win.manual_close_ch1_play.clicked.connect(lambda: self.manual_close_contours('ch1'))
        self.main_win.manual_close_ch2_play.clicked.connect(lambda: self.manual_close_contours('ch2'))
        self.main_win.manual_close_ch3_play.clicked.connect(lambda: self.manual_close_contours('ch3'))
        self.main_win.manual_close_ch4_play.clicked.connect(lambda: self.manual_close_contours('ch4'))

        # Run tuples
        self.main_win.slc_tuple_ch1_play.clicked.connect(lambda: mA.close_slcs_tuple(controller = self, ch_name='ch1'))
        self.main_win.slc_tuple_ch2_play.clicked.connect(lambda: mA.close_slcs_tuple(controller = self, ch_name='ch2'))
        self.main_win.slc_tuple_ch3_play.clicked.connect(lambda: mA.close_slcs_tuple(controller = self, ch_name='ch3'))
        self.main_win.slc_tuple_ch4_play.clicked.connect(lambda: mA.close_slcs_tuple(controller = self, ch_name='ch4'))

        # Slices and tuples
        self.main_win.next_slice_ch1.clicked.connect(lambda: mA.next_prev_slice_in_tuple(next=True, controller = self, ch_name ='ch1'))
        self.main_win.next_slice_ch2.clicked.connect(lambda: mA.next_prev_slice_in_tuple(next=True, controller = self, ch_name ='ch2'))
        self.main_win.next_slice_ch3.clicked.connect(lambda: mA.next_prev_slice_in_tuple(next=True, controller = self, ch_name ='ch3'))
        self.main_win.next_slice_ch4.clicked.connect(lambda: mA.next_prev_slice_in_tuple(next=True, controller = self, ch_name ='ch4'))

        self.main_win.prev_slice_ch1.clicked.connect(lambda: mA.next_prev_slice_in_tuple(next=False, controller = self, ch_name ='ch1'))
        self.main_win.prev_slice_ch2.clicked.connect(lambda: mA.next_prev_slice_in_tuple(next=False, controller = self, ch_name ='ch2'))
        self.main_win.prev_slice_ch3.clicked.connect(lambda: mA.next_prev_slice_in_tuple(next=False, controller = self, ch_name ='ch3'))
        self.main_win.prev_slice_ch4.clicked.connect(lambda: mA.next_prev_slice_in_tuple(next=False, controller = self, ch_name ='ch4'))

        self.main_win.next_tuple_ch1.clicked.connect(lambda: mA.next_prev_tuple_to_manually_close(next=True, controller=self, ch_name='ch1'))
        self.main_win.next_tuple_ch2.clicked.connect(lambda: mA.next_prev_tuple_to_manually_close(next=True, controller=self, ch_name='ch2'))
        self.main_win.next_tuple_ch3.clicked.connect(lambda: mA.next_prev_tuple_to_manually_close(next=True, controller=self, ch_name='ch3'))
        self.main_win.next_tuple_ch4.clicked.connect(lambda: mA.next_prev_tuple_to_manually_close(next=True, controller=self, ch_name='ch4'))

        self.main_win.prev_tuple_ch1.clicked.connect(lambda: mA.next_prev_tuple_to_manually_close(next=False, controller=self, ch_name='ch1'))
        self.main_win.prev_tuple_ch2.clicked.connect(lambda: mA.next_prev_tuple_to_manually_close(next=False, controller=self, ch_name='ch2'))
        self.main_win.prev_tuple_ch3.clicked.connect(lambda: mA.next_prev_tuple_to_manually_close(next=False, controller=self, ch_name='ch3'))
        self.main_win.prev_tuple_ch4.clicked.connect(lambda: mA.next_prev_tuple_to_manually_close(next=False, controller=self, ch_name='ch4'))

        #SELECTING CONTOURS
        self.main_win.select_contours_ch1_play.clicked.connect(lambda: self.select_contours('ch1'))
        self.main_win.select_contours_ch2_play.clicked.connect(lambda: self.select_contours('ch2'))
        self.main_win.select_contours_ch3_play.clicked.connect(lambda: self.select_contours('ch3'))
        self.main_win.select_contours_ch4_play.clicked.connect(lambda: self.select_contours('ch4'))

        #Run tuples
        self.main_win.slc_tuple_select_ch1_play.clicked.connect(lambda: mA.select_slcs_tuple(controller = self, ch_name='ch1'))
        self.main_win.slc_tuple_select_ch2_play.clicked.connect(lambda: mA.select_slcs_tuple(controller = self, ch_name='ch2'))
        self.main_win.slc_tuple_select_ch3_play.clicked.connect(lambda: mA.select_slcs_tuple(controller = self, ch_name='ch3'))
        self.main_win.slc_tuple_select_ch4_play.clicked.connect(lambda: mA.select_slcs_tuple(controller = self, ch_name='ch4'))
        
        #Next tuple
        self.main_win.next_group_ch1.clicked.connect(lambda: mA.next_tuple_select(controller=self, ch_name='ch1'))
        self.main_win.next_group_ch2.clicked.connect(lambda: mA.next_tuple_select(controller=self, ch_name='ch2'))
        self.main_win.next_group_ch3.clicked.connect(lambda: mA.next_tuple_select(controller=self, ch_name='ch3'))
        self.main_win.next_group_ch4.clicked.connect(lambda: mA.next_tuple_select(controller=self, ch_name='ch4'))

        self.main_win.select_manually_slc_ch1_play.clicked.connect(lambda: mA.modify_selected_contours(controller = self, ch_name='ch1'))
        self.main_win.select_manually_slc_ch2_play.clicked.connect(lambda: mA.modify_selected_contours(controller = self, ch_name='ch2'))
        self.main_win.select_manually_slc_ch3_play.clicked.connect(lambda: mA.modify_selected_contours(controller = self, ch_name='ch3'))
        self.main_win.select_manually_slc_ch4_play.clicked.connect(lambda: mA.modify_selected_contours(controller = self, ch_name='ch4'))

        #Run slc
        self.main_win.slc_manually_select_ch1_play.clicked.connect(lambda: mA.select_slc(controller=self, ch_name='ch1'))
        self.main_win.slc_manually_select_ch2_play.clicked.connect(lambda: mA.select_slc(controller=self, ch_name='ch2'))
        self.main_win.slc_manually_select_ch3_play.clicked.connect(lambda: mA.select_slc(controller=self, ch_name='ch3'))
        self.main_win.slc_manually_select_ch4_play.clicked.connect(lambda: mA.select_slc(controller=self, ch_name='ch4'))

        #Next slc 
        self.main_win.next_slc_select_ch1.clicked.connect(lambda: mA.next_slc_select(controller=self, ch_name='ch1'))
        self.main_win.next_slc_select_ch2.clicked.connect(lambda: mA.next_slc_select(controller=self, ch_name='ch2'))
        self.main_win.next_slc_select_ch3.clicked.connect(lambda: mA.next_slc_select(controller=self, ch_name='ch3'))
        self.main_win.next_slc_select_ch4.clicked.connect(lambda: mA.next_slc_select(controller=self, ch_name='ch4'))

    def init_morphoHeart_tab(self): 
        #Process and Analyse Tab
        self.main_win.keeplargest_play.clicked.connect(lambda: self.run_keeplargest())
        self.main_win.cleanup_play.clicked.connect(lambda: self.run_cleanup())
        self.main_win.trimming_play.clicked.connect(lambda: self.run_trimming())
        self.main_win.orientation_play.clicked.connect(lambda: self.run_axis_orientation())
        self.main_win.chNS_play.clicked.connect(lambda: self.run_chNS())
        self.main_win.measure_wholeAll_play.clicked.connect(lambda: self.run_measure_whole())

        self.main_win.centreline_clean_play.clicked.connect(lambda: self.run_centreline_clean())
        self.main_win.centreline_ML_play.clicked.connect(lambda: self.run_centreline_ML())
        self.main_win.centreline_vmtk_play.clicked.connect(lambda: self.run_centreline_vmtk())
        self.main_win.centreline_select.clicked.connect(lambda: self.run_centreline_select())
        # self.main_win.centreline_play.clicked.connect(lambda: self.run_centreline())

        #HEATMAPS
        # self.main_win.heatmaps3D_play.clicked.connect(lambda: self.run_heatmaps3D())
        #Heatmap Indiv Play buttons
        self.main_win.hm_play1.clicked.connect(lambda: self.run_heatmaps3D(btn=1))
        self.main_win.hm_play2.clicked.connect(lambda: self.run_heatmaps3D(btn=2))
        self.main_win.hm_play3.clicked.connect(lambda: self.run_heatmaps3D(btn=3))
        self.main_win.hm_play4.clicked.connect(lambda: self.run_heatmaps3D(btn=4))
        self.main_win.hm_play5.clicked.connect(lambda: self.run_heatmaps3D(btn=5))
        self.main_win.hm_play6.clicked.connect(lambda: self.run_heatmaps3D(btn=6))
        self.main_win.hm_play7.clicked.connect(lambda: self.run_heatmaps3D(btn=7))
        self.main_win.hm_play8.clicked.connect(lambda: self.run_heatmaps3D(btn=8))
        self.main_win.hm_play9.clicked.connect(lambda: self.run_heatmaps3D(btn=9))
        self.main_win.hm_play10.clicked.connect(lambda: self.run_heatmaps3D(btn=10))
        self.main_win.hm_play11.clicked.connect(lambda: self.run_heatmaps3D(btn=11))
        self.main_win.hm_play12.clicked.connect(lambda: self.run_heatmaps3D(btn=12))

        self.main_win.hm2d_play1.clicked.connect(lambda: self.run_heatmaps2D(btn=1))
        self.main_win.hm2d_play2.clicked.connect(lambda: self.run_heatmaps2D(btn=2))
        self.main_win.hm2d_play3.clicked.connect(lambda: self.run_heatmaps2D(btn=3))
        self.main_win.hm2d_play4.clicked.connect(lambda: self.run_heatmaps2D(btn=4))
        self.main_win.hm2d_play5.clicked.connect(lambda: self.run_heatmaps2D(btn=5))
        self.main_win.hm2d_play6.clicked.connect(lambda: self.run_heatmaps2D(btn=6))
        self.main_win.hm2d_play7.clicked.connect(lambda: self.run_heatmaps2D(btn=7))
        self.main_win.hm2d_play8.clicked.connect(lambda: self.run_heatmaps2D(btn=8))
        self.main_win.hm2d_play9.clicked.connect(lambda: self.run_heatmaps2D(btn=9))
        self.main_win.hm2d_play10.clicked.connect(lambda: self.run_heatmaps2D(btn=10))
        self.main_win.hm2d_play11.clicked.connect(lambda: self.run_heatmaps2D(btn=11))
        self.main_win.hm2d_play12.clicked.connect(lambda: self.run_heatmaps2D(btn=12))

        # SEGMENTS
        # self.main_win.segments_play.clicked.connect(lambda: self.run_segments())
        #Segments Indiv Play buttons
        #Cut 1
        self.main_win.cut1_play_segm1.clicked.connect(lambda: self.run_segments(btn='Cut1_1'))
        self.main_win.cut1_play_segm2.clicked.connect(lambda: self.run_segments(btn='Cut1_2'))
        self.main_win.cut1_play_segm3.clicked.connect(lambda: self.run_segments(btn='Cut1_3'))
        self.main_win.cut1_play_segm4.clicked.connect(lambda: self.run_segments(btn='Cut1_4'))
        self.main_win.cut1_play_segm5.clicked.connect(lambda: self.run_segments(btn='Cut1_5'))
        self.main_win.cut1_play_segm6.clicked.connect(lambda: self.run_segments(btn='Cut1_6'))
        self.main_win.cut1_play_segm7.clicked.connect(lambda: self.run_segments(btn='Cut1_7'))
        self.main_win.cut1_play_segm8.clicked.connect(lambda: self.run_segments(btn='Cut1_8'))
        self.main_win.cut1_play_segm9.clicked.connect(lambda: self.run_segments(btn='Cut1_9'))
        self.main_win.cut1_play_segm10.clicked.connect(lambda: self.run_segments(btn='Cut1_10'))
        self.main_win.cut1_play_segm11.clicked.connect(lambda: self.run_segments(btn='Cut1_11'))
        self.main_win.cut1_play_segm12.clicked.connect(lambda: self.run_segments(btn='Cut1_12'))
        #Cut 2
        self.main_win.cut2_play_segm1.clicked.connect(lambda: self.run_segments(btn='Cut2_1'))
        self.main_win.cut2_play_segm2.clicked.connect(lambda: self.run_segments(btn='Cut2_2'))
        self.main_win.cut2_play_segm3.clicked.connect(lambda: self.run_segments(btn='Cut2_3'))
        self.main_win.cut2_play_segm4.clicked.connect(lambda: self.run_segments(btn='Cut2_4'))
        self.main_win.cut2_play_segm5.clicked.connect(lambda: self.run_segments(btn='Cut2_5'))
        self.main_win.cut2_play_segm6.clicked.connect(lambda: self.run_segments(btn='Cut2_6'))
        self.main_win.cut2_play_segm7.clicked.connect(lambda: self.run_segments(btn='Cut2_7'))
        self.main_win.cut2_play_segm8.clicked.connect(lambda: self.run_segments(btn='Cut2_8'))
        self.main_win.cut2_play_segm9.clicked.connect(lambda: self.run_segments(btn='Cut2_9'))
        self.main_win.cut2_play_segm10.clicked.connect(lambda: self.run_segments(btn='Cut2_10'))
        self.main_win.cut2_play_segm11.clicked.connect(lambda: self.run_segments(btn='Cut2_11'))
        self.main_win.cut2_play_segm12.clicked.connect(lambda: self.run_segments(btn='Cut2_12'))

        #Angles segments
        self.main_win.play_angle_dir1_cut1.clicked.connect(lambda: self.run_angles(btn= 'dir1_cut1'))
        self.main_win.play_angle_dir2_cut1.clicked.connect(lambda: self.run_angles(btn= 'dir2_cut1'))
        self.main_win.play_angle_dir3_cut1.clicked.connect(lambda: self.run_angles(btn= 'dir3_cut1'))

        self.main_win.play_angle_dir1_cut2.clicked.connect(lambda: self.run_angles(btn= 'dir1_cut2'))
        self.main_win.play_angle_dir2_cut2.clicked.connect(lambda: self.run_angles(btn= 'dir2_cut2'))
        self.main_win.play_angle_dir3_cut2.clicked.connect(lambda: self.run_angles(btn= 'dir3_cut2'))

        # Ellipsoids
        self.main_win.ellipsoids_play.clicked.connect(lambda: self.run_ellipsoids())

        # SECTIONS
        # self.main_win.sections_play.clicked.connect(lambda: self.run_sections())
        #Cut 1
        self.main_win.cut1_play_sect1.clicked.connect(lambda: self.run_sections(btn='Cut1_1'))
        self.main_win.cut1_play_sect2.clicked.connect(lambda: self.run_sections(btn='Cut1_2'))
        self.main_win.cut1_play_sect3.clicked.connect(lambda: self.run_sections(btn='Cut1_3'))
        self.main_win.cut1_play_sect4.clicked.connect(lambda: self.run_sections(btn='Cut1_4'))
        self.main_win.cut1_play_sect5.clicked.connect(lambda: self.run_sections(btn='Cut1_5'))
        self.main_win.cut1_play_sect6.clicked.connect(lambda: self.run_sections(btn='Cut1_6'))
        self.main_win.cut1_play_sect7.clicked.connect(lambda: self.run_sections(btn='Cut1_7'))
        self.main_win.cut1_play_sect8.clicked.connect(lambda: self.run_sections(btn='Cut1_8'))
        self.main_win.cut1_play_sect9.clicked.connect(lambda: self.run_sections(btn='Cut1_9'))
        self.main_win.cut1_play_sect10.clicked.connect(lambda: self.run_sections(btn='Cut1_10'))
        self.main_win.cut1_play_sect11.clicked.connect(lambda: self.run_sections(btn='Cut1_11'))
        self.main_win.cut1_play_sect12.clicked.connect(lambda: self.run_sections(btn='Cut1_12'))
        #Cut 2
        self.main_win.cut2_play_sect1.clicked.connect(lambda: self.run_sections(btn='Cut2_1'))
        self.main_win.cut2_play_sect2.clicked.connect(lambda: self.run_sections(btn='Cut2_2'))
        self.main_win.cut2_play_sect3.clicked.connect(lambda: self.run_sections(btn='Cut2_3'))
        self.main_win.cut2_play_sect4.clicked.connect(lambda: self.run_sections(btn='Cut2_4'))
        self.main_win.cut2_play_sect5.clicked.connect(lambda: self.run_sections(btn='Cut2_5'))
        self.main_win.cut2_play_sect6.clicked.connect(lambda: self.run_sections(btn='Cut2_6'))
        self.main_win.cut2_play_sect7.clicked.connect(lambda: self.run_sections(btn='Cut2_7'))
        self.main_win.cut2_play_sect8.clicked.connect(lambda: self.run_sections(btn='Cut2_8'))
        self.main_win.cut2_play_sect9.clicked.connect(lambda: self.run_sections(btn='Cut2_9'))
        self.main_win.cut2_play_sect10.clicked.connect(lambda: self.run_sections(btn='Cut2_10'))
        self.main_win.cut2_play_sect11.clicked.connect(lambda: self.run_sections(btn='Cut2_11'))
        self.main_win.cut2_play_sect12.clicked.connect(lambda: self.run_sections(btn='Cut2_12'))

        # SEGMENTS-SECTIONS
        # self.main_win.segm_sect_play.clicked.connect(lambda: self.run_segm_sect())
        #sCut1 - Cut1
        self.main_win.scut1_cut1_play_sect1.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut1_1'))
        self.main_win.scut1_cut1_play_sect2.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut1_2'))
        self.main_win.scut1_cut1_play_sect3.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut1_3'))
        self.main_win.scut1_cut1_play_sect4.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut1_4'))
        self.main_win.scut1_cut1_play_sect5.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut1_5'))
        self.main_win.scut1_cut1_play_sect6.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut1_6'))
        self.main_win.scut1_cut1_play_sect7.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut1_7'))
        self.main_win.scut1_cut1_play_sect8.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut1_8'))
        self.main_win.scut1_cut1_play_sect9.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut1_9'))
        self.main_win.scut1_cut1_play_sect10.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut1_10'))
        self.main_win.scut1_cut1_play_sect11.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut1_11'))
        self.main_win.scut1_cut1_play_sect12.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut1_12'))

        #sCut1 - Cut2
        self.main_win.scut1_cut2_play_sect1.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut2_1'))
        self.main_win.scut1_cut2_play_sect2.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut2_2'))
        self.main_win.scut1_cut2_play_sect3.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut2_3'))
        self.main_win.scut1_cut2_play_sect4.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut2_4'))
        self.main_win.scut1_cut2_play_sect5.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut2_5'))
        self.main_win.scut1_cut2_play_sect6.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut2_6'))
        self.main_win.scut1_cut2_play_sect7.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut2_7'))
        self.main_win.scut1_cut2_play_sect8.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut2_8'))
        self.main_win.scut1_cut2_play_sect9.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut2_9'))
        self.main_win.scut1_cut2_play_sect10.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut2_10'))
        self.main_win.scut1_cut2_play_sect11.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut2_11'))
        self.main_win.scut1_cut2_play_sect12.clicked.connect(lambda: self.run_segm_sect(btn='sCut1_o_Cut2_12'))

        #sCut2 - Cut1
        self.main_win.scut2_cut1_play_sect1.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut1_1'))
        self.main_win.scut2_cut1_play_sect2.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut1_2'))
        self.main_win.scut2_cut1_play_sect3.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut1_3'))
        self.main_win.scut2_cut1_play_sect4.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut1_4'))
        self.main_win.scut2_cut1_play_sect5.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut1_5'))
        self.main_win.scut2_cut1_play_sect6.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut1_6'))
        self.main_win.scut2_cut1_play_sect7.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut1_7'))
        self.main_win.scut2_cut1_play_sect8.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut1_8'))
        self.main_win.scut2_cut1_play_sect9.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut1_9'))
        self.main_win.scut2_cut1_play_sect10.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut1_10'))
        self.main_win.scut2_cut1_play_sect11.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut1_11'))
        self.main_win.scut2_cut1_play_sect12.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut1_12'))

        #sCut2 - Cut2
        self.main_win.scut2_cut2_play_sect1.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut2_1'))
        self.main_win.scut2_cut2_play_sect2.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut2_2'))
        self.main_win.scut2_cut2_play_sect3.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut2_3'))
        self.main_win.scut2_cut2_play_sect4.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut2_4'))
        self.main_win.scut2_cut2_play_sect5.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut2_5'))
        self.main_win.scut2_cut2_play_sect6.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut2_6'))
        self.main_win.scut2_cut2_play_sect7.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut2_7'))
        self.main_win.scut2_cut2_play_sect8.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut2_8'))
        self.main_win.scut2_cut2_play_sect9.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut2_9'))
        self.main_win.scut2_cut2_play_sect10.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut2_10'))
        self.main_win.scut2_cut2_play_sect11.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut2_11'))
        self.main_win.scut2_cut2_play_sect12.clicked.connect(lambda: self.run_segm_sect(btn='sCut2_o_Cut2_12'))

    def init_morphoCell_tab(self): 
        #Isosurface
        self.main_win.chB_play.clicked.connect(lambda: self.run_isosurface('chB'))
        self.main_win.chC_play.clicked.connect(lambda: self.run_isosurface('chC'))
        self.main_win.chD_play.clicked.connect(lambda: self.run_isosurface('chD'))

        self.main_win.remove_cells_play.clicked.connect(lambda: self.run_remove_cells())

        self.main_win.cut1_segm_cell_play.clicked.connect(lambda: self.run_segments_mC('Cut1'))
        self.main_win.cut2_segm_cell_play.clicked.connect(lambda: self.run_segments_mC('Cut2'))

        self.main_win.ind_segments_play.clicked.connect(lambda: self.run_IND_segm())

        self.main_win.zone1_cell_play.clicked.connect(lambda: self.run_zones('Zone1'))
        self.main_win.zone2_cell_play.clicked.connect(lambda: self.run_zones('Zone2'))
        self.main_win.zone3_cell_play.clicked.connect(lambda: self.run_zones('Zone3'))
                
    def init_multip_analysis_win(self): 

        self.multip_analysis_win.button_see_proj_settings.clicked.connect(lambda: self.show_proj_settings(parent_win=self.main_win))

        #Action buttons
        
    #Functions related to API  
    # Project Related  
    def set_proj_meas_param(self):
        if self.meas_param_win.validate_params() == True: 
            self.meas_param_win.get_final_parameters(self)  
            for stype in ['segm', 'sect']: 
                ck_type = getattr(self.new_proj_win, 'tick_'+stype)
                if ck_type.isChecked():
                    getattr(self.new_proj_win, 'button_set_'+stype).setEnabled(True)
        else: 
            return 
        
    def new_proj(self):
        if self.new_proj_win.validate_set_all():
            self.new_proj_win.win_msg("Creating and saving new project...")
            temp_dir = self.new_proj_win.check_template()
            if temp_dir != False: 
                self.new_proj_win.button_new_proj.setChecked(True)

                proj_dict = {'name': self.new_proj_win.lineEdit_proj_name.text().strip(), 
                            'notes' : self.new_proj_win.textEdit_ref_notes.toPlainText().strip(),
                            'date' : str(self.new_proj_win.dateEdit.date().toPyDate()),
                            'analysis' : self.new_proj_win.checked_analysis, 
                            'dir_proj' : self.new_proj_win.proj_dir, 
                            'heart_default': self.new_proj_win.heart_analysis.isChecked()}
                
                self.proj = mHC.Project(proj_dict, new=True)
                self.new_proj_win.mH_settings['chs_all'] = self.ch_all
                self.new_proj_win.mH_settings['params'] = self.mH_params

                self.proj.set_settings(settings={'mH': {'settings':self.new_proj_win.mH_settings, 
                                                        'params': self.new_proj_win.mH_user_params},
                                                'mC': {'settings': self.new_proj_win.mC_settings,
                                                    'params': self.new_proj_win.mC_user_params}})
                
                self.proj.set_workflow()
                self.proj.create_proj_dir()
        
                self.proj.save_project(temp_dir = temp_dir)
                self.new_proj_win.button_add_organ.setEnabled(True)
                print('\n>>> New Project: ',self.proj.__dict__.keys())
                self.new_proj_win.win_msg("New project '"+self.new_proj_win.lineEdit_proj_name.text()+"' has been created and saved! Continue by creating an organ as part of this project. ")
        
    def load_proj(self, parent_win):
        cwd = Path().absolute().home()
        path_folder = QFileDialog.getExistingDirectory(self.load_proj_win, directory=str(cwd), 
                                                       caption="Select the Project's directory")
        proj_name = str(Path(path_folder).name)[2:] #Removing the R_
        proj_name_us = proj_name.replace(' ', '_')
        json_name = 'mH_'+proj_name_us+'_project.json'
        proj_settings_path = Path(path_folder) / 'settings' / json_name
        win = getattr(self, parent_win)
        if proj_settings_path.is_file(): 
            proj_dict = {'name': proj_name, 
                         'dir': path_folder}
            self.proj = mHC.Project(proj_dict, new=False)
            print('Loaded project:',self.proj.__dict__.keys())
            print('Project[organs]:',self.proj.organs)
            win.proj = self.proj
            #Fill window with project info
            win.fill_proj_info(proj = self.proj)
            win.init_tables(load=True)
        else: 
            win.button_browse_proj.setChecked(False)
            win.init_tables(load=False)
            win.win_msg('*There is no settings file for a project within the selected directory. Please select a new directory.')

    def select_template(self): 
        selection = self.new_proj_win_from_temp.cB_temp_names.currentText()
        if selection != '--select--': 
            if hasattr(self, 'proj_temp'): 
                pass
            else: 
                jsonDict_name = selection+'.json'
                json2open_dir = mH_config.path_templates / jsonDict_name
                if json2open_dir.is_file():
                    with open(json2open_dir, "r") as read_file:
                        print(">> "+jsonDict_name+": Opening JSON encoded data")
                        load_dict = json.load(read_file)
                    
                load_dict = make_Paths(load_dict)
                self.proj_temp = mHC.TempProj(proj_dict = load_dict)
            
            #Open new window  filled with information of template
            self.new_proj_win_from_temp.close()
            #Create new proj window and show
            if self.new_proj_win == None: 
                self.new_proj_win = CreateNewProj(controller=self, template=True)
                self.init_create_new_proj()
            self.new_proj_win.init_template(proj = self.proj_temp)
            self.new_proj_win.show()
        else: 
            self.new_proj_win_from_temp.win_msg('*Select a valid template to be use for the new project...', self.new_proj_win_from_temp.button_select_temp)
            return

    #Organ related
    def new_organ(self): 
        if self.new_organ_win.validate_organ(self.proj): 
            self.new_organ_win.win_msg('Creating and saving new organ...')
            if self.new_organ_win.check_selection(self.proj):
                if self.new_organ_win.check_shapes(self.proj): 
                    self.new_organ_win.button_create_new_organ.setChecked(True)
                    self.new_organ_win.win_msg('Creating organ "'+self.new_organ_win.lineEdit_organ_name.text()+'"')
                    name = self.new_organ_win.lineEdit_organ_name.text().strip()
                    notes = self.new_organ_win.textEdit_ref_notes.toPlainText().strip()
                    strain = self.new_organ_win.cB_strain.currentText()
                    stage = self.new_organ_win.cB_stage.currentText()
                    genotype = self.new_organ_win.cB_genotype.currentText()
                    manipulation = self.new_organ_win.cB_manipulation.currentText()
                    im_or = self.new_organ_win.cB_stack_orient.currentText()
                    custom_angle = self.new_organ_win.cust_angle.text()
                    res_units = self.new_organ_win.resolution
                    resolution = [res_units[axis]['scaling'] for axis in ['x','y','z']]
                    units = [res_units[axis]['units'] for axis in ['x','y','z']]
                    date = str(self.new_organ_win.dateEdit.date().toPyDate())

                    organ_settings = {'project': {'user': self.proj.user_projName,
                                                'mH': self.proj.mH_projName,
                                                'dict_dir_info': self.proj.dir_info},
                                        'user_organName': name,
                                        'user_organNotes': notes,
                                        'im_orientation': im_or,
                                        'custom_angle': custom_angle,
                                        'resolution': resolution,
                                        'im_res_units': units,
                                        'stage': stage, 
                                        'strain': strain, 
                                        'genotype': genotype,
                                        'manipulation': manipulation, 
                                        'date_created': date
                                            }
                    organ_dict = {'settings': organ_settings, 
                                'img_dirs': self.new_organ_win.img_dirs}

                    self.organ = mHC.Organ(project=self.proj, organ_dict=organ_dict, new = True)
                    self.new_organ_win.lab_filled_organ_dir.setText(str(self.organ.dir_res()))

                    self.proj.add_organ(self.organ)
                    self.organ.save_organ()
                    self.new_organ_win.win_msg('New organ "'+name+'" has been created as part of "'+self.proj.user_projName+'" project.')
                    self.new_organ_win.go_to_main_window.setEnabled(True)
                else: 
                    self.new_organ_win.button_create_new_organ.setChecked(False)
                    return
            else: 
                self.new_organ_win.button_create_new_organ.setChecked(False)
                return 
        else:
            self.new_organ_win.button_create_new_organ.setChecked(False)
            return 
    
    def load_organ(self, proj, organ_to_load, single_organ=True):
        loaded_organ = proj.load_organ(organ_to_load = organ_to_load)
        if not hasattr(loaded_organ, 'obj_temp'):
                loaded_organ.obj_temp = {}
        if single_organ: 
            self.organ = loaded_organ
            print('-------------Loaded Organ:-------------')
            print('organ.workflow: ', self.organ.workflow)
            print('self.organ.obj_temp: ',self.organ.obj_temp)
            print('self.organ.mH_settings: ',self.organ.mH_settings)
            print('self.organ.mH_settings[wf_info]: ',self.organ.mH_settings['wf_info'])
            print('self.organ.submeshes: ', self.organ.submeshes)
        else: 
            return loaded_organ

    def load_multip_proj_and_organs(self, proj_org, single_proj):
        
        #Transform the list of dictionaries into a dataframe
        df_pando = pd.DataFrame(proj_org) 
        #Get all proj and save them in dict of proj
        proj_info = df_pando[['user_projName', 'proj_path']]
        unique_proj_info = proj_info.drop_duplicates()
        dict_projs = {}
        for nn, row in unique_proj_info.iterrows(): 
            proj_name = row['user_projName']
            proj_path = row['proj_path']
            proj_dict = {'name': proj_name, 
                         'dir': str(proj_path)}
            proj = mHC.Project(proj_dict, new=False)
            dict_projs[nn] = {'proj_name': proj_name, 
                             'proj_path': proj_path, 
                             'proj': proj}
        # print('dict_projs:', dict_projs)

        dict_organs = {}
        proj_num = []
        for index, row in df_pando.iterrows():
            # print(index, row)
            #Get values from organ
            org_proj_name = row['user_projName']
            org_proj_path = row['proj_path']
            for nk in dict_projs.keys(): 
                pj_proj_name = dict_projs[nk]['proj_name']
                pj_proj_path = dict_projs[nk]['proj_path']
                if org_proj_name == pj_proj_name and str(org_proj_path) == str(pj_proj_path):
                    proj_num.append(nk)
                    break
            org_name = row['user_organName']
            organ = self.load_organ(proj = dict_projs[nk]['proj'], organ_to_load = org_name, single_organ=False)
            dict_organs[index] = {'organ_name': org_name, 
                                  'organ': organ}
        df_pando['proj_num'] = proj_num

        return df_pando, dict_projs, dict_organs

    #Channel segmentation related
    def mask_ch(self, ch_name): 
        mA.mask_channel(controller=self, ch_name=ch_name)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def mask_npy(self, ch_name): 
        mA.mask_npy_channel(controller=self, ch_name=ch_name)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def autom_close_contours(self, ch_name):
        mA.autom_close_contours(controller=self, ch_name=ch_name)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def manual_close_contours(self, ch_name):
        mA.manual_close_contours(controller=self, ch_name=ch_name)
    
    def select_contours(self, ch_name): 
        mA.select_contours(controller=self, ch_name=ch_name)

    #Analysis related
    def run_keeplargest(self):
        mA.run_keeplargest(controller=self)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_cleanup(self):
        mA.run_cleanup(controller=self)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)
    
    def run_trimming(self):
        mA.run_trimming(controller=self)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_axis_orientation(self):
        mA.run_axis_orientation(controller=self)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_chNS(self):
        mA.run_chNS(controller=self)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_measure_whole(self): 
        mA.run_measure(controller=self)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_centreline_clean(self):
        mA.run_centreline_clean(controller=self)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_centreline_ML(self):
        mA.run_centreline_ML(controller=self)

    def run_centreline_vmtk(self): 
        mA.run_centreline_vmtk(controller=self)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)
    
    def run_centreline_select(self): 
        mA.run_centreline_select(controller=self)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    # def run_centreline(self):
    #     mA.run_centreline(controller=self)

    def run_heatmaps3D(self, btn=None):
        mA.run_heatmaps3D(controller=self, btn = btn)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_heatmaps2D(self, btn=None):
        mA.run_heatmaps2D(controller=self, btn = btn)
        # if not mH_config.dev:
        #     self.main_win.save_project_and_organ_pressed(alert_on = False)

    # def run_heatmaps2D(self):
    #     mA.run_heatmaps2D(controller=self)

    def run_segments(self, btn=None):
        mA.run_segments(controller=self, btn=btn)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_angles(self, btn=None): 
        done = mA.run_angles(controller=self, btn=btn)
        if not mH_config.dev and done:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_ellipsoids(self): 
        mA.run_ellipsoids(controller=self)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False) 
    
    def run_sections(self, btn=None):
        mA.run_sections(controller=self, btn=btn)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)
    
    def run_segm_sect(self, btn=None): 
        mA.run_segm_sect(controller=self, btn=btn)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_isosurface(self, btn): 
        mA.run_isosurface(controller=self, btn=btn)

    def run_remove_cells(self): 
        mA.run_remove_cells(controller=self)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_segments_mC(self, btn): 
        mA.run_segments_mC(controller=self, btn=btn)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_IND_segm(self): 
        mA.run_IND_segm(controller=self)
        if not mH_config.dev:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    def run_zones(self, zone):
        res = mA.run_zones(controller=self, zone=zone)
        if not mH_config.dev and res != None:
            self.main_win.save_project_and_organ_pressed(alert_on = False)

    #Actions Main Win
    def open_new_organ_and_project(self): 
        #Close welcome window
        self.load_proj_win = None
        self.main_win.close()
        if self.main_win.prompt.output in ['Discard', 'Save All']: 
            self.main_win = None
            self.load_proj_win = LoadProj() 
            self.init_load_proj()
            self.load_proj_win.button_go_back.clicked.connect(lambda: self.clear_win_show_welcome(parent = 'load_proj_win'))
            self.load_proj_win.show()

    def create_new_project(self): 
        self.new_proj_win = None
        self.main_win.close()
        if self.main_win.prompt.output in ['Discard', 'Save All']: 
            self.main_win = None
            self.new_proj_win = CreateNewProj(controller=self)
            self.init_create_new_proj()
            self.new_proj_win.button_go_back.clicked.connect(lambda: self.clear_win_show_welcome(parent = 'new_proj_win'))
            self.new_proj_win.show()

    def open_another_organ_same_project(self):
        self.main_win.close()
        if self.main_win.prompt.output in ['Discard', 'Save All']: 
            self.main_win = None
            self.load_proj_win.button_go_back.clicked.connect(lambda: self.clear_win_show_welcome(parent = 'load_proj_win'))
            self.load_proj_win.go_to_main_window.setChecked(False)
            self.load_proj_win.organ_selected = None
            self.load_proj_win.organs_2del = None
            self.load_proj_win.multi_organ_checkboxes = None
            self.load_proj_win.multi_organs_added = []
            self.load_proj_win.added_organs_checkboxes = None
            self.load_proj_win.load_proj_organs(proj = self.load_proj_win.proj)
            self.load_proj_win.show()
        
    def create_new_organ_same_project(self): 
        if hasattr(self, 'new_organ_win'):
            delattr(self, 'new_organ_win')
        self.main_win.close()
        if self.main_win.prompt.output in ['Discard', 'Save All']: 
            self.main_win = None
            self.new_organ_win = NewOrgan(proj = self.proj)
            self.init_new_organ_win()
            self.new_organ_win.button_go_back.clicked.connect(lambda: self.clear_win_show_welcome(parent = 'new_organ_win'))
            self.new_organ_win.show()
    
    def clear_win_show_welcome(self, parent):
        if parent != None: 
            getattr(self, parent).close()
        else: 
            for win in self.wins: 
                if hasattr(self, win): 
                    getattr(self, win).close()
        
        for win in self.wins: 
            setattr(self, win, None)

        self.show_welcome()

def main():
    app = QtWidgets.QApplication(sys.argv)
    screen = app.primaryScreen()
    #print('screen.name', screen.name())
    size = screen.size()
    #print('screen.size: ', size.width(), size.height())
    controller = Controller()
    controller.show_welcome()
    try: 
        sys.exit(app.exec())
    except: 
        print('Exiting')

if __name__ == '__main__':
    main()
