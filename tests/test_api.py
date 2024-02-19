'''
morphoHeart - test_api.py

@author: Juliana Sanchez-Posada

'''
import pytest
from datetime import datetime


from src.modules.mH_classes import Project, Organ

@pytest.fixture(scope='class')
def fix_project():
    user_projName= 'My Project'
    proj_notes= 'This is my project'
    proj1 = Project(user_projName, proj_notes)
    return proj1 


class TestProject:
    def test_create_project(self):
        proj1 = Project()
        now_str = datetime.now().strftime('%Y%m%d')

        assert 'mH_Proj-' in proj1.mH_projName 
        assert  now_str in proj1.mH_projName

    # def test_create_gWF (self):
    #     user_projName= 'My Project'
    #     proj_notes= 'This is my project'
    #     proj1 = Project(user_projName, proj_notes)

    #     assert proj1.user_projName == user_projName.replace(' ', '_')
    #     assert proj1.mH_projName is not None
    #     assert isinstance(proj1.dict_projInfo, dict)  
    #     assert proj1.dict_projInfo['proj_notes'] == proj_notes
    #     assert isinstance(proj1.dict_Workflow, dict) 
    #     assert isinstance(proj1.dict_UserPipeline, dict) 
    #     assert isinstance(proj1.organ, Organ)

    # def test_create_folders(self):
    #     pass

@pytest.mark.usefixtures('fix_project')
class TestOrgan():
        organ1 = Organ(project=fix_project)
        assert isinstance(organ1, Organ)
        assert organ1.parent_project == fix_project
        
    def test_create_folders(self):
        organ1 = Organ(project=fix_project)
        assert isinstance(organ1, Organ)
        pass
    def test_create_organ(self):
