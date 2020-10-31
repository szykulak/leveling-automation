from abaqusGui import *
from abaqusConstants import ALL
import osutils, os


###########################################################################
# Class definition
###########################################################################

class LevellingGUI1_plugin(AFXForm):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, owner):
        
        # Construct the base class.
        #
        AFXForm.__init__(self, owner)
        self.radioButtonGroups = {}

        self.cmd = AFXGuiCommand(mode=self, method='levellingFunction',
            objectName='plugin', registerQuery=False)
        pickedDefault = ''
        self.initial_flatnessKw = AFXFloatKeyword(self.cmd, 'initial_flatness', True)
        self.plate_thicknessKw = AFXIntKeyword(self.cmd, 'plate_thickness', True)
        self.roll_diameterKw = AFXFloatKeyword(self.cmd, 'roll_diameter', True)
        self.distance_between_rollsKw = AFXFloatKeyword(self.cmd, 'distance_between_rolls', True)
        self.upper_rolls_countKw = AFXIntKeyword(self.cmd, 'upper_rolls_count', True)
        self.lower_rolls_countKw = AFXIntKeyword(self.cmd, 'lower_rolls_count', True)
        self.roll_velocityKw = AFXFloatKeyword(self.cmd, 'roll_velocity', True)
        self.friction_coefficient1Kw = AFXFloatKeyword(self.cmd, 'friction_coefficient1', True)
        self.friction_coefficient2Kw = AFXFloatKeyword(self.cmd, 'friction_coefficient2', True)
        self.mesh_densityKw = AFXFloatKeyword(self.cmd, 'mesh_density', True)
        self.calibration_inputKw = AFXFloatKeyword(self.cmd, 'calibration_input', True)
        self.calibration_outputKw = AFXFloatKeyword(self.cmd, 'calibration_output', True)
        self.material_nameKw = AFXStringKeyword(self.cmd, 'material_name', True)
	self.material_libraryKw = AFXStringKeyword(self.cmd, 'material_library', True, '')

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstDialog(self):

        import levellingGUI1DB
        return levellingGUI1DB.LevellingGUI1DB(self)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def doCustomChecks(self):

        # Try to set the appropriate radio button on. If the user did
        # not specify any buttons to be on, do nothing.
        #
        for kw1,kw2,d in self.radioButtonGroups.values():
            try:
                value = d[ kw1.getValue() ]
                kw2.setValue(value)
            except:
                pass
        return True

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def okToCancel(self):

        # No need to close the dialog when a file operation (such
        # as New or Open) or model change is executed.
        #
        return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Register the plug-in
#
thisPath = os.path.abspath(__file__)
thisDir = os.path.dirname(thisPath)

toolset = getAFXApp().getAFXMainWindow().getPluginToolset()
toolset.registerGuiMenuButton(
    buttonText='Plate levelling', 
    object=LevellingGUI1_plugin(toolset),
    messageId=AFXMode.ID_ACTIVATE,
    icon=None,
    kernelInitString='import plugin',
    applicableModules=ALL,
    version='N/A',
    author='N/A',
    description='N/A',
    helpUrl='N/A'
)
