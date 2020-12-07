from abaqusConstants import *
from abaqusGui import *
from kernelAccess import mdb, session
import os

thisPath = os.path.abspath(__file__)
thisDir = os.path.dirname(thisPath)


###########################################################################
# Class definition
###########################################################################

class LevellingGUI1DB(AFXDataDialog):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form):

        # Construct the base class.
        #

        AFXDataDialog.__init__(self, form, 'Plate leveling model generator',
            self.OK|self.CANCEL, DIALOG_ACTIONS_SEPARATOR)
            

        okBtn = self.getActionButton(self.ID_CLICKED_OK)
        okBtn.setText('OK')
        
        VFrame_1 = FXVerticalFrame(p=self, opts=0, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        fileName = os.path.join(thisDir, 'final_reference.png')
        icon = afxCreatePNGIcon(fileName)
        FXLabel(p=VFrame_1, text='', ic=icon)
		
	'''
        VAligner_14 = AFXVerticalAligner(p=self, opts=0, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        AFXTextField(p=VAligner_14, ncols=12, labelText='Initial plate flatness - p [mm]:', tgt=form.initial_flatnessKw, sel=0)
        AFXTextField(p=VAligner_14, ncols=12, labelText='Plate thickness [mm]:', tgt=form.plate_thicknessKw, sel=0)
        AFXTextField(p=VAligner_14, ncols=12, labelText='Roller diameter - d [mm]:', tgt=form.roll_diameterKw, sel=0)
        AFXTextField(p=VAligner_14, ncols=12, labelText='Distance between rollers - s [mm]:', tgt=form.distance_between_rollsKw, sel=0)
        AFXTextField(p=VAligner_14, ncols=12, labelText='Number of upper rollers: ', tgt=form.upper_rolls_countKw, sel=0)
        AFXTextField(p=VAligner_14, ncols=12, labelText='Number of lower rollers:', tgt=form.lower_rolls_countKw, sel=0)
        AFXTextField(p=VAligner_14, ncols=12, labelText='Roller rotation velocity [mm/s]:', tgt=form.roll_velocityKw, sel=0)
        AFXTextField(p=VAligner_14, ncols=12, labelText='Friction coefficient - f1:', tgt=form.friction_coefficient1Kw, sel=0)
        AFXTextField(p=VAligner_14, ncols=12, labelText='Friction coefficient - f2: ', tgt=form.friction_coefficient2Kw, sel=0)
        AFXTextField(p=VAligner_14, ncols=12, labelText='Mesh density: ', tgt=form.mesh_densityKw, sel=0)
        AFXTextField(p=VAligner_14, ncols=12, labelText='Calibration input value:', tgt=form.calibration_inputKw, sel=0)
        AFXTextField(p=VAligner_14, ncols=12, labelText='Calibration output value:', tgt=form.calibration_outputKw, sel=0)
        ComboBox_4 = AFXComboBox(p=VAligner_14, ncols=0, nvis=1, text='Material:       ', tgt=form.material_nameKw, sel=0)
        ComboBox_4.setMaxVisible(10)
        ComboBox_4.appendItem(text='Steel S235JRH')
        ComboBox_4.appendItem(text='Steel S355J2H')
        ComboBox_4.appendItem(text='Steel S700')
	ComboBox_4.appendItem(text='Custom') #the indentation is messed up in notepad ++
	# current_item = ComboBox_4.getCurrentItem()
	# if current_item == 3:
	fileHandler = LevellingGUI1DBFileHandler(form, 'material_library', '*.lib')
	fileTextHf = FXHorizontalFrame(p=VAligner_14, opts=0, x=0, y=0, w=0, h=0,
	pl=0, pr=0, pt=0, pb=0, hs=DEFAULT_SPACING, vs=DEFAULT_SPACING)
	# Note: Set the selector to indicate that this widget should not be
	#       colored differently from its parent when the 'Color layout managers'
	#       button is checked in the RSG Dialog Builder dialog.
	'''
	HFrame_5 = FXHorizontalFrame(p=self, opts=LAYOUT_FILL_X, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        AFXTextField(p=HFrame_5, ncols=12, labelText='Initial plate flatness - p [mm]:         ', tgt=form.initial_flatnessKw, sel=0)
        AFXTextField(p=HFrame_5, ncols=12, labelText='Number of upper rollers:', tgt=form.upper_rolls_countKw, sel=0)
        AFXTextField(p=HFrame_5, ncols=12, labelText='Roller rotation velocity [mm/s]:', tgt=form.roll_velocityKw, sel=0)
        HFrame_4 = FXHorizontalFrame(p=self, opts=LAYOUT_FILL_X, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        AFXTextField(p=HFrame_4, ncols=12, labelText='Plate thickness [mm]:                       ', tgt=form.plate_thicknessKw, sel=0)
        AFXTextField(p=HFrame_4, ncols=12, labelText='Number of lower rollers:', tgt=form.lower_rolls_countKw, sel=0)
        AFXTextField(p=HFrame_4, ncols=12, labelText='Mesh density:                               ', tgt=form.mesh_densityKw, sel=0)
        HFrame_3 = FXHorizontalFrame(p=self, opts=LAYOUT_FILL_X, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        AFXTextField(p=HFrame_3, ncols=12, labelText='Distance between rollers - s [mm]: ', tgt=form.distance_between_rollsKw, sel=0)
        AFXTextField(p=HFrame_3, ncols=12, labelText='Friction coefficient - f1:  ', tgt=form.friction_coefficient1Kw, sel=0)
        AFXTextField(p=HFrame_3, ncols=12, labelText='Calibration input value:              ', tgt=form.calibration_inputKw, sel=0)
        HFrame_2 = FXHorizontalFrame(p=self, opts=LAYOUT_FILL_X, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
        AFXTextField(p=HFrame_2, ncols=12, labelText='Roller diameter - d [mm]:                ', tgt=form.roll_diameterKw, sel=0)
        AFXTextField(p=HFrame_2, ncols=12, labelText='Friction coefficient - f2:  ', tgt=form.friction_coefficient2Kw, sel=0)
        AFXTextField(p=HFrame_2, ncols=12, labelText='Calibration output value:            ', tgt=form.calibration_outputKw, sel=0)
        HFrame_1 = FXHorizontalFrame(p=self, opts=LAYOUT_FILL_X, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0)
	AFXTextField(p=HFrame_1, ncols=12, labelText='Time step [s]:                                      ', tgt=form.rolls_right_time_periodKw, sel=0)
        ComboBox_1 = AFXComboBox(p=HFrame_1, ncols=0, nvis=1, text='Material: ', tgt=form.material_nameKw, sel=0)
        ComboBox_1.setMaxVisible(10)
	'''
        ComboBox_1.appendItem(text='Item 1')
        ComboBox_1.appendItem(text='Item 2')
        ComboBox_1.appendItem(text='Item 3')
		'''
	ComboBox_1.appendItem(text='Steel S235JRH')
        ComboBox_1.appendItem(text='Steel S355J2H')
        ComboBox_1.appendItem(text='Steel S700')
	ComboBox_1.appendItem(text='Custom') 
		
        fileHandler = LevellingGUI1DBFileHandler(form, 'material_library', '*.lib')
        fileTextHf = FXHorizontalFrame(p=HFrame_1, opts=0, x=0, y=0, w=0, h=0,
            pl=0, pr=0, pt=0, pb=0, hs=DEFAULT_SPACING, vs=DEFAULT_SPACING)
	fileTextHf.setSelector(99)
	AFXTextField(p=fileTextHf, ncols=12, labelText='Custom material (if "Custom" chosen above ):', tgt=form.material_libraryKw, sel=0,
	opts=AFXTEXTFIELD_STRING|LAYOUT_CENTER_Y)
	icon = afxGetIcon('fileOpen', AFX_ICON_SMALL )
	FXButton(p=fileTextHf, text='	Select File\nFrom Dialog', ic=icon, tgt=fileHandler, sel=AFXMode.ID_ACTIVATE,
	opts=BUTTON_NORMAL|LAYOUT_CENTER_Y, x=0, y=0, w=0, h=0, pl=1, pr=1, pt=1, pb=1)




###########################################################################
# Class definition
###########################################################################

class LevellingGUI1DBFileHandler(FXObject):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form, keyword, patterns='*'):

        self.form = form
        self.patterns = patterns
        self.patternTgt = AFXIntTarget(0)
        exec('self.fileNameKw = form.%sKw' % keyword)
        self.readOnlyKw = AFXBoolKeyword(None, 'readOnly', AFXBoolKeyword.TRUE_FALSE)
        FXObject.__init__(self)
        FXMAPFUNC(self, SEL_COMMAND, AFXMode.ID_ACTIVATE, LevellingGUI1DBFileHandler.activate)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def activate(self, sender, sel, ptr):

       fileDb = AFXFileSelectorDialog(getAFXApp().getAFXMainWindow(), 'Select a File',
           self.fileNameKw, self.readOnlyKw,
           AFXSELECTFILE_ANY, self.patterns, self.patternTgt)
       fileDb.setReadOnlyPatterns('*.odb')
       fileDb.create()
       fileDb.showModal()

