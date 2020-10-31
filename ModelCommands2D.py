from abaqus import *
from abaqusConstants import *
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import sys
from os.path import expanduser
import gauss_elimination
import pickle
import os
from material import createMaterialFromDataString

	
class ModelCommands2D(object):
	def __init__(self, roll_diameter, plate_thickness, initial_flatness, mesh_density, upper_rolls_count, lower_rolls_count, distance_between_rolls,
					  roll_velocity, material_name, calibration_input, calibration_output, friction_coefficient1, friction_coefficient2, material_library):
		self.roll_diameter = roll_diameter 
		self.plate_thickness = plate_thickness 
		self.initial_flatness = initial_flatness 
		self.mesh_density = mesh_density
		self.upper_rolls_count = upper_rolls_count 
		self.lower_rolls_count = lower_rolls_count 
		self.distance_between_rolls = distance_between_rolls
		self.roll_velocity = roll_velocity 
		self.material_name = material_name
		self.calibration_input = calibration_input 
		self.calibration_output = calibration_output 
		self.friction_coefficient1 = friction_coefficient1
		self.friction_coefficient2 = friction_coefficient2
 		self.input_plate_length = (self.lower_rolls_count * self.distance_between_rolls) + self.distance_between_rolls
		self.model_name = 'Leveling_D{}G{}_T{}_D{}c_S{}c_F{}_{}'.format(str(self.lower_rolls_count), str(self.upper_rolls_count), 
		str(self.plate_thickness), str(self.roll_diameter), str(self.distance_between_rolls), str(self.initial_flatness), str(self.material_name).replace(" ","_") )
		self.material_library = material_library
	
	def create_model(self):
		mdb.Model(name=self.model_name, absoluteZero=0,
        stefanBoltzmann=5.6697E-11, modelType=STANDARD_EXPLICIT, 
        copyInteractions=ON, copyConnectors=ON, copyConstraints=ON)
	
	def create_upper_and_lower_roll(self):
		s = mdb.models[self.model_name].ConstrainedSketch(
        name='__profile__', sheetSize=200.0)
		g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
		s.setPrimaryObject(option=STANDALONE)
		s.ConstructionLine(point1=(0.0, 0.0), point2=(7.77105331420898, 0.0))
		s.HorizontalConstraint(entity=g[2], addUndoState=False)
		s.ConstructionLine(point1=(0.0, 0.0), point2=(0.0, -7.69877910614014))
		s.VerticalConstraint(entity=g[3], addUndoState=False)
		s.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, self.roll_diameter/2.0), point2=(self.roll_diameter/2.0, 0.0), 
			direction=CLOCKWISE)
		s.copyMirror(mirrorLine=g[2], objectList=(g[4], ))
		s.copyMirror(mirrorLine=g[3], objectList=(g[4], g[5]))
		p = mdb.models[self.model_name].Part(name='Lowerroll', 
			dimensionality=TWO_D_PLANAR, type=ANALYTIC_RIGID_SURFACE)
		p = mdb.models[self.model_name].parts['Lowerroll']
		p.AnalyticRigidSurf2DPlanar(sketch=s)
		v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
		p.ReferencePoint(point=p.InterestingPoint(edge=e[0], rule=CENTER))
		s = p.edges
		side2Edges = s.getSequenceFromMask(mask=('[#f ]', ), )
		p.Surface(side2Edges=side2Edges, name='surface')
		p = mdb.models[self.model_name].parts['Lowerroll']
		session.viewports['Viewport: 1'].setValues(displayedObject=p)
		del mdb.models[self.model_name].sketches['__profile__']
		p = mdb.models[self.model_name].Part(
			name='Upper roll', 
			objectToCopy=mdb.models[self.model_name].parts['Lowerroll'])
		pass 
	
	def create_input_plate(self): 
		x1 = (float(self.input_plate_length)/2.0) * -1.0
		x2 = float(self.input_plate_length)/2.0
		s = mdb.models[self.model_name].ConstrainedSketch(
        name='__profile__', sheetSize=5000.0)
		g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
		s.setPrimaryObject(option=STANDALONE)
		s.rectangle(point1=(x1, self.plate_thickness), point2=(x2, 0.0)) 
		p = mdb.models[self.model_name].Part(
        name='input_plate', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
		p = mdb.models[self.model_name].parts['input_plate']
		p.BaseShell(sketch=s)
		e = p.edges
		edges = e.getSequenceFromMask(mask=('[#8 ]', ), )
		p.Set(edges=edges, name='connect')
		s = p.edges
		side1Edges = s.getSequenceFromMask(mask=('[#f ]', ), )
		p.Surface(side1Edges=side1Edges, name='all')
		side1Edges = s.getSequenceFromMask(mask=('[#8 ]', ), )
		p.Surface(side1Edges=side1Edges, name='connection')
		side1Edges = s.getSequenceFromMask(mask=('[#c ]', ), )
		p.Surface(side1Edges=side1Edges, name='lower')
		side1Edges = s.getSequenceFromMask(mask=('[#9 ]', ), )
		p.Surface(side1Edges=side1Edges, name='upper')
		pass
	
	
	def create_plate_with_edge_waves(self):
		s1 = mdb.models[self.model_name].ConstrainedSketch(
        name='__profile__', sheetSize=5000.0)
		g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
		s1.setPrimaryObject(option=STANDALONE)
		
		# top 
		x1 = 0
		y1 = self.plate_thickness
		# bottom 
		x2 = 0
		y2 = 0 
		top_start = (x1, y1)
		bottom_start = (x2, y2)
		
		top_spots = []
		bottom_spots = []
		
		top_spots.append(top_start)
		bottom_spots.append(bottom_start)
		
		while x1 < 3500 and x2 < 3500:
			x1+=500
			x2+=500
			
			if (x1 == 500 or x1 == 1500 or x1 == 2500 or x1 == 3500):
				y1 = self.initial_flatness/2.0
				y2 = (self.initial_flatness/2.0) - self.plate_thickness
			elif (x1 == 1000 or x1 == 2000 or x1 == 3000):
				y1 = (self.initial_flatness/2.0) - self.initial_flatness
				y2 = -(self.initial_flatness/2.0) - self.plate_thickness
			top_spots.append((x1,y1))
			bottom_spots.append((x2,y2))
		
		s1.Line(point1=top_start, point2=bottom_start)
		s1.VerticalConstraint(entity=g[2], addUndoState=False)
		s1.Line(point1=top_spots[-1], point2=bottom_spots[-1])
		s1.VerticalConstraint(entity=g[3], addUndoState=False)
		s1.Spline(points=tuple(top_spots))
		s1.Spline(points=tuple(bottom_spots))
	
		p = mdb.models[self.model_name].Part(name='plate', 
        dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
		p.BaseShell(sketch=s1)
		
		# datum planes
		p = mdb.models[self.model_name].parts['plate']
		offset = 0.0
		for i in range(0,8):
			p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=offset)
			offset+=500
		
		s = p.edges
		side1Edges = s.getSequenceFromMask(mask=('[#2 ]', ), )
		p.Surface(side1Edges=side1Edges, name='Left')
		side1Edges = s.getSequenceFromMask(mask=('[#1 ]', ), )
		p.Surface(side1Edges=side1Edges, name='Upper')
		side1Edges = s.getSequenceFromMask(mask=('[#4 ]', ), )
		p.Surface(side1Edges=side1Edges, name='Lower')
		side1Edges = s.getSequenceFromMask(mask=('[#8 ]', ), )
		p.Surface(side1Edges=side1Edges, name='Right')
		
		
		f = p.faces
		faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
		p.Set(faces=faces, name='all')
		
		pass
	
		
	def create_levelling_table(self):
		s = mdb.models[self.model_name].ConstrainedSketch(
        name='__profile__', sheetSize=8000.0)
		g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
		s.setPrimaryObject(option=STANDALONE)
		s.Line(point1=(0.0, 0.0), point2=(8000.0, 0.0))
		s.HorizontalConstraint(entity=g[2], addUndoState=False)
		p = mdb.models[self.model_name].Part(name='table', 
			dimensionality=TWO_D_PLANAR, type=ANALYTIC_RIGID_SURFACE)
		p.AnalyticRigidSurf2DPlanar(sketch=s)
		v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
		p.ReferencePoint(point=v1[0])
		pass 
		
	
	def create_example_materials(self):
		if self.material_name == 'Steel S235JRH':
			mdb.models[self.model_name].Material(
			name='Steel_S235JRH')
			mdb.models[self.model_name].materials['Steel_S235JRH'].Density(
			table=((7.85e-09, ), ))
			mdb.models[self.model_name].materials['Steel_S235JRH'].Elastic(
			table=((205000.0, 0.3), ))
			mdb.models[self.model_name].materials['Steel_S235JRH'].Plastic(
			hardening=COMBINED, dataType=PARAMETERS, table=((282.0, 122.4, 430.0), 
			))
			mdb.models[self.model_name].materials['Steel_S235JRH'].plastic.CyclicHardening(
			parameters=ON, table=((282.0, 35.0, 1.12), ))
		if self.material_name == 'Steel S355J2H':
			mdb.models[self.model_name].Material(
			name='Steel_S355J2H')
			mdb.models[self.model_name].materials['Steel_S355J2H'].Density(
			table=((7.85e-09, ), ))
			mdb.models[self.model_name].materials['Steel_S355J2H'].Elastic(
			table=((215000.0, 0.3), ))
			mdb.models[self.model_name].materials['Steel_S355J2H'].Plastic(
			hardening=COMBINED, dataType=PARAMETERS, table=((461.0, 23.5, 139.0), 
			))
			mdb.models[self.model_name].materials['Steel_S355J2H'].plastic.CyclicHardening(
			parameters=ON, table=((461.0, 55.0, 2.38), ))
		if self.material_name == 'Steel S700':
			mdb.models[self.model_name].Material(
			name='Steel_S700')
			mdb.models[self.model_name].materials['Steel_S700'].Density(
			table=((7.85e-09, ), ))
			mdb.models[self.model_name].materials['Steel_S700'].Elastic(
			table=((190000.0, 0.3), ))
			mdb.models[self.model_name].materials['Steel_S700'].Plastic(
			hardening=COMBINED, dataType=PARAMETERS, numBackstresses=4, table=((
			650.0, 23300.0, 500.0, 11330.0, 340.0, 8000.0, 180.0, 550.0, 20.0), ))
			mdb.models[self.model_name].materials['Steel_S700'].plastic.CyclicHardening(
			parameters=ON, table=((650.0, -60.0, 100.0), ))
		if self.material_name == 'Custom':
			# lib = open(self.material_library, "r")
			# todo parse the material
			# print(destination[1][2]) #material name
			# print(destination[1][4]['Data']) #material data
			material_data = self.read_material_data_from_file(self.material_library)
			material_name = material_data[1][2]
			material_data_string = material_data[1][4]['Data']
			createMaterialFromDataString(self.model_name, material_name, version, material_data_string)
			self.material_name = material_name
			prev_model_name = self.model_name
			# self.model_name = '{}{}'.format(self.model_name, material_name)
			self.model_name = self.model_name.replace('Custom', material_name)
			mdb.models.changeKey(fromName=prev_model_name, 
				toName=self.model_name)
		
	
	def mesh_parts(self):
		parts_to_mesh = ['input_plate', 'plate']
		for index, part in enumerate(parts_to_mesh):
			p = mdb.models[self.model_name].parts[part]
			p.seedPart(size=self.mesh_density, deviationFactor=0.1, minSizeFactor=0.1)
			elemType1 = mesh.ElemType(elemCode=CPE4R, elemLibrary=STANDARD, 
			secondOrderAccuracy=OFF, hourglassControl=DEFAULT, 
			distortionControl=DEFAULT)
			elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD, 
			secondOrderAccuracy=OFF, distortionControl=DEFAULT)
			f = p.faces
			faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
			pickedRegions =(faces, )
			p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
			p.generateMesh()
			p.deleteMesh(regions=pickedRegions)
			f = p.faces
			pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
			p.setMeshControls(regions=pickedRegions, technique=STRUCTURED)
			p.generateMesh()
		pass
		
	def create_and_assign_sections(self):
		material_name = self.material_name.replace(" ", "_")
		
		parts_to_assign = ['input_plate', 'plate']
		
		mdb.models[self.model_name].HomogeneousSolidSection(
			name='Section-1', material=material_name, thickness=2000.0)
			
		for index, part in enumerate(parts_to_assign):
			p = mdb.models[self.model_name].parts[part]
			f = p.faces
			faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
			region = p.Set(faces=faces, name='Set-1')
			p = mdb.models[self.model_name].parts[part]
			p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
				offsetType=MIDDLE_SURFACE, offsetField='', 
				thicknessAssignment=FROM_SECTION)
			del mdb.models[self.model_name].parts[part].sets['Set-1']
			f = p.faces
			faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
			region=regionToolset.Region(faces=faces)
			mdb.models[self.model_name].parts[part].sectionAssignments[0].setValues(
			region=region)
		pass
	
	def create_mesh_sets(self):
		p = mdb.models[self.model_name].parts['plate']
		
		# element sets
		e = p.elements
		elements_per_width = int(self.plate_thickness/self.mesh_density)
		elements_per_length = int(len(e)/elements_per_width)
		
		upper_elements = e[0:elements_per_length]
		p.Set(elements=upper_elements, name='Upper_Elements')
		
		inside_start = int(elements_per_width/2)*elements_per_length
		inside_end = inside_start + elements_per_length
		inside_elements = e[inside_start:inside_end]
		p.Set(elements=inside_elements, name='Inside_Elements')

		lower_start = (elements_per_width-1)*elements_per_length
		lower_end = lower_start + elements_per_length
		lower_elements = e[lower_start:lower_end]
		p.Set(elements=lower_elements, name='Lower_Elements')
		
		
		# peeq500 - peeq3000
		index = int(elements_per_length/7)*6
		set_number = 500
		
		for i in range(6,0,-1):
			peeq_set = e[index:index+1]
			for j in range(1, elements_per_width): # od 1 bo pierwszy juz mamy
				index-=elements_per_length
				peeq_set+=e[index:index+1]
			set_name = 'PEEQ'+str(set_number)
			p.Set(elements=peeq_set, name=set_name)
			set_number +=500
			if i != 1:
				index=int(elements_per_length/7)*(i-1)
		
		# node sets
		n = p.nodes
		nodes_per_width = int(self.plate_thickness/self.mesh_density)+1
		nodes_per_length = int(len(e)/elements_per_width)+1
		upper_nodes = n[0:nodes_per_length]
		p.Set(nodes=upper_nodes, name='Upper_nodes')
		lower_nodes_start = (nodes_per_width-1)*nodes_per_length
		lower_nodes_end = lower_nodes_start + nodes_per_length
		lower_nodes = n[lower_nodes_start:lower_nodes_end]
		p.Set(nodes=lower_nodes, name='Lower_nodes')
				
		pass 
	
	
	
	def prepare_assembly(self):
		a = mdb.models[self.model_name].rootAssembly
		a.DatumCsysByDefault(CARTESIAN)
		lower_ref_points = tuple()
		
		# create instances  
		p = mdb.models[self.model_name].parts['input_plate']
		a.Instance(name='input_plate-1', part=p, dependent=ON)
		plate_vector = (-1)*(float(self.input_plate_length)/2.0)
		a.translate(instanceList=('input_plate-1', ), vector=(plate_vector, 0.0, 0.0))

		
		p = mdb.models[self.model_name].parts['Lowerroll']
		lower_aux_roll_x = (float(self.roll_diameter)/2.0)*(-1.0)
		r = p.referencePoints
		refPoints=(r[2], )
		p.Set(referencePoints=refPoints, name='RP')
		a.Instance(name='Roller-2', part=p, dependent=ON)
		a.translate(instanceList=('Roller-2', ), vector=(500.0, lower_aux_roll_x, 0.0))
	
		for i in range(0, self.lower_rolls_count):
			lower_roll_x = -24.0 - float(i*self.distance_between_rolls) 
			lower_roll_y = (float(self.roll_diameter)/2.0)*(-1.0)
			lower_roll_z = 0.0
			instance_name = 'Lowerroll-'+str(i+1)
			r = p.referencePoints
			refPoints=(r[2], )
			p.Set(referencePoints=refPoints, name='RP')
			a.Instance(name=instance_name, part=p, dependent=ON)
			a.translate(instanceList=(instance_name, ), vector=(lower_roll_x, lower_roll_y, lower_roll_z))
			rp = a.instances[instance_name].referencePoints
			lower_ref_points+=(rp[2], )
			
		a.Set(referencePoints=lower_ref_points, name='All lower rolls')
				
		
		
		p = mdb.models[self.model_name].parts['Upper roll']
		lower_aux_roll_y = float(self.roll_diameter + self.plate_thickness)
		r = p.referencePoints
		refPoints=(r[2], )
		p.Set(referencePoints=refPoints, name='RP')
		a.Instance(name='Roller-1', part=p, dependent=ON)
		a.translate(instanceList=('Roller-1', ), vector=(500.0, lower_aux_roll_y , 0.0))
		
		for i in range(0, self.upper_rolls_count):
			upper_roll_starting_x = -24.0 - float(self.distance_between_rolls/2.0)
			upper_roll_x = upper_roll_starting_x - float(i*self.distance_between_rolls)
			upper_roll_y = (float(self.roll_diameter)/2.0)+float(self.plate_thickness)
			upper_roll_z = 0.0
			instance_name = 'Upper roll-'+str(i+1)
			r = p.referencePoints
			refPoints=(r[2], )
			p.Set(referencePoints=refPoints, name='RP')
			a.Instance(name=instance_name, part=p, dependent=ON)
			a.translate(instanceList=(instance_name, ), vector=(upper_roll_x, upper_roll_y, upper_roll_z))
			
		p = mdb.models[self.model_name].parts['table']
		
		for i in range(0, 4):
			instance_name = 'table-'+str(i+1)
			a.Instance(name=instance_name, part=p, dependent=ON)
			s = a.instances[instance_name].edges
			side2Edges = s.getSequenceFromMask(mask=('[#1 ]', ), )
			a.Surface(side2Edges=side2Edges, name=instance_name)
			
		table_x1 = (8000.0 + float(self.input_plate_length)/2.0)*(-1)
		table_x2 = -770
		table_y1 = float(2.0*self.plate_thickness)
		table_y2 = float(self.plate_thickness)*(-1.0)
		table_y3 = float(self.plate_thickness)+float(self.initial_flatness/2.0)
		table_y4 = (-1.0)*table_y3
		
		a.translate(instanceList=('table-4', ), vector=(table_x1, table_y1, 0.0))
		r = p.referencePoints
		refPoints=(r[2], )
		p.Set(referencePoints=refPoints, name='RP')
		a.translate(instanceList=('table-3', ), vector=(table_x1, table_y2, 0.0))
		r = p.referencePoints
		refPoints=(r[2], )
		p.Set(referencePoints=refPoints, name='RP')
		a.translate(instanceList=('table-2', ), vector=(table_x2, table_y3, 0.0))
		r = p.referencePoints
		refPoints=(r[2], )
		p.Set(referencePoints=refPoints, name='RP')
		a.translate(instanceList=('table-1', ), vector=(table_x2, table_y4, 0.0))
		r = p.referencePoints
		refPoints=(r[2], )
		p.Set(referencePoints=refPoints, name='RP')
		
		p = mdb.models[self.model_name].parts['plate']
		a.Instance(name='plate-1', part=p, dependent=ON)
		
		e1 = a.instances['input_plate-1'].edges
		edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
		a.Set(edges=edges1, name='begin')
		a.Set(edges=edges1, name='symmetry')
		
		s1 = a.instances['input_plate-1'].edges
		side1Edges1 = s1.getSequenceFromMask(mask=('[#7 ]', ), )
		s2 = a.instances['plate-1'].edges
		side1Edges2 = s2.getSequenceFromMask(mask=('[#d ]', ), )
		a.Surface(side1Edges=side1Edges1+side1Edges2, name='whole plate')
		
		a = mdb.models[self.model_name].rootAssembly
		region1=a.instances['plate-1'].surfaces['Left']
		a = mdb.models[self.model_name].rootAssembly
		region2=a.instances['input_plate-1'].surfaces['connection']
		mdb.models[self.model_name].Tie(name='Constraint-1', 
        master=region1, slave=region2, positionToleranceMethod=COMPUTED, 
        adjust=ON, tieRotations=ON, thickness=ON)
			
		pass
	
	def create_surface_to_surface_contact(self, master_region, slave_region, interaction_name, friction_coefficient):
		property_name = self.create_interaction_properties(friction_coefficient)
		mdb.models[self.model_name].SurfaceToSurfaceContactStd(
			name=interaction_name, createStepName='Initial', master=master_region, slave=slave_region, 
			sliding=FINITE, thickness=ON, contactTracking=ONE_CONFIG, 
			interactionProperty=property_name, adjustMethod=NONE, 
			initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

		pass
		
	def create_interaction_properties(self, friction_coefficient): 
		property_name = 'Frictionless'+str(friction_coefficient).replace(".", "")
		mdb.models[self.model_name].ContactProperty(
		property_name)
		mdb.models[self.model_name].interactionProperties[property_name].TangentialBehavior(
		formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
		pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
		table=((friction_coefficient, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
		fraction=0.005, elasticSlipStiffness=None)
		return property_name
		
	def create_interactions(self):
		a = mdb.models[self.model_name].rootAssembly
		
		for i in range (0,2):
			friction_coefficient = self.friction_coefficient1
			instance_name = 'Roller-'+str(i+1)
			region1=a.instances[instance_name].surfaces['surface'] # master
			region2=a.surfaces['whole plate'] # slave
			interaction_name = 'R'+str(i+1)
			self.create_surface_to_surface_contact(region1, region2, interaction_name, friction_coefficient) 
		
		# upper roll surf to surf interactions for initial step 
		for i in range(0, self.lower_rolls_count):
			instance_name = 'Lowerroll-'+str(i+1)
			region1=a.instances[instance_name].surfaces['surface']
			region2=a.surfaces['whole plate']
			interaction_name = 'L'+str(i+1)
			if i < 3:
				friction_coefficient = self.friction_coefficient2
			else: 
				friction_coefficient = self.friction_coefficient1
			self.create_surface_to_surface_contact(region1, region2, interaction_name, friction_coefficient)
		
		# upper roll surf to surf interactions for initial step 
		for i in range(0, self.upper_rolls_count):
			instance_name = 'Upper roll-'+str(i+1)
			region1=a.instances[instance_name].surfaces['surface']
			region2=a.surfaces['whole plate']
			interaction_name = 'U'+str(i+1)
			if i < 2:
				friction_coefficient = self.friction_coefficient2
			else: 
				friction_coefficient = self.friction_coefficient1
			self.create_surface_to_surface_contact(region1, region2, interaction_name, friction_coefficient)
		
		for i in range(0,4):
			friction_coefficient = self.friction_coefficient1
			master_surface_name = 'table-'+str(i+1)
			region1=a.surfaces[master_surface_name]
			region2=a.surfaces['whole plate']
			interaction_name = 'T'+str(i+1)
			self.create_surface_to_surface_contact(region1, region2, interaction_name, friction_coefficient) 
			
		property = self.create_interaction_properties(self.friction_coefficient1)	
		mdb.models[self.model_name].ContactStd(name='General', 
        createStepName='Initial')
		mdb.models[self.model_name].interactions['General'].includedPairs.setValuesInStep(
        stepName='Initial', useAllstar=ON)
		mdb.models[self.model_name].interactions['General'].contactPropertyAssignments.appendInStep(
        stepName='Initial', assignments=((GLOBAL, SELF, property), ))
			
		pass

	def create_pinned_bc_input_plate(self):
		a = mdb.models[self.model_name].rootAssembly
		e1 = a.instances['input_plate-1'].edges 
		edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
		region = regionToolset.Region(edges=edges1)
		mdb.models[self.model_name].PinnedBC(name='Begin', 
			createStepName='Initial', region=region, localCsys=None)
		pass
	
	def create_displacement_bc(self, region, bc_name):
		mdb.models[self.model_name].DisplacementBC(
        name=bc_name, createStepName='Initial', region=region, u1=SET, 
        u2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, 
        fieldName='', localCsys=None)
		pass
	
	def create_encastre_bc(self,region, bc_name):
		mdb.models[self.model_name].EncastreBC(name=bc_name, 
				createStepName='Initial', region=region, localCsys=None)
		pass 
	
	def apply_boundary_conditions(self): # bcs for initial step
		a = mdb.models[self.model_name].rootAssembly
		self.create_pinned_bc_input_plate()
		region = a.sets['All lower rolls']
		self.create_displacement_bc(region, 'Bottom_rolls')
		
		# auxillary rolls 
		for i in range(0,2):
			instance_name = 'Roller-'+str(i+1)
			bc_name = 'R'+str(i+1)
			region = a.instances[instance_name].sets['RP']
			self.create_displacement_bc(region, bc_name)
		
		# table 
		for i in range(0,4):
			set_name = 'table-'+str(i+1)
			region = a.instances[set_name].sets['RP']
			bc_name = 'T'+str(i+1)
			
			if i == 0 or i == 1: 
				self.create_encastre_bc(region, bc_name)
			else:
				self.create_displacement_bc(region, bc_name)
				
				
		#upper rolls
		for i in range(0,self.upper_rolls_count):
			instance_name = 'Upper roll-'+str(i+1)
			region = a.instances[instance_name].sets['RP']
			bc_name = 'Upper-'+str(i+1)
			self.create_displacement_bc(region, bc_name)

		pass
	
	def create_further_steps(self):
		self.create_static_step('Rolls-down', 'Initial', 10000, 1e-05, 1e-20, 0.05, 1)
		self.create_static_step('Stabilize', 'Rolls-down', 10000, 0.001, 1e-10, 0.05, 1)
		self.create_static_step('Rolls-right','Rolls-down', 100000, 1e-05, 1e-10, 1.0, 6.13)
		self.set_gravity_in_step('Rolls-down', -9810.0)
		
		pass
	
	def set_gravity_in_step(self, step_name, value):
		mdb.models[self.model_name].Gravity(name='Load-1', 
        createStepName=step_name, comp2=value, distributionType=UNIFORM, 
        field='')
		pass
		
	def create_static_step(self, step_name, previous_step_name, max_num_inc, initial_inc, min_inc, max_inc, time_period):
		if time_period != 1: 
			mdb.models[self.model_name].StaticStep(
			name=step_name, previous=previous_step_name, timePeriod = time_period, maxNumInc=max_num_inc, 
			stabilizationMagnitude=0.0002, stabilizationMethod=DAMPING_FACTOR, 
			continueDampingFactors=False, adaptiveDampingRatio=None, 
			initialInc=initial_inc, minInc=min_inc, maxInc=max_inc, nlgeom=ON)
		else: 
			mdb.models[self.model_name].StaticStep(
			name=step_name, previous=previous_step_name, maxNumInc=max_num_inc, 
			stabilizationMagnitude=0.0002, stabilizationMethod=DAMPING_FACTOR, 
			continueDampingFactors=False, adaptiveDampingRatio=None, 
			initialInc=initial_inc, minInc=min_inc, maxInc=max_inc, nlgeom=ON)
		pass
	
	def calibrate_upper_rolls(self):
		# calculate coefficients of line equation
		matrix = [[2.0, 1.0, self.calibration_input], [self.upper_rolls_count*2.0, 1.0, self.calibration_output]]
		coeffs = gauss_elimination.gauss_elimination(matrix, 2)
		
		for i in range (0, self.upper_rolls_count):
			y = (coeffs[0] * (2.0*(i+1))) + coeffs[1] # y = ax + b 
			if self.plate_thickness - y != 0:
				u2 = (self.plate_thickness - y) * -1.0
			else:
				u2 = self.plate_thickness - y
			mdb.models[self.model_name].boundaryConditions['Upper-'+str(i+1)].setValuesInStep(
			stepName='Rolls-down', u2=u2) 
			
		r2_u2 = - 1.0 * (self.roll_diameter/2.0) + 2.0
		mdb.models[self.model_name].boundaryConditions['R1'].setValuesInStep(
        stepName='Rolls-down', u2=r2_u2)
		pass 
	
	def set_displacements_and_rotation_values(self): # todo najpierw przeliczyć roll velocity na rad/s 
		radius = float(self.roll_diameter)/2.0
		angular_velocity = float(self.roll_velocity)/radius
		upper_ur3 = -1.0 * angular_velocity 
		mdb.models[self.model_name].boundaryConditions['Bottom_rolls'].setValuesInStep(
			stepName='Rolls-right', u1=4600.0)
		mdb.models[self.model_name].boundaryConditions['R1'].setValuesInStep(
			stepName='Rolls-right', u1=4600.0, ur3=upper_ur3)
		mdb.models[self.model_name].boundaryConditions['R2'].setValuesInStep(
			stepName='Rolls-right', u1=4600.0, ur3=angular_velocity)
		mdb.models[self.model_name].boundaryConditions['T3'].setValuesInStep(
			stepName='Rolls-right', u1=4600.0)
		mdb.models[self.model_name].boundaryConditions['T4'].setValuesInStep(
			stepName='Rolls-right', u1=4600.0)
		
		for i in range(0, self.upper_rolls_count):
			mdb.models[self.model_name].boundaryConditions['Upper-'+str(i+1)].setValuesInStep(
				stepName='Rolls-right', u1=4600.0, ur3=upper_ur3)
		pass
	
	
	
	def create_job(self):
		mdb.Job(name=self.model_name,
        model=self.model_name, description='', 
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numGPUs=0)
		
		pass
	
	def set_field_output_requests(self):
		mdb.models[self.model_name].fieldOutputRequests['F-Output-1'].setValues(
        variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'RM'), 
        timeInterval=0.2)
		pass
	
	def set_history_output_requests(self):
		regionDef=mdb.models[self.model_name].rootAssembly.sets['All lower rolls']
		mdb.models[self.model_name].historyOutputRequests['H-Output-1'].setValues(
			variables=('RF1', 'RF2', 'RM1', 'RM2'), region=regionDef, 
			sectionPoints=DEFAULT, rebar=EXCLUDE)
		roll_name = 'Upper roll-'+str(self.upper_rolls_count)
		regionDef=mdb.models[self.model_name].rootAssembly.allInstances[roll_name].sets['RP']
		mdb.models[self.model_name].HistoryOutputRequest(
			name='RF_U1', createStepName='Rolls-down', variables=('RF1', 'RF2', 
			'RM1', 'RM2'), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
		regionDef=mdb.models[self.model_name].rootAssembly.sets['All lower rolls']
		mdb.models[self.model_name].HistoryOutputRequest(
			name='RF_D', createStepName='Rolls-down', variables=('RF1', 'RF2', 
			'RM1', 'RM2'), frequency=1, region=regionDef, sectionPoints=DEFAULT, 
			rebar=EXCLUDE)
		del mdb.models[self.model_name].historyOutputRequests['H-Output-1']
		pass
	
	def save_mdb(self):
		home = expanduser('~/'+self.model_name+'.cae')
		mdb.saveAs(home) 
		pass
	
	def read_material_data_from_file(self, material_library):
			
		# lib = open("C:\\SIMULIA\\Commands\\abaqus_plugins\\Steel235jrh1.lib", "rb")
		# material = pickle.load(lib)
		# for keys in material:
		#     print(keys, '=>', material[keys])
		# lib.close()
		
		# empt_dict = {} 
		# path = "C:\\SIMULIA\\Commands\\abaqus_plugins\\material.lib"
		# original = path
		# destination = "C:\\SIMULIA\\Commands\\abaqus_plugins\\materialunix.lib" #todo usuwac po dodaniu materialu
		
		origin = self.material_library
		destination = expanduser('~/'+'tmpunix.lib')

		content = ''
		outsize = 0
		with open(origin, 'rb') as infile:
			content = infile.read()
		with open(destination, 'wb') as output:
			for line in content.splitlines():
				outsize += len(line) + 1
				output.write(line + str.encode('\n'))

		if os.path.getsize(destination) > 0:
			with open(destination, "rb") as f:
				unpickler = pickle.Unpickler(f)
				# if file is not empty scores will be equal
				# to the value unpickled
				destination = unpickler.load()
				# os.remove(expanduser('~/'+'tmpunix.lib'))
				return destination
				# print(destination)
				# print(destination[1][2]) #material name
				# print(destination[1][4]['Data']) #material data
				# for keys in destination:
				#     print(keys, '=>', destination[keys])


			
