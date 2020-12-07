import ModelCommands2D

def levellingFunction(roll_diameter, plate_thickness, initial_flatness, mesh_density, upper_rolls_count, lower_rolls_count, distance_between_rolls,
					  roll_velocity, material_name, calibration_input, calibration_output, friction_coefficient1, friction_coefficient2, material_library, rolls_right_time_period): 
	modelCommands = ModelCommands2D.ModelCommands2D(roll_diameter, plate_thickness, initial_flatness, mesh_density, upper_rolls_count, lower_rolls_count, distance_between_rolls,
					  roll_velocity, material_name, calibration_input, calibration_output, friction_coefficient1, friction_coefficient2, material_library, rolls_right_time_period) 
	modelCommands.create_model()
	modelCommands.create_upper_and_lower_roll()
	modelCommands.create_input_plate()
	modelCommands.create_plate_with_edge_waves()
	modelCommands.create_levelling_table()
	modelCommands.create_example_materials()
	modelCommands.mesh_parts()
	modelCommands.create_and_assign_sections()
	modelCommands.create_mesh_sets()
	modelCommands.prepare_assembly()
	# modelCommands.create_interaction_properties()
	modelCommands.create_interactions()
	modelCommands.apply_boundary_conditions()
	modelCommands.create_further_steps()
	modelCommands.calibrate_upper_rolls()
	modelCommands.set_displacements_and_rotation_values()
	modelCommands.set_field_output_requests()
	modelCommands.set_history_output_requests()
	modelCommands.create_job()
	modelCommands.save_mdb()