package main

import "core:log"
import "core:mem/virtual"
import "core:os"
import "core:slice"

import "../fem"
import "../fem/infra"
import "../fem/serialization"

import "conduction"

main :: proc() {
	context.logger = log.create_console_logger()

	arena: virtual.Arena

	_ = virtual.arena_init_growing(&arena)

	context.allocator = virtual.arena_allocator(&arena)

	mesh_path := os.args[1]

	mesh, warn, err := serialization.gmsh_parse(mesh_path)
	if err != nil {log.panic(err)}
	if warn != nil {log.warn(warn)}

	fem.setup_default_rules()

	params: conduction.Model_Parameters

	params.isothermal_bcs = make(map[fem.Boundary_ID]conduction.Isothermal_Int)
	params.materials = make(map[fem.Section_ID]conduction.Material_Int)

	id := mesh.boundary_names["constraint"] or_else 0
	//id2 := mesh.boundary_names["right"] or_else 1

	params.isothermal_bcs[id] = {
		procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr) -> conduction.Isothermal_BC {
			return 100
		},
	}

	// params.isothermal_bcs[id2] = {
	// 	procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr) -> conduction.Isothermal_BC {
	// 		return 15
	// 	},
	// }

	params.materials[mesh.elements[0].section] = {
		procedure = proc(
			mapped: fem.Mapped_Element,
			time: f64,
			data: rawptr,
			temperature: []f64,
			out: conduction.Material,
		) {
			slice.fill(out.k, 1)
			slice.fill(out.cp, 2)
			slice.fill(out.rho, 1)
		},
	}


	conduction.solve(params, mesh, .Linear, &arena)
}

// main :: proc() {
// 	context.logger = log.create_console_logger()

// 	ca := infra.Checkpoint_Allocator{}
// 	infra.ca_init(&ca)
// 	defer infra.ca_deinit(&ca)

// 	mesh_path := os.args[1]

// 	mesh, warn, err := serialization.gmsh_parse(mesh_path)
// 	if err != nil {log.panic(err)}
// 	if warn != nil {log.warn(warn)}

// 	fem.setup_default_rules()

// 	// system setup

// 	sys_desc := fem.System_Description{}

// 	phi_handle := fem.description_add_variable(
// 		&sys_desc,
// 		{bd = {.Lagrange, .Linear}, components = 1},
// 	)

// 	fem.description_couple(&sys_desc, phi_handle, phi_handle)

// 	system := fem.system_from_description(mesh, sys_desc, context.allocator)

// 	ics := fem.system_vector(system)
// 	cm := fem.system_constraint_mask(system)


// 	ni := fem.nli_create(system, ics, 1e-4, 5)

// 	for name, id in mesh.boundary_names {
// 		if name == "free" {continue}

// 		fem.system_mark_boundary(system, phi_handle, id, cm)

// 		bdi := fem.boundary_dof_iter_create(system, phi_handle, id)
// 		for bdof in fem.boundary_dof_iter_next(&bdi) {
// 			fem.nli_current(ni)[bdof.gdof] = 100
// 		}
// 	}

// 	context.allocator = infra.ca_allocator(&ca)

// 	for iter in fem.nli_step(&ni) {
// 		infra.ca_check(&ca); defer infra.ca_rewind(&ca)

// 		defer fem.nli_update(&ni)

// 		for element, id in mesh.elements {
// 			infra.ca_check(&ca); defer infra.ca_rewind(&ca)

// 			ls := fem.system_local_problem(system, fem.Entity_ID(id))
// 			defer fem.system_scatter(
// 				system,
// 				fem.Entity_ID(id),
// 				fem.nli_jacobian(ni),
// 				fem.nli_residual(ni),
// 				cm,
// 				ls,
// 			)
// 			// weak form
// 			mapped := fem.map_quadrature(element, {.Interior, .Quad_3, 0})
// 			space := fem.basis_grad_space(mapped, fem.system_var_bd(system, phi_handle), .Scalar)

// 			u_coeffs := fem.system_gather_var_coeffs(
// 				system,
// 				phi_handle,
// 				fem.Entity_ID(id),
// 				fem.nli_current(ni),
// 			)

// 			//log.info(u_coeffs)

// 			for point in 0 ..< fem.space_points(space) {
// 				for test in 0 ..< fem.space_arity(space) {
// 					grad_v := fem.space_gradient(space, mapped, point, test)
// 					for trial in 0 ..< fem.space_arity(space) {
// 						grad_u := fem.space_gradient(space, mapped, point, trial)

// 						k_ij := fem.dot(grad_u, grad_v) * fem.dV(mapped, point)

// 						fem.local_system_mat_add(ls, phi_handle, test, phi_handle, trial, k_ij)
// 						fem.local_system_rhs_add(ls, phi_handle, test, -k_ij * u_coeffs[trial])
// 					}
// 				}
// 			}
// 		}

// 		fem.system_finalize_constraints(system, fem.nli_jacobian(ni), fem.nli_residual(ni), cm)

// 		fem.nli_should_continue(&ni) or_break

// 		fem.sparse_solve(
// 			fem.nli_jacobian(ni),
// 			fem.nli_du(ni),
// 			fem.nli_residual(ni),
// 			tol = fem.nli_linear_tolerance(ni),
// 		)
// 	}

// 	od: serialization.Output_Variable_Data = {system, phi_handle, fem.nli_current(ni)}
// 	of := serialization.output_field_from_system_variable(
// 		fem.Grad_Space(.Scalar),
// 		&od,
// 		"phi_wahoo",
// 	)

// 	viz_mesh := serialization.vtk_create_visualization_mesh(mesh, .Linear)

// 	serialization.write_vtu("output.vtu", mesh, viz_mesh, {of})
// }
