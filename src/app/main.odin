// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
package main

import "core:log"
import "core:math/linalg"
import "core:mem/virtual"
import "core:time"

import "../fem"
import fio "../fem/serialization"
import "../la"
import "../models"

main :: proc() {
	context.logger = log.create_console_logger()

	arena: virtual.Arena
	if err := virtual.arena_init_growing(&arena); err != nil {log.panic("Allocation failure")}
	defer virtual.arena_destroy(&arena)
	context.allocator = virtual.arena_allocator(&arena)

	log.info("loading mesh...")
	mesh, warn, err := fio.gmsh_parse("./meshes/part.msh")

	if err != nil {log.error(err); return}
	if warn != nil {log.warn(warn)}

	fem.populate_basis_tables()

	viz_mesh := fio.vtk_create_visualization_mesh(mesh, .Visualization_Linear)


	iso_wall := mesh.boundary_names["left"]
	flux_wall := mesh.boundary_names["right"]

	are_constrained := make(map[fem.Boundary_ID]struct{})
	are_constrained[iso_wall] = {}

	// config start
	temperature := models.lagrange_scalar_field(mesh, .Quadratic, are_constrained)
	log.infof("solving %d dofs...", len(temperature.coeffs))

	constraints := make(map[fem.Boundary_ID]models.Constraint)
	weak := make(map[fem.Boundary_ID]models.Variational_BC)

	constraints[iso_wall] = models.isothermal(50)
  	weak[flux_wall] = models.fixed_flux(10)
  	models.field_set_constraints(mesh, temperature, constraints)

  	stiffness_sparsity := fem.create_sparsity_from_layout(mesh, temperature, temperature)
	stiffness := la.sparse_mat_from_sparsity(stiffness_sparsity)
	load := make(la.Vec, len(temperature.coeffs))

	for element, id in mesh.elements {
		defer free_all(context.temp_allocator)
		k, f := models.conduction(mesh, fem.Entity_ID(id), temperature, 10.0, weak)
		fem.scatter_local(fem.Entity_ID(id), stiffness, load, k, f, temperature, temperature)
	}

	fem.system_finalize_constraints(stiffness, load, temperature.layout)

	solve_time := time.now()
	iter, resid, converged := la.sparse_cg_solve(stiffness, temperature.coeffs, load)
	log.infof("solve took: %v", time.since(solve_time))

	log.infof("Sparse solve: converged? %v, residual: %v, iterations: %v", converged, resid, iter)

	fio.write_vtu("output.vtu", mesh, viz_mesh, {models.lagrange_scalar_field_visualizer(&temperature, "temperature")})
}
