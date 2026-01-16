// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
package main

import "core:log"
import "core:math/linalg"
import "core:mem/virtual"

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

	mesh, warn, err := fio.gmsh_parse("./meshes/box.msh")

	if err != nil {log.error(err); return}
	if warn != nil {log.warn(warn)}

	fem.populate_basis_tables()

	viz_mesh := fio.vtk_create_visualization_mesh(mesh, .Visualization_Quadratic)

	// config start
	temperature := models.lagrange_scalar_field(mesh, .Linear)

	constraints := make(map[fem.Boundary_ID]models.Constraint)
	weak := make(map[fem.Boundary_ID]models.Variational_BC)
	iso_wall := mesh.boundary_names["left"]
	free := mesh.boundary_names["bottom"]
	sd := mesh.boundary_names["top"]
	constraints[iso_wall] = models.isothermal(273)
	weak[free] = models.fixed_flux(10)
	weak[sd] = models.fixed_flux(-10)
	models.field_discover_constraints(mesh, temperature, constraints)

	models.conduction(mesh, temperature, 1.0, weak)

	fio.write_vtu("output.vtu", mesh, viz_mesh, {models.lagrange_scalar_field_visualizer(&temperature, "temperature")})
}
