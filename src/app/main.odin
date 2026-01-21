// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
package main

import "core:log"
import "core:math/linalg"
import "core:mem/virtual"
import "core:time"
import "core:encoding/json"
import "core:os/os2"

import "../fem"
import fio "../fem/serialization"
import "../la"
import "../models"

Property :: f64

Schema :: struct {
    mesh: string,
    model: enum{Conduction, Elastic_Beam},
    output: struct {
    	path: string,
    },
    sections: map[string]Section_Info,
    fields: map[string]Field_Info,

}

Section_Info :: map[string]Property

Field_Info :: struct {
    order: fem.Order,
    boundaries: map[string]Boundary_Info,
}

Boundary_Info :: struct{
    type: string,
    value: Property,
}


main :: proc() {
	context.logger = log.create_console_logger()

	arena: virtual.Arena
	if err := virtual.arena_init_growing(&arena); err != nil {log.panic("Allocation failure")}
	defer virtual.arena_destroy(&arena)
	context.allocator = virtual.arena_allocator(&arena)

	fem.populate_basis_tables()

	if len(os2.args) < 2 {
		log.error("No config file provided.")
		os2.exit(1)
	}

	config_path := os2.args[1]

	json_data, open_err := os2.read_entire_file_from_path(config_path, context.allocator)

	if open_err != nil {
		log.error(os2.error_string(open_err))
		os2.exit(1)
	}

	config: Schema
	json.unmarshal(json_data, &config, .MJSON)

	log.info("Loading config...")
	mesh, warn, err := fio.gmsh_parse(config.mesh)

	if err != nil {log.error(err); return}
	if warn != nil {log.warn(warn)}

	switch config.model {
	case .Conduction: run_conduction(mesh, config)
	case .Elastic_Beam: unimplemented("beams")
	case: log.error("Unknown model encountered.")
	}
}

run_conduction :: proc(mesh: fem.Mesh, c: Schema) {
	temperature_config, exists := c.fields["temperature"]
	if !exists {
		log.error("Temperature field was not defined on config.")
		return
	}

	constraints := make(map[fem.Boundary_ID]models.Constraint)
	weak := make(map[fem.Boundary_ID]models.Variational_BC)
	mat_props := make(map[fem.Section_ID]models.Material_Property)

	for boundary_name, boundary_id in mesh.boundary_names {
		boundary, exists := temperature_config.boundaries[boundary_name]
		if !exists {
			log.errorf("Boundary %s was not defined on temperature", boundary_name)
			return
		}

		switch boundary.type {
		case "isothermal": constraints[boundary_id] = models.isothermal(boundary.value)
		case "adiabatic": weak[boundary_id] = models.fixed_flux(0)
		case "flux": weak[boundary_id] = models.fixed_flux(boundary.value)
		case:
			log.errorf("%s is not a valid boundary for conduction", boundary.type)
		}
	}

	for section_name, section_id in mesh.section_names {
		section, exists := c.sections[section_name]
		if !exists {
			log.errorf("Section %s was not defined in config", section_name)
			return
		}
		if conductivity, exists := section["thermal_conductivity"]; exists {
			mat_props[section_id] = models.property_conductivity(conductivity)
		}else{
			log.errorf("Expected thermal_conductivity defined in section %s", section_name)
			return
		}
	}
	temperature := models.lagrange_scalar_field(mesh, temperature_config.order, constraints)

	log.infof("solving %d degrees of freedom", len(temperature.coeffs))

	stiffness_sparsity := fem.create_sparsity_from_layout(mesh, temperature, temperature)
	stiffness := la.sparse_mat_from_sparsity(stiffness_sparsity)
	load := make(la.Vec, len(temperature.coeffs))

	for element, id in mesh.elements {
		defer free_all(context.temp_allocator)
		k, f := models.conduction(mesh, fem.Entity_ID(id), temperature, mat_props[element.section], weak, context.temp_allocator)
		fem.scatter_local(fem.Entity_ID(id), stiffness, load, k, f, temperature, temperature)
	}

	solve_time := time.now()
	iter, resid, converged := la.sparse_cg_solve(stiffness, temperature.coeffs, load)
	log.infof("solve took: %v", time.since(solve_time))

	viz_mesh := fio.vtk_create_visualization_mesh(mesh, .Visualization_Linear)

	fio.write_vtu(c.output.path, mesh, viz_mesh, {models.lagrange_scalar_field_visualizer(&temperature, "temperature")})
}