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
