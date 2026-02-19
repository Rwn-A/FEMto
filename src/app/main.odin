package main

import "core:log"
import "core:mem/virtual"

import "core:math/linalg"

import fem "../fe_core"
import fio "../serialization"
import models "../models"
import "../la"

main :: proc() {
	context.logger = log.create_console_logger()

	arena: virtual.Arena
	if err := virtual.arena_init_growing(&arena); err != nil {log.panic("Allocation failure")}
	defer virtual.arena_destroy(&arena)
	context.allocator = virtual.arena_allocator(&arena)

	fem.setup_quadrature_rules()
	fem.setup_subcell_rules()

	mesh, warn, err := fio.gmsh_parse("./meshes/line.msh")

	isotherm := mesh.boundary_names["Gamma_left"]
	iso2 := mesh.boundary_names["Gamma_right"]

	constrained := make(map[fem.Boundary_ID]f64)

	constrained[isotherm] = 273
	constrained[iso2] = 400


    model_params := models.Conduction_Params{ soln_order = .Linear, isothermal_bnds = constrained }

    model := models.model_conduction(&model_params)

    solve_config := models.Solver_Config{
        tc = {},
        output = {
            output_dir = "./output",
            frequency = 1,
            prefix = "test_sim",
            output_rule = .Split_Linear,
        },

        linsolve_rel_tol = 1e-7,
        linsolve_max_iter = 1000,
        non_linear_rel_tol = 1e-5,
        non_linear_max_iter = 12,
    }

    res := models.solve_model(solve_config, &model, mesh)

    if res.converged == false {
        log.error(res)
    }else{
        log.info("Simulation complete!")
    }

}
