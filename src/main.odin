package main

import "core:log"
import "core:os"
import "core:mem/virtual"

import "fem"
import "cfg"

import "fem/infra"
import "models/conduction"
import "models/small_strain"
import "models/inc_flow"

main :: proc() {
    context.logger = log.create_console_logger()

	arena: virtual.Arena
	if err := virtual.arena_init_growing(&arena); err != nil {log.panicf("Allocation Failed. %v", err)}
	defer virtual.arena_destroy(&arena)

	context.allocator = virtual.arena_allocator(&arena)

	if !run() {
	   os.exit(1)
	 }
}

run :: proc() -> bool {
	cores := os.get_processor_core_count()

	prt: infra.Parallel_Runtime
	infra.parallel_runtime_init(&prt, cores - 2) // cores - 1 total threads.
	defer infra.parallel_runtime_shutdown(&prt)

	if len(os.args) < 2 {
		log.errorf("Expected a path to config file.\nusage: %s my_config.json.", os.args[0])
		return false
	}

	config_path := os.args[1]

	fem.setup_default_rules()

	log.info("Loading configuration...")

	schema := cfg.load_root(config_path) or_return
    mesh := cfg.load_mesh(schema.mesh) or_return
    out_h := cfg.output_handler_create(schema.sim_name, schema.output, mesh) or_return
    plug_libs, plugins := cfg.load_plugins(schema.plugins) or_return
    defer cfg.unload_plugins(plug_libs)

    switch schema.model {
        case "conduction":
            cd := conduction.configure_driver(mesh, &out_h, &plugins, schema) or_return
            conduction.drive(cd, mesh, &prt)
        case "small_strain":
            ssd := small_strain.configure_driver(mesh, &out_h, &plugins, schema) or_return
            small_strain.drive(ssd, mesh, &prt)
        case "inc_flow":
            driver := inc_flow.configure_driver(mesh, &out_h, &plugins, schema) or_return
            inc_flow.drive(driver, mesh, &prt)
        case:
            log.errorf("Unknown model %s", schema.model); return false
    }


    return true

}