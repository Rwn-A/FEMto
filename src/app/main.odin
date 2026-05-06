package main

import "core:log"
import "core:mem/virtual"
import "core:os"
import "core:path/filepath"
import "core:slice"

import "../../vendor/toml"
import "../fem"
import "../fem/infra"
import fio "../fem/serialization"

import "cfg"
import "conduction"
import le "linear_elasticity"

main :: proc() {
	context.logger = log.create_console_logger()

	// all setup & config memory in here,
	// also passed to models to allocate there top level state
	arena: virtual.Arena
	if err := virtual.arena_init_growing(&arena); err != nil {log.panicf("Allocation Failed. %v", err)}
	defer virtual.arena_destroy(&arena)

	context.allocator = virtual.arena_allocator(&arena)

	cores := os.get_processor_core_count()

	prt: infra.Parallel_Runtime
	infra.parallel_runtime_init(&prt, cores - 2) // cores - 1 total threads.
	defer infra.parallel_runtime_shutdown(&prt)

	if len(os.args) < 2 {
		log.errorf("Expected a path to config file.\nusage: %s my_config.json.", os.args[0])
		os.exit(1)
	}

	config_path := os.args[1]


	fem.setup_default_rules()

	log.info("Loading configuration...")

	schema, schema_ok := cfg.load_schema(config_path)
	if !schema_ok {log.info("Failed to parse configuration."); os.exit(1)}

	general_params, outputter, load_ok := cfg.load_general_params(schema)

	if !load_ok {log.info("Failed to parse configuration."); os.exit(1)}

	switch schema.model {
	case "conduction":
		model_cfg, model_cfg_ok := conduction.load_model_config(schema, general_params)
		if !model_cfg_ok {log.info("Failed to parse configuration."); os.exit(1)}
		if !conduction.run_simulation(general_params, &outputter, model_cfg, &arena, &prt) {os.exit(1)}
	case "elasticity":
		model_cfg, model_cfg_ok := le.load_model_config(schema, general_params)
		if !model_cfg_ok {log.info("Failed to parse configuration."); os.exit(1)}
		if !le.run_simulation(general_params, &outputter, model_cfg, &arena, &prt) {os.exit(1)}
	case:
		log.errorf("config: Unknown model %s", schema.model); os.exit(1)
	}
}
