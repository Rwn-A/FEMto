package main

import "core:log"
import "core:mem/virtual"
import "core:os"
import "core:path/filepath"
import "core:slice"

import "../../vendor/toml"
import "../fem"
import "../fem/infra"
import "../fem/serialization"

import "conduction"
import le"linear_elasticity"
import "config"

main :: proc() {
	context.logger = log.create_console_logger()

	arena: virtual.Arena

	if err := virtual.arena_init_growing(&arena); err != nil {
		log.panicf("Allocation Failed. %v", err)
	}
	defer virtual.arena_destroy(&arena)

	context.allocator = virtual.arena_allocator(&arena)

	if len(os.args) < 2 {
		log.errorf("Expected a path to config file.\nusage: %s my_config.json.", os.args[0])
		os.exit(1)
	}

	config_path := os.args[1]

	fem.setup_default_rules()

	log.info("Loading configuration...")

	data, err := os.read_entire_file_from_path(config_path, context.allocator)

	if err != nil {
		log.fatalf(os.error_string(err))
		os.exit(1)

	}

	config_dir := filepath.dir(config_path)
	os.chdir(config_dir)

	root, toml_err := toml.parse(string(data), config_path, context.temp_allocator)
	defer free_all(context.temp_allocator)

	if toml_err.type != .None {
		log.fatal(toml.format_error(&toml_err))
		os.exit(1)
	}
	schema := config.DEFAULT_SCHEMA
	schema.sources = new(toml.List)

	if !config.unmarshal(root, &schema, context.temp_allocator) {
		os.exit(1)
	}

	cfg, parse_ok := config.load_general_config(schema, &arena)

	if !parse_ok {
		log.info("Failed to parse configuration.")
		os.exit(1)
	}


	cores := os.get_processor_core_count()

	prt: infra.Parallel_Runtime
	infra.parallel_runtime_init(&prt, cores - 2) // with main thread that means we uses cores - 1 total threads.
	defer infra.parallel_runtime_shutdown(&prt)

	switch schema.model {
	case "conduction":
		model_cfg, model_cfg_ok := conduction.load_model_config(schema, cfg)
		if !model_cfg_ok {
			log.info("Failed to parse configuration.")
			os.exit(1)
		}
		//free_all(context.temp_allocator)
		if !conduction.run_simulation(cfg, model_cfg, &arena, &prt) {
			os.exit(1)
		}
	case "elasticity":
		model_cfg, model_cfg_ok := le.load_model_config(schema, cfg)
		if !model_cfg_ok {
			log.info("Failed to parse configuration.")
			os.exit(1)
		}
		free_all(context.temp_allocator)
		if !le.run_simulation(cfg, model_cfg, &arena, &prt) {
			os.exit(1)
		}
	case:
		log.errorf("Unknown model %s", schema.model)
		os.exit(1)
	}


	// //schema, schema_ok := config.read_config(config_path, &arena)
	// schema, schema_ok := config.read_toml_config(config_path, &arena)

	// if !schema_ok {
	// 	log.info("Failed to load configuration file.")
	// 	os.exit(1)
	// }

	// cfg, parse_ok := config.load_general_config(schema, &arena)

	// if !parse_ok {
	// 	log.info("Failed to parse configuration.")
	// 	os.exit(1)
	// }

	// if schema.model == "" {
	// 	log.error("Expected a key `model` in configuration.")
	// 	os.exit(1)
	// }

	// prt: infra.Parallel_Runtime
	// infra.parallel_runtime_init(&prt, 4)
	// defer infra.parallel_runtime_shutdown(&prt)

	// model_params := config.Model_Schema{schema.linear_solver, schema.fields, schema.sections}

	// switch schema.model {
	// case "conduction":
	// 	if !conduction.run_simulation(cfg, model_params, &arena, &prt) {
	// 		os.exit(1)
	// 	}
	// case:
	// 	log.errorf("Unknown model %s", schema.model)
	// 	os.exit(1)
	// }


}
