package conduction

import "core:log"
import "core:math"
import "core:slice"

import "../../fem"
import "../cfg"

Model_Config :: struct {
	params:      Model_Parameters,
	bd:          fem.Basis_Descriptor,
	time_scheme: fem.Time_Scheme,
	ics:         map[fem.Section_ID]Initial_Condition_Int,
}

SUPPORTED_FAMILIES := bit_set[fem.Basis_Family]{.Lagrange}
SUPPORTED_TIME_SCHEMES := bit_set[fem.Time_Scheme]{.BE, .BDF2}

Plug_Context :: cfg.Plugin_Context(
	Material_Int,
	Source_Int,
	BC_Int,
	Initial_Condition_Int,
)

load_model_config :: proc(
	schema: cfg.Schema,
	gcfg: cfg.General_Params,
	allocator := context.allocator
) -> (model_cfg: Model_Config, ok: bool) {
	context.allocator = allocator

	//load plugins
	if len(schema.plugins) > 0 { log.info("Loading user plugins...") }

	plug_ctx: Plug_Context
	cfg.init_plugin_context(&plug_ctx)

	builtin_register(&plug_ctx)

	for plugin in schema.plugins{
		// not solving pure conduction so plugins model must be conduction or we skip it.
		if schema.model != "conduction" && plugin.model != "conduction" {continue}
		if !cfg.load_plugin(plugin.path, &plug_ctx ) { return {}, false }
	}

	if len(schema.plugins) > 0 { log.info("Plugins loaded...") }


	var, iso := cfg.load_bcs(BC_Int, Isothermal_Int, Variational_Int, schema, gcfg.mesh, plug_ctx.bcs) or_return
	model_cfg.params.variational_bcs = var
	model_cfg.params.isothermal_bcs = iso
	model_cfg.params.materials = cfg.load_material(Material_Int, schema, gcfg.mesh, plug_ctx.materials) or_return
	model_cfg.params.sources = cfg.load_sources(Source_Int, schema, gcfg.mesh, plug_ctx.sources) or_return

	field_schema := struct {
		order:              fem.Basis_Order,
		basis_family:       fem.Basis_Family,
		time_scheme:        fem.Time_Scheme,
		initial_conditions: []^cfg.Table,
	} {
		order        = .Linear,
		basis_family = .Lagrange,
		time_scheme  = .BE,
	}

	// optional field config
	temperature_table, exists := cfg.table_get_opt(^cfg.Table, schema.fields, "temperature")

	if exists {
		cfg.unmarshal(temperature_table, &field_schema) or_return

		cfg.valid_enum_option(SUPPORTED_FAMILIES, field_schema.basis_family) or_return
		cfg.valid_enum_option(SUPPORTED_TIME_SCHEMES, field_schema.time_scheme) or_return

		model_cfg.ics = cfg.load_ics(Initial_Condition_Int, field_schema.initial_conditions, gcfg.mesh, plug_ctx.ics) or_return

	}

	model_cfg.bd.order = field_schema.order
	model_cfg.bd.family = field_schema.basis_family
	model_cfg.time_scheme = field_schema.time_scheme


	return model_cfg, true
}
