// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
package serialization

import fem"../fe_core"

MAX_OUTPUT_FIELD_COMPONENTS :: 9

Output_Proc :: #type proc(ctx: fem.Element_Context, point: int, data: rawptr) -> [MAX_OUTPUT_FIELD_COMPONENTS]f64

Output_Field :: struct {
	friendly_name:  string,
	components:     int,
	data:           rawptr,
	value_provider: Output_Proc,
	geometry_required: fem.Geometry_Options,
}
