// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
/*
 Element context is the primary mechanism for efficient basis and geometry evaluation.
 An element context consists of a known reference context, and a computed geometry based on the physical element.

 Element context is used for quadrature, post-processing, visualization etc.
*/
package fem

import "core:slice"

Geometry_Option :: enum {
	Physical_Point,
	Jac_Inv,
	Jac_Det,
	Facet_Jac_Inv,
	Facet_Jac_Det,
	Facet_Normal,
}

Geometry_Options :: bit_set[Geometry_Option]

@(private)
Facet_Options :: bit_set[Geometry_Option]{.Facet_Jac_Inv, .Facet_Jac_Det, .Facet_Normal}

Point_Geometry :: struct {
	physical_point:                         Vec3,
	jac, jac_inv, facet_jac, facet_jac_inv: Mat3,
	jac_det, facet_jac_det:                 f64,
	facet_normal:                           Vec3,
}

// defined at startup
Element_Ref_Context :: struct {
	points:               []Vec3,
	facet_indices:        []int,
	basis_lagrange_table: [Order]Lagrange_Context_Table,
}

MAX_REF_CONTEXTS :: 512 // should be plenty
Ref_Context_ID :: distinct uint
REFERENCE_CONTEXTS := [MAX_REF_CONTEXTS]Element_Ref_Context{}

//created per actual element, transient structure.
Element_Context :: struct {
	using ref_context: Element_Ref_Context,
	element:           Mesh_Element,
	geo:               #soa[]Point_Geometry,
	geo_opt:           Geometry_Options,
}

@(private)
register_reference_context :: proc(
	element_type: Element_Type,
	points: []Vec3,
	facet_indices: []int = {},
	allocator := context.allocator,
) -> (
	context_index: Ref_Context_ID,
) {
	@(static) next_id: Ref_Context_ID
	defer next_id += 1

	if next_id == MAX_REF_CONTEXTS {panic("bug: Too many reference contexts.")}

	new_ctx := Element_Ref_Context {
		points               = points,
		facet_indices        = facet_indices,
		basis_lagrange_table = populate_lagrange_context_table(element_type, points, allocator),
	}

	REFERENCE_CONTEXTS[next_id] = new_ctx

	return next_id
}

// Only geometry specified in `geo_opt` will be computed.
element_create_context :: proc(
	element: Mesh_Element,
	ref_context_id: Ref_Context_ID,
	geo_opt: Geometry_Options,
	allocator := context.allocator,
) -> (
	ctx: Element_Context,
) {
	ctx.ref_context = REFERENCE_CONTEXTS[ref_context_id]
	ctx.element = element
	ctx.geo_opt = geo_opt
	ctx.geo = make(#soa[]Point_Geometry, len(ctx.points), allocator)

	if Facet_Options & ctx.geo_opt !=
	   {} {assert(len(ctx.facet_indices) != 0, "Cannot get surface geometry for a reference context defined on interior.")}

	if element.affine && Facet_Options & ctx.geo_opt == {} {
		geo := compute_geometry(ctx, 0)
		for &entry in ctx.geo {entry = geo}
	} else {
		for point in 0 ..< len(ctx.points) {
			geo := compute_geometry(ctx, point)
			ctx.geo[point] = geo
		}
	}

	compute_geometry :: proc(ctx: Element_Context, point: int) -> (geo: Point_Geometry) {
		if ctx.geo_opt == {} {return {}}

		geo.jac = compute_jacobian_context(ctx, point)

		if Facet_Options & ctx.geo_opt != {} {
			geo.facet_jac = compute_facet_jacobian(ctx.element, geo.jac, ctx.facet_indices[point])
		}

		if .Physical_Point in ctx.geo_opt {
			geo.physical_point = compute_physical_point_context(ctx, point)
		}

		if .Jac_Inv in ctx.geo_opt {
			geo.jac_inv = inverse_from_jacobian(element_dim(ctx.element.type), geo.jac)
		}

		if .Jac_Det in ctx.geo_opt {
			geo.jac_det = determinant_from_jacobian(element_dim(ctx.element.type), geo.jac)
		}

		if .Facet_Jac_Inv in ctx.geo_opt {
			geo.facet_jac_inv = inverse_from_jacobian(element_facet_dim(ctx.element.type), geo.facet_jac)
		}

		if .Facet_Jac_Det in ctx.geo_opt {
			geo.facet_jac_det = determinant_from_jacobian(element_facet_dim(ctx.element.type), geo.facet_jac)
		}

		if .Facet_Normal in ctx.geo_opt {
			geo.facet_normal = facet_normal_from_jacobian(
				element_facet_type(ctx.element.type, ctx.facet_indices[point]),
				geo.facet_jac,
				geo.jac,
			)
		}


		return geo
	}

	return ctx
}


ctx_destroy :: proc(ctx: Element_Context, allocator := context.allocator) {
	delete(ctx.geo, allocator)
}

// returns the facet index and a boolean if it was a boundary.
ctx_boundary_point :: proc(ctx: Element_Context, point: int) -> (Boundary_ID, int, bool) {
	assert(len(ctx.facet_indices) != 0, "Only facet-defined reference contexts can be queried for boundary.")
	id, is_bnd := ctx.element.boundaries[ctx.facet_indices[point]].?
	return id, ctx.facet_indices[point], is_bnd
}

ctx_physical_point :: #force_inline proc(ctx: Element_Context, point: int) -> Vec3 {
	assert(.Physical_Point in ctx.geo_opt)
	return ctx.geo.physical_point[point]
}

ctx_jacobian_determinant :: #force_inline proc(ctx: Element_Context, point: int) -> f64 {
	assert(.Jac_Det in ctx.geo_opt)
	return ctx.geo.jac_det[point]
}

ctx_inverse_jacobian :: #force_inline proc(ctx: Element_Context, point: int) -> Mat3 {
	assert(.Jac_Inv in ctx.geo_opt)
	return ctx.geo.jac_inv[point]
}

ctx_facet_inverse_jacobian :: #force_inline proc(ctx: Element_Context, point: int) -> Mat3 {
	assert(.Facet_Jac_Inv in ctx.geo_opt)
	return ctx.geo.facet_jac_inv[point]
}

ctx_facet_jacobian_determinant :: #force_inline proc(ctx: Element_Context, point: int) -> f64 {
	assert(.Facet_Jac_Det in ctx.geo_opt)
	return ctx.geo.facet_jac_det[point]
}

ctx_facet_normal :: #force_inline proc(ctx: Element_Context, point: int) -> Vec3 {
	assert(.Facet_Normal in ctx.geo_opt)
	return ctx.geo.facet_normal[point]
}
