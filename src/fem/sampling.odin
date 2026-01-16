// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
/*
  Allows for iteration and lazy geometry computation over some known sets of reference points.
  A set of points may be defined on the facet of an element, special facet-aware iteration is in place.
  A set of points may be weighted, this is used in quadrature but possibly in other sets down the road.

  At each specific point, called a sample, a sample context is bound this lazily computes geometry.

  A sample object can also be used to query other data that is precomputed/cached at sample points.
  Perhaps its basis values, or connectivity information for output meshes, or whatever else.
*/
package fem

import "core:math"

Sample_Group :: enum {
	Integration_1,
	Integration_3,
	Facet_Integration_1,
	Facet_Integration_3,
	Visualization_Linear,
	Visualization_Quadratic,
}

Sample :: struct {
	group: Sample_Group,
	index: int,
}

Facet_Sample :: struct {
	using interior_sample: Sample,
	facet_index:           int,
}

Facet_Groups :: bit_set[Sample_Group]{.Facet_Integration_1, .Facet_Integration_3}


@(private = "file")
Interior_Geometry :: struct {
	ready:              bit_set[enum {
		Jacobian,
		Scale,
		Pullback,
	}],
	point:              Vec3,
	jacobian, pullback: Mat3,
	scale:              f64,
}

@(private = "file")
Facet_Geometry :: struct {
	ready:              bit_set[enum {
		Jacobian,
		Scale,
		Normal,
		Pullback,
	}],
	normal:             Vec3,
	jacobian, pullback: Mat3,
	scale:              f64,
}

Sample_Context :: struct {
	element:       Element,
	element_id:    Entity_ID,
	interior_geom: Interior_Geometry,
	facet_geom:    [MAX_FACETS]Facet_Geometry,
	s:             Sample,
	fs:            Facet_Sample,
	facet_valid:   bool,
}

Sample_Iterator :: struct {
	current_index: int,
	current_facet: int,
}

Facet_Filter :: bit_set[0 ..< MAX_FACETS]

// Iterating sample groups

sample_context_create :: proc(mesh: Mesh, element_id: Entity_ID) -> (Sample_Iterator, Sample_Context) {
	return {}, Sample_Context{element = mesh.elements[element_id], element_id = element_id}
}

iterate_sample_group :: proc(si: ^Sample_Iterator, ctx: ^Sample_Context, group: Sample_Group) -> (more: bool) {
	if si.current_index >= len(SAMPLE_POINTS[ctx.element.type][group]) {return false}

	s := Sample {
		group = group,
		index = si.current_index,
	}

	ctx.s = s
	ctx.facet_valid = false

	if !ctx.element.affine {
		ctx.interior_geom = {}
		ctx.facet_geom = {}
	}

	si.current_index += 1

	return true
}

// filters OUT interface facets, leaving only boundary facets.
create_interface_facet_filter :: proc(element: Mesh_Element) -> (f: Facet_Filter) {
	for adj, facet_index in element.adjacency {
		if _, is_bnd := adj.(Boundary_ID); is_bnd {f += {facet_index}}
	}
	return ~f
}

iterate_facet_sample_group :: proc(
	si: ^Sample_Iterator,
	ctx: ^Sample_Context,
	group: Sample_Group,
	facet_filter: Facet_Filter = {},
) -> (
	more: bool,
) {
	assert(group in Facet_Groups)

	for si.current_index < len(SAMPLE_POINTS[ctx.element.type][group]) {
		facet_index := SAMPLE_FACET_INDICES[ctx.element.type][group][si.current_index]

		if facet_index != si.current_facet {si.current_facet = facet_index}
		if si.current_facet not_in facet_filter {break}

		si.current_index += 1
	}

	// completely done
	if si.current_index >= len(SAMPLE_POINTS[ctx.element.type][group]) {return false}

	fs := Facet_Sample {
		group       = group,
		index       = si.current_index,
		facet_index = si.current_facet,
	}

	ctx.fs = fs
	ctx.s = ctx.fs.interior_sample
	ctx.facet_valid = true

	if !ctx.element.affine {
		ctx.interior_geom = {}
		ctx.facet_geom = {}
	}

	si.current_index += 1

	return true
}

sample_iterator_reset :: proc(si: ^Sample_Iterator) {
	si.current_index = 0
	si.current_facet = 0
}

// Returned point is in reference space
sample_point :: #force_inline proc(element_type: Element_Type, s: Sample) -> Vec3 {
	return SAMPLE_POINTS[element_type][s.group][s.index]
}


// Lazy geometry computation

sample_physical_point :: proc(elg: ^Sample_Context) -> Vec3 {
	return compute_physical_point_sample(elg.element, elg.s)
}

sample_jacobian :: proc(ctx: ^Sample_Context) -> Mat3 {
	geom := &ctx.interior_geom
	if .Jacobian in geom.ready {return geom.jacobian}
	geom.ready += {.Jacobian}
	geom.jacobian = compute_jacobian_sample(ctx.element, ctx.s)
	return geom.jacobian
}

sample_facet_jacobian :: proc(ctx: ^Sample_Context) -> Mat3 {
	assert(ctx.facet_valid)
	geom := &ctx.facet_geom[ctx.fs.facet_index]
	if .Jacobian in geom.ready {return geom.jacobian}
	geom.ready += {.Jacobian}
	geom.jacobian = compute_facet_jacobian(ctx.element, sample_jacobian(ctx), ctx.fs.facet_index)
	return geom.jacobian
}

sample_scale :: proc(ctx: ^Sample_Context) -> f64 {
	geom := &ctx.interior_geom
	if .Scale in geom.ready {return geom.scale}
	geom.ready += {.Scale}
	geom.scale = measure_from_jacobian(element_dim(ctx.element.type), sample_jacobian(ctx))
	return geom.scale
}

sample_facet_scale :: proc(ctx: ^Sample_Context) -> f64 {
	assert(ctx.facet_valid)
	geom := &ctx.facet_geom[ctx.fs.facet_index]
	if .Scale in geom.ready {return geom.scale}
	geom.ready += {.Scale}
	geom.scale = measure_from_jacobian(element_facet_dim(ctx.element.type), sample_facet_jacobian(ctx))
	return geom.scale
}

sample_pullback :: proc(ctx: ^Sample_Context) -> Mat3 {
	geom := &ctx.interior_geom
	if .Pullback in geom.ready {return geom.pullback}
	geom.ready += {.Pullback}
	geom.pullback = pullback_from_jacobian(element_dim(ctx.element.type), sample_jacobian(ctx))
	return geom.pullback
}

sample_facet_pullback :: proc(ctx: ^Sample_Context) -> Mat3 {
	assert(ctx.facet_valid)
	geom := &ctx.facet_geom[ctx.fs.facet_index]
	if .Pullback in geom.ready {return geom.pullback}
	geom.ready += {.Pullback}
	geom.pullback = pullback_from_jacobian(element_facet_dim(ctx.element.type), sample_facet_jacobian(ctx))
	return geom.pullback
}

sample_facet_normal :: proc(ctx: ^Sample_Context) -> Vec3 {
	assert(ctx.facet_valid)
	geom := &ctx.facet_geom[ctx.fs.facet_index]
	if .Normal in geom.ready {return geom.normal}
	geom.ready += {.Normal}
	geom.normal = compute_facet_normal(
		element_facet_type(ctx.element.type, ctx.fs.facet_index),
		sample_facet_jacobian(ctx),
		sample_jacobian(ctx),
	)
	return geom.normal
}


// Tables

@(rodata)
SAMPLE_POINTS := [Element_Type][Sample_Group][]Vec3 {
	.Point = {
		.Integration_1 = POINT_QUADRATURE[.Integration_1],
		.Integration_3 = POINT_QUADRATURE[.Integration_3],
		.Facet_Integration_1 = POINT_QUADRATURE[.Facet_Integration_1],
		.Facet_Integration_3 = POINT_QUADRATURE[.Facet_Integration_3],
		.Visualization_Linear = POINT_VISUALIZATION[.Visualization_Linear],
		.Visualization_Quadratic = POINT_VISUALIZATION[.Visualization_Quadratic],
	},
	.Line = {
		.Integration_1 = LINE_QUADRATURE[.Integration_1],
		.Integration_3 = LINE_QUADRATURE[.Integration_3],
		.Facet_Integration_1 = LINE_QUADRATURE[.Facet_Integration_1],
		.Facet_Integration_3 = LINE_QUADRATURE[.Facet_Integration_3],
		.Visualization_Linear = LINE_VISUALIZATION[.Visualization_Linear],
		.Visualization_Quadratic = LINE_VISUALIZATION[.Visualization_Quadratic],
	},
	.Triangle = {
		.Integration_1 = TRIANGLE_QUADRATURE[.Integration_1],
		.Integration_3 = TRIANGLE_QUADRATURE[.Integration_3],
		.Facet_Integration_1 = TRIANGLE_QUADRATURE[.Facet_Integration_1],
		.Facet_Integration_3 = TRIANGLE_QUADRATURE[.Facet_Integration_3],
		.Visualization_Linear = TRIANGLE_VISUALIZATION[.Visualization_Linear],
		.Visualization_Quadratic = TRIANGLE_VISUALIZATION[.Visualization_Quadratic],
	},
	.Quadrilateral = {
		.Integration_1 = QUADRILATERAL_QUADRATURE[.Integration_1],
		.Integration_3 = QUADRILATERAL_QUADRATURE[.Integration_3],
		.Facet_Integration_1 = QUADRILATERAL_QUADRATURE[.Facet_Integration_1],
		.Facet_Integration_3 = QUADRILATERAL_QUADRATURE[.Facet_Integration_3],
		.Visualization_Linear = QUADRILATERAL_VISUALIZATION[.Visualization_Linear],
		.Visualization_Quadratic = QUADRILATERAL_VISUALIZATION[.Visualization_Quadratic],
	},
	.Tetrahedron = {
		.Integration_1 = TETRAHEDRON_QUADRATURE[.Integration_1],
		.Integration_3 = TETRAHEDRON_QUADRATURE[.Integration_3],
		.Facet_Integration_1 = TETRAHEDRON_QUADRATURE[.Facet_Integration_1],
		.Facet_Integration_3 = TETRAHEDRON_QUADRATURE[.Facet_Integration_3],
		.Visualization_Linear = TETRAHEDRON_VISUALIZATION[.Visualization_Linear],
		.Visualization_Quadratic = TETRAHEDRON_VISUALIZATION[.Visualization_Quadratic],
	},
}

// valid for facet groups only
@(rodata)
SAMPLE_FACET_INDICES := #partial [Element_Type][Sample_Group][]int {
	.Point = {},
	.Line = #partial{.Facet_Integration_1 = {0, 1}, .Facet_Integration_3 = {0, 1}},
	.Quadrilateral = #partial{.Facet_Integration_1 = {0, 1, 2, 3}, .Facet_Integration_3 = {0, 0, 1, 1, 2, 2, 3, 3}},
	.Triangle = #partial{.Facet_Integration_1 = {0, 1, 2}, .Facet_Integration_3 = {0, 0, 1, 1, 2, 2}},
	.Tetrahedron = #partial{.Facet_Integration_1 = {0, 1, 2, 3}, .Facet_Integration_3 = {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3}},
}
