// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
package models

import "../fem"
import "../la"

beam :: proc(mesh: fem.Mesh, allocator := context.allocator) {
	assert(mesh.dim == .D1)
	context.allocator = allocator
}

membrane :: proc(mesh: fem.Mesh, allocator := context.allocator) {
	assert(mesh.dim == .D2)
	context.allocator = allocator
}

solid :: proc(mesh: fem.Mesh, allocator := context.allocator) {
	assert(mesh.dim == .D3)
	context.allocator = allocator
}
