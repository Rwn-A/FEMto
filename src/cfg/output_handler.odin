/*
 Consumes default output schema, and then can have output fields added to it as config is parsed.
*/
package cfg

import "core:mem"
import "core:os"
import "core:log"
import "core:slice"
import "core:fmt"

import "../fem"
import fio"../fem/serialization"

MAX_OUTPUT_FIELDS :: 32

Output_Handler :: struct {
    path_base: string, // output_directory/sim_name
    formats: bit_set[Output_Format],
    viz_mesh: fio.VTK_Raw_Mesh,
    frequency: int,

    output_fields: [dynamic; MAX_OUTPUT_FIELDS]fio.Output_Field,

    pvd_paths: [dynamic]string,
    pvd_times: [dynamic]f64,

    // used for allocating filenames to go into the pvd's, context allocator in hot loop is usually mass freed.
    allocator: mem.Allocator
}

output_handler_create :: proc(sim_name: string, output: Output_Schema, mesh: fem.Mesh, allocator := context.allocator) -> (Output_Handler, bool) {
    context.allocator = allocator

    if !os.exists(output.directory) {
        if err := os.make_directory(output.directory); err != nil {
    		log.error(os.error_string(err))
    		return {}, false
	   }
    }

    path_base, err := os.join_path({output.directory, sim_name}, allocator)
    if err != nil {
        log.errorf(os.error_string(err)); return {}, false
    }

    return {
        path_base = path_base,
        viz_mesh =  fio.vtk_create_visualization_mesh(mesh, output.order),
        frequency = output.frequency if output.frequency != 0 else 1,
        formats = slice.enum_slice_to_bitset(output.formats, bit_set[Output_Format]),
        pvd_paths = make([dynamic]string),
        pvd_times = make([dynamic]f64),
        allocator = allocator,
    }, true
}

output_handler_destroy :: proc(out_h: ^Output_Handler) {
    context.allocator = out_h.allocator
    fio.vtk_destroy_visualization_mesh(out_h.viz_mesh)
    for &name in out_h.pvd_paths {delete(name)}
    delete(out_h.path_base)
    delete(out_h.pvd_paths)
    delete(out_h.pvd_times)
}

output_handler_add_field :: proc(out_h: ^Output_Handler, field: fio.Output_Field) {
    if append(&out_h.output_fields, field) == 0 {
        log.infof("Too many output fields, increase MAX_OUTPUT_FIELDS current value %d", MAX_OUTPUT_FIELDS)
    }
}

output_handler_step :: proc(
	out_h: ^Output_Handler,
	mesh: fem.Mesh,
	step: fem.Timestep,
	force: bool = false, // useful for ics
) {
	should_output := (step.step % out_h.frequency == 0) || force

	if !should_output {return}

	for format in out_h.formats {
		switch format {
		case .CSV:
			unimplemented("csv")
		case .VTU:
			filename := fmt.aprintf("%s_%d.vtu", out_h.path_base, step.step, allocator = out_h.allocator)
			fio.write_vtu(filename, mesh, out_h.viz_mesh, step.time, out_h.output_fields[:])
			append(&out_h.pvd_paths, filename)
			append(&out_h.pvd_times, step.time)
		}
	}
}


// writes the current pvd, does not clear paths/times
output_handler_flush_pvd :: proc(out_h: ^Output_Handler) {
    context.allocator = out_h.allocator

    if .VTU not_in out_h.formats {return}

    pvd_path := fmt.aprintf("%s.pvd", out_h.path_base)
    defer delete(pvd_path)

    fio.write_pvd(pvd_path, out_h.pvd_paths[:], out_h.pvd_times[:])
}