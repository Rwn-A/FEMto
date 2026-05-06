package cfg

import "core:os"
import "core:log"
import "core:slice"
import "core:mem"
import "core:fmt"
import "core:mem/virtual"

import "../../fem"
import fio"../../fem/serialization"

Outputter :: struct {
    viz_mesh:  fio.VTK_Raw_Mesh,
	sim_name, directory:  string,
	frequency: int,
	is_transient: bool,
	formats: bit_set[Output_Format],
	pvd_times: [dynamic]f64,
	pvd_names: [dynamic]string,
	allocator: mem.Allocator,
}

outputter_create :: proc(schema: Schema, mesh: fem.Mesh, allocator := context.allocator) -> (output: Outputter, ok: bool) {
    context.allocator = allocator
    output.allocator = allocator

	output.viz_mesh = fio.vtk_create_visualization_mesh(mesh, schema.output.order)

	output.directory = schema.output.directory
	output.frequency = schema.output.frequency
	output.sim_name = schema.sim_name

	for format in schema.output.formats {output.formats += {format}}
	output.formats = slice.enum_slice_to_bitset(schema.output.formats, bit_set[Output_Format])

	output.is_transient = !(schema.time_control == {})

	output.pvd_names = make([dynamic]string)
	output.pvd_times = make([dynamic]f64)

	if os.exists(output.directory) {return output, true}

	if err := os.make_directory(output.directory); err != nil {
		log.error(os.error_string(err))
		return {}, false
	}

	return output, true
}

output_destroy :: proc(outputter: ^Outputter) {
    delete(outputter.pvd_times)
    delete(outputter.pvd_names)
}

output_step :: proc(outputter: ^Outputter, mesh: fem.Mesh, step: fem.Timestep, ofs: []fio.Output_Field, force: bool = false) {
    should_output := (step.step % outputter.frequency == 0) || !outputter.is_transient || force

    if !should_output {return}

    for format in outputter.formats{
        switch format{
            case .CSV: unimplemented("csv")
            case .VTU:
                 filename := get_filename(outputter^, step.step, "vtu")
                 fio.write_vtu(filename, mesh, outputter.viz_mesh, step.time, ofs)
                 if outputter.is_transient {
                    append(&outputter.pvd_names, filename)
		            append(&outputter.pvd_times, step.time)
                 }
        }
    }
}


// Writes the pvd file in the case of a transient simulation with VTU format.
// Clears the internal pvd arrays.
output_flush :: proc(outputter: ^Outputter) {
    defer {
        clear(&outputter.pvd_names)
        clear(&outputter.pvd_times)
    }

    if .VTU in outputter.formats && outputter.is_transient {
        filename := get_filename(outputter^, 0, "pvd")
        fio.write_pvd(filename, outputter.pvd_names[:], outputter.pvd_times[:])
    }
}


@(private)
get_filename :: proc(outputter: Outputter, step: int, extension: string) -> string {
    context.allocator = outputter.allocator

    filename: string
    switch extension {
        case "pvd": filename = fmt.aprintf("%s.pvd", outputter.sim_name)
        case:  filename = fmt.aprintf("%s_%d.%s", outputter.sim_name, step, extension)
    }

	path, _ := os.join_path({outputter.directory, filename}, context.allocator)
	return path
}
