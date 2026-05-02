/*
 Simple parallel for implementation.

 Does not use work stealing so ideal for similar work loads between threads.
 Main thread does participiate.
*/
package infra

import "base:intrinsics"
import "base:runtime"
import "core:mem/virtual"
import "core:sync"
import "core:thread"

import "core:log"
import "core:testing"
import "core:time"

// Disable multi-threading, no matter how many threads you ask for it always just uses 1.
USE_THREADS :: #config(USE_THREADS, true)

// Correct for most systems, not sure how to get it programtically.
CACHE_LINE_SIZE :: 64

// [start, end)
Range :: struct {
	start, end: int,
}

// Do not copy after initialization
Parallel_Runtime :: struct {
	allocator:       runtime.Allocator,
	spawned_threads: []^thread.Thread, // the threads internal id `tid` is the index into this array, main thread id is 1 + the last index here.
	owner_thread_id: int, // OS thread id, not the internal id.
	total_threads:   int,
	// sync stuff, any below fields are sync primitives or atomic.
	start_barrier:   sync.Barrier,
	done_barrier:    sync.Barrier,
	should_shutdown: bool,
	// state updated per call to parallel_for_*
	work_ranges:     []Range, //indexed by tid
	user_data:       rawptr,
	user_proc:       proc(data: rawptr, range: Range),
}

// Set on initialization, value from 0 -> total_threads - 1, useful for accessing thread specific resources from the data parameter.
// modifying this is a big no-no.
@(thread_local)
tid: int

// Sets up the runtime for subsequent calls, spawns `helper_threads` total threads, cannot copy runtime after this.
parallel_runtime_init :: proc(
	rt: ^Parallel_Runtime,
	helper_threads: int,
	allocator := context.allocator,
	loc := #caller_location,
) {
	assert(helper_threads >= 1, "Must use 1 or more helper threads.", loc)
	context.allocator = allocator

	rt.allocator = allocator
	rt.total_threads = helper_threads + 1 when USE_THREADS else 1
	rt.owner_thread_id = sync.current_thread_id()

	rt.work_ranges = make([]Range, rt.total_threads)

	sync.barrier_init(&rt.start_barrier, rt.total_threads)
	sync.barrier_init(&rt.done_barrier, rt.total_threads)
	when USE_THREADS {
		rt.spawned_threads = make([]^thread.Thread, helper_threads)
		for &spawned_thread, i in rt.spawned_threads {
			spawned_thread = thread.create_and_start_with_poly_data2(rt, i + 1, thread_proc) // tid 0 is main hence + 1
		}
	}

	//main thread, 0 tid, allows for multiple runtimes wiht the same main.
	tid = 0

	thread_proc :: proc(rt: ^Parallel_Runtime, thread_index: int) {
		tid = thread_index
		for {
			sync.barrier_wait(&rt.start_barrier)

			//we might have been woken up solely to exit this loop, hence check first.
			if sync.atomic_load(&rt.should_shutdown) {break}

			//do said work
			work_range := rt.work_ranges[tid]

			//make sure we actually have some work to do.
			if work_range.start < work_range.end {
				rt.user_proc(rt.user_data, work_range)
			}

			//signal we are done
			sync.barrier_wait(&rt.done_barrier)
		}
	}
}

// Returns the part of a range a specifc thread would be tasked with if run with `parallel_for_*.`
parallel_runtime_partition :: proc(prt: ^Parallel_Runtime, range: Range, tid: int, loc := #caller_location) -> Range {
	assert(tid >= 0 && tid < prt.total_threads, "tid is out of bounds.", loc)
	assert(prt_is_initialized(prt), "Runtime must be initialized.", loc)

	work := range.end - range.start
	chunk_size := work / prt.total_threads
	chunk_start := range.start + tid * chunk_size
	chunk_end := range.start + (tid + 1) * chunk_size

	// last tid handles last chunk + remainder.
	if tid == prt.total_threads - 1 {chunk_end = range.end}

	return {chunk_start, chunk_end}
}

parallel_runtime_shutdown :: proc(prt: ^Parallel_Runtime, loc := #caller_location) {
	assert(prt_is_initialized(prt), "Runtime must be initialized.", loc)
	assert(sync.current_thread_id() == prt.owner_thread_id, "Must call from thread that created runtime.", loc)

	context.allocator = prt.allocator
	sync.atomic_store(&prt.should_shutdown, true)

	sync.barrier_wait(&prt.start_barrier)

	thread.join_multiple(..prt.spawned_threads)
	for &helper in prt.spawned_threads {thread.destroy(helper)}

	delete(prt.work_ranges)
	delete(prt.spawned_threads)
}


// Runtime polymorphic parallel for.
// The body procedure should process all entries in the range it's given.
parallel_for_untyped :: proc(
	prt: ^Parallel_Runtime,
	range: Range,
	data: rawptr,
	body: proc(data: rawptr, range: Range),
	loc := #caller_location,
) {
	assert(sync.current_thread_id() == prt.owner_thread_id, "Must call from thread that created runtime.", loc)
	assert(prt_is_initialized(prt), "Runtime must be initialized.", loc)
	prt.user_data = data
	prt.user_proc = body
	for &work_range, tid in prt.work_ranges {work_range = parallel_runtime_partition(prt, range, tid)}

	sync.barrier_wait(&prt.start_barrier)

	main_range := prt.work_ranges[tid]
	if main_range.start < main_range.end {prt.user_proc(prt.user_data, main_range)}

	// signal main is finished and wait for everyone else
	sync.barrier_wait(&prt.done_barrier)
}

// generic wrapper over untyped parallel for, see `parallel_for_untyped` for more docs.
parallel_for :: proc(
	prt: ^Parallel_Runtime,
	range: Range,
	data: $T,
	body: proc(data: T, range: Range),
	loc := #caller_location,
) {
	Wrapped_Data :: struct {
		data: T,
		body: proc(data: T, range: Range),
	}

	wrapped_data := Wrapped_Data{data, body}

	untyped_body :: proc(data: rawptr, range: Range) {
		wrapped_data := cast(^Wrapped_Data)(data)
		wrapped_data.body(wrapped_data.data, range)
	}

	parallel_for_untyped(prt, range, &wrapped_data, untyped_body, loc)
}

in_range :: #force_inline proc(val: int, range: Range) -> bool {return val >= range.start && val < range.end}
prt_is_initialized :: #force_inline proc(prt: ^Parallel_Runtime) -> bool {return prt != nil && prt.total_threads >= 1}
